# Job Cancellation Documentation

## Overview

The Kanad backend API now supports comprehensive job cancellation for experiments running locally or on cloud quantum backends (IBM Quantum, BlueQubit). This feature allows users to gracefully stop long-running experiments and free up computational resources.

## Features

- **Local Job Cancellation**: Stop VQE/SQD optimization loops running on local classical simulators
- **IBM Quantum Cancellation**: Cancel jobs submitted to IBM Quantum hardware/cloud
- **BlueQubit Cancellation**: Cancel jobs running on BlueQubit GPU simulators
- **Graceful Shutdown**: Partial results are saved with "cancelled" status
- **Thread-Safe**: Uses threading events for safe cancellation across workers
- **Status Tracking**: Experiments track cancellation timestamps and status

## API Endpoint

### Cancel Experiment

**Endpoint**: `PATCH /api/v1/experiments/{experiment_id}/cancel`

**Description**: Cancel a running or queued experiment.

**Path Parameters**:
- `experiment_id` (integer): ID of the experiment to cancel

**Response** (200 OK):
```json
{
  "status": "cancelled",
  "experiment_id": 123,
  "local_cancelled": true,
  "cloud_cancelled": false,
  "cloud_job_id": null,
  "message": "Experiment cancellation requested successfully",
  "warning": null
}
```

**Error Responses**:
- `404 Not Found`: Experiment ID does not exist
- `400 Bad Request`: Experiment is not in a cancellable state (already completed, failed, or cancelled)

**Example cURL**:
```bash
# Cancel experiment with ID 5
curl -X PATCH http://localhost:8000/api/v1/experiments/5/cancel

# Pretty print response
curl -X PATCH http://localhost:8000/api/v1/experiments/5/cancel | jq
```

**Example Python**:
```python
import requests

# Cancel experiment
response = requests.patch('http://localhost:8000/api/v1/experiments/5/cancel')
result = response.json()

print(f"Status: {result['status']}")
print(f"Local cancelled: {result['local_cancelled']}")
print(f"Cloud cancelled: {result['cloud_cancelled']}")
```

## Cancellation Behavior

### For Local Jobs (VQE/SQD on Classical Simulators)

1. **During Queue Wait**: If experiment is queued but not yet running, it's marked as cancelled immediately when removed from queue
2. **During Optimization**: Cancellation flag is checked at every VQE/SQD iteration
3. **Graceful Stop**: Current iteration completes, then optimization stops
4. **Partial Results**: Convergence history up to cancellation point is saved
5. **Status Update**: Experiment status changed to "cancelled" with timestamp

**Typical Cancellation Time**: Within 1-2 optimization iterations (< 10 seconds for most molecules)

### For IBM Quantum Jobs

1. **Job Lookup**: Uses stored `cloud_job_id` to identify IBM Runtime job
2. **Cancellation Request**: Calls `QiskitRuntimeService.job(job_id).cancel()`
3. **IBM Processing**: IBM Quantum backend processes cancellation
4. **Status Sync**: Experiment marked as cancelled in Kanad database

**Notes**:
- IBM jobs that are already running on hardware may complete despite cancellation
- Jobs in queue are cancelled immediately
- Requires valid IBM API token and credentials

### For BlueQubit Jobs

1. **Job Lookup**: Uses stored `cloud_job_id` to identify BlueQubit job
2. **Cancellation Request**: Calls BlueQubit SDK `cancel_job(job_id)`
3. **GPU Release**: BlueQubit frees GPU resources immediately
4. **Status Sync**: Experiment marked as cancelled

**Notes**:
- BlueQubit jobs cancel very quickly (usually < 1 second)
- Requires valid BlueQubit API token

## Database Schema

### New Fields in `experiments` Table

| Column | Type | Description |
|--------|------|-------------|
| `cloud_job_id` | VARCHAR(255) | Cloud provider job ID (IBM/BlueQubit), NULL for local jobs |
| `cloud_backend` | VARCHAR(50) | Backend type: 'ibm', 'bluequbit', or NULL for local |
| `cancelled_at` | TIMESTAMP | UTC timestamp when cancellation was requested |

### Migration

Run the migration script to add these columns:

```bash
cd /home/mk/deeprealm/kanad
python api/migrations/add_cancellation_fields.py
```

Or manually apply:

```sql
ALTER TABLE experiments ADD COLUMN cloud_job_id VARCHAR(255);
ALTER TABLE experiments ADD COLUMN cloud_backend VARCHAR(50);
ALTER TABLE experiments ADD COLUMN cancelled_at TIMESTAMP;
```

## Implementation Details

### Architecture

```
User Request (PATCH /cancel)
    │
    ├─> Check experiment status (queued/running)
    │
    ├─> [If cloud job exists]
    │   ├─> IBM: IBMBackend.cancel_job()
    │   └─> BlueQubit: BlueQubitBackend.cancel_job()
    │
    ├─> [For local jobs]
    │   └─> JobQueue.cancel_job() → Set threading.Event()
    │
    └─> Update database
        ├─> status = "cancelled"
        ├─> cancelled_at = now()
        └─> completed_at = now()
```

### Thread-Safe Cancellation

The `JobQueue` class maintains a dictionary of `threading.Event` objects:

```python
self.cancellation_flags = {}  # experiment_id -> threading.Event()
```

**Workflow**:
1. When job is added: `cancellation_flags[exp_id] = threading.Event()`
2. When cancellation requested: `cancellation_flags[exp_id].set()`
3. In VQE/SQD callback: Check `if cancellation_flags[exp_id].is_set()` → raise exception
4. After job completes/cancelled: `del cancellation_flags[exp_id]`

### VQE Solver Integration

The progress callback in VQE solver checks for cancellation:

```python
def progress_callback(iteration: int, energy: float, params):
    # Check for cancellation
    if job_queue.is_cancelled(experiment_id):
        raise ExperimentCancelledException(experiment_id)

    # Normal progress tracking
    convergence_data.append({'iteration': iteration, 'energy': energy})
```

### Error Handling

Cancellation is treated as a controlled exception, not a failure:

```python
except ExperimentCancelledException:
    experiment.status = "cancelled"
    experiment.cancelled_at = datetime.now()
    experiment.error_message = "Cancelled by user"
```

## Testing

### Test Local VQE Cancellation

```python
# Start a long-running H2O experiment
import requests

experiment = {
    "name": "H2O VQE Test",
    "molecule": {
        "smiles": "O",
        "basis": "sto-3g"
    },
    "configuration": {
        "method": "VQE",
        "ansatz": "ucc",
        "max_iterations": 1000
    },
    "execute_immediately": True
}

response = requests.post('http://localhost:8000/api/v1/experiments', json=experiment)
exp_id = response.json()['id']

# Wait 5 seconds
import time
time.sleep(5)

# Cancel it
cancel_response = requests.patch(f'http://localhost:8000/api/v1/experiments/{exp_id}/cancel')
print(cancel_response.json())

# Check status
status = requests.get(f'http://localhost:8000/api/v1/experiments/{exp_id}/status')
print(f"Status: {status.json()['status']}")  # Should be "cancelled"
```

### Test IBM Quantum Cancellation

```python
# Submit to IBM Quantum
experiment = {
    "name": "H2 IBM Test",
    "molecule": {
        "smiles": "[H][H]",
        "basis": "sto-3g"
    },
    "configuration": {
        "method": "VQE",
        "backend": "ibm",
        "backend_name": "ibm_brisbane",
        "ibm_token": "your_token_here"
    },
    "execute_immediately": True
}

response = requests.post('http://localhost:8000/api/v1/experiments', json=experiment)
exp_id = response.json()['id']

# Cancel immediately (while in IBM queue)
cancel_response = requests.patch(f'http://localhost:8000/api/v1/experiments/{exp_id}/cancel')
print(cancel_response.json())
# Output: {"cloud_cancelled": true, "cloud_job_id": "c1234567890..."}
```

### Test BlueQubit Cancellation

```python
# Submit to BlueQubit
experiment = {
    "name": "H2 BlueQubit Test",
    "molecule": {
        "smiles": "[H][H]",
        "basis": "sto-3g"
    },
    "configuration": {
        "method": "VQE",
        "backend": "bluequbit",
        "device": "gpu",
        "bluequbit_token": "your_token_here"
    },
    "execute_immediately": True
}

response = requests.post('http://localhost:8000/api/v1/experiments', json=experiment)
exp_id = response.json()['id']

# Cancel
cancel_response = requests.patch(f'http://localhost:8000/api/v1/experiments/{exp_id}/cancel')
print(cancel_response.json())
```

## Limitations & Edge Cases

### Limitations

1. **Hardware Jobs**: IBM Quantum jobs running on real hardware may complete despite cancellation request
2. **Iteration Granularity**: Local jobs cancel at iteration boundaries (not mid-iteration)
3. **Cloud API Delays**: Cloud cancellation depends on provider API responsiveness
4. **No Undo**: Once cancelled, experiment cannot be resumed from checkpoint

### Edge Cases

1. **Already Completed**: Cannot cancel completed experiments (returns 400 error)
2. **Missing Credentials**: If cloud credentials are unavailable, cloud cancellation fails but local cancellation succeeds
3. **Network Issues**: Cloud cancellation may timeout; local cancellation always works
4. **Race Conditions**: If experiment completes naturally at same moment as cancellation, status may be "completed" instead of "cancelled"

## Future Enhancements

Potential improvements for future versions:

1. **Checkpointing**: Save intermediate parameters for resume capability
2. **Batch Cancellation**: Cancel multiple experiments in one request
3. **Auto-Cancel**: Cancel experiments exceeding time/cost budgets
4. **Cancel Hooks**: Custom callbacks when cancellation occurs
5. **WebSocket Notifications**: Real-time cancellation status updates

## Support

For issues or questions:
- GitHub Issues: [kanad/issues](https://github.com/your-org/kanad/issues)
- Email: support@example.com
- Docs: [Full API Documentation](https://docs.example.com)

---

**Last Updated**: 2025-10-08
**Version**: 1.0.0
**Author**: Kanad Backend Team
