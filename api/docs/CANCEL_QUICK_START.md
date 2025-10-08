# Job Cancellation Quick Start Guide

## TL;DR

Cancel any running or queued experiment:

```bash
curl -X PATCH http://localhost:8000/api/v1/experiments/{experiment_id}/cancel
```

## Basic Usage

### 1. Cancel a Running Experiment

```bash
# Get list of running experiments
curl http://localhost:8000/api/v1/experiments?status=running

# Cancel experiment with ID 42
curl -X PATCH http://localhost:8000/api/v1/experiments/42/cancel
```

**Response**:
```json
{
  "status": "cancelled",
  "experiment_id": 42,
  "local_cancelled": true,
  "cloud_cancelled": false,
  "message": "Experiment cancellation requested successfully"
}
```

### 2. Check Cancellation Status

```bash
# Get experiment details
curl http://localhost:8000/api/v1/experiments/42

# Check just the status
curl http://localhost:8000/api/v1/experiments/42/status
```

**Response**:
```json
{
  "id": 42,
  "status": "cancelled",
  "cancelled_at": "2025-10-08T15:30:45.123456",
  "convergence_data": [...],  // Partial results saved
  "error_message": "Cancelled by user"
}
```

## Python Examples

### Cancel a Local VQE Job

```python
import requests
import time

# Start experiment
response = requests.post('http://localhost:8000/api/v1/experiments', json={
    "name": "H2 Test",
    "molecule": {"smiles": "[H][H]", "basis": "sto-3g"},
    "configuration": {"method": "VQE", "max_iterations": 1000},
    "execute_immediately": True
})

exp_id = response.json()['id']
print(f"Experiment {exp_id} started")

# Wait a bit
time.sleep(5)

# Cancel it
cancel_response = requests.patch(f'http://localhost:8000/api/v1/experiments/{exp_id}/cancel')
print(f"Cancelled: {cancel_response.json()['status']}")

# Get partial results
results = requests.get(f'http://localhost:8000/api/v1/experiments/{exp_id}').json()
print(f"Partial iterations: {len(results['convergence_data'])}")
```

### Cancel IBM Quantum Job

```python
import requests

# Submit to IBM Quantum
response = requests.post('http://localhost:8000/api/v1/experiments', json={
    "name": "H2 IBM Test",
    "molecule": {"smiles": "[H][H]", "basis": "sto-3g"},
    "configuration": {
        "method": "VQE",
        "backend": "ibm",
        "backend_name": "ibm_brisbane",
        "ibm_token": "your_token_here"
    },
    "execute_immediately": True
})

exp_id = response.json()['id']
cloud_job_id = response.json().get('cloud_job_id')

print(f"IBM job submitted: {cloud_job_id}")

# Cancel immediately (while in queue)
cancel_response = requests.patch(f'http://localhost:8000/api/v1/experiments/{exp_id}/cancel')
result = cancel_response.json()

print(f"Cloud cancelled: {result['cloud_cancelled']}")
print(f"Cloud job ID: {result['cloud_job_id']}")
```

### Cancel BlueQubit Job

```python
import requests

# Submit to BlueQubit
response = requests.post('http://localhost:8000/api/v1/experiments', json={
    "name": "H2 BlueQubit Test",
    "molecule": {"smiles": "[H][H]", "basis": "sto-3g"},
    "configuration": {
        "method": "VQE",
        "backend": "bluequbit",
        "device": "gpu",
        "bluequbit_token": "your_token_here"
    },
    "execute_immediately": True
})

exp_id = response.json()['id']

# Cancel
cancel_response = requests.patch(f'http://localhost:8000/api/v1/experiments/{exp_id}/cancel')
print(f"BlueQubit job cancelled: {cancel_response.json()}")
```

## Frontend Integration

### React Example

```typescript
// Cancel button component
const CancelButton = ({ experimentId }) => {
  const [cancelling, setCancelling] = useState(false);

  const handleCancel = async () => {
    setCancelling(true);
    try {
      const response = await fetch(
        `http://localhost:8000/api/v1/experiments/${experimentId}/cancel`,
        { method: 'PATCH' }
      );

      if (response.ok) {
        const result = await response.json();
        console.log('Cancelled:', result);
        // Update UI
      } else {
        const error = await response.json();
        console.error('Cancellation failed:', error.detail);
      }
    } finally {
      setCancelling(false);
    }
  };

  return (
    <button
      onClick={handleCancel}
      disabled={cancelling}
      className="btn-cancel"
    >
      {cancelling ? 'Cancelling...' : 'Cancel Experiment'}
    </button>
  );
};
```

### Vue.js Example

```vue
<template>
  <button @click="cancelExperiment" :disabled="cancelling">
    {{ cancelling ? 'Cancelling...' : 'Cancel' }}
  </button>
</template>

<script>
export default {
  props: ['experimentId'],
  data() {
    return {
      cancelling: false
    };
  },
  methods: {
    async cancelExperiment() {
      this.cancelling = true;
      try {
        const response = await fetch(
          `http://localhost:8000/api/v1/experiments/${this.experimentId}/cancel`,
          { method: 'PATCH' }
        );

        if (response.ok) {
          const result = await response.json();
          this.$emit('cancelled', result);
        } else {
          const error = await response.json();
          this.$emit('error', error.detail);
        }
      } finally {
        this.cancelling = false;
      }
    }
  }
};
</script>
```

## Common Scenarios

### Scenario 1: Cancel Long-Running Calculation

```bash
# User starts an experiment that's taking too long
# They realize they used wrong parameters

# Find the experiment
curl http://localhost:8000/api/v1/experiments?status=running | jq '.experiments[0]'

# Cancel it
curl -X PATCH http://localhost:8000/api/v1/experiments/123/cancel

# Start a new one with correct parameters
curl -X POST http://localhost:8000/api/v1/experiments -d '...'
```

### Scenario 2: Free Up Resources

```bash
# Multiple experiments running, need to free resources

# List running experiments
curl http://localhost:8000/api/v1/experiments?status=running | jq '.experiments[].id'

# Cancel low-priority ones
curl -X PATCH http://localhost:8000/api/v1/experiments/45/cancel
curl -X PATCH http://localhost:8000/api/v1/experiments/46/cancel
```

### Scenario 3: Emergency Stop All Jobs

```bash
# Emergency: need to stop all running jobs
# (requires scripting)

# Get all running experiment IDs
for id in $(curl -s http://localhost:8000/api/v1/experiments?status=running | jq -r '.experiments[].id'); do
  curl -X PATCH http://localhost:8000/api/v1/experiments/$id/cancel
  echo "Cancelled experiment $id"
done
```

## Error Handling

### Handle "Cannot Cancel Completed" Error

```python
import requests

def safe_cancel(experiment_id):
    response = requests.patch(f'http://localhost:8000/api/v1/experiments/{experiment_id}/cancel')

    if response.status_code == 200:
        print(f"Experiment {experiment_id} cancelled successfully")
        return True
    elif response.status_code == 400:
        error = response.json()
        if "Cannot cancel" in error['detail']:
            print(f"Experiment {experiment_id} already finished")
            return False
    elif response.status_code == 404:
        print(f"Experiment {experiment_id} not found")
        return False
    else:
        print(f"Unexpected error: {response.status_code}")
        return False

# Usage
safe_cancel(42)
```

### Check if Experiment is Cancellable

```python
def is_cancellable(experiment_id):
    response = requests.get(f'http://localhost:8000/api/v1/experiments/{experiment_id}/status')

    if response.status_code == 200:
        status = response.json()['status']
        return status in ['queued', 'running']

    return False

# Usage
if is_cancellable(42):
    requests.patch(f'http://localhost:8000/api/v1/experiments/42/cancel')
else:
    print("Experiment cannot be cancelled")
```

## FAQs

### Q: How long does cancellation take?

**A**:
- Local jobs: Within 1-2 VQE/SQD iterations (1-10 seconds typical)
- IBM Quantum: Immediate for queued jobs, may complete for running jobs
- BlueQubit: < 1 second

### Q: Are partial results saved?

**A**: Yes! All convergence data up to the cancellation point is preserved in the `convergence_data` field.

### Q: Can I resume a cancelled experiment?

**A**: Not currently. This is a planned feature for v2.0 (checkpointing).

### Q: What happens to cloud costs if I cancel?

**A**:
- IBM: You're charged for queue time and any execution that occurred
- BlueQubit: Free tier is not charged; paid tier charged for actual execution time

### Q: Can I cancel someone else's experiment?

**A**: Currently no authorization checks (planned enhancement). Any user can cancel any experiment.

### Q: What if cancellation fails?

**A**: The local job will still be cancelled. Cloud cancellation failure is logged as a warning but doesn't prevent local cancellation.

## Troubleshooting

### Problem: Cancellation not working

**Check**:
1. Experiment status: `curl http://localhost:8000/api/v1/experiments/42/status`
2. API server logs: `tail -f api/logs/api.log`
3. Job queue status: Check if workers are running

### Problem: Cloud job not cancelled

**Check**:
1. API token validity (IBM/BlueQubit)
2. Network connectivity to cloud provider
3. Job may have already completed

### Problem: Partial results not saved

**Check**:
1. Database connection
2. Convergence callback was being called (check logs)
3. Cancellation happened before first iteration

## Next Steps

- **Full Documentation**: See `/api/docs/JOB_CANCELLATION.md`
- **Implementation Details**: See `CANCELLATION_IMPLEMENTATION_REPORT.md`
- **Run Tests**: `python api/tests/test_cancellation.py`
- **Report Issues**: Create GitHub issue with label `cancellation`

---

**Quick Links**:
- [Full API Documentation](./JOB_CANCELLATION.md)
- [Implementation Report](../../CANCELLATION_IMPLEMENTATION_REPORT.md)
- [Test Suite](../tests/test_cancellation.py)
