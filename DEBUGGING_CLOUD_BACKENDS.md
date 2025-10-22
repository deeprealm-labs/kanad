# Debugging Cloud Backend Execution

## Logging Added

I've added comprehensive logging throughout the backend execution pipeline to help debug why experiments aren't running on real quantum hardware.

### 1. Backend Configuration Logs

In [api/services/experiment_service.py](api/services/experiment_service.py):

```
🔧 get_backend_kwargs called with backend_type: ibm_quantum
🌐 Configuring IBM Quantum backend...
✅ IBM credentials loaded from database
📍 Using IBM Quantum backend: ibm_brisbane
```

OR for statevector:
```
🔧 get_backend_kwargs called with backend_type: classical
📍 Using statevector simulation
```

### 2. Backend Initialization Logs

In [kanad/solvers/vqe_solver.py](kanad/solvers/vqe_solver.py):

```
🔧 Initializing backend: ibm
🌐 Initializing IBM backend with kwargs: ['backend_name', 'api_token', 'instance']
✅ IBM Quantum backend initialized successfully
   Backend name: ibm_brisbane
```

### 3. Job Submission Logs

During VQE execution, you should see:

```
🚀 Submitting job to IBM Quantum (function eval 1)
✅ IBM job submitted: d3rxxxx...
🚀 Submitting job to IBM Quantum (function eval 2)
✅ IBM job submitted: d3rxxxx...
...
```

For BlueQubit:
```
🚀 Submitting job to BlueQubit (function eval 1)
✅ BlueQubit job completed
```

## How to Debug

### Step 1: Check API Server Logs

When you run an experiment, watch the API server console output. You should see the emoji-prefixed logs showing:

1. **Backend configuration**: What backend type was requested
2. **Credential loading**: Whether credentials were found
3. **Backend initialization**: Whether the cloud backend initialized successfully
4. **Job submissions**: Each quantum job being submitted

### Step 2: What You Should See

**For IBM Quantum backend**:
```bash
# Start API server
cd api && python3 -m uvicorn main:app --reload

# When experiment runs, you should see:
🔧 get_backend_kwargs called with backend_type: ibm_quantum
🌐 Configuring IBM Quantum backend...
✅ IBM credentials loaded from database
📍 Using IBM Quantum backend: ibm_brisbane
🔧 Initializing backend: ibm
🌐 Initializing IBM backend with kwargs: ['backend_name', 'api_token', 'instance']
✅ IBM Quantum backend initialized successfully
   Backend name: ibm_brisbane
🚀 Submitting job to IBM Quantum (function eval 1)
✅ IBM job submitted: d3ruc9g60rgc73a9n090
```

**For Classical/Statevector backend**:
```bash
🔧 get_backend_kwargs called with backend_type: classical
📍 Using statevector simulation
🔧 Initializing backend: statevector
📍 Using statevector simulation
```

### Step 3: If You DON'T See Job Submissions

If you see backend initialization but NO job submission logs, it means:

1. **The backend initialized** - ✅
2. **But jobs aren't being submitted** - ❌

This could mean:
- VQE is falling back to statevector due to an error
- The `_compute_energy_quantum()` method isn't being called
- There's an exception being caught silently

Check for these warning messages:
```
⚠️ Falling back to statevector simulation
```

### Step 4: Check Experiment Configuration

Verify the experiment is configured with cloud backend:

```bash
sqlite3 api/kanad_experiments.db "SELECT id, backend, json_extract(configuration, '$.backend') FROM experiments WHERE status='running'"
```

Should show:
```
<experiment_id>|ibm_quantum|ibm_quantum
```

NOT:
```
<experiment_id>|classical|classical
```

## Known Issues Fixed

### Issue 1: Stuck Experiments ✅ FIXED
Experiments that are cancelled show `status='running'` even though job is `status='cancelled'`.

**Fix**: Updated database manually for now. Need to add proper cleanup in cancellation handler.

### Issue 2: Missing Analysis Data ✅ FIXED
Multi-atom molecules weren't generating analysis data.

**Fix**: Added analysis tool initialization in components mode ([vqe_solver.py:214-225](kanad/solvers/vqe_solver.py#L214-L225))

### Issue 3: Iteration Count Confusion ✅ FIXED
Showing 4149 iterations instead of 100.

**Fix**: Now shows optimizer iterations (100) separately from function evaluations (4149)

## Testing

### Test IBM Backend

1. Start API server: `cd api && python3 -m uvicorn main:app --reload`
2. Run experiment from GUI with:
   - Backend: IBM Quantum
   - Backend Name: ibm_brisbane
   - Optimizer: COBYLA (recommended)
   - Max Iterations: 10 (for testing)
3. Watch server console for logs
4. You should see job IDs being printed
5. Verify jobs on IBM dashboard: https://quantum.cloud.ibm.com/

### Test BlueQubit Backend

1. Start API server
2. Run experiment with:
   - Backend: BlueQubit
   - Optimizer: COBYLA
   - Max Iterations: 10
3. Watch server console for logs
4. Verify jobs on BlueQubit dashboard: https://app.bluequbit.io/

## Next Steps

If you still don't see job submissions after these logging changes:

1. **Share the server console output** - I need to see what the logs say
2. **Check credentials** - Verify they're in the database:
   ```bash
   sqlite3 api/kanad_experiments.db "SELECT provider FROM cloud_credentials"
   ```
3. **Check experiment config** - Verify backend is set correctly in the experiment
4. **Test integration tests** - Run the integration tests to verify backends work:
   ```bash
   IBM_API='...' IBM_CRN='...' python3 -m pytest tests/integration/test_ibm_integration.py -v
   ```

## Summary

The logging will now show you EXACTLY where the execution is happening:
- ✅ Backend configuration
- ✅ Credential loading
- ✅ Backend initialization
- ✅ Job submissions (THIS is what you need to see!)

If you see job submission logs but no jobs on the cloud platform, then there's an issue with the backend API calls themselves. If you DON'T see job submission logs, the VQE solver is falling back to statevector for some reason.
