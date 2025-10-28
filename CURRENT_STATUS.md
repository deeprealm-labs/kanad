# Current Status - Excited States Progress Bar & Config Issues

## Investigation Summary (2025-10-28)

### What We Found

The server IS running properly and code changes ARE loading. Investigation of the logs revealed:

#### 1. Config IS Being Sent Correctly
From logs (experiment fcd5f162):
```
🔍 Backend config received: {
  'method': 'EXCITED_STATES',
  'excited_method': 'vqe',
  'excited_n_states': 3,  ← USER SELECTED 3 STATES
  ...
}
```

#### 2. But n_states Ends Up Wrong
From the same experiment:
```
🔍 DEBUG Excited States Config:
   config.get('excited_method'): vqe
   config.get('excitedMethod'): None
   FINAL method: vqe
   FINAL n_states: 5  ← SHOULD BE 3!
```

### Root Cause

The problem is in [api/services/experiment_service.py:683](api/services/experiment_service.py#L683):

```python
# OLD BUGGY CODE (before my fix):
n_states = config.get('excited_n_states') or config.get('excitedNStates') or config.get('n_states', 5)
```

The bug: When `excited_n_states` exists and equals `3`, this SHOULD work. But the logs show it's not. This suggests there's something else going on - perhaps the config dict structure is nested differently than expected.

### My Latest Fix

I changed the code to be more explicit (lines 683-707):

```python
# NEW CODE (should be loading now):
_excited_n_states = config.get('excited_n_states')
_excitedNStates = config.get('excitedNStates')
_n_states = config.get('n_states')

print(f"🔍 DEBUG Excited States Config:")
print(f"   config.get('excited_method'): {config.get('excited_method')}")
print(f"   config.get('excitedMethod'): {config.get('excitedMethod')}")
print(f"   config.get('excited_n_states'): {_excited_n_states} (type: {type(_excited_n_states)})")
print(f"   config.get('excitedNStates'): {_excitedNStates} (type: {type(_excitedNStates)})")
print(f"   config.get('n_states'): {_n_states} (type: {type(_n_states)})")

# Use the first non-None value
if _excited_n_states is not None:
    n_states = _excited_n_states
elif _excitedNStates is not None:
    n_states = _excitedNStates
elif _n_states is not None:
    n_states = _n_states
else:
    n_states = 5
```

This will show us:
- The exact values being read
- The Python types (int, str, NoneType, etc.)
- Which branch is being taken

### Next Steps for User

**Please run a NEW experiment** (not reuse old ones) and check the terminal logs for the new debug output. It should now show:

```
🔍 DEBUG Excited States Config:
   config.get('excited_method'): vqe
   config.get('excitedMethod'): None
   config.get('excited_n_states'): 3 (type: <class 'int'>)  ← Should show type
   config.get('excitedNStates'): None (type: <class 'NoneType'>)
   config.get('n_states'): None (type: <class 'NoneType'>)
   FINAL method: vqe
   FINAL n_states: 3  ← Should be 3 now!
```

If it STILL shows 5, the debug output will tell us why (which variable is None, what type it is, etc.).

### Other Issues Found

#### VQE with Quantum Backend Still Crashes
Error from logs:
```
AttributeError: 'str' object has no attribute 'get'
  File "api/services/experiment_service.py", line 643
    backend_kwargs = get_backend_kwargs(backend_type, config)
```

BUT: Looking at line 643 in my current code, it should be:
```python
backend_type, backend_kwargs = get_backend_kwargs(config, experiment_id)
```

This means the error is from OLD code still in cache. Need to investigate further.

#### Progress Bar & Convergence Graph
Changes to add these were made in:
- [kanad/utils/vqe_solver.py:959-976](kanad/utils/vqe_solver.py#L959-L976) - Progress broadcasting
- [kanad/solvers/excited_states_solver.py](kanad/solvers/excited_states_solver.py) - State-by-state logs

BUT: These files are in `kanad/` directory which is outside the `api/` directory that uvicorn's `--reload` watches. This means **these changes won't auto-reload**.

### To Test Everything

1. **Kill and restart server manually** to load kanad/ changes:
   ```bash
   pkill -9 -f uvicorn
   find /home/mk/deeprealm/kanad -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
   cd /home/mk/deeprealm/kanad
   . env/bin/activate
   cd api
   python3 -m uvicorn main:app --reload --port 8000
   ```

2. **Create a NEW excited states experiment** (don't reuse old ones)

3. **Check terminal for new debug output** - Should show types and actual values

4. **Check if progress bar updates** - Should go from 0% → 20% → 80% → 100%

5. **Check if convergence graph shows orange bars** - Should appear in real-time

### Files Modified in This Session

1. [api/services/experiment_service.py:683-707](api/services/experiment_service.py#L683-L707)
   - Added explicit debug logging with types
   - Changed n_states reading to be more explicit

### Status

🔄 **NEEDS TESTING** - User must run new experiment to see if fixes work

The debug output will definitively show us whether the config is structured correctly or if there's another issue.
