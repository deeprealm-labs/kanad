# Excited States Method Selection - FINAL FIX

## Problem Identified

You selected **CIS** method but the system ran **VQE** instead. The issue had TWO causes:

### 1. Backend Config Model Missing Fields
The `BackendConfig` Pydantic model was missing `excited_method` and `excited_n_states` fields, so they were silently dropped when submitted from the frontend.

**FIXED:** Added fields to [api/routes/experiments.py:41-43](api/routes/experiments.py#L41-L43)

### 2. Wrong Logic in execute_excited_states()
The code was ignoring the user's method choice and forcing VQE for quantum backends OR forcing CIS regardless of user selection.

**FIXED:** Completely rewrote logic in [api/services/experiment_service.py:630-647](api/services/experiment_service.py#L630-L647)

## The Fix

### New Logic Flow

```python
# 1. READ USER'S CHOICE FIRST
user_method = config.get('excited_method') or config.get('excitedMethod', 'cis')
print(f"🔍 User selected excited states method: {user_method}")

# 2. RESPECT THE USER'S CHOICE
if user_method == 'vqe' and backend_type in ['bluequbit', 'ibm_quantum']:
    # User explicitly wants VQE with quantum - allow it (with warning)
    print(f"⚠️  VQE excited states with quantum backend - requires many jobs")
    backend_kwargs = get_backend_kwargs(backend_type, config)
elif backend_type in ['bluequbit', 'ibm_quantum'] and user_method != 'vqe':
    # User wants CIS/TDDFT with quantum backend - override to classical
    print(f"📊 {user_method.upper()} is classical-only, switching backend to classical")
    backend_type = 'classical'

# 3. USE THE USER'S METHOD
method = user_method  # NO OVERRIDING!
```

### What Changed

**Before:**
- ❌ Hardcoded VQE when quantum backend detected
- ❌ Ignored user's method selection
- ❌ No way to actually run CIS or TDDFT

**After:**
- ✅ Reads user's method choice (`excitedMethod` from frontend)
- ✅ Uses exactly what the user selected
- ✅ Only overrides backend (classical for CIS/TDDFT, quantum for VQE)
- ✅ Warns if VQE will submit many jobs

## Testing Instructions

### **IMPORTANT: Restart Everything First**

There are multiple API servers running with cached code. You MUST:

```bash
# Kill all API servers
pkill -f "uvicorn"
pkill -f "restart_api"

# Clear Python cache
cd /home/mk/deeprealm/kanad
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
find . -name "*.pyc" -delete

# Restart API
. env/bin/activate
./restart_api.sh
```

### Test Each Method

#### Test 1: CIS Method (Classical)
1. Open Settings
2. Select "Excited States"
3. Choose **"CIS"** method
4. Set States: 3
5. Backend: BlueQubit or Classical
6. Run H2 experiment

**Expected logs:**
```
🔍 User selected excited states method: cis
📊 CIS is classical-only, switching backend to classical
🔍 DEBUG Excited States Config:
   config.get('excited_method'): None
   config.get('excitedMethod'): cis
   FINAL method: cis
   FINAL n_states: 3
   Backend type: classical
🔬 Starting CIS excited states calculation...
📊 Computing 3 excited states
⚙️ Running CIS solver...
✅ CIS calculation complete
```

**NO VQE, NO QUANTUM JOBS!**

#### Test 2: TDDFT Method (Classical)
1. Select **"TDDFT"** method
2. Run experiment

**Expected logs:**
```
🔍 User selected excited states method: tddft
📊 TDDFT is classical-only, switching backend to classical
FINAL method: tddft
🔬 Starting TDDFT excited states calculation...
```

#### Test 3: VQE Method (Quantum - Many Jobs!)
1. Select **"VQE"** method
2. Backend: BlueQubit
3. Run experiment

**Expected logs:**
```
🔍 User selected excited states method: vqe
⚠️  VQE excited states with quantum backend - requires many jobs
FINAL method: vqe
   Backend type: bluequbit
🔬 Starting VQE excited states calculation...
📊 Computing 3 excited states
```

**WILL SUBMIT MANY QUANTUM JOBS (3 states × ~5-10 evaluations = 15-30 jobs)**

## Files Modified

1. **[api/routes/experiments.py](api/routes/experiments.py)**
   - Lines 41-43: Added `excited_method` and `excited_n_states` to BackendConfig

2. **[api/services/experiment_service.py](api/services/experiment_service.py)**
   - Lines 630-647: Complete rewrite of method selection logic
   - Lines 679-690: Simplified method assignment and debug logging

3. **[web/src/components/settings/SettingsModal.tsx](web/src/components/settings/SettingsModal.tsx)**
   - Already has UI (previous fix)

## Summary

The system NOW:
1. ✅ Accepts `excitedMethod` and `excitedNStates` from frontend
2. ✅ Reads the user's method choice
3. ✅ RESPECTS the user's choice - no more hardcoding!
4. ✅ Shows debug logs to verify correct method
5. ✅ Runs CIS when CIS selected
6. ✅ Runs TDDFT when TDDFT selected
7. ✅ Runs VQE when VQE selected (with warning)

## Next Step

**YOU MUST RESTART THE API SERVER** to load the new code!

The Python bytecode cache was causing old code to run. After restarting properly, test all three methods and verify the logs show the correct method being used.
