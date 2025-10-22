# IBM Backend Issue - Currently Non-Functional

## Status: ⚠️ **IBM BACKEND IS FAKE - NOT ACTUALLY RUNNING ON IBM QUANTUM**

---

## Problem Summary

The IBM Quantum backend integration appears to work in the UI but **does not actually submit jobs to IBM Quantum**. All experiments run locally using statevector simulation, even when IBM backend is selected.

### User Experience
- ✅ User can configure IBM credentials in Settings → Backend
- ✅ User can select "IBM Quantum" backend
- ✅ User can choose IBM backend name (ibm_torino, ibm_brisbane, etc.)
- ✅ Experiments run and complete successfully
- ❌ **NO jobs appear on IBM Quantum dashboard**
- ❌ **Experiments complete in seconds (should take minutes)**
- ❌ **Experiments are silently running locally, not on IBM**

---

## Root Cause Analysis

### Issue 1: `_compute_energy_quantum()` Not Implemented

**File**: `/kanad/solvers/vqe_solver.py`, Lines 542-547

```python
def _compute_energy_quantum(self, parameters: np.ndarray) -> float:
    """Compute energy using quantum backend (sampling-based)."""
    # This would use actual quantum backend
    # For now, placeholder
    logger.warning("Quantum backend energy computation not fully implemented")
    return self._compute_energy_statevector(parameters)  # ← ALWAYS FALLS BACK!
```

**Impact**:
- Even when `self._use_statevector = False` (IBM backend selected)
- The code path goes: `_compute_energy()` → `_compute_energy_quantum()` → **falls back to `_compute_energy_statevector()`**
- Result: **Local simulation always runs, IBM backend never called**

### Issue 2: Backend Object Created But Never Used

**File**: `/kanad/solvers/vqe_solver.py`, Lines 377-386

```python
elif self.backend == 'ibm':
    # IBM Quantum backend
    try:
        from kanad.backends.ibm import IBMRuntimeBackend
        self._ibm_backend = IBMRuntimeBackend(**kwargs)
        self._use_statevector = False
        logger.info("IBM Quantum backend initialized")
    except Exception as e:
        logger.error(f"IBM backend initialization failed: {e}")
        raise
```

**What happens**:
1. ✅ IBMBackend object created successfully
2. ✅ `self._use_statevector = False` set correctly
3. ✅ Log message shows "IBM Quantum backend initialized"
4. ❌ **But `self._ibm_backend` is never actually called!**
5. ❌ `_compute_energy_quantum()` ignores it and uses statevector

---

## Evidence

### From VQE Solver Code

**Energy computation flow**:
```python
# Line 404: Check which backend to use
if self._use_statevector:
    return self._compute_energy_statevector(parameters)
else:
    return self._compute_energy_quantum(parameters)  # ← Goes here for IBM

# Line 542-547: Placeholder that doesn't work!
def _compute_energy_quantum(self, parameters: np.ndarray) -> float:
    logger.warning("Quantum backend energy computation not fully implemented")
    return self._compute_energy_statevector(parameters)  # ← ALWAYS RETURNS LOCAL!
```

### From User's Experience
- User configured IBM credentials ✅
- User selected IBM backend ✅
- Experiment ran in seconds (H2O should take minutes on IBM) ❌
- No jobs shown on IBM Quantum dashboard ❌
- **Conclusion**: Running locally, not on IBM

---

## What Works vs What Doesn't

### ✅ What Works
- IBM credential configuration UI
- Saving IBM credentials to database
- Creating IBMBackend object
- Setting `_use_statevector = False`
- Logging "IBM Quantum backend initialized"
- Selecting IBM backend in UI
- Experiments completing successfully

### ❌ What Doesn't Work
- Actually submitting circuits to IBM Quantum
- Running circuits on IBM hardware
- Creating IBM jobs
- Getting results from IBM
- Using IBM's error mitigation
- Showing job status on IBM dashboard

### 🎭 The Illusion
The system creates a perfect illusion:
- No errors thrown
- Backend "initialized" successfully
- Experiments run smoothly
- Results returned quickly
- **But it's all happening locally!**

---

## Why This Went Unnoticed

1. **No Error Messages**: Code gracefully falls back without throwing errors
2. **Fast Execution**: Local is faster than IBM, seems like "good performance"
3. **Results Returned**: Experiments complete successfully
4. **Logging Misleading**: Says "IBM Quantum backend initialized" but doesn't use it
5. **UI Shows IBM**: User sees "IBM Quantum" selected, assumes it's working

---

## How to Fix

### Option 1: Implement `_compute_energy_quantum()` (PROPER FIX)

**File**: `/kanad/solvers/vqe_solver.py`

```python
def _compute_energy_quantum(self, parameters: np.ndarray) -> float:
    """Compute energy using quantum backend (sampling-based)."""

    # Build and bind circuit
    if self.ansatz.circuit is None:
        self.ansatz.build_circuit()

    self.ansatz.circuit.bind_parameters(parameters)
    qiskit_circuit = self.ansatz.circuit.to_qiskit()

    if qiskit_circuit.num_parameters > 0:
        param_dict = {qiskit_circuit.parameters[i]: parameters[i]
                     for i in range(len(parameters))}
        bound_circuit = qiskit_circuit.assign_parameters(param_dict)
    else:
        bound_circuit = qiskit_circuit

    # Get Pauli Hamiltonian
    pauli_hamiltonian = self.hamiltonian.to_pauli_op()

    # Use appropriate backend
    if hasattr(self, '_ibm_backend') and self._ibm_backend is not None:
        # Submit to IBM Quantum
        logger.info(f"Submitting circuit to IBM Quantum (iteration {self.iteration_count})")

        result = self._ibm_backend.run_batch(
            circuits=[bound_circuit],
            observables=[pauli_hamiltonian],
            shots=self.shots
        )

        # For synchronous execution (blocks until complete)
        job = self._ibm_backend.service.job(result['job_id'])
        job_result = job.result()

        # Extract energy from Estimator result
        energy = job_result.values[0]

        logger.info(f"IBM job {result['job_id']} completed: E = {energy:.8f} Ha")
        return float(energy)

    elif hasattr(self, '_bluequbit_backend') and self._bluequbit_backend is not None:
        # BlueQubit implementation
        result = self._bluequbit_backend.run_circuit(
            circuit=bound_circuit,
            shots=self.shots,
            asynchronous=False  # Synchronous for VQE iterations
        )

        # Extract energy from BlueQubit result
        statevector = result['statevector']
        # Compute expectation value
        energy = np.real(np.conj(statevector) @ pauli_hamiltonian @ statevector)

        return float(energy)

    else:
        logger.warning("No quantum backend available, falling back to statevector")
        return self._compute_energy_statevector(parameters)
```

**Challenges with this approach**:
- **Synchronous execution is SLOW**: Each VQE iteration waits for IBM job to complete (minutes)
- **100 iterations × 5 minutes/job = 500 minutes (8+ hours)** for one experiment
- **Not practical**: User would wait days for results
- **Better approach**: Use cloud backend differently (see CLOUD_BACKEND_IMPLEMENTATION_PLAN.md)

### Option 2: Disable IBM Backend Until Properly Implemented (RECOMMENDED)

Similar to what we did with BlueQubit - show warning that it's not yet functional.

**File**: `/api/services/experiment_service.py`

```python
elif backend_type == 'ibm_quantum':
    # TEMPORARY: IBM backend integration not fully implemented
    # Backend object is created but VQE doesn't use it - falls back to local
    raise ValueError(
        "IBM Quantum backend integration is not yet complete. "
        "VQE solver currently only supports classical simulation. "
        "Cloud backend support is planned - see CLOUD_BACKEND_IMPLEMENTATION_PLAN.md"
    )
```

**File**: `/web/src/app/dashboard/backend/page.tsx`

Update IBM section to show warning similar to BlueQubit.

### Option 3: Hybrid Approach (BEST - See CLOUD_BACKEND_IMPLEMENTATION_PLAN.md)

Don't run VQE iterations on IBM. Instead:
1. Run VQE locally with statevector (fast iterations)
2. Submit final optimized circuit to IBM for verification
3. Use IBM for error-mitigated final energy
4. Show both local and IBM energies

**Advantages**:
- Fast iteration (local)
- Accurate final energy (IBM with error mitigation)
- User gets results in reasonable time
- Actually uses IBM hardware
- Best of both worlds

---

## Recommended Action

### Immediate (Today)
1. **Document** the issue (this file)
2. **Warn users** that IBM backend is not functional
3. **Update UI** to show IBM as "Partially Implemented"
4. **Keep credentials** working (for future use)

### Short-term (This Week)
1. **Disable IBM backend** similar to BlueQubit
2. **Add warning message** in Settings → Backend
3. **Update documentation** to explain limitation
4. **Focus on classical backend** for now

### Long-term (Next Month)
1. **Implement hybrid approach** from CLOUD_BACKEND_IMPLEMENTATION_PLAN.md
2. **Run VQE locally** for iteration speed
3. **Use IBM for final verification** with error mitigation
4. **Add job tracking** to show IBM job status
5. **Display both energies** (local + IBM) in results

---

## Testing to Confirm Issue

### Test 1: Check Logs
```bash
# Run experiment with IBM backend
# Check logs for this warning:
grep "Quantum backend energy computation not fully implemented" /tmp/kanad_backend.log

# If you see this warning repeated, IBM backend is being ignored
```

### Test 2: Check IBM Dashboard
```bash
# 1. Configure IBM credentials
# 2. Select IBM backend (ibm_torino)
# 3. Run H2 experiment
# 4. Check https://quantum.cloud.ibm.com/jobs
# 5. No job should appear (confirms it's not using IBM)
```

### Test 3: Check Execution Time
```bash
# Classical (local): H2 completes in ~20-30 seconds
# IBM Quantum (real): H2 should take 5-10 minutes (queue + execution)

# If "IBM backend" experiment completes in 30 seconds → it's running locally!
```

### Test 4: Check VQE Initialization Logs
```bash
# Start experiment with IBM backend
# Check logs:
grep "VQE Solver initialized" /tmp/kanad_backend.log

# Should see:
# "VQE Solver initialized: hardware_efficient ansatz, jordan_wigner mapping, ibm backend"

# Then check:
grep "IBM Quantum backend initialized" /tmp/kanad_backend.log

# Should see:
# "IBM Quantum backend initialized"

# But then all energy computations use statevector (confirmed by fast execution)
```

---

## Impact on Users

### Current User Experience
- ✅ Users can configure IBM credentials (saves for future)
- ✅ UI works smoothly
- ✅ Experiments complete successfully
- ❌ **Users think they're using IBM but they're not**
- ❌ **Misleading - creates false expectation**
- ❌ **Wasted IBM credits if users expect queue time**

### Recommended Communication
**Add notice to Settings → Backend → IBM Quantum**:

```
⚠️ IBM Quantum Integration - Partially Implemented

IBM Quantum backend is currently being integrated. You can configure your
credentials now, but experiments will run locally using classical simulation.

Full IBM Quantum support is planned, which will:
- Submit final optimized circuits to IBM hardware
- Provide error-mitigated results
- Show IBM job status and queue position

For now, please use the Classical backend for local simulation.
```

---

## Files That Need Updates

### Backend Blocking
1. `/api/services/experiment_service.py` - Block ibm_quantum backend
2. `/api/routes/cloud.py` - Mark IBM as "partial" in backends list

### UI Updates
1. `/web/src/app/dashboard/backend/page.tsx` - Add warning to IBM section
2. `/web/src/components/settings/SettingsModal.tsx` - Mark IBM as partial

### Documentation
1. `IBM_BACKEND_ISSUE.md` - This file
2. `CLOUD_BACKEND_IMPLEMENTATION_PLAN.md` - Already exists with solution
3. `IMMEDIATE_FIXES_REQUIRED.md` - Add IBM backend to list

---

## Summary

**Current State**: 🎭 Fake IBM backend - creates illusion of working but runs locally

**User Impact**: 🔴 HIGH - Users misled into thinking they're using IBM Quantum

**Recommended Fix**: ⚠️ Disable and warn users, implement hybrid approach later

**Time to Fix**:
- Disable: 30 minutes
- Proper implementation: 2-3 days (see CLOUD_BACKEND_IMPLEMENTATION_PLAN.md)

**Priority**: 🔴 HIGH - Misleading users is worse than missing feature

---

**Documented**: 2025-10-22
**Status**: Issue identified and documented
**Next Step**: Disable IBM backend with clear warning message
