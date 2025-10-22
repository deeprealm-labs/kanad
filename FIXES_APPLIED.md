# Fixes Applied for SQD/EXCITED_STATES Issues

## Issues Reported by User
1. ❌ **Iteration counter showing 9 instead of 6** - Shows stale data from previous experiment
2. ❌ **Analysis data not displaying** - Even though backend generates it
3. ❌ **Backend selection ignored** - SQD doesn't actually use quantum circuits
4. ⚠️  **EXCITED_STATES redirects to SQD** - Confusing UX

## Fixes Applied

### 1. Iteration Counter - FIXED ✅
**File**: `web/src/components/simulation/ExperimentMonitor.tsx`

**Problem**: Iteration counter accumulated from previous experiments, showing stale values.

**Solution A - State Reset** (lines 65-89):
```typescript
// Reset state when experiment ID changes (new experiment)
useEffect(() => {
  console.log("🔄 Resetting state for new experiment:", experimentId);

  // Reset all state
  setCurrentIteration(0);  // <-- RESET TO 0
  setConvergenceData([]);
  setProgress(0);
  setLogs([]);
  setStatus("queued");
  setResults(null);
  setExperiment(null);

  // Close any existing WebSocket
  if (wsRef.current) {
    wsRef.current.close();
    wsRef.current = null;
  }

  // Stop any existing polling
  if (pollingRef.current) {
    pollingRef.current();
    pollingRef.current = null;
  }
}, [experimentId]);  // <-- Triggers when experimentId changes
```

**Solution B - Correct Final Count** (lines 236-249):
```typescript
} else if (exp.status === "completed") {
  setProgress(100);
  // Set final iteration count based on convergence data from results
  if (exp.results) {
    let finalIterations = 0;
    if (exp.results.convergence_history && Array.isArray(exp.results.convergence_history)) {
      finalIterations = exp.results.convergence_history.length;
    } else if (exp.results.energy_history && Array.isArray(exp.results.energy_history)) {
      finalIterations = exp.results.energy_history.length;
    }
    if (finalIterations > 0) {
      setCurrentIteration(finalIterations);  // <-- SET CORRECT VALUE
    }
  }
```

**Expected Behavior**:
- When new experiment starts: Counter resets to 0
- When experiment completes: Counter shows actual number of iterations (e.g., 6 for SQD)
- No more stale values from previous experiments

**How to Verify**:
1. Run an SQD experiment → Should show "Current Iteration: 6" when done
2. Run another experiment → Should reset to 0, not accumulate from 6

---

### 2. Analysis Data Display - NEEDS BROWSER CONSOLE CHECK ❓
**Files**:
- `web/src/components/simulation/ExperimentMonitor.tsx` (lines 192-196)
- `web/src/components/simulation/AnalysisResults.tsx` (lines 8-16)

**Problem**: Analysis not displaying even though backend generates and stores it.

**Debugging Added**:
```typescript
// In ExperimentMonitor.tsx (lines 192-196)
if (exp.results) {
  console.log("📊 Received results:", exp.results);
  console.log("📊 Has analysis in results?", !!exp.results.analysis);
  if (exp.results.analysis) {
    console.log("📊 Analysis keys:", Object.keys(exp.results.analysis));
  }
  setResults(exp.results);
}

// In AnalysisResults.tsx (lines 8-16)
console.log("AnalysisResults received results:", results);
console.log("Has analysis?", !!results?.analysis);
if (results?.analysis) {
  console.log("Analysis keys:", Object.keys(results.analysis));
  console.log("Energy components:", results.analysis.energy_components);
  console.log("Bonding:", results.analysis.bonding);
  console.log("Properties:", results.analysis.properties);
}
```

**Status**:
- ✅ Backend generates analysis (confirmed by backend logs)
- ✅ Database stores analysis (confirmed by SQL query)
- ❓ Frontend receives analysis (NEEDS VERIFICATION via browser console)
- ❓ Component renders analysis (NEEDS VERIFICATION via browser console)

**Next Step**:
See `BROWSER_CONSOLE_CHECK.md` for detailed instructions on what to check in browser console.

**Possible Root Causes**:
1. API not returning analysis field
2. JSON serialization/deserialization issue
3. Component render condition not met (`status === "completed" && results`)
4. Data structure mismatch (e.g., `results.analysis` vs `results.analysisData`)

---

### 3. Method-Specific Configuration Display - FIXED ✅
**File**: `web/src/components/simulation/ExperimentMonitor.tsx` (lines 471-517)

**Problem**: SQD experiments showed VQE fields (Ansatz) instead of SQD fields (Subspace, States).

**Solution**:
```typescript
{/* VQE-specific fields */}
{experimentConfig?.backendSettings?.method === "VQE" && (
  <div className="flex justify-between">
    <span className="text-muted-foreground">Ansatz:</span>
    <span className="font-quando font-medium">
      {experimentConfig?.backendSettings?.ansatz || "HEA"}
    </span>
  </div>
)}

{/* SQD-specific fields */}
{experimentConfig?.backendSettings?.method === "SQD" && (
  <>
    <div className="flex justify-between">
      <span className="text-muted-foreground">Subspace:</span>
      <span className="font-quando font-medium">
        {experimentConfig?.backendSettings?.subspaceDim || 10}
      </span>
    </div>
    <div className="flex justify-between">
      <span className="text-muted-foreground">States:</span>
      <span className="font-quando font-medium">
        {experimentConfig?.backendSettings?.nStates || 3}
      </span>
    </div>
  </>
)}

{/* Excited States-specific fields */}
{experimentConfig?.backendSettings?.method === "EXCITED_STATES" && (
  <>
    <div className="flex justify-between">
      <span className="text-muted-foreground">ES Method:</span>
      <span className="font-quando font-medium">
        {experimentConfig?.backendSettings?.excited_method?.toUpperCase() || "CIS"}
      </span>
    </div>
    <div className="flex justify-between">
      <span className="text-muted-foreground">States:</span>
      <span className="font-quando font-medium">
        {experimentConfig?.backendSettings?.nStates || 5}
      </span>
    </div>
  </>
)}
```

**Expected Behavior**:
- VQE shows: Ansatz, Optimizer
- SQD shows: Subspace, States
- EXCITED_STATES shows: ES Method, States

---

### 4. Method-Aware Max Iterations Display - FIXED ✅
**File**: `web/src/components/simulation/ExperimentMonitor.tsx` (lines 544-565)

**Problem**: Hardcoded "/42" in iteration display for all methods.

**Solution**:
```typescript
{(() => {
  const method = experimentConfig?.backendSettings?.method;
  let maxIterations;

  if (method === "VQE") {
    maxIterations = experimentConfig?.backendSettings?.maxIterations || 100;
  } else if (method === "SQD") {
    maxIterations = 6; // SQD has 7 stages (0-6)
  } else if (method === "EXCITED_STATES") {
    // If using quantum backend, EXCITED_STATES redirects to SQD (7 stages)
    const backend = experimentConfig?.backendSettings?.backend;
    if (backend === "bluequbit" || backend === "ibm_quantum") {
      maxIterations = 6; // Redirected to SQD
    } else {
      maxIterations = experimentConfig?.backendSettings?.nStates || 5;
    }
  } else {
    maxIterations = 100;
  }

  return `Iteration ${currentIteration} / ${maxIterations}`;
})()}
```

**Expected Behavior**:
- VQE: Shows "Iteration X / 100" (or configured maxIterations)
- SQD: Shows "Iteration X / 6"
- EXCITED_STATES with quantum: Shows "Iteration X / 6"
- EXCITED_STATES with classical: Shows "Iteration X / 5"

---

## Issues NOT Yet Fixed

### 1. SQD Doesn't Use Quantum Circuits ❌
**File**: `kanad/solvers/sqd_solver.py`

**Problem**:
- `qiskit-addon-sqd` is installed and detected (line 92-98)
- But `_has_sqd_addon` flag is **never checked** in `solve()` method
- All computation is classical numpy (lines 307-414)
- No quantum circuits are actually generated or executed

**Impact**:
- Backend selection (BlueQubit/IBM) is initialized but never used
- User selects quantum backend but gets classical simulation
- This is the "huge carelessness in backend" user complained about

**To Fix** (requires significant refactoring):
1. Implement quantum circuit generation using qiskit-addon-sqd
2. Check `_has_sqd_addon` flag in `solve()` method
3. Actually use `self._bluequbit_backend` or `self._ibm_backend` for execution
4. Submit circuits to chosen backend instead of classical numpy operations

**Code Location**:
```python
# kanad/solvers/sqd_solver.py:92-98
try:
    import qiskit_addon_sqd
    self._has_sqd_addon = True  # <-- SET BUT NEVER CHECKED!
    logger.info("qiskit-addon-sqd available")
except ImportError:
    self._has_sqd_addon = False
    logger.warning("qiskit-addon-sqd not installed. Using simplified implementation.")

# kanad/solvers/sqd_solver.py:316-414
def solve(self, n_states: int = 3, callback=None) -> Dict[str, Any]:
    # ... entire method uses numpy, no quantum circuits!
    H_sub = self._project_hamiltonian(basis)  # Classical matrix operations
    eigenvalues, eigenvectors = np.linalg.eigh(H_sub)  # Classical diagonalization
    # Should use qiskit-addon-sqd to generate circuits and submit to backend!
```

---

### 2. EXCITED_STATES Redirect is Confusing UX ⚠️
**File**: `api/services/experiment_service.py` (lines 585-594)

**Problem**:
- User selects "EXCITED_STATES" method
- User selects "bluequbit" backend
- Backend silently redirects to SQD execution
- Frontend shows "Method: SQD" instead of "EXCITED_STATES"
- No clear indication to user that their choice was overridden

**Current Implementation**:
```python
# api/services/experiment_service.py:585-594
# Check if using quantum backend
backend_type = config.get('backend', 'classical')

# For quantum backends, use SQD instead of classical excited states
if backend_type in ['bluequbit', 'ibm_quantum']:
    print(f"⚠️  Excited States with quantum backend detected - using SQD instead")
    print(f"📊 Quantum excited states computed via Subspace Quantum Diagonalization")

    # Redirect to SQD execution
    return execute_sqd(molecule, config, job_id, experiment_id)
```

**Better UX Options**:
1. Show clear warning/info message in frontend
2. Display as "EXCITED_STATES (via SQD)" instead of just "SQD"
3. Disable quantum backend selection for EXCITED_STATES in settings
4. Add tooltip explaining EXCITED_STATES only supports classical methods

---

## Summary

### ✅ Fixed in This Session
1. Iteration counter reset on new experiment
2. Iteration counter set correctly on completion
3. Method-specific configuration display (VQE/SQD/EXCITED_STATES)
4. Method-aware max iterations display
5. Debug logging for analysis data flow

### ❓ Needs Browser Console Verification
1. Analysis data display (data exists in backend and database, but frontend display unknown)

### ❌ Known Issues (Not Fixed)
1. SQD doesn't actually use quantum circuits (always runs classically)
2. EXCITED_STATES redirect to SQD is confusing UX
3. qiskit-addon-sqd installed but not integrated

### 📋 Next Steps
1. **Immediate**: Run new experiment and check browser console logs (see `BROWSER_CONSOLE_CHECK.md`)
2. **Short-term**: Fix analysis display issue once root cause identified
3. **Long-term**: Implement actual quantum SQD using qiskit-addon-sqd
4. **Long-term**: Improve EXCITED_STATES UX (warning message or disable quantum backend option)

---

## Files Modified
- `web/src/components/simulation/ExperimentMonitor.tsx`
- `web/src/components/experiment/ExperimentReport.tsx` (previous session)
- `web/src/components/settings/SettingsModal.tsx` (previous session)
- `web/src/lib/types.ts` (previous session)
- `api/services/experiment_service.py` (debug logging added)
- `kanad/solvers/sqd_solver.py` (callback support added)

## New Documentation
- `BROWSER_CONSOLE_CHECK.md` - Browser debugging instructions
- `FIXES_APPLIED.md` - This file
