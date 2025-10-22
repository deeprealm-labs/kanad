# Cloud Backend Issues - Fixed

## Problems Identified

### 1. Max Iterations Not Respected ✅ FIXED
**Symptom**: User set `maxIterations=100` but experiment ran 4150 iterations

**Root Cause**:
- Frontend `BackendSettings` type didn't include `maxIterations` field
- Dashboard state didn't include `maxIterations`
- Function evaluations were being counted instead of optimizer iterations

**Fixes Applied**:
1. Added `maxIterations?: number` to [BackendSettings](web/src/lib/types.ts#L28) and [UserSettings](web/src/lib/types.ts#L159)
2. Added default `maxIterations: 100` to [dashboard state](web/src/app/dashboard/page.tsx#L29)
3. Fixed iteration counting in [VQE solver](kanad/solvers/vqe_solver.py#L825) to use `result.nit` (optimizer iterations) instead of function evaluations

### 2. Excessive Quantum Job Submissions ⚠️ WARNING ADDED
**Symptom**: 4150 function evaluations when only 100 iterations expected

**Root Cause**:
- Gradient-based optimizers (SLSQP, L-BFGS-B) call objective function ~40x per iteration
- Each function evaluation = 1 quantum job submission to IBM/BlueQubit
- With `maxiter=100`, this creates ~4000 quantum jobs!
- **This is extremely expensive and slow for cloud backends**

**Fix Applied**:
- Added warning in [VQE solver](kanad/solvers/vqe_solver.py#L811-L817) when using gradient-based optimizers with cloud backends
- Recommends using COBYLA or POWELL optimizers for fewer job submissions

**Example**:
```
WARNING: Using gradient-based optimizer 'SLSQP' with cloud backend 'ibm'.
This will submit many quantum jobs (~4000 function evaluations).
Consider using COBYLA or POWELL for fewer job submissions.
```

### 3. No Analysis Data Available ✅ FIXED
**Symptom**: "No analysis data available" despite `enable_analysis=True`

**Root Cause**:
- For multi-atom molecules, VQE uses "components mode" API
- Components mode didn't initialize analysis tools (`EnergyAnalyzer`, `BondingAnalyzer`, etc.)
- When `_add_analysis_to_results()` tried to access analyzers, they didn't exist

**Fix Applied**:
- Added analysis tool initialization in [components mode](kanad/solvers/vqe_solver.py#L214-L225)
- Now initializes analyzers when `enable_analysis=True` and molecule is available

## Technical Details

### Iteration Counting Fix
**Before**:
```python
'iterations': self.iteration_count  # Function evaluations (4150)
```

**After**:
```python
'iterations': result.nit if hasattr(result, 'nit') else self.iteration_count,  # Optimizer iterations (~100)
'function_evaluations': self.iteration_count,  # Separate field for function evals
```

### Analysis Tools Initialization
**Added** to `_init_from_components_mode()`:
```python
if self.enable_analysis and self.molecule is not None:
    try:
        from kanad.analysis import EnergyAnalyzer, BondingAnalyzer, MolecularPropertyAnalyzer
        self.energy_analyzer = EnergyAnalyzer(self.hamiltonian)
        self.bonding_analyzer = BondingAnalyzer(self.molecule)
        self.property_analyzer = MolecularPropertyAnalyzer(self.molecule)
        self.atoms = self.molecule.atoms
        logger.info("Analysis tools initialized in components mode")
    except Exception as e:
        logger.warning(f"Failed to initialize analysis tools: {e}")
        self.enable_analysis = False
```

## Recommendations for Cloud Backends

### Optimizer Selection
For IBM Quantum and BlueQubit backends:

| Optimizer | Function Evals/Iteration | Recommended for Cloud? |
|-----------|-------------------------|------------------------|
| **COBYLA** | ~1-2 | ✅ **Best choice** |
| **POWELL** | ~2-5 | ✅ **Good choice** |
| SLSQP | ~40 | ❌ Too many jobs |
| L-BFGS-B | ~50 | ❌ Too many jobs |

### Max Iterations
For cloud backends, use lower `maxIterations` values:
- **Local classical**: 100-1000 iterations OK
- **IBM Quantum**: 10-50 iterations (each = 40+ jobs with SLSQP!)
- **BlueQubit**: 20-100 iterations (faster than IBM but still costs credits)

### Cost Estimation
With `maxIterations=100` and SLSQP optimizer:
- **Function evaluations**: ~4000
- **IBM jobs submitted**: ~4000 (queue time: hours to days!)
- **BlueQubit credits**: ~4000 credits

With `maxIterations=20` and COBYLA optimizer:
- **Function evaluations**: ~40
- **IBM jobs submitted**: ~40 (queue time: 30min - 2hrs)
- **BlueQubit credits**: ~40 credits

## Files Modified

1. [web/src/lib/types.ts](web/src/lib/types.ts#L28) - Added `maxIterations` field
2. [web/src/app/dashboard/page.tsx](web/src/app/dashboard/page.tsx#L29) - Added default `maxIterations: 100`
3. [kanad/solvers/vqe_solver.py:825](kanad/solvers/vqe_solver.py#L825) - Fixed iteration counting
4. [kanad/solvers/vqe_solver.py:811-817](kanad/solvers/vqe_solver.py#L811-L817) - Added optimizer warning
5. [kanad/solvers/vqe_solver.py:214-225](kanad/solvers/vqe_solver.py#L214-L225) - Initialize analysis tools

## Testing

Run with cloud backend:
```bash
# Should now show:
# - Correct iteration count (100 optimizer iterations, not 4000)
# - Warning about SLSQP + cloud backend
# - Analysis data populated (dipole moment, bond analysis, etc.)
```

## Next Steps

Consider adding to GUI:
1. Optimizer selector with recommendations for cloud backends
2. Warning when selecting gradient-based optimizer + cloud backend
3. Cost estimator showing expected job submissions
4. Backend-specific defaults (COBYLA for cloud, SLSQP for classical)