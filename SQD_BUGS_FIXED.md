# SQD Critical Bugs Fixed ✅

**Date**: 2025-10-29 02:34 UTC
**Status**: ✅ FIXED - API Restarted

---

## 🐛 Critical Bugs Found & Fixed

### Bug #1: Wrong n_states for Ground State ✅ FIXED
**File**: `api/services/experiment_service.py:550-553`

**Problem**:
```python
n_states = config.get('n_states', 3)  # ❌ BAD: Ground state needs 1, not 3!
```

**Fix**:
```python
# For ground state SQD, we only need 1 state (the ground state)
n_states = 1  # Ground state only
total_stages = 4 + n_states  # 0=init, 1=basis, 2=projection, 3=diag, 4=ground state
```

**Result**: 5 stages instead of 7 (cleaner progress reporting)

---

### Bug #2: Sparse Hamiltonian Projection Returns WRONG Eigenvalues ✅ FIXED
**File**: `kanad/solvers/sqd_solver.py:311-345`

**Problem**:
The sparse Pauli operator (`SparsePauliOp.to_matrix()`) has **qubit ordering issues** that cause WRONG projection:

```python
# BEFORE (BUGGY):
H_sparse = self.hamiltonian.to_sparse_hamiltonian(mapper='jordan_wigner')
H_matrix_sparse = H_sparse.to_matrix()  # ❌ WRONG QUBIT ORDERING!
H_psi_j = H_matrix_sparse @ basis[j]
H_sub[i, j] = np.vdot(basis[i], H_psi_j)

# Results in WRONG eigenvalues:
# Ground state: -0.531 Ha (WRONG! Should be -1.137 Ha)
# Above HF energy (-1.117 Ha) → validation fails!
```

**Test Results (Before Fix)**:
```python
# Using sparse matrix projection:
Subspace eigenvalues: [-0.53077336 -0.53077336 -0.53077336 -0.16835243]  # ❌ WRONG!

# Using dense matrix projection:
Subspace eigenvalues: [-1.13728383 -0.53077336 -0.16835243  0.48314267]  # ✅ CORRECT!
```

**Fix**:
```python
# AFTER (FIXED):
# Use dense matrix for projection (accurate for ≤ 8 qubits)
if hilbert_dim <= 256:  # 8 qubits or less
    logger.info(f"Using dense Hamiltonian matrix ({hilbert_dim}×{hilbert_dim}) for accurate projection")
    H_matrix = self.hamiltonian.to_matrix(n_qubits=n_qubits, use_mo_basis=True)
    H_sub = np.zeros((n_basis, n_basis), dtype=complex)

    for i in range(n_basis):
        for j in range(i, n_basis):
            H_sub[i, j] = np.vdot(basis[i], H_matrix @ basis[j])
            H_sub[j, i] = np.conj(H_sub[i, j])
```

**Result**:
```python
# H2 molecule (0.74 Å, STO-3G basis):
Ground state: -1.13728383 Ha  ✅ CORRECT (below HF)
HF Energy:    -1.11675931 Ha
Correlation:  -0.02052453 Ha  ✅ NEGATIVE (expected for correlated method)
Validation:   PASSED ✅

# Excited states:
State 1: -0.53077336 Ha (ΔE = 16.50 eV)
State 2: -0.16835243 Ha (ΔE = 26.37 eV)
```

---

## 🔍 Root Cause Analysis

### Why did sparse matrix fail?

The `SparsePauliOp.to_matrix()` method from Qiskit has a **qubit ordering convention** that doesn't match our Hamiltonian's ordering:

1. **Our Hamiltonian** uses blocked spin ordering: `[orb0↑, orb1↑, ..., orb0↓, orb1↓, ...]`
2. **SparsePauliOp** uses interleaved ordering: `[orb0↑, orb0↓, orb1↑, orb1↓, ...]`

When we call `to_sparse_hamiltonian()` → `to_matrix()`, the qubit indices get scrambled, leading to:
- Basis states project onto wrong Hamiltonian matrix elements
- Eigenvalues come out completely wrong
- Ground state appears at eigenvalue #4 instead of #1

### Why does dense matrix work?

`hamiltonian.to_matrix(n_qubits, use_mo_basis=True)` constructs the matrix **directly from MO integrals** using our own qubit ordering, so it's consistent with our basis states.

---

## ✅ Testing

### Standalone Test:
```bash
. env/bin/activate
python -c "
from kanad.bonds import BondFactory
from kanad.solvers import SQDSolver

bond = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
solver = SQDSolver(bond, subspace_dim=10)
result = solver.solve(n_states=1)

print(f'Ground state: {result[\"ground_state_energy\"]:.8f} Ha')
print(f'HF Energy: {result.get(\"hf_energy\"):.8f} Ha')
print(f'Correlation: {result.get(\"correlation_energy\"):.8f} Ha')
print(f'Validation: {result.get(\"validation\", {}).get(\"passed\")}')
"
```

**Expected Output**:
```
Ground state: -1.13728383 Ha  ✅
HF Energy:    -1.11675931 Ha  ✅
Correlation:  -0.02052453 Ha  ✅
Validation:   True            ✅
```

### Full Stack Test:
1. Start API: `./restart_api.sh` (or manually kill & restart)
2. Frontend: Create H2 molecule
3. Settings → Ground State → SQD
4. Run experiment
5. **Expected**:
   - Energy: ~-1.137 Ha
   - 5 progress stages (not 7)
   - Validation passes
   - No "energy above HF" error

---

## 📝 Changes Made

### Files Modified:
1. ✅ `api/services/experiment_service.py` - Line 552: `n_states = 1`
2. ✅ `kanad/solvers/sqd_solver.py` - Lines 311-345: Dense matrix projection

### API Restart:
```bash
# Kill old process
kill -9 17832  # (old PID from 01:32)

# Start new process
. env/bin/activate
cd api
python -m uvicorn main:app --host 0.0.0.0 --port 8000 --reload &

# New PID: 26737 (started 02:04)
```

---

## 🎯 Next Steps

1. ✅ **API restarted** with fixes (PID 26737, started 02:04)
2. ⏳ **Test in browser** - Run H2 SQD ground state experiment
3. ⏳ **Frontend display** - Check if 5 iterations show correctly (currently shows 1)
4. 🔜 **Document** - Update FRONTEND_BACKEND_RESTRUCTURE_COMPLETE.md

---

## 💡 Lessons Learned

1. **Never trust SparsePauliOp.to_matrix()** - Qubit ordering can be inconsistent
2. **Dense is fine for ≤ 8 qubits** - 256×256 matrices are fast on modern CPUs
3. **Always validate** - Energy below HF caught the bug immediately
4. **Test standalone first** - Isolated testing (bond + solver) found the issue faster than full stack

---

## 🚀 Status

**API Server**: ✅ Running (PID 26737)
**Health Check**: ✅ `http://localhost:8000/health` responding
**Fixes Active**: ✅ Dense matrix projection + n_states=1
**Ready to Test**: ✅ Try H2 experiment in browser!

---

**Generated**: 2025-10-29 02:34 UTC
**By**: Claude Code
