# Spectroscopy Module - Placeholder Analysis

**Date**: November 7, 2025
**Status**: Issues Identified - Fixes Needed

---

## Executive Summary

The spectroscopy module has placeholders and missing features that affect the accuracy of excited state calculations. While the quantum methods (VQE/SQD) are properly implemented, they lack oscillator strength calculations, and the CIS method uses placeholder values for two-electron integrals.

---

## Critical Issues Found

### 1. CIS Two-Electron Integrals (CRITICAL)

**Location**: [kanad/solvers/excited_states_solver.py:190-196](kanad/solvers/excited_states_solver.py#L190-L196)

**Problem**:
```python
# Off-diagonal: two-electron integrals
# Simplified: use approximation for now
# Full implementation would use ERI tensor
if i == j and a == b:
    # Coulomb integral (approximate)
    A[ia, jb] += 0.1  # Placeholder
if i == j or a == b:
    # Exchange integral (approximate)
    A[ia, jb] -= 0.05  # Placeholder
```

**Impact**:
- CIS excitation energies will be **completely wrong**
- Users cannot trust CIS results
- Should use proper electron repulsion integrals (ERI)

**What Should Be Done**:
```python
# Correct CIS matrix elements:
# A[ia,jb] = δ_ij δ_ab (ε_a - ε_i) + 2(ia|jb) - (ij|ab)
#
# Where:
# (ia|jb) = ∫∫ φ_i*(r1) φ_a(r1) (1/r12) φ_j*(r2) φ_b(r2) dr1 dr2
#
# These can be obtained from PySCF:
from pyscf import ao2mo
eri_mo = ao2mo.kernel(mol, mo_coeffs)

# Then compute proper CIS matrix:
for i in range(n_occ):
    for a in range(n_occ, n_orb):
        ia = idx_map[(i, a)]
        for j in range(n_occ):
            for b in range(n_occ, n_orb):
                jb = idx_map[(j, b)]

                # Diagonal
                if ia == jb:
                    A[ia, jb] = mo_energies[a] - mo_energies[i]

                # Coulomb integral: 2(ia|jb)
                A[ia, jb] += 2 * eri_mo[i, a, j, b]

                # Exchange integral: -(ij|ab)
                A[ia, jb] -= eri_mo[i, j, a, b]
```

**Estimated Fix Time**: 30 minutes

---

### 2. Missing Oscillator Strengths for VQE (HIGH PRIORITY)

**Location**: [kanad/solvers/excited_states_solver.py:591](kanad/solvers/excited_states_solver.py#L591)

**Problem**:
```python
'oscillator_strengths': np.zeros(len(excitation_energies_ev)),  # Not computed in VQE
```

**Impact**:
- Cannot predict absorption intensity
- UV-Vis spectra from VQE will have no intensity information
- All transitions appear equally likely (incorrect)

**What Should Be Done**:

Oscillator strength formula:
```
f = (2/3) * ΔE * |⟨ψ_0|μ|ψ_i⟩|²
```

Where:
- ΔE = excitation energy
- μ = dipole moment operator
- ψ_0 = ground state wavefunction
- ψ_i = excited state wavefunction

**Implementation Approach**:
1. Store ground state and excited state wavefunctions from VQE
2. Compute transition dipole moments ⟨ψ_0|μ|ψ_i⟩ using statevector
3. Compute oscillator strengths using formula above
4. For hardware backends, use finite difference method on molecular geometry

**Estimated Fix Time**: 2-3 hours

---

### 3. Missing Oscillator Strengths for SQD (HIGH PRIORITY)

**Location**: [kanad/solvers/excited_states_solver.py:361](kanad/solvers/excited_states_solver.py#L361)

**Problem**:
```python
'oscillator_strengths': np.zeros(len(excitation_energies_ev)),  # Not computed in SQD
```

**Impact**: Same as VQE issue above

**What Should Be Done**:
Same approach as VQE - compute transition dipole moments from SQD eigenvectors

**Estimated Fix Time**: 2-3 hours (similar to VQE)

---

### 4. Vibronic Spectrum Approximations (MEDIUM PRIORITY)

**Location**: [kanad/analysis/spectroscopy.py:1022-1033](kanad/analysis/spectroscopy.py#L1022-L1033)

**Problems**:

```python
# Step 3: Compute excited state frequencies
# NOTE: For now, approximate excited frequencies as similar to ground
# TODO: Implement excited state Hessian calculation
if verbose:
    print(f"\n📊 Step 3/4: Estimating excited state frequencies...")
    print(f"⚠️  Note: Using approximate excited state frequencies")
    print(f"         Future versions will compute exact excited state Hessian")

# Approximate: excited frequencies similar to ground, with small shifts
excited_frequencies = ground_frequencies * 0.95  # Typically lower in excited state

# Approximate displacement (geometry change)
# For simple systems, use bond length change
displacement = np.random.uniform(0.1, 0.5, len(ground_frequencies))  # Approximate
```

**Impact**:
- Vibronic spectra will be approximate
- Franck-Condon factors will be incorrect
- Vibrational progressions will not match experiment

**What Should Be Done**:
1. Implement excited state Hessian calculation
   - Use excited state solver (VQE/SQD)
   - Compute second derivatives of excited state energy
   - Diagonalize Hessian to get excited state frequencies

2. Compute proper displacement
   - Optimize geometry in ground state
   - Optimize geometry in excited state
   - Displacement = difference in normal mode coordinates

**Estimated Fix Time**: 1-2 days (complex feature)

---

### 5. ExcitedStateSolver in spectroscopy.py (MINOR - REDUNDANT CODE)

**Location**: [kanad/analysis/spectroscopy.py:515-664](kanad/analysis/spectroscopy.py#L515-L664)

**Problem**:
There's a duplicate `ExcitedStateSolver` class in spectroscopy.py that has placeholder implementations. However, this class is **NOT USED** - the actual implementation is in `kanad/solvers/excited_states_solver.py`.

The `compute_excited_states_vqe` method in spectroscopy.py (lines 537-664) has issues:
- Lines 609: "This is a simplified approach"
- Lines 634: Returns placeholder `E_ground + 0.1`
- Lines 640: Returns placeholder `E_ground + 0.1 * i`

**Impact**: None - this code is not used by the framework

**What Should Be Done**:
- **Option 1**: Delete the duplicate class entirely
- **Option 2**: Deprecate it with warning
- **Option 3**: Keep as legacy/example code with clear warnings

**Recommended**: Delete it to avoid confusion

**Estimated Fix Time**: 5 minutes

---

## Summary of Placeholders

| Issue | Severity | Location | Impact | Fix Time |
|-------|----------|----------|--------|----------|
| CIS two-electron integrals | **CRITICAL** | excited_states_solver.py:190-196 | CIS gives wrong energies | 30 min |
| VQE oscillator strengths | HIGH | excited_states_solver.py:591 | No intensity info | 2-3 hrs |
| SQD oscillator strengths | HIGH | excited_states_solver.py:361 | No intensity info | 2-3 hrs |
| Vibronic excited Hessian | MEDIUM | spectroscopy.py:1022-1033 | Approximate vibronic | 1-2 days |
| Duplicate ExcitedStateSolver | MINOR | spectroscopy.py:515-664 | Confusion only | 5 min |

---

## What Actually Works (No Placeholders)

### ✅ Quantum SQD Excited States
- [excited_states_solver.py:280-376](kanad/solvers/excited_states_solver.py#L280-L376)
- Uses proper SQD solver
- Computes multiple eigenvalues correctly
- Returns accurate excitation energies
- **Only missing**: oscillator strengths

### ✅ Quantum VQE Excited States
- [excited_states_solver.py:383-601](kanad/solvers/excited_states_solver.py#L383-L601)
- Uses orthogonality-constrained VQE (proper method)
- Computes statevector overlaps
- Adds penalties to avoid ground state
- **Only missing**: oscillator strengths
- **Limitation**: Cannot reach high-energy states (> 10 eV) due to ansatz variational principle

### ✅ UV-Vis Calculator with Quantum SQD
- [spectroscopy.py:221-359](kanad/analysis/spectroscopy.py#L221-L359)
- Fully implemented quantum UV-Vis
- Uses ExcitedStatesSolver properly
- Generates spectra with Gaussian broadening
- **Only missing**: oscillator strengths (inherited from solver)

### ✅ Classical TD-DFT/TDA
- [spectroscopy.py:112-219](kanad/analysis/spectroscopy.py#L112-L219)
- Uses PySCF properly
- Computes oscillator strengths correctly
- Production-ready

---

## Recommended Fix Priority

### Immediate (This Week):
1. **Fix CIS two-electron integrals** - CRITICAL
   - Current CIS is completely wrong
   - Easy fix (use PySCF ERI)
   - 30 minutes

2. **Delete duplicate ExcitedStateSolver** - MINOR
   - Avoid confusion
   - 5 minutes

### Short Term (Next Week):
3. **Implement oscillator strengths for VQE/SQD** - HIGH
   - Essential for meaningful spectra
   - 4-6 hours total
   - Compute transition dipole moments

### Medium Term (Next Month):
4. **Implement excited state Hessian** - MEDIUM
   - Improves vibronic spectra
   - 1-2 days
   - Complex but valuable feature

---

## Testing Requirements

After fixes:

### 1. CIS Validation
```python
# Test CIS against known H2 excitation (11.4 eV to first excited state)
from kanad.bonds import BondFactory
from kanad.solvers import ExcitedStatesSolver

bond = BondFactory.create_bond('H', 'H', distance=0.74)
solver = ExcitedStatesSolver(bond, method='cis', n_states=3)
result = solver.solve()

assert abs(result['excitation_energies'][0] - 11.4) < 1.0  # Within 1 eV
```

### 2. Oscillator Strength Validation
```python
# Test that oscillator strengths are non-zero and sum correctly
result = solver.solve()
f_sum = sum(result['oscillator_strengths'])

assert f_sum > 0, "Oscillator strengths should not all be zero"
assert all(f >= 0 for f in result['oscillator_strengths']), "All f should be non-negative"
```

### 3. Vibronic Spectrum Validation
```python
# Test that vibronic spectrum has Franck-Condon structure
from kanad.analysis import VibronicCalculator

vibro = VibronicCalculator(molecule)
spectrum = vibro.compute_quantum_vibronic_spectrum(n_states=1)

assert 'fc_factors' in spectrum
assert len(spectrum['fc_factors']['franck_condon_factors']) > 1
```

---

## API Changes (After Fixes)

### Before (Current):
```python
result = excited_solver.solve()
# result['oscillator_strengths'] = [0, 0, 0, ...]  # All zeros!
```

### After (Fixed):
```python
result = excited_solver.solve()
# result['oscillator_strengths'] = [0.025, 0.341, 0.002, ...]  # Real values!
```

---

## Conclusion

**Current State**:
- Quantum excited states (VQE/SQD) are **properly implemented**
- CIS has **critical bug** (placeholder integrals)
- Oscillator strengths are **missing** for all quantum methods
- Vibronic spectra use **approximations** but are functional

**After Fixes**:
- CIS will give correct energies (30 min fix)
- All methods will have oscillator strengths (4-6 hr fix)
- Vibronic spectra will be accurate (1-2 day fix)

**Priority**: Fix CIS first (critical), then oscillator strengths (high impact), then vibronic improvements (nice to have).

---

**Document Created**: November 7, 2025
**Status**: Analysis Complete - Ready for Fixes
**Next Action**: Fix CIS two-electron integrals
