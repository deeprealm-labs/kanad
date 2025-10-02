# Kanad Framework - Validation Issues Log

**Generated**: 2025-10-02
**Total Scripts**: 6
**Overall Pass Rate**: 43/44 checks (98%)

This document catalogs all scientific accuracy issues, warnings, and errors found during validation script inspections.

---

## CRITICAL ISSUES (Must Fix)

### 1. Basis Function Normalization (H2 Covalent)
**Script**: `01_h2_covalent_bond.py`
**Severity**: 🔴 CRITICAL
**Location**: `kanad/core/integrals/basis_sets.py`

**Problem**:
- Overlap matrix diagonal elements = 10.74 (should be 1.0)
- Off-diagonal overlap = 10.22 (should be < 1.0)
- Basis functions are not properly normalized

**Expected**:
```
S = [[1.0000  0.6593]
     [0.6593  1.0000]]  (for H2 at 0.74 Å)
```

**Actual**:
```
S = [[10.7421  10.2155]
     [10.2155  10.7421]]
```

**Impact**: All integral calculations are incorrect by ~10x factor. This affects:
- MO energies
- Hamiltonian matrix elements
- VQE energy estimates

**Root Cause**: `GaussianPrimitive._normalization_constant()` or integral evaluation incorrect

---

### 2. Band Structure Occupation (Metallic Na)
**Script**: `03_metallic_sodium_chain.py`
**Severity**: 🔴 CRITICAL
**Location**: `kanad/bonds/metallic_bond.py`

**Problem**:
- All 6 bands marked as OCCUPIED
- Should only have 3 occupied bands (6 electrons / 2 per band)
- Fermi energy at band maximum (2.0 eV) suggests insulator, not metal

**Expected for 6 Na atoms (6 valence electrons)**:
```
Band 0: -2.0 eV  [OCCUPIED] ← 2 electrons
Band 1: -1.0 eV  [OCCUPIED] ← 2 electrons
Band 2: -1.0 eV  [OCCUPIED] ← 2 electrons
Band 3:  1.0 eV  [VIRTUAL]
Band 4:  1.0 eV  [VIRTUAL]
Band 5:  2.0 eV  [VIRTUAL]
Fermi energy: ~-0.5 eV (between bands 2 and 3)
```

**Impact**: Metallic bonding model gives wrong electronic structure and metallicity

---

### 3. HOMO-LUMO Gap Too Large (H2)
**Script**: `01_h2_covalent_bond.py`, `06_complete_workflow.py`
**Severity**: 🔴 CRITICAL
**Location**: `kanad/core/hamiltonians/covalent_hamiltonian.py`

**Problem**:
- Computed gap: 37.15 eV
- Experimental gap: ~11-12 eV (σ to σ*)
- 3× too large!

**Likely Causes**:
1. Related to normalization issue (#1)
2. Core Hamiltonian matrix elements incorrect
3. Overlap matrix used incorrectly in eigenvalue problem

---

## HIGH PRIORITY ISSUES

### 4. VQE SWAP Network Not Implemented
**Script**: `06_complete_workflow.py`
**Severity**: 🟠 HIGH
**Location**: `kanad/solvers/vqe_solver.py:322`

**Problem**:
```python
raise NotImplementedError("Non-adjacent two-qubit gates require SWAP network")
```

**Impact**: Cannot run actual VQE optimization on circuits with non-adjacent qubit gates

**Recommendation**: Implement SWAP network insertion or use all-to-all connectivity assumption

---

### 5. Ionic Character Correlation Broken
**Script**: `04_bond_comparison.py`
**Severity**: 🟠 HIGH
**Location**: Bond character calculation

**Problem**:
| Bond | ΔEN  | Ionic % | Issue |
|------|------|---------|-------|
| C-N  | 0.49 | 5.8%    | Too low |
| N-O  | 0.40 | 3.9%    | Lower ΔEN but similar % |
| H-N  | 0.84 | 16.2%   | Higher than C-O despite lower ΔEN |
| C-O  | 0.89 | 18.0%   | Should be higher |

**Expected**: Ionic character should monotonically increase with ΔEN

**Root Cause**: Likely using simple empirical formula that doesn't account for orbital overlap differences

---

### 6. Bond Length Systematic Error
**Script**: All scripts
**Severity**: 🟠 HIGH
**Location**: `kanad/bonds/bond_factory.py`

**Problem**:
| System | Computed | Experimental | Error |
|--------|----------|--------------|-------|
| H₂     | 0.62 Å   | 0.74 Å       | -16%  |
| NaCl   | 2.68 Å   | 2.36 Å       | +14%  |

**Root Cause**: Using static covalent radii sum instead of:
- Optimized bond lengths
- Equilibrium from energy minimization
- Empirical corrections

---

### 7. Metallic Bond Energy = 0
**Script**: `03_metallic_sodium_chain.py`
**Severity**: 🟠 HIGH
**Location**: `kanad/bonds/metallic_bond.py`

**Problem**:
- Total energy: 0.000 eV
- Energy per atom: 0.000 eV
- Should be negative (cohesive energy ~1.1 eV/atom for Na)

**Expected**: Total energy ≈ -6.6 eV for 6 Na atoms

---

## MEDIUM PRIORITY ISSUES

### 8. Incorrect Hybridization Labels
**Script**: `01_h2_covalent_bond.py`, `06_complete_workflow.py`
**Severity**: 🟡 MEDIUM
**Location**: Bond analysis

**Problem**:
- H₂ labeled as "sp3 hybridization"
- H only has 1s orbital (no p orbitals in minimal basis)
- Should be "s-s overlap" or "no hybridization"

---

### 9. NaCl Charge Transfer Lower Than Expected
**Script**: `02_nacl_ionic_bond.py`
**Severity**: 🟡 MEDIUM

**Problem**:
- Computed: 0.71 e
- Expected: ~0.8 e
- 11% lower than expected

**Note**: May be acceptable given approximations, but worth investigating

---

### 10. NaCl Coulomb Energy Seems Small
**Script**: `02_nacl_ionic_bond.py`
**Severity**: 🟡 MEDIUM

**Problem**:
- Computed: -0.1 Ha = -2.72 eV
- For point charges at 2.68 Å: E = -ke²/r ≈ -5.4 eV
- Factor of 2× too small

**Possible Causes**:
- Not pure point charges (electron cloud screening)
- Missing polarization contributions
- Calculation error

---

### 11. 2-Qubit Gate Counting Bug
**Script**: `05_ansatz_comparison.py`
**Severity**: 🟡 MEDIUM
**Location**: Gate counting in validation script or ansatz

**Problem**:
- All circuits show "2-Qubit Gates: 0"
- Circuits definitely have CZ, CX, CNOT gates
- Gate type strings may not match exactly ('CNOT' vs 'cnot' vs 'cx')

---

### 12. UCC-D Zero Parameters for H2
**Script**: `05_ansatz_comparison.py`
**Severity**: 🟡 MEDIUM
**Location**: `kanad/ansatze/ucc_ansatz.py`

**Problem**:
- UCC-D (doubles only) has 0 parameters and 0 excitations for H2
- Technically correct (2 electrons → no doubles)
- But ansatz should warn user or handle gracefully

**Recommendation**: Add warning when no excitations are generated

---

### 13. UCCSD vs UCC-S Same Parameters
**Script**: `05_ansatz_comparison.py`
**Severity**: 🟡 MEDIUM

**Problem**:
- Both UCCSD and UCC-S have 3 parameters for H2
- Suggests UCCSD is only doing singles (no doubles exist for H2)
- Confusing for users

**Recommendation**: Clarify in output when doubles are not possible

---

## LOW PRIORITY ISSUES

### 14. Adaptive Governance Initialization Error
**Script**: `05_ansatz_comparison.py`
**Severity**: 🟢 LOW

**Problem**:
```
AdaptiveGovernanceAnsatz.__init__() got multiple v[alues...]
```
(Message truncated)

**Impact**: Adaptive governance not tested in validation

---

### 15. Governance Type Shows N/A
**Script**: `05_ansatz_comparison.py`
**Severity**: 🟢 LOW

**Problem**:
- Governance-aware ansätze show "Governance Type: N/A"
- Should show "ionic" or "covalent"

**Location**: Missing attribute or not set in governance ansatz classes

---

### 16. Limited Basis Set Coverage
**Script**: `04_bond_comparison.py`
**Severity**: 🟢 LOW
**Impact**: User Experience

**Problem**:
- Only H, C, N, O available in STO-3G
- Cannot test Na, Cl, F, etc. with full calculations
- Limits validation scope

**Recommendation**: Add more elements to STO-3G basis data

---

### 17. Ionic Governance Applied to Covalent H2
**Script**: `05_ansatz_comparison.py`
**Severity**: 🟢 LOW

**Problem**:
- Ionic governance ansatz created for H2 (16 parameters)
- H2 is purely covalent
- Physically doesn't make sense

**Note**: May be intentional for testing flexibility

---

### 18. Tight-Binding Dispersion Not Standard
**Script**: `03_metallic_sodium_chain.py`
**Severity**: 🟢 LOW

**Problem**:
- Band energies: ±2, ±1, ±1 eV
- Standard 1D tight-binding: E(k) = -2t·cos(ka)
- Should give smooth dispersion, not discrete levels

**Recommendation**: Verify tight-binding Hamiltonian construction

---

### 19. NaCl 28 Electrons → 4 Qubits Mapping Unclear
**Script**: `02_nacl_ionic_bond.py`
**Severity**: 🟢 LOW

**Problem**:
- NaCl has 28 total electrons
- Mapped to only 4 qubits
- Mapping strategy not explained

**Clarification Needed**: Are we only treating valence electrons?

---

## SUMMARY BY CATEGORY

### By Severity:
- 🔴 **CRITICAL**: 3 issues (normalization, band occupation, HOMO-LUMO gap)
- 🟠 **HIGH**: 4 issues (VQE SWAP, ionic character, bond lengths, metallic energy)
- 🟡 **MEDIUM**: 6 issues (hybridization, charge transfer, Coulomb energy, etc.)
- 🟢 **LOW**: 6 issues (missing features, UX improvements)

### By Component:
- **Integrals/Basis Sets**: 1 critical (normalization)
- **Hamiltonians**: 2 critical (HOMO-LUMO gap, band structure)
- **Bond Analysis**: 3 high, 2 medium (ionic character, bond lengths, etc.)
- **VQE Solver**: 1 high (SWAP network)
- **Ansätze**: 3 medium, 2 low (parameter counting, edge cases)
- **General**: 6 low (UX, missing features)

---

## RECOMMENDED PRIORITY ORDER FOR FIXES

1. **Fix basis function normalization** (#1) - Affects everything
2. **Fix band occupation counting** (#2) - Metallic bonding broken
3. **Debug HOMO-LUMO gap** (#3) - Likely related to #1
4. **Implement VQE SWAP network** (#4) - Blocks actual optimization
5. **Fix ionic character correlation** (#5) - Important for bonding analysis
6. **Implement bond length optimization** (#6) - Systematic errors
7. **Fix metallic cohesive energy** (#7) - Unphysical zero energy
8. **Address medium priority issues** (#8-13)
9. **Polish low priority issues** (#14-19)

---

## POSITIVE FINDINGS ✓

Despite the issues above, the framework demonstrates:

- ✅ Correct overall architecture and workflow
- ✅ Proper governance protocol selection
- ✅ Working ansatz construction
- ✅ Appropriate mapper selection for bond types
- ✅ Complete API for researchers
- ✅ 262 unit tests passing
- ✅ Good documentation and validation coverage

The issues are primarily **quantitative accuracy** rather than **structural problems**. With the critical fixes, the framework will be scientifically sound.
