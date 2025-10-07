# Kanad Framework Validation Suite

Comprehensive validation scripts for testing all solvers, mappers, Hamiltonians, and ansatze.

## Validation Scripts

### 1. VQE Solver Validation (`01_vqe_solver_validation.py`)

Tests the Variational Quantum Eigensolver with multiple configurations:

**Configurations Tested:**
- **Ansatze**: Governance, Real-Amplitudes, Hardware-Efficient
- **Mappers**: Jordan-Wigner, Bravyi-Kitaev
- **Optimizers**: SLSQP, Powell, COBYLA
- **Hamiltonians**: Covalent (H2), Ionic (LiH)

**Validation Criteria:**
- Energy within 50 mHa of reference for small molecules
- Governance ansatz achieves best accuracy
- Different optimizers converge to same energy

**Expected Results:**
- H2 with Governance + SLSQP: 0-5 mHa error
- H2 with Real-Amplitudes: < 50 mHa error
- LiH with Governance: < 100 mHa error

---

### 2. SQD Solver Validation (`02_sqd_solver_validation.py`)

Tests Subspace Quantum Diagonalization for ground and excited states:

**Configurations Tested:**
- Different subspace dimensions (10, 15, 20)
- Different circuit depths (3, 4)
- Comparison with exact diagonalization

**Validation Criteria:**
- Ground state energy matches exact within 50 mHa
- Multiple eigenvalues are properly ordered
- Larger subspace gives better accuracy

**Expected Results:**
- SQD with dim=20 should match exact < 10 mHa
- All excited states > ground state
- Convergence improves with subspace size

---

### 3. Excited States Solver Validation (`03_excited_states_validation.py`)

Tests excited state calculations with multiple methods:

**Methods Tested:**
- CIS (Configuration Interaction Singles)
- QPE (Quantum Phase Estimation)
- Comparison with exact diagonalization

**Validation Criteria:**
- All excitation energies are positive
- States are properly ordered (E₀ < E₁ < E₂ < ...)
- Methods agree with each other

**Expected Results:**
- CIS ground state < 100 mHa error
- Excited states properly ordered
- QPE matches CIS within reasonable tolerance

---

### 4. Mapper Comparison (`04_mapper_comparison.py`)

Tests that different qubit mappers produce consistent results:

**Mappers Tested:**
- Jordan-Wigner (JW)
- Bravyi-Kitaev (BK)
- Parity (if available)

**Validation Criteria:**
- All mappers produce same ground state energy (< 1 mHa difference)
- Qubit requirements are correct
- Works with different Hamiltonians (covalent, ionic)

**Expected Results:**
- JW and BK energies differ < 1 mHa
- Both work on H2 and LiH
- Hamiltonian structure is correct

---

### 5. Hamiltonian Comparison (`05_hamiltonian_comparison.py`)

Tests different Hamiltonian types and governance protocols:

**Hamiltonians Tested:**
- CovalentHamiltonian (H2)
- IonicHamiltonian (LiH)
- MetallicHamiltonian (Na₂)

**Validation Criteria:**
- Governance protocols are correctly assigned
- Energy values are physically reasonable
- Hamiltonian properties match bond type
- Basis sets propagate correctly

**Expected Results:**
- Each Hamiltonian has correct protocol
- LiH needs more qubits than H2
- All three protocols have different entanglement strategies
- Basis propagates from bond to Hamiltonian

---

## Running Validations

### Run All Validations

```bash
cd tests/validation
python run_all_validations.py
```

This will:
1. Run all validation scripts sequentially
2. Capture output and results
3. Generate comprehensive summary report
4. Exit with code 0 if all pass, 1 if any fail

### Run Individual Validation

```bash
cd tests/validation
python 01_vqe_solver_validation.py
```

Each script is standalone and reports:
- ✓/✗ status for each test
- Energy errors in mHa
- Detailed comparison with references
- Overall pass/fail status

---

## Interpreting Results

### Success Criteria

**VQE Solver:**
- ✓ 5/5 tests pass
- H2 energies within 50 mHa
- All configurations work

**SQD Solver:**
- ✓ 3/4 tests pass (50% acceptable)
- Ground state within 50 mHa
- Eigenvalues properly ordered

**Excited States:**
- ✓ 4/6 tests pass (50% acceptable)
- Excitation energies positive
- Methods agree

**Mappers:**
- ✓ All tests pass
- JW and BK differ < 1 mHa
- Both work on all molecules

**Hamiltonians:**
- ✓ 7/8 tests pass (70% acceptable)
- All protocols working
- Governance applied correctly

### Failure Investigation

If a validation fails:

1. **Check error messages** - Look for specific exceptions
2. **Compare energies** - Are they close but outside tolerance?
3. **Check qubits** - Do qubit counts make sense?
4. **Test dependencies** - Is Qiskit/NumPy working?
5. **Run with traceback** - Add `-v` flag for verbose output

---

## Reference Energies

### H2 (0.74 Å, sto-3g)
- **Exact**: -1.137284 Ha
- **HF**: -1.116759 Ha
- **Correlation**: -0.020525 Ha

### LiH (1.56 Å, sto-3g)
- **Exact**: -7.882324 Ha
- **HF**: -7.861865 Ha
- **Correlation**: -0.020459 Ha

### Energy Tolerances

- **Excellent**: < 1 mHa (chemical accuracy)
- **Good**: < 10 mHa
- **Acceptable**: < 50 mHa (small molecules)
- **Acceptable**: < 100 mHa (larger molecules)

---

## Adding New Validations

To add a new validation script:

1. Create `0X_your_validation.py` in this directory
2. Follow the existing format:
   - Print clear section headers
   - Validate against known references
   - Track results in list
   - Print summary at end
   - Include "PASSED" or "FAILED" in output

3. Add to `validation_scripts` list in `run_all_validations.py`

4. Document in this README

---

## Expected Runtime

- **VQE Validation**: ~60-120s
- **SQD Validation**: ~30-60s
- **Excited States**: ~20-40s
- **Mapper Comparison**: ~60-90s
- **Hamiltonian Comparison**: ~40-60s

**Total**: ~3-6 minutes for full validation suite

---

## Continuous Integration

These validations can be run in CI/CD:

```yaml
- name: Run Kanad Validations
  run: |
    cd tests/validation
    python run_all_validations.py
```

Exit code 0 = all pass, 1 = at least one failure.
