# Phase 1: Fix Failing Tests & Expand Coverage
## Immediate Action Plan - Week 1

**Status:** 🚀 Ready to Execute
**Current:** 321/344 tests passing (93%)
**Target:** 344/344 tests passing (100%)
**Timeline:** 5-7 days

---

## Day 1-2: Fix VQE API Tests (12 tests)

### Issue
VQESolver API changed from accepting `hamiltonian` directly to accepting `bond` object.

### Old API (Failing Tests)
```python
solver = VQESolver(
    hamiltonian=hamiltonian,
    ansatz=ansatz,
    mapper=mapper,
    backend='classical'
)
```

### New API (Required)
```python
from kanad.bonds import BondFactory

bond = BondFactory.create_bond('H', 'H', distance=0.74)
solver = VQESolver(
    bond=bond,
    ansatz_type='ucc',
    mapper_type='jordan_wigner',
    backend='statevector'
)
```

### Tests to Fix

1. `tests/unit/test_vqe.py::TestVQESolver::test_vqe_solver_creation`
2. `tests/unit/test_vqe.py::TestVQESolver::test_vqe_energy_computation`
3. `tests/unit/test_vqe.py::TestVQESolver::test_vqe_optimization`
4. `tests/unit/test_vqe.py::TestVQESolver::test_vqe_with_ucc_ansatz`
5. `tests/unit/test_vqe.py::TestVQESolver::test_vqe_callback`
6. `tests/unit/test_vqe.py::TestVQESolver::test_vqe_energy_history`
7. `tests/unit/test_vqe.py::TestVQESolver::test_vqe_variance_computation`
8. `tests/unit/test_vqe.py::TestVQEIntegration::test_vqe_with_different_optimizers`
9. `tests/unit/test_vqe.py::TestVQEIntegration::test_vqe_parameter_initialization`
10. `tests/unit/test_vqe.py::TestVQEIntegration::test_vqe_reproducibility`
11. `tests/unit/test_vqe.py::TestVQEIntegration::test_vqe_converges_below_initial`
12. `tests/unit/test_qiskit_integration.py::TestQiskitVQE::test_vqe_classical_vs_qiskit`

### Action Items
- [ ] Read current `test_vqe.py` file
- [ ] Update all VQESolver instantiations to use new API
- [ ] Verify tests pass
- [ ] Document API changes in test comments

---

## Day 3-4: Fix Qiskit Integration Tests (7 tests)

### Issue
Import paths changed for Qiskit backend modules.

### Tests to Fix

1. `tests/unit/test_qiskit_integration.py::TestQiskitBackend::test_aer_simulator_init`
2. `tests/unit/test_qiskit_integration.py::TestQiskitBackend::test_aer_statevector_init`
3. `tests/unit/test_qiskit_integration.py::TestQiskitBackend::test_get_estimator`
4. `tests/unit/test_qiskit_integration.py::TestQiskitBackend::test_backend_info`
5. `tests/unit/test_qiskit_integration.py::TestQiskitVQE::test_vqe_classical_backend`
6. `tests/unit/test_qiskit_integration.py::TestQiskitVQE::test_vqe_aer_backend`
7. `tests/unit/test_qiskit_integration.py::TestQiskitVQE::test_vqe_classical_vs_qiskit`

### Possible Issues
- QiskitBackend class moved or renamed
- Import path changes
- API changes in Qiskit backend

### Action Items
- [ ] Investigate current Qiskit backend structure
- [ ] Update import statements
- [ ] Fix any API changes
- [ ] Verify Qiskit compatibility

---

## Day 5: Fix Import Tests (3 tests)

### Tests to Fix

1. `tests/unit/test_imports.py::test_main_import`
2. `tests/unit/test_imports.py::test_solvers_import`
3. `tests/unit/test_imports.py::test_backends_import`
4. `tests/unit/test_imports.py::test_ibm_solvers_import`

### Issue
Some modules were removed or reorganized, causing import assertions to fail.

### Action Items
- [ ] Review current module structure
- [ ] Update import assertions to match reality
- [ ] Remove tests for deprecated modules
- [ ] Add tests for new modules

---

## Day 6-7: Expand Test Coverage

### Priority Areas

#### Core Components
```
tests/unit/test_molecule_extended.py - NEW
tests/unit/test_basis_sets_comprehensive.py - NEW
tests/unit/test_spin_states.py - NEW
tests/unit/test_charged_molecules.py - NEW
```

#### Hamiltonians
```
tests/unit/test_ionic_hamiltonian_extended.py - NEW
tests/unit/test_metallic_hamiltonian_extended.py - NEW
tests/unit/test_hamiltonian_properties.py - NEW (symmetry, hermiticity)
```

#### Mappers
```
tests/unit/test_parity_mapper.py - NEW
tests/unit/test_mapper_equivalence.py - NEW
tests/unit/test_qubit_reduction.py - NEW
```

---

## Tools & Scripts to Create

### 1. Test Runner Script
```bash
#!/bin/bash
# run_failing_tests.sh - Run only failing tests

pytest tests/unit/test_vqe.py::TestVQESolver -v
pytest tests/unit/test_qiskit_integration.py -v
pytest tests/unit/test_imports.py -v
```

### 2. Coverage Reporter
```bash
#!/bin/bash
# check_coverage.sh - Generate coverage report

pytest --cov=kanad --cov-report=html --cov-report=term-missing tests/unit/
open htmlcov/index.html
```

### 3. API Migration Helper
```python
# scripts/migrate_vqe_tests.py - Automated test migration

import re
import sys

def migrate_vqe_test(test_content):
    """Convert old VQE API to new API."""
    # Pattern matching and replacement
    # ...
    return updated_content
```

---

## Success Metrics

### Day 1-2
- [ ] 12 VQE tests fixed
- [ ] Tests passing: 333/344 (97%)

### Day 3-4
- [ ] 7 Qiskit tests fixed
- [ ] Tests passing: 340/344 (99%)

### Day 5
- [ ] 4 import tests fixed
- [ ] Tests passing: 344/344 (100%) ✅

### Day 6-7
- [ ] 20+ new tests added
- [ ] Code coverage > 90%
- [ ] All critical paths tested

---

## Continuous Integration Setup

### GitHub Actions Workflow
```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: '3.9'
      - name: Install dependencies
        run: |
          pip install -e .
          pip install pytest pytest-cov
      - name: Run tests
        run: pytest --cov=kanad tests/unit/
      - name: Upload coverage
        uses: codecov/codecov-action@v2
```

---

## Documentation Updates

### Test Documentation
- [ ] Update CONTRIBUTING.md with testing guidelines
- [ ] Add test writing guide
- [ ] Document test structure
- [ ] Add examples of good tests

### API Documentation
- [ ] Document VQESolver new API
- [ ] Add migration guide from old to new API
- [ ] Update all examples in docs
- [ ] Add API reference

---

## Risk Mitigation

### Potential Issues

**Issue 1: Tests fail after fix**
- *Mitigation:* Run tests individually first
- *Backup:* Keep original test file as `.bak`

**Issue 2: API changes break examples**
- *Mitigation:* Update examples immediately after test fixes
- *Verification:* Run example validation script

**Issue 3: Coverage tools slow down tests**
- *Mitigation:* Run coverage only on CI, not locally
- *Solution:* Use pytest-xdist for parallel testing

---

## Next Phase Preview

After Phase 1 completion, we immediately move to:

**Phase 2: Scientific Validation**
- Set up PySCF comparison framework
- Create reference energy database
- Implement automated validation tests
- Document accuracy metrics

---

## Ready to Start?

**Commands to begin:**

```bash
# 1. Create working branch
git checkout -b fix/phase1-failing-tests

# 2. Run failing tests to confirm issues
pytest tests/unit/test_vqe.py::TestVQESolver -v

# 3. Start fixing (Day 1)
# Open tests/unit/test_vqe.py and begin updates

# 4. Verify fixes
pytest tests/unit/test_vqe.py::TestVQESolver::test_vqe_solver_creation -v

# 5. Commit progress
git add tests/unit/test_vqe.py
git commit -m "Fix: Update VQE tests for new API"
```

**Let's begin Phase 1! Should I start fixing the VQE tests now?**
