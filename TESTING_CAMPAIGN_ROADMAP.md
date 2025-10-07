# Kanad Testing Campaign Roadmap
## Path to Commercial Release

**Status:** Phase 1 - In Progress
**Goal:** Production-grade quantum chemistry framework with scientific accuracy validation
**Timeline:** 4-6 weeks for complete validation campaign

---

## Overview

Before commercial release with GUI, we must ensure:
1. ✅ **100% unit test coverage** - All components tested
2. ✅ **Scientific accuracy** - Match or exceed reference implementations
3. ✅ **Benchmark validation** - Competitive with PySCF, Qiskit Nature, PennyLane
4. ✅ **Real-world validation** - Examples work correctly on real hardware
5. ✅ **Performance metrics** - Documented speed and accuracy

---

## Phase 1: Fix Failing Tests & Expand Coverage
**Duration:** 1 week
**Current Status:** 321/344 tests passing (93%)
**Target:** 100% tests passing, 95%+ code coverage

### 1.1 Fix 22 Failing Tests

**VQE API Tests (12 tests)**
- [ ] Update `test_vqe.py::test_vqe_solver_creation` - New API with `bond` parameter
- [ ] Update `test_vqe.py::test_vqe_energy_computation` - New API
- [ ] Update `test_vqe.py::test_vqe_optimization` - New API
- [ ] Update `test_vqe.py::test_vqe_with_ucc_ansatz` - New API
- [ ] Update `test_vqe.py::test_vqe_callback` - New API
- [ ] Update `test_vqe.py::test_vqe_energy_history` - New API
- [ ] Update `test_vqe.py::test_vqe_variance_computation` - New API
- [ ] Update `test_vqe.py::test_vqe_with_different_optimizers` - New API
- [ ] Update `test_vqe.py::test_vqe_parameter_initialization` - New API
- [ ] Update `test_vqe.py::test_vqe_reproducibility` - New API
- [ ] Update `test_vqe.py::test_vqe_converges_below_initial` - New API
- [ ] Update `test_qiskit_integration.py::test_vqe_classical_vs_qiskit` - New API

**Qiskit Integration Tests (7 tests)**
- [ ] Fix `test_qiskit_integration.py::test_aer_simulator_init` - Import paths
- [ ] Fix `test_qiskit_integration.py::test_aer_statevector_init` - Import paths
- [ ] Fix `test_qiskit_integration.py::test_get_estimator` - Import paths
- [ ] Fix `test_qiskit_integration.py::test_backend_info` - Import paths
- [ ] Fix `test_qiskit_integration.py::test_vqe_classical_backend` - API + imports
- [ ] Fix `test_qiskit_integration.py::test_vqe_aer_backend` - API + imports
- [ ] Fix `test_qiskit_integration.py::test_vqe_classical_vs_qiskit` - API + imports

**Import Tests (3 tests)**
- [ ] Fix `test_imports.py::test_main_import` - Update for removed modules
- [ ] Fix `test_imports.py::test_solvers_import` - Update import paths
- [ ] Fix `test_imports.py::test_backends_import` - Update for new structure
- [ ] Fix `test_imports.py::test_ibm_solvers_import` - Update paths

### 1.2 Expand Test Coverage

**Missing Coverage Areas:**

**Core Components**
- [ ] Add tests for `Molecule` edge cases (empty, single atom, large molecules)
- [ ] Add tests for all basis sets (STO-3G, 6-31G, cc-pVDZ, def2-SVP)
- [ ] Add tests for spin states (singlet, doublet, triplet, quartet)
- [ ] Add tests for charged molecules (+1, -1, +2, -2)
- [ ] Add tests for excited states

**Hamiltonians**
- [ ] Add comprehensive ionic Hamiltonian tests
- [ ] Add comprehensive metallic Hamiltonian tests
- [ ] Add periodic Hamiltonian tests (if applicable)
- [ ] Add Hamiltonian symmetry tests
- [ ] Add Hamiltonian hermiticity tests

**Mappers**
- [ ] Add Parity mapper comprehensive tests
- [ ] Add Hybrid Orbital mapper edge cases
- [ ] Add mapper equivalence tests (same energy different mappings)
- [ ] Add qubit reduction tests

**Ansatze**
- [ ] Add governance-aware ansatz comprehensive tests
- [ ] Add adaptive ansatz tests
- [ ] Add parameter gradient tests
- [ ] Add circuit depth optimization tests

**Solvers**
- [ ] Add SQD solver comprehensive tests
- [ ] Add QPE solver tests
- [ ] Add Excited states solver tests
- [ ] Add solver convergence tests
- [ ] Add solver noise resilience tests

**I/O**
- [ ] Add SMILES parser edge cases (radicals, charged species, aromatics)
- [ ] Add XYZ file format validation tests
- [ ] Add crystal builder comprehensive tests
- [ ] Add PDB file I/O tests (if implemented)

**Backends**
- [ ] Add IBM backend error handling tests
- [ ] Add BlueQubit backend error handling tests
- [ ] Add backend timeout tests
- [ ] Add job cancellation tests
- [ ] Add batch job tests

**Governance**
- [ ] Add protocol violation tests
- [ ] Add constraint enforcement tests
- [ ] Add operator validation tests
- [ ] Add circuit topology validation tests

### 1.3 Test Organization

Create systematic test structure:
```
tests/
├── unit/                      # Unit tests (existing)
│   ├── test_*.py
│   └── io/
├── integration/               # Integration tests
│   ├── test_bond_to_energy.py
│   ├── test_full_vqe_workflow.py
│   └── test_cloud_backends.py
├── validation/               # Scientific validation (NEW)
│   ├── test_h2_dissociation.py
│   ├── test_reference_energies.py
│   ├── test_molecular_properties.py
│   └── benchmarks/
└── performance/              # Performance tests (NEW)
    ├── test_speed_benchmarks.py
    └── test_memory_usage.py
```

---

## Phase 2: Scientific Validation Campaign
**Duration:** 2 weeks
**Target:** Match reference implementations to <1 mHa accuracy

### 2.1 Reference Energy Validation

Compare against established quantum chemistry codes:

**Test Set 1: Small Molecules (Exact Reference)**
| Molecule | Geometry | Basis | Reference (Ha) | Source | Status |
|----------|----------|-------|----------------|--------|--------|
| H2 | 0.74 Å | STO-3G | -1.117 | PySCF | ⏳ |
| H2 | 0.74 Å | 6-31G | -1.131 | PySCF | ⏳ |
| H2+ | 1.06 Å | STO-3G | -0.587 | Exact | ⏳ |
| HeH+ | 0.772 Å | STO-3G | -2.908 | PySCF | ⏳ |
| LiH | 1.596 Å | STO-3G | -7.863 | PySCF | ⏳ |
| LiH | 1.596 Å | 6-31G | -7.987 | PySCF | ⏳ |
| H2O | Optimized | STO-3G | -74.964 | PySCF | ⏳ |
| H2O | Optimized | 6-31G | -76.010 | PySCF | ⏳ |
| NH3 | Optimized | STO-3G | -55.454 | PySCF | ⏳ |
| CH4 | Optimized | STO-3G | -39.727 | PySCF | ⏳ |

**Test Set 2: Diatomics (FCI Reference)**
| Molecule | Geometry | Basis | HF Energy | FCI Energy | Source | Status |
|----------|----------|-------|-----------|------------|--------|--------|
| N2 | 1.098 Å | STO-3G | -107.496 | -107.654 | Literature | ⏳ |
| CO | 1.128 Å | STO-3G | -111.226 | -111.358 | Literature | ⏳ |
| F2 | 1.412 Å | STO-3G | -197.068 | -197.213 | Literature | ⏳ |
| BH | 1.232 Å | STO-3G | -24.963 | -25.012 | Literature | ⏳ |
| OH | 0.970 Å | STO-3G | -74.805 | -74.891 | Literature | ⏳ |

**Test Set 3: Polyatomics**
| Molecule | Description | Basis | HF Energy | MP2 Energy | Status |
|----------|-------------|-------|-----------|------------|--------|
| H2O2 | Hydrogen peroxide | STO-3G | -149.582 | -149.721 | ⏳ |
| N2H4 | Hydrazine | STO-3G | -110.462 | -110.683 | ⏳ |
| C2H6 | Ethane | STO-3G | -78.306 | -78.489 | ⏳ |
| C2H4 | Ethylene | STO-3G | -77.073 | -77.234 | ⏳ |
| C2H2 | Acetylene | STO-3G | -75.856 | -75.992 | ⏳ |

### 2.2 Molecular Properties Validation

**Dipole Moments**
| Molecule | Basis | Reference (Debye) | Kanad | Error | Status |
|----------|-------|-------------------|-------|-------|--------|
| H2O | STO-3G | 2.19 | ? | ? | ⏳ |
| H2O | 6-31G | 2.36 | ? | ? | ⏳ |
| NH3 | STO-3G | 1.89 | ? | ? | ⏳ |
| HF | STO-3G | 2.11 | ? | ? | ⏳ |
| CO | STO-3G | 0.25 | ? | ? | ⏳ |

**Bond Dissociation Curves**
- [ ] H2 dissociation curve (0.5 - 3.0 Å)
- [ ] LiH dissociation curve
- [ ] N2 dissociation curve
- [ ] Verify no unphysical behavior

**Vibrational Frequencies**
| Molecule | Mode | Reference (cm⁻¹) | Kanad | Error | Status |
|----------|------|------------------|-------|-------|--------|
| H2O | ν1 (symmetric) | 3832 | ? | ? | ⏳ |
| H2O | ν2 (bend) | 1649 | ? | ? | ⏳ |
| H2O | ν3 (antisym) | 3943 | ? | ? | ⏳ |
| NH3 | ν1 | 3506 | ? | ? | ⏳ |

### 2.3 Quantum Algorithm Validation

**VQE Convergence**
- [ ] Test convergence to HF energy (without correlation)
- [ ] Test convergence below HF (with correlation)
- [ ] Compare UCC vs Hardware-Efficient ansatze
- [ ] Verify parameter gradients are correct

**QPE Accuracy**
- [ ] Test phase estimation accuracy
- [ ] Compare with exact eigenvalues
- [ ] Test on small molecules

**SQD Accuracy**
- [ ] Test subspace diagonalization
- [ ] Compare with full diagonalization
- [ ] Verify orthogonality

---

## Phase 3: Benchmark Against Existing Frameworks
**Duration:** 2 weeks
**Target:** Competitive or better performance

### 3.1 Frameworks to Compare Against

**Primary Comparisons:**
1. **PySCF** - Classical quantum chemistry (HF, MP2, CCSD)
2. **Qiskit Nature** - Quantum algorithms for chemistry
3. **PennyLane** - Differentiable quantum chemistry
4. **OpenFermion** - Fermionic-qubit mappings
5. **Psi4** - Ab initio quantum chemistry

### 3.2 Benchmark Tests

**Test 1: Energy Accuracy**
```
Molecules: H2, LiH, H2O, NH3, CH4
Basis: STO-3G, 6-31G
Methods: HF, VQE-UCC, VQE-HEA
Metric: Energy error vs reference
```

**Test 2: Speed Performance**
```
Molecule sizes: 2-20 atoms
Metric: Time to HF convergence, Time to VQE convergence
Compare: Kanad vs PySCF vs Qiskit Nature
```

**Test 3: Qubit Efficiency**
```
Molecules: H2, LiH, H2O, BeH2
Metric: Number of qubits required
Compare different mappers: JW, BK, Parity
```

**Test 4: Circuit Depth**
```
Molecules: H2, LiH, H2O
Ansatze: UCC, HEA, Governance-aware
Metric: Circuit depth, gate count
```

**Test 5: Noise Resilience**
```
Test on: IBM simulators with noise models
Molecules: H2, LiH
Metric: Energy accuracy with noise
Compare: Error mitigation strategies
```

### 3.3 Benchmark Targets

| Metric | Target | Rationale |
|--------|--------|-----------|
| Energy Accuracy | <1 mHa error vs PySCF | Chemical accuracy |
| HF Speed | Within 2x of PySCF | Acceptable overhead |
| VQE Convergence | <100 iterations | Efficient optimization |
| Qubit Count | Match or better than Qiskit | Efficient mappings |
| Circuit Depth | 20% reduction vs standard | Governance optimization |

---

## Phase 4: Real-Life Examples Validation
**Duration:** 1 week
**Target:** All examples work correctly on real hardware

### 4.1 Cloud Backend Testing

**IBM Quantum Hardware**
- [ ] Test on `ibm_brisbane` (127 qubits)
- [ ] Test on `ibm_kyoto` (127 qubits)
- [ ] Test error mitigation (resilience levels 0, 1, 2)
- [ ] Test batch mode submission
- [ ] Verify job status tracking
- [ ] Test job cancellation

**BlueQubit Simulators**
- [ ] Test GPU simulator (36 qubits)
- [ ] Test CPU simulator (34 qubits)
- [ ] Test MPS simulator (40+ qubits)
- [ ] Verify statevector accuracy
- [ ] Test shot-based simulation

### 4.2 Example Validation

**Drug Molecules** (`cloud_vqe_drug_molecules.py`)
- [ ] Aspirin fragment - verify energy is reasonable
- [ ] Caffeine fragment - verify energy
- [ ] Ibuprofen fragment - verify energy
- [ ] Paracetamol - verify energy
- [ ] Compare with literature values

**Materials Science** (`cloud_materials_science.py`)
- [ ] LiH battery material - verify dissociation curve
- [ ] FeH catalyst - verify spin states
- [ ] H2 storage - verify bond lengths

**SMILES Library** (`cloud_smiles_molecules.py`)
- [ ] Test all 30+ molecules
- [ ] Verify no failures
- [ ] Check energy trends make sense
- [ ] Validate against PubChem data

### 4.3 Documentation

- [ ] Update all examples with correct outputs
- [ ] Add expected results to README
- [ ] Create troubleshooting guide
- [ ] Add performance tips

---

## Phase 5: GUI Development Preparation
**Duration:** Planning phase

### 5.1 Requirements Gathering

**User Personas:**
1. **Research Chemist** - Needs molecule builder, property calculator
2. **Quantum Developer** - Needs circuit customization, mapper selection
3. **Student** - Needs educational interface, visualizations
4. **Industrial User** - Needs batch processing, API access

**Core Features:**
- [ ] Molecule builder (draw structures, SMILES input)
- [ ] Calculation setup (basis, method, backend selection)
- [ ] Job submission and monitoring
- [ ] Results visualization (energies, orbitals, properties)
- [ ] Batch processing interface
- [ ] Export/import functionality

### 5.2 Technology Stack

**Frontend Options:**
- Electron + React (desktop app)
- Web interface (Flask/FastAPI + React)
- Jupyter widgets (for researchers)

**Visualization:**
- py3Dmol for molecular structures
- Plotly for graphs and charts
- Matplotlib for publication-quality plots

### 5.3 API Design

Clean API for GUI to interact with:
```python
# High-level API
from kanad.api import KanadAPI

api = KanadAPI()
job = api.calculate(
    molecule="CCO",  # SMILES or file
    method="vqe",
    basis="sto-3g",
    backend="bluequbit_gpu"
)
results = job.wait()
```

---

## Testing Campaign Execution Plan

### Week 1: Fix Failing Tests
- Days 1-2: Fix VQE API tests
- Days 3-4: Fix Qiskit integration tests
- Day 5: Fix import tests
- Days 6-7: Expand test coverage (core components)

### Week 2: Expand Coverage
- Days 1-2: Hamiltonian tests
- Days 3-4: Mapper and ansatz tests
- Days 5-6: Solver tests
- Day 7: I/O and backend tests

### Week 3: Scientific Validation
- Days 1-3: Reference energy validation (small molecules)
- Days 4-5: Diatomic and polyatomic validation
- Days 6-7: Molecular properties validation

### Week 4: Continued Validation
- Days 1-2: Bond dissociation curves
- Days 3-4: Quantum algorithm validation
- Days 5-7: Documentation and fixes

### Week 5: Benchmarking
- Days 1-2: Setup benchmark framework
- Days 3-4: Run energy and speed benchmarks
- Days 5-7: Analyze results and optimize

### Week 6: Real-Life Examples
- Days 1-2: Cloud backend testing
- Days 3-4: Example validation
- Days 5-7: Documentation and final testing

---

## Success Criteria

### Must Have (Required for Commercial Release)
- [ ] 100% unit tests passing
- [ ] >95% code coverage
- [ ] All reference energies within 1 mHa
- [ ] All molecular properties within 5% of reference
- [ ] Successful execution on IBM and BlueQubit
- [ ] All examples work correctly
- [ ] Complete documentation

### Should Have (Competitive Edge)
- [ ] Speed within 2x of PySCF
- [ ] Qubit efficiency equal or better than Qiskit
- [ ] Circuit depth 20% better than standard
- [ ] Noise resilience better than baseline

### Nice to Have (Future Enhancements)
- [ ] GPU acceleration
- [ ] Active space automation
- [ ] Error mitigation strategies
- [ ] Distributed computing

---

## Deliverables

1. **Test Suite** - Comprehensive, automated, CI/CD ready
2. **Validation Report** - Scientific accuracy documentation
3. **Benchmark Report** - Performance comparison document
4. **Example Suite** - Working real-world examples
5. **API Documentation** - Clean, well-documented API
6. **User Guide** - Complete usage documentation
7. **Developer Guide** - Contributing and extending guide

---

## Next Steps

**Immediate Actions:**
1. Create test scripts for fixing 22 failing tests
2. Set up validation test framework
3. Establish PySCF comparison pipeline
4. Create benchmark automation scripts

**This Week:**
- Fix all 22 failing tests
- Expand core component test coverage
- Set up CI/CD pipeline

Would you like me to start with Phase 1 immediately?
