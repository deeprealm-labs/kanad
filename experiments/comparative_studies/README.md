# Kanad Framework - Comparative Studies Campaign

## Overview

This campaign explores the Kanad framework through systematic comparative experiments focusing on **bonds**, **solvers**, and **backends** modules.

**Objective**: Understand how different configurations affect computational results, performance, and accuracy.

## Campaign Structure

### Experiment 1: Bond Type Comparison
**File**: `bond_types/01_bond_type_comparison.py`

**Purpose**: Compare how different bond types are represented and computed

**Test Cases**:
- H₂ (pure covalent bond)
- LiH (polar covalent/ionic bond)
- NaCl (ionic bond)
- Auto-detection vs manual bond type specification

**Metrics Collected**:
- Energy accuracy (Hartree)
- Qubit count
- Circuit depth and gate count
- Computation time
- Representation type used
- Charge transfer and bond order

**Expected Insights**:
- How bond type affects qubit requirements
- Accuracy of auto-detection mechanism
- Performance differences between ionic/covalent/metallic representations

---

### Experiment 2: Ansatz Benchmarking
**File**: `ansatz_benchmark/02_ansatz_benchmark.py`

**Purpose**: Compare different variational ansatz types on the same molecule

**Test Cases**:
- UCC (Unitary Coupled Cluster)
- Hardware Efficient Ansatz
- Governance-Aware Ansatz (physics-driven)

**Test Molecule**: H₂ at 0.74 Å (well-studied reference)

**Metrics Collected**:
- Energy accuracy vs reference
- Circuit depth and complexity
- Number of parameters
- Convergence speed (iterations)
- Computation time
- Gate count

**Expected Insights**:
- Accuracy vs circuit complexity trade-offs
- How governance-aware ansätze perform vs standard approaches
- Which ansatz converges fastest

---

### Experiment 3: Mapper Efficiency Tests
**File**: `mapper_efficiency/03_mapper_efficiency.py`

**Purpose**: Compare different fermion-to-qubit mapping schemes

**Mappers Tested**:
- Jordan-Wigner (sequential, local operators)
- Bravyi-Kitaev (tree-based, logarithmic weight)
- Hybrid Orbital (MO pair encoding, covalent-specific)

**Test Molecules**:
- H₂ (covalent bond)
- LiH (ionic/polar bond)

**Metrics Collected**:
- Qubit count
- Operator weight/complexity
- Energy accuracy
- Circuit depth and gates
- Computation time

**Expected Insights**:
- Mapper suitability for different bond types
- Performance characteristics (qubit count vs gate complexity)
- Accuracy differences between mapping schemes

---

## How to Run

### Run Individual Experiments

```bash
# Experiment 1: Bond Types
cd experiments/comparative_studies
python bond_types/01_bond_type_comparison.py

# Experiment 2: Ansatz Benchmark
python ansatz_benchmark/02_ansatz_benchmark.py

# Experiment 3: Mapper Efficiency
python mapper_efficiency/03_mapper_efficiency.py
```

### Run Full Campaign

```bash
cd experiments/comparative_studies
python run_all_experiments.py
```

This will execute all experiments sequentially and provide a comprehensive summary.

---

## Results

Results are saved as JSON files for each experiment:

1. **Bond Type Comparison**: `bond_types/results_bond_comparison.json`
2. **Ansatz Benchmark**: `ansatz_benchmark/results_ansatz_benchmark.json`
3. **Mapper Efficiency**: `mapper_efficiency/results_mapper_efficiency.json`

### Result Structure

Each result file contains an array of experiment runs with:

```json
{
  "molecule": "H2",
  "bond_type": "covalent",
  "representation_type": "LCAORepresentation",
  "num_qubits": 4,
  "energy_hartree": -1.174476,
  "computation_time": 2.34,
  "circuit_depth": 42,
  "gate_count": 156,
  "converged": true,
  "analysis": {
    "charge_transfer": 0.02,
    "bond_order": 0.98,
    "correlation_energy": -0.04
  },
  "status": "success"
}
```

---

## Framework Components Tested

### Bonds Module
- `BondFactory.create_bond()` - Automatic bond creation
- `BondFactory.quick_bond_info()` - Bond type prediction
- `IonicBond`, `CovalentBond` - Different bond implementations

### Solvers Module
- `VQESolver` - Variational Quantum Eigensolver
- Ansatz types: `ucc`, `hardware_efficient`, `governance`
- Mapper types: `jordan_wigner`, `bravyi_kitaev`, `hybrid_orbital`

### Backends Module
- `statevector` - Exact state vector simulation (used for accuracy)
- Can extend to `qasm`, `ibm`, `bluequbit` for noise/hardware studies

---

## Analysis Guidelines

### Energy Accuracy
- Compare computed energies to reference values
- Error in mHa (milli-Hartree): 1 mHa = 0.001 Ha
- Chemical accuracy target: ~1 mHa

### Performance Metrics
- **Qubit count**: Lower is better (hardware requirements)
- **Circuit depth**: Lower is better (reduces decoherence)
- **Gate count**: Lower is better (execution time)
- **Convergence speed**: Fewer iterations is better

### Trade-offs to Explore
1. Accuracy vs Circuit Complexity
2. Computation Time vs Precision
3. Representation Choice vs Qubit Count
4. Ansatz Expressivity vs Trainability

---

## Next Steps

After completing this campaign, potential follow-up studies:

1. **Backend Comparison**: Test QASM simulator vs real hardware (IBM/BlueQubit)
2. **Basis Set Studies**: Compare STO-3G, 6-31G, cc-pVDZ
3. **Molecule Complexity**: Scale to H₂O, NH₃, larger molecules
4. **Custom Governance**: Design new governance protocols
5. **Error Mitigation**: Test readout error mitigation, ZNE

---

## Requirements

- Python 3.8+
- Kanad framework installed
- Dependencies: numpy, scipy, qiskit, pyscf
- Recommended: 8GB RAM minimum

## Estimated Runtime

- Individual experiments: 5-15 minutes each
- Full campaign: 15-30 minutes (depends on system)

---

## Notes

- All experiments use `statevector` backend for exact, deterministic results
- Results are reproducible (no stochastic sampling)
- Experiments are designed to complete on modest hardware
- For production use, consider upgrading to cloud backends for larger molecules

---

## Citation

If using these experiments for research, please cite the Kanad framework appropriately.

---

**Created**: 2025-10-21
**Framework Version**: Kanad (current)
**Campaign**: Option A - Comparative Studies
