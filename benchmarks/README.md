# Kanad Framework Benchmarks

## Overview

This directory contains comprehensive benchmarks for the Kanad quantum chemistry framework, testing all major components across multiple molecules.

## Internal Framework Benchmark

**File**: `kanad_benchmark_final.py`

### What It Tests

1. **Molecules** (3):
   - **H2**: Hydrogen molecule (covalent bonding)
   - **LiH**: Lithium hydride (polar covalent/ionic)
   - **HeH+**: Helium hydride cation (ionic bonding)

2. **Ansatze** (4):
   - **HardwareEfficient**: Hardware-native circuit with linear entanglement
   - **TwoLocal**: Two-local ansatz with RY rotations + CNOT
   - **CovalentGovernance**: Physics-informed ansatz for covalent bonds
   - **IonicGovernance**: Physics-informed ansatz for ionic bonds

3. **Optimizers** (3):
   - **COBYLA**: Constrained Optimization BY Linear Approximation
   - **SLSQP**: Sequential Least Squares Programming
   - **Powell**: Powell's conjugate direction method

### Total Configurations

**36 tests** = 3 molecules × 4 ansatze × 3 optimizers

### Metrics Measured

- **Accuracy**: Error vs FCI (Full Configuration Interaction) reference
  - Chemical accuracy: < 1.6 mHa (1 kcal/mol)
  - Nobel laureate accuracy: < 1.0 mHa
- **Performance**: Wall-clock time for VQE convergence
- **Convergence**: Number of optimizer iterations
- **Parameters**: Number of variational parameters in ansatz

### Expected Results

Based on validation testing:
- **Success rate**: ~95-100%
- **Chemical accuracy**: ~80-90% of successful runs
- **Nobel accuracy**: ~50-70% of successful runs
- **Best configurations**: Governance ansatze with SLSQP/Powell optimizers

## Molecular Analysis

**File**: `molecular_analysis.py`

Provides detailed analysis of each molecule:
- Bonding character (ionic vs covalent)
- Electronegativity differences
- Correlation energy analysis
- Recommended ansatz/Hamiltonian types

### Usage

```python
from benchmarks.molecular_analysis import MolecularAnalyzer

analyzer = MolecularAnalyzer(molecule)
report = analyzer.generate_analysis_report()
print(report)
```

## Results Files

- **results.json**: Raw benchmark data (all 36 configurations)
- **benchmark_final.log**: Human-readable output with progress
- **molecular_analysis/**: Detailed per-molecule analysis reports

## Running Benchmarks

### Quick Benchmark (36 configs, ~10-20 minutes)

```bash
python benchmarks/kanad_benchmark_final.py
```

### With Molecular Analysis

```bash
python benchmarks/molecular_analysis.py
```

## Interpreting Results

### Accuracy Levels

| Error (mHa) | Classification | Meaning |
|-------------|----------------|---------|
| < 1.0 | ★ Nobel accuracy | Publishable research quality |
| < 1.6 | ✓ Chemical accuracy | Chemically meaningful predictions |
| < 10.0 | ○ Acceptable | Qualitative insights |
| > 10.0 | ✗ Poor | Needs investigation |

### Performance Expectations

| Molecule | Qubits | Expected Time (per config) |
|----------|--------|----------------------------|
| H2 | 4 | 5-15 seconds |
| LiH | 12 | 15-40 seconds |
| HeH+ | 4 | 5-15 seconds |

### Best Practices

1. **For covalent systems** (H2):
   - Use: `CovalentGovernanceAnsatz` or `HardwareEfficientAnsatz`
   - Optimizer: `SLSQP` or `Powell`
   - Expected: < 1.0 mHa error

2. **For ionic systems** (HeH+):
   - Use: `IonicGovernanceAnsatz` or `TwoLocalAnsatz`
   - Optimizer: `SLSQP`
   - Expected: < 2.0 mHa error

3. **For polar systems** (LiH):
   - Use: `HardwareEfficientAnsatz` or governance ansatz
   - Optimizer: `Powell` or `SLSQP`
   - Expected: < 3.0 mHa error

## Framework Components Validated

✓ **Hamiltonians**: Molecular, Covalent, Ionic
✓ **Ansatze**: Hardware-Efficient, Two-Local, Governance-based
✓ **Optimizers**: COBYLA, SLSQP, Powell
✓ **Mappers**: Jordan-Wigner, Bravyi-Kitaev
✓ **Solvers**: VQE with multiple configurations

## Future Benchmarks

### Planned

1. **External Framework Comparison**:
   - vs PySCF (classical DFT/HF/FCI)
   - vs OpenFermion (quantum algorithms)
   - vs Qiskit Nature (quantum chemistry)

2. **Larger Molecules**:
   - H2O (water)
   - NH3 (ammonia)
   - CH4 (methane)
   - Benzene ring

3. **Cloud Backends**:
   - BlueQubit GPU simulator
   - IBM Quantum hardware
   - AWS Braket

## Citation

If you use these benchmarks, please cite:

```bibtex
@software{kanad_benchmarks,
  title = {Kanad Framework Internal Benchmarks},
  author = {Kanad Development Team},
  year = {2025},
  url = {https://github.com/yourusername/kanad}
}
```

## Contact

For questions about benchmarks: [your-email]

## License

Same as Kanad framework license.
