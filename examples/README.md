# Kanad Examples

Example scripts demonstrating Kanad's capabilities.

## Directory Structure

```
examples/
├── basic/                    # Basic usage examples (coming soon)
├── advanced/                 # Advanced features and solvers
└── qiskit_integration/       # Qiskit backend examples
```

## Advanced Solvers

### [advanced_solvers_demo.py](advanced/advanced_solvers_demo.py)

Comprehensive demonstration of governance-based advanced solvers:

- **Excited States Solver**: TD-VQE for electronic excitations
- **Vibrational Structure**: Harmonic frequencies and normal modes
- **Protein Folding**: Covalent governance for peptide bonds
- **Alloy Formation**: Metallic governance for metallurgy
- **Active Space Selection**: Qubit reduction techniques

```bash
python examples/advanced/advanced_solvers_demo.py
```

## Qiskit Integration

### [01_h2_vqe_comparison.py](qiskit_integration/01_h2_vqe_comparison.py)
Compare Kanad's classical solver with Qiskit VQE for H₂ molecule.

### [02_hardware_execution_template.py](qiskit_integration/02_hardware_execution_template.py)
Template for running on IBM Quantum hardware.

### [03_cuquantum_gpu_acceleration.py](qiskit_integration/03_cuquantum_gpu_acceleration.py)
GPU-accelerated simulations using NVIDIA cuQuantum.

```bash
python examples/qiskit_integration/01_h2_vqe_comparison.py
```

## Requirements

```bash
pip install -r requirements.txt
```

For GPU support:
```bash
pip install cuquantum-python
```

## More Examples

Check the [tests/validation/](../tests/validation/) directory for scientific validation examples demonstrating:
- Various molecular systems
- Metallic bonding and band structures
- Temperature-dependent properties
- Complete VQE workflows
