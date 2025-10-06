# Kanad Cloud Testing Examples

Real-world quantum chemistry simulations running on cloud quantum backends (IBM Quantum and BlueQubit).

## Overview

These scripts demonstrate how to use Kanad to run quantum chemistry calculations on real cloud hardware and simulators. They showcase:

- **Complex molecular input** via SMILES strings
- **Custom Hamiltonian configurations** for different chemistry problems
- **Multiple mapper strategies** (Jordan-Wigner, Bravyi-Kitaev)
- **Different ansatz types** (UCC, hardware-efficient)
- **Real cloud execution** (IBM Quantum hardware, BlueQubit GPU simulators)

## Available Tests

### 1. Drug Molecule VQE (`cloud_vqe_drug_molecules.py`)

Tests VQE energy calculations on pharmaceutical molecules:
- Aspirin fragments
- Caffeine fragments
- Ibuprofen fragments
- Paracetamol

**Features:**
- SMILES-based molecular input
- Automatic 3D geometry generation
- UCC and hardware-efficient ansatze
- Cloud backend execution

**Usage:**
```bash
# BlueQubit GPU backend (fast, free)
python examples/cloud_vqe_drug_molecules.py --backend bluequbit

# IBM Quantum backend (real hardware)
python examples/cloud_vqe_drug_molecules.py --backend ibm
```

### 2. Materials Science (`cloud_materials_science.py`)

Quantum simulations for energy materials:
- **LiH** - Battery electrode materials (Li-ion batteries)
- **FeH** - Transition metal catalysts (Fischer-Tropsch synthesis)
- **H2** - Hydrogen storage materials (fuel cells)

**Features:**
- Multiple basis sets comparison (STO-3G, 6-31G)
- Different mapper strategies
- Dissociation curve calculations
- Custom Hamiltonian configurations

**Usage:**
```bash
# Test battery materials on BlueQubit
python examples/cloud_materials_science.py --backend bluequbit

# Test on IBM Quantum
python examples/cloud_materials_science.py --backend ibm
```

### 3. Custom Workflow (`cloud_custom_workflow.py`)

Advanced comparative study across multiple configurations:
- **Multiple basis sets** (STO-3G, 6-31G)
- **Multiple mappers** (Jordan-Wigner, Bravyi-Kitaev)
- **Multiple ansatze** (UCC, HW-efficient shallow/deep)
- **Multiple backends** (all available)

**Features:**
- Comprehensive configuration matrix
- Performance comparison across backends
- JSON result export
- Automated analysis and reporting

**Usage:**
```bash
python examples/cloud_custom_workflow.py
```

## Setup

### 1. Install Dependencies

```bash
pip install -e .
pip install rdkit  # For SMILES parsing
```

### 2. Get API Tokens

#### BlueQubit (Free GPU simulators)

1. Sign up at https://app.bluequbit.io
2. Get your API token
3. Set environment variable:
   ```bash
   export BLUE_TOKEN=your_token_here
   ```

#### IBM Quantum

1. Create account at https://quantum.ibm.com
2. Get your API token from the dashboard
3. Set environment variables:
   ```bash
   export IBM_API=your_token_here
   # Optional for IBM Cloud:
   export IBM_CRN=your_cloud_resource_name
   ```

### 3. Verify Setup

```bash
# Check environment variables
echo $BLUE_TOKEN
echo $IBM_API

# Test import
python -c "from kanad.backends.bluequbit import BlueQubitBackend; print('✓ BlueQubit ready')"
python -c "from kanad.backends.ibm import IBMBackend; print('✓ IBM ready')"
```

## Cloud Backend Comparison

### BlueQubit

**Pros:**
- ✓ Free GPU simulators (36 qubits)
- ✓ Very fast (GPU-accelerated)
- ✓ Statevector simulation (exact results)
- ✓ No queue times
- ✓ MPS tensor networks (40+ qubits)

**Cons:**
- ✗ Simulators only (no real quantum hardware)

**Best for:**
- Rapid prototyping
- Algorithm development
- Large-scale classical simulation

### IBM Quantum

**Pros:**
- ✓ Real quantum hardware (127+ qubits)
- ✓ Multiple backend options
- ✓ Error mitigation built-in
- ✓ Qiskit Runtime primitives

**Cons:**
- ✗ Queue times (especially for premium hardware)
- ✗ Limited free tier access
- ✗ Noise in real hardware results

**Best for:**
- Real quantum hardware testing
- NISQ algorithm validation
- Publication-quality results

## Example Output

```
================================================================================
Running VQE on: Aspirin Fragment (Phenol)
SMILES: c1ccc(cc1)O
================================================================================
Converting SMILES to molecule...
  Atoms: 7
  Electrons: 40
Building Hamiltonian with sto-3g basis...
  Orbitals: 15
  Qubits needed: 30
Building ucc ansatz...
  Parameters: 45
Building quantum circuit...
  Circuit depth: 127
  Circuit gates: 342

Submitting to bluequbit backend...
Job completed successfully!

================================================================================
RESULTS for Aspirin Fragment (Phenol):
  HF Energy:          -305.284567 Hartree
  VQE Energy:         -305.298234 Hartree
  Correlation:        -0.013667 Hartree
================================================================================
```

## Advanced Usage

### Custom Molecule from SMILES

```python
from kanad.io import SMILESConverter
from kanad.backends.bluequbit import BlueQubitBackend

# Convert SMILES to molecule
converter = SMILESConverter()
molecule = converter.smiles_to_molecule('CCO', name='ethanol')

# Build Hamiltonian
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
hamiltonian = MolecularHamiltonian(molecule, basis_name='sto-3g')

# Run on cloud
backend = BlueQubitBackend(device='gpu')
# ... build circuit and run
```

### Custom Ansatz Configuration

```python
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz

# UCC ansatz (chemically-inspired)
ucc = UCCAnsatz(
    n_qubits=10,
    n_electrons=6,
    include_singles=True,
    include_doubles=True
)

# Hardware-efficient (shallow for NISQ)
hw_eff = HardwareEfficientAnsatz(
    n_qubits=10,
    n_electrons=6,
    n_layers=2,  # Shallow for noisy hardware
    entanglement='linear'
)
```

### Multiple Mappers Comparison

```python
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.mappers.bravyi_kitaev_mapper import BravyiKitaevMapper

# Jordan-Wigner (simple, local)
jw = JordanWignerMapper()
qubit_ham_jw = jw.map_hamiltonian(hamiltonian)

# Bravyi-Kitaev (efficient for some systems)
bk = BravyiKitaevMapper()
qubit_ham_bk = bk.map_hamiltonian(hamiltonian)

# Compare Pauli term counts
print(f"JW terms: {len(qubit_ham_jw.terms)}")
print(f"BK terms: {len(qubit_ham_bk.terms)}")
```

## Performance Tips

### 1. Start Small
- Test with small molecules (H2, H2O) first
- Use STO-3G basis initially
- Increase complexity gradually

### 2. Optimize Circuits
- Use shallow ansatze for NISQ hardware
- Apply circuit optimization passes
- Consider active space reduction for large molecules

### 3. Backend Selection
- **Development/Testing**: Use BlueQubit GPU (fast, free)
- **Algorithm Validation**: Use IBM simulators
- **Final Results**: Use IBM quantum hardware

### 4. Resource Management
- Monitor qubit counts (stay under backend limits)
- Track circuit depth (shallow is better for NISQ)
- Use batch mode for IBM to minimize queue time

## Troubleshooting

### "API token not found"
```bash
export BLUE_TOKEN=your_token
export IBM_API=your_token
```

### "Too many qubits"
- Reduce basis set (6-31G → STO-3G)
- Use active space selection
- Try frozen core approximation

### "Circuit too deep"
- Use shallower ansatz (fewer layers)
- Reduce UCC excitation orders
- Use hardware-efficient ansatz

### "IBM queue too long"
- Use simulator backends during development
- Switch to BlueQubit for rapid iteration
- Use batch mode to minimize jobs

## Citation

If you use these examples in your research:

```bibtex
@software{kanad2025,
  title = {Kanad: Governance-Driven Quantum Chemistry Framework},
  author = {Kanad Framework Team},
  year = {2025},
  url = {https://github.com/yourusername/kanad}
}
```

## Support

- **Issues**: https://github.com/yourusername/kanad/issues
- **Docs**: See main README.md
- **Examples**: This directory

## License

MIT License - see LICENSE file
