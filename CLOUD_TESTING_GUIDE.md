# Kanad Cloud Testing Guide

Complete guide for running real quantum chemistry calculations on cloud quantum backends.

## Quick Start

### 1. Set up API tokens

```bash
# For BlueQubit (free GPU simulators)
export BLUE_TOKEN=your_bluequbit_token

# For IBM Quantum (real hardware)
export IBM_API=your_ibm_token
```

### 2. Run tests

```bash
# Run all tests on BlueQubit
./examples/run_cloud_tests.sh --backend bluequbit

# Run specific test on IBM
./examples/run_cloud_tests.sh --backend ibm --test drug_molecules

# Run drug molecules test
python examples/cloud_vqe_drug_molecules.py --backend bluequbit

# Run materials science test
python examples/cloud_materials_science.py --backend bluequbit

# Run custom workflow
python examples/cloud_custom_workflow.py

# Run SMILES library test
python examples/cloud_smiles_molecules.py --backend bluequbit
```

## Available Test Scripts

### 1. Drug Molecules VQE (`cloud_vqe_drug_molecules.py`)

Tests pharmaceutical molecules:
- Aspirin fragments (phenol core)
- Caffeine fragments (xanthine core)
- Ibuprofen fragments
- Paracetamol
- Benzene (aromatic reference)

**What it demonstrates:**
- SMILES-based molecular input
- Automatic 3D geometry generation
- UCC ansatz for accurate chemistry
- Cloud backend execution
- Result analysis and reporting

**Run:**
```bash
python examples/cloud_vqe_drug_molecules.py --backend bluequbit
```

### 2. Materials Science (`cloud_materials_science.py`)

Tests energy materials:
- **LiH** - Battery electrodes (Li-ion batteries)
- **FeH** - Catalysts (transition metals)
- **H2** - Hydrogen storage (dissociation curves)

**What it demonstrates:**
- Multiple basis sets (STO-3G, 6-31G)
- Different mappers (Jordan-Wigner, Bravyi-Kitaev)
- Multiple ansatze (UCC, hardware-efficient)
- Material-specific configurations
- Property calculations (bond dissociation)

**Run:**
```bash
python examples/cloud_materials_science.py --backend bluequbit
```

### 3. Custom Workflow (`cloud_custom_workflow.py`)

Comprehensive configuration comparison:
- Tests multiple molecules
- All available mappers
- All ansatz types
- All backends simultaneously
- Generates comparison reports

**What it demonstrates:**
- Advanced workflow patterns
- Parallel backend execution
- Configuration matrix testing
- Result aggregation
- JSON export for analysis

**Run:**
```bash
python examples/cloud_custom_workflow.py
```

### 4. SMILES Library (`cloud_smiles_molecules.py`)

High-throughput molecule testing:
- 30+ molecules organized by category
- Diatomics, small molecules, organics, aromatics, drugs, ionic
- Batch processing
- Systematic testing

**What it demonstrates:**
- Large-scale molecule processing
- Category-based organization
- Automated screening workflows
- Result persistence

**Run:**
```bash
# All categories
python examples/cloud_smiles_molecules.py --backend bluequbit

# Specific categories
python examples/cloud_smiles_molecules.py --category aromatics drug_fragments

# Custom qubit limit
python examples/cloud_smiles_molecules.py --max-qubits 15
```

## Framework Components Used

### I/O Module
```python
from kanad.io import SMILESConverter

# Convert SMILES to molecule
converter = SMILESConverter(optimize_geometry=True)
molecule = converter.smiles_to_molecule('CCO', charge=0, multiplicity=1)
```

### Hamiltonians
```python
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian

# Build Hamiltonian
hamiltonian = MolecularHamiltonian(molecule, basis_name='sto-3g')
```

### Mappers
```python
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.mappers.bravyi_kitaev_mapper import BravyiKitaevMapper

# Create mapper
mapper = JordanWignerMapper()
qubit_ham = mapper.map_hamiltonian(hamiltonian)
```

### Ansatze
```python
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz

# UCC ansatz (chemically-inspired)
ansatz = UCCAnsatz(
    n_qubits=n_qubits,
    n_electrons=n_electrons,
    include_singles=True,
    include_doubles=True
)

# Hardware-efficient (NISQ-friendly)
ansatz = HardwareEfficientAnsatz(
    n_qubits=n_qubits,
    n_electrons=n_electrons,
    n_layers=2
)
```

### Cloud Backends

#### BlueQubit
```python
from kanad.backends.bluequbit import BlueQubitBackend

backend = BlueQubitBackend(
    device='gpu',  # 'cpu', 'gpu', 'mps.gpu'
    api_token=os.getenv('BLUE_TOKEN')
)

result = backend.run_circuit(circuit, shots=1024)
```

#### IBM Quantum
```python
from kanad.backends.ibm import IBMBackend

backend = IBMBackend(
    backend_name='ibm_brisbane',  # or 'ibmq_qasm_simulator'
    api_token=os.getenv('IBM_API')
)

result = backend.run_batch(
    circuits=[circuit],
    observables=[observable],
    shots=1024,
    optimization_level=1,
    resilience_level=1
)
```

## Custom Configuration Examples

### Example 1: Custom Molecule from Scratch

```python
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian

# Create atoms
h1 = Atom('H', [0.0, 0.0, 0.0])
h2 = Atom('H', [0.0, 0.0, 0.74])

# Create molecule
molecule = Molecule([h1, h2], charge=0, spin=0)

# Build Hamiltonian
hamiltonian = MolecularHamiltonian(molecule, basis_name='6-31g')
```

### Example 2: Custom Mapper + Ansatz Combination

```python
from kanad.core.mappers.bravyi_kitaev_mapper import BravyiKitaevMapper
from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz

# Use Bravyi-Kitaev for efficient mapping
mapper = BravyiKitaevMapper()
qubit_ham = mapper.map_hamiltonian(hamiltonian)

# Shallow circuit for NISQ hardware
ansatz = HardwareEfficientAnsatz(
    n_qubits=n_qubits,
    n_electrons=n_electrons,
    n_layers=2,  # Shallow for noise resilience
    entanglement='linear'
)
```

### Example 3: SMILES Input with Custom Properties

```python
from kanad.io import SMILESConverter

converter = SMILESConverter(optimize_geometry=True)

# Create molecule from SMILES
molecule = converter.smiles_to_molecule(
    'c1ccc(cc1)O',  # Phenol
    charge=0,
    multiplicity=1,
    name='phenol'
)

# Common molecules library
common = converter.common_molecules()
h2o = converter.smiles_to_molecule(common['H2O'])
benzene = converter.smiles_to_molecule(common['benzene'])
```

### Example 4: Parallel Backend Execution

```python
from concurrent.futures import ThreadPoolExecutor

backends = {
    'bluequbit_gpu': BlueQubitBackend(device='gpu'),
    'bluequbit_cpu': BlueQubitBackend(device='cpu'),
    'ibm_sim': IBMBackend(backend_name='ibmq_qasm_simulator')
}

def run_on_backend(name, backend, circuit):
    return backend.run_circuit(circuit)

# Run in parallel
with ThreadPoolExecutor(max_workers=3) as executor:
    futures = {
        executor.submit(run_on_backend, name, backend, circuit): name
        for name, backend in backends.items()
    }

    for future in as_completed(futures):
        name = futures[future]
        result = future.result()
        print(f"{name}: {result}")
```

## Performance Optimization

### 1. Active Space Selection

For large molecules, use active space to reduce qubits:

```python
from kanad.solvers.active_space import ActiveSpaceSelector

selector = ActiveSpaceSelector(molecule, hamiltonian)
active_orbitals, active_electrons = selector.select_active_space(
    n_orbitals=6,  # Keep only 6 most important orbitals
    method='natural_orbitals'
)

# Build reduced Hamiltonian
reduced_ham = hamiltonian.get_active_space_hamiltonian(active_orbitals)
```

### 2. Circuit Optimization

```python
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Optimize1qGatesDecomposition, CXCancellation

# Optimize circuit before submission
pm = PassManager([
    Optimize1qGatesDecomposition(),
    CXCancellation()
])

optimized_circuit = pm.run(circuit)
```

### 3. Basis Set Selection

| Basis Set | Accuracy | Qubits | Recommended For |
|-----------|----------|--------|-----------------|
| STO-3G    | Low      | Few    | Testing, large molecules |
| 6-31G     | Medium   | More   | Production, medium molecules |
| cc-pVDZ   | High     | Many   | Benchmark, small molecules |

## Troubleshooting

### Issue: "Too many qubits"

**Solutions:**
1. Use STO-3G basis instead of 6-31G
2. Apply active space selection
3. Use frozen core approximation
4. Try smaller molecule fragments

### Issue: "Circuit too deep"

**Solutions:**
1. Use hardware-efficient ansatz with fewer layers
2. Reduce UCC excitation orders
3. Apply circuit optimization passes
4. Use shallower ansatz for NISQ hardware

### Issue: "API token not found"

```bash
# Check if set
echo $BLUE_TOKEN
echo $IBM_API

# Set properly
export BLUE_TOKEN=your_token
export IBM_API=your_token
```

### Issue: "Job queue timeout"

**Solutions:**
1. Use simulator backends during development
2. Switch to BlueQubit for rapid iteration
3. Use IBM batch mode to reduce queue position
4. Run during off-peak hours

## Results Interpretation

### Energy Values

```python
result = {
    'hf_energy': -1.117,      # Hartree-Fock reference
    'vqe_energy': -1.136,     # VQE result
    'correlation_energy': -0.019  # Electron correlation captured
}
```

- **HF Energy**: Mean-field approximation (baseline)
- **VQE Energy**: Should be lower (more negative) than HF
- **Correlation Energy**: Difference = electron correlation effects

### Quality Metrics

- **Good VQE**: Correlation energy is negative and significant
- **Excellent VQE**: Within 1 mHa of FCI (full CI)
- **Poor VQE**: VQE energy higher than HF (optimization failed)

## Citation

```bibtex
@software{kanad2025,
  title = {Kanad: Governance-Driven Quantum Chemistry Framework},
  author = {Kanad Framework Team},
  year = {2025},
  url = {https://github.com/yourusername/kanad}
}
```

## Support

- **Documentation**: See `examples/README_CLOUD_TESTS.md`
- **Issues**: GitHub Issues
- **Examples**: All scripts in `examples/` directory

## License

MIT License
