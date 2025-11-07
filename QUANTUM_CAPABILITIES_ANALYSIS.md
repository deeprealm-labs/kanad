# Kanad Framework: Quantum Execution Capabilities Analysis

## Executive Summary

The Kanad framework implements a **3-tier quantum execution architecture**:
1. **Quantum-Native**: Algorithms that execute on real quantum hardware
2. **Quantum-Enhanced**: Classical algorithms using quantum-computed properties
3. **Classical**: Pure classical simulations for comparison/fallback

---

## 1. QUANTUM-NATIVE CAPABILITIES
### Can Execute on Real Quantum Hardware (IBM Quantum, BlueQubit)

### 1.1 Solvers (Core Quantum Algorithms)

#### ✅ VQE (Variational Quantum Eigensolver)
- **File**: `/home/user/kanad/kanad/solvers/vqe_solver.py`
- **Quantum Execution Path**:
  ```python
  # Lines 121-151: IBM Backend
  if hasattr(self, '_ibm_backend') and self._ibm_backend is not None:
      result = self._ibm_backend.run_batch(
          circuits=[bound_circuit],
          observables=[pauli_hamiltonian],
          shots=self.shots
      )
  
  # Lines 922-1027: BlueQubit Backend
  elif hasattr(self, '_bluequbit_backend') and self._bluequbit_backend is not None:
      job_info = self._bluequbit_backend.run_circuit(
          circuit=bound_circuit,
          shots=use_shots,
          asynchronous=True
      )
  ```
- **Hardware Modes**:
  - IBM: Real quantum processors (127+ qubits), cloud simulators
  - BlueQubit: GPU/CPU simulators (36 qubits), MPS tensor (40+ qubits)
- **Usage**:
  ```python
  solver = VQESolver(bond, backend='ibm', backend_name='ibm_brisbane')
  solver = VQESolver(bond, backend='bluequbit', device='gpu')
  ```
- **Auto-Optimization**: Automatically switches to SPSA optimizer for cloud backends (20x efficiency gain)

#### ✅ SQD (Subspace Quantum Diagonalization)
- **File**: `/home/user/kanad/kanad/solvers/sqd_solver.py`
- **Quantum Execution**: Generates quantum subspace basis using short-depth circuits
- **Lines 492-580**: Quantum statevector generation with governance filtering
- **Hardware Execution**:
  ```python
  solver = SQDSolver(bond, backend='ibm', subspace_dim=10)
  solver = SQDSolver(bond, backend='bluequbit', subspace_dim=15)
  ```
- **Advantage**: Lower circuit depth than VQE, more noise-resistant

#### ✅ Excited States Solver (Quantum Methods)
- **File**: `/home/user/kanad/kanad/solvers/excited_states_solver.py`
- **Quantum Methods**: QPE, VQE (state-averaged), SQD
- **Usage**:
  ```python
  solver = ExcitedStatesSolver(bond, method='sqd', backend='ibm', n_states=5)
  solver = ExcitedStatesSolver(bond, method='vqe', backend='bluequbit')
  ```

#### ✅ Hi-VQE (Hierarchical VQE)
- **File**: `/home/user/kanad/kanad/utils/hivqe_solver_mixin.py`
- **Quantum Execution**: Iterative subspace expansion using quantum circuits
- **Hardware Support**: IBM, BlueQubit
- **Usage**:
  ```python
  solver = VQESolver(bond, mode='hivqe', backend='ibm', hivqe_max_iterations=10)
  ```

### 1.2 Backend Infrastructure

#### ✅ IBM Quantum Backend
- **File**: `/home/user/kanad/kanad/backends/ibm/backend.py`
- **Capabilities**:
  - Real quantum hardware (127+ qubit processors)
  - Cloud simulators
  - Batch mode (parallel jobs)
  - Session mode (reserved hardware for iterative algorithms)
  - Qiskit Runtime primitives (SamplerV2, EstimatorV2)
  - Auto-configured error mitigation (ZNE, dynamical decoupling)
- **Supported Backends**: ibm_brisbane, ibm_torino, ibm_kyoto, etc.
- **Code Example**:
  ```python
  from kanad.backends.ibm import IBMBackend
  backend = IBMBackend(backend_name='ibm_brisbane')
  result = backend.run_batch(circuits, observables, shots=1024)
  ```

#### ✅ BlueQubit Backend
- **File**: `/home/user/kanad/kanad/backends/bluequbit/backend.py`
- **Capabilities**:
  - GPU simulators (free, 36 qubits, fast)
  - CPU simulators (34 qubits)
  - MPS tensor network (40+ qubits)
  - Statevector and sampling modes
- **Devices**: 'cpu', 'gpu', 'mps.cpu', 'mps.gpu', 'pauli-path'
- **Code Example**:
  ```python
  from kanad.backends.bluequbit import BlueQubitBackend
  backend = BlueQubitBackend(device='gpu')
  result = backend.run_circuit(circuit, shots=1024)
  ```

---

## 2. QUANTUM-ENHANCED CAPABILITIES
### Classical Algorithms Using Quantum Results

### 2.1 Analysis Modules (Use Quantum Density Matrices)

All analysis modules can operate in two modes:
1. **Classical Mode**: Use Hartree-Fock density matrix
2. **Quantum-Enhanced Mode**: Use quantum density from VQE/SQD

#### 🔄 Property Calculator
- **File**: `/home/user/kanad/kanad/analysis/property_calculator.py`
- **Quantum Enhancement**:
  - Lines 81-94: Automatically detects and uses quantum density matrix
  ```python
  if hasattr(self.hamiltonian, 'get_density_matrix'):
      density_matrix = self.hamiltonian.get_density_matrix()  # Quantum if available
  ```
- **Quantum-Enhanced Properties**:
  - Dipole moment (lines 45-173)
  - Polarizability (lines 673+)
  - Quadrupole moment
- **Usage**:
  ```python
  # Run quantum solver first
  vqe_result = VQESolver(bond, backend='ibm').solve()
  # Then compute properties with quantum density
  props = PropertyCalculator(hamiltonian).compute_dipole_moment()
  ```

#### 🔄 UV-Vis Spectroscopy
- **File**: `/home/user/kanad/kanad/analysis/spectroscopy.py`
- **Hybrid Modes**:
  1. **Classical**: TDA, TDDFT, CIS (lines 112-220)
  2. **Quantum**: `quantum_sqd` method (lines 221+)
- **First Production Quantum UV-Vis Calculator**:
  ```python
  result = uv_calc.compute_excitations(
      n_states=5, 
      method='quantum_sqd',
      backend='ibm',  # Real quantum hardware!
      subspace_dim=15
  )
  ```

#### 🔄 DOS Calculator
- **File**: `/home/user/kanad/kanad/analysis/dos_calculator.py`
- **Quantum Method**: `compute_quantum_dos()` (lines 473+)
- **World's First**: Governance-aware quantum DOS
- **Usage**:
  ```python
  dos = DOSCalculator().compute_quantum_dos(
      bond=h2_bond, 
      backend='ibm',
      subspace_dim=20
  )
  ```

#### 🔄 NMR Calculator
- **File**: `/home/user/kanad/kanad/analysis/nmr_calculator.py`
- **Quantum Enhancement**: `compute_quantum_chemical_shifts()` (lines 357+)
- **Hybrid Approach**: Classical magnetic integrals + quantum density

#### 🔄 Raman/IR Spectroscopy
- **File**: `/home/user/kanad/kanad/analysis/raman_calculator.py`
- **Quantum Polarizability**: `_compute_quantum_polarizability()` (lines 394+)
- **Can run on quantum hardware** for improved accuracy

#### 🔄 Thermochemistry
- **File**: `/home/user/kanad/kanad/analysis/thermochemistry.py`
- **Quantum Method**: `compute_quantum_thermochemistry()` (lines 514+)
- **Enhancement**: Uses quantum vibrational frequencies

#### 🔄 Vibrational Analysis
- **File**: `/home/user/kanad/kanad/analysis/vibrational_analysis.py`
- **Quantum Spectroscopy**: `compute_quantum_vibronic_spectrum()` (lines 890+)
- **Hardware Support**: IBM, BlueQubit

### 2.2 Application Modules (Orchestrate Quantum Solvers)

#### 🔄 Drug Discovery Platform
- **File**: `/home/user/kanad/kanad/applications/drug_discovery.py`
- **Quantum Integration**:
  - Binding affinity calculations using VQE/SQD
  - pH-dependent properties with quantum protonation states
  - Metabolite prediction with quantum transition states
- **Competitive Advantage**: <1 kcal/mol error vs 3 kcal/mol for classical
- **Usage**: Orchestrates quantum solvers internally

#### 🔄 Materials Scout
- **File**: `/home/user/kanad/kanad/applications/materials_scout.py`
- **Quantum Features**:
  - Band structure from quantum DOS
  - Quantum-enhanced property predictions

#### 🔄 Catalyst Optimizer
- **File**: `/home/user/kanad/kanad/applications/catalyst_optimizer.py`
- **Quantum Integration**: Reaction barriers from quantum solvers

#### 🔄 Alloy Designer
- **File**: `/home/user/kanad/kanad/applications/alloy_designer.py`
- **Quantum Features**: Electronic structure from quantum calculations

---

## 3. PURELY CLASSICAL CAPABILITIES
### No Quantum Hardware Execution

### 3.1 Classical Simulations
- **Statevector Simulation**: Exact classical simulation (exponential cost)
  ```python
  solver = VQESolver(bond, backend='statevector')  # Classical simulation
  ```
- **Hartree-Fock**: SCF calculations (all Hamiltonians)
- **Classical TD-DFT**: Excited states without quantum hardware
- **Classical NMR**: Using classical density matrices

### 3.2 Molecular Construction
- **Files**: `kanad/io/`, `kanad/core/molecule.py`
- **Pure Classical**: Geometry parsing, coordinate transformations

### 3.3 Governance Protocols
- **Files**: `kanad/governance/protocols/`
- **Pure Classical**: Configuration validation, excitation filtering

---

## 4. BACKEND FLAG IMPLEMENTATION

### 4.1 Architecture Pattern

All solvers follow this pattern:
```python
class SomeSolver(BaseSolver):
    def __init__(self, bond, backend='statevector', **kwargs):
        self.backend = backend
        self._init_backend(**kwargs)
    
    def _init_backend(self, **kwargs):
        if backend == 'statevector':
            self._use_statevector = True
        elif backend == 'bluequbit':
            from kanad.backends.bluequbit import BlueQubitBackend
            self._bluequbit_backend = BlueQubitBackend(**kwargs)
            self._use_statevector = False
        elif backend == 'ibm':
            from kanad.backends.ibm import IBMBackend
            self._ibm_backend = IBMBackend(**kwargs)
            self._use_statevector = False
```

### 4.2 Execution Routing

**VQE Example** (`vqe_solver.py`, lines 467-1041):
```python
def _compute_energy(self, parameters):
    if self._use_statevector:
        return self._compute_energy_statevector(parameters)
    else:
        return self._compute_energy_quantum(parameters)

def _compute_energy_quantum(self, parameters):
    # Route to correct backend
    if hasattr(self, '_ibm_backend') and self._ibm_backend is not None:
        # IBM quantum hardware path
        result = self._ibm_backend.run_batch(...)
        job_id = result['job_id']
        # Wait for hardware execution
        job_result = job.result()
        return energy
    
    elif hasattr(self, '_bluequbit_backend') and self._bluequbit_backend is not None:
        # BlueQubit cloud path
        job_info = self._bluequbit_backend.run_circuit(...)
        result = self._bluequbit_backend.wait_for_job(job_id)
        return energy
```

### 4.3 Backend Configuration

**API Layer** (`api/routes/experiments.py`):
```python
class BackendConfig(BaseModel):
    method: str = "VQE"
    backend: str = "classical"  # classical, ibm_quantum, bluequbit
    backend_name: Optional[str] = None  # For IBM: ibm_brisbane, etc.
    bluequbit_device: Optional[str] = "gpu"  # For BlueQubit
```

**Experiment Service** maps this to solver initialization:
- `backend='classical'` → `backend='statevector'` (local simulation)
- `backend='ibm_quantum'` → `backend='ibm'` + IBM credentials
- `backend='bluequbit'` → `backend='bluequbit'` + BlueQubit credentials

---

## 5. QUANTUM EXECUTION FLOW EXAMPLES

### 5.1 VQE on IBM Quantum
```python
from kanad.bonds import BondFactory
from kanad.solvers import VQESolver

# Create bond
bond = BondFactory.create_bond('H', 'H', distance=0.74)

# Initialize VQE with IBM backend
solver = VQESolver(
    bond=bond,
    ansatz_type='ucc',
    mapper_type='jordan_wigner',
    optimizer='SPSA',  # Auto-selected for cloud
    max_iterations=50,
    backend='ibm',
    backend_name='ibm_brisbane',
    api_token='YOUR_IBM_TOKEN'
)

# Solve (submits to real quantum hardware)
result = solver.solve()
# Output:
# 🔧 Initializing backend: ibm
# ✅ Connected to IBM Quantum: ibm_brisbane
# 🚀 Submitting job to IBM Quantum (function eval 1)
# ✅ IBM job submitted: cx7g3...
# 🔗 Track job at: https://quantum.ibm.com/jobs/cx7g3...
```

### 5.2 SQD on BlueQubit
```python
from kanad.solvers import SQDSolver

solver = SQDSolver(
    bond=bond,
    subspace_dim=10,
    backend='bluequbit',
    device='gpu',
    api_token='YOUR_BLUEQUBIT_TOKEN'
)

result = solver.solve()
# Output:
# 🔧 BlueQubit device selected: gpu
# ✅ BlueQubit job submitted: bq_job_123
# 🔗 Track job at: https://app.bluequbit.io/jobs/bq_job_123
```

### 5.3 Quantum-Enhanced Analysis
```python
# Step 1: Run quantum solver
vqe_result = VQESolver(bond, backend='ibm').solve()
# This populates hamiltonian.quantum_rdm1

# Step 2: Analysis automatically uses quantum density
from kanad.analysis import PropertyCalculator
calc = PropertyCalculator(bond.hamiltonian)

# Uses quantum density if available, falls back to HF
dipole = calc.compute_dipole_moment()
print(f"Quantum dipole: {dipole['dipole_magnitude']:.3f} D")

# Explicit quantum method
from kanad.analysis import UVVisCalculator
uv_calc = UVVisCalculator(bond.molecule)
spectrum = uv_calc.compute_excitations(
    n_states=5,
    method='quantum_sqd',  # Runs on quantum hardware
    backend='ibm'
)
```

---

## 6. QUANTUM VS CLASSICAL DECISION TREE

```
User Request
    │
    ├─ backend='statevector'
    │   └─> Classical simulation (exact, exponential cost)
    │
    ├─ backend='ibm'
    │   ├─> Initialize IBMBackend
    │   ├─> Transpile circuits for hardware
    │   ├─> Submit to IBM Quantum
    │   ├─> Apply error mitigation (ZNE, DD)
    │   └─> Return quantum result
    │
    ├─ backend='bluequbit'
    │   ├─> Initialize BlueQubitBackend
    │   ├─> Select device (gpu/cpu/mps)
    │   ├─> Submit to BlueQubit cloud
    │   └─> Return quantum result
    │
    └─ backend='classical' (API alias)
        └─> Maps to 'statevector'
```

---

## 7. KEY IMPLEMENTATION DETAILS

### 7.1 Quantum Density Matrix Propagation

VQE/SQD compute quantum density and store it in the Hamiltonian:
```python
# vqe_solver.py, lines 1613-1615
if hasattr(self.hamiltonian, 'set_quantum_density_matrix'):
    self.hamiltonian.set_quantum_density_matrix(quantum_density)
```

Analysis modules retrieve it:
```python
# property_calculator.py, lines 85-86
if hasattr(self.hamiltonian, 'get_density_matrix'):
    density_matrix = self.hamiltonian.get_density_matrix()  # Quantum!
```

### 7.2 Auto-Optimization for Cloud Backends

VQE automatically switches optimizers for efficiency:
```python
# vqe_solver.py, lines 1473-1500
if self.backend in ['ibm', 'bluequbit']:
    if self.optimizer_method not in ['SPSA', 'COBYLA']:
        print("Auto-switching to SPSA")
        print("Expected speedup: 20x fewer quantum jobs")
        self.optimizer_method = 'SPSA'
```

### 7.3 Cancellation Support

Both backends support job cancellation:
```python
# vqe_solver.py, lines 855-875 (IBM)
if self.experiment_id and self.job_id:
    check_cancellation(self.experiment_id, self.job_id)
    # If cancelled, cancel IBM job and raise
    self._ibm_backend.cancel_job(job_id)
```

### 7.4 Cloud Job Tracking

Results include cloud provider information:
```python
# vqe_solver.py, lines 1716-1733
if self.cloud_provider:
    self.results['cloud_provider'] = self.cloud_provider
    self.results['cloud_job_ids'] = self.cloud_job_ids
    if self.cloud_provider == 'ibm':
        self.results['cloud_job_urls'] = [
            f"https://quantum.ibm.com/jobs/{jid}" 
            for jid in self.cloud_job_ids
        ]
```

---

## 8. SUMMARY TABLE

| Component | Quantum Hardware | Quantum-Enhanced | Classical |
|-----------|-----------------|------------------|-----------|
| **VQESolver** | ✅ IBM, BlueQubit | - | ✅ Statevector |
| **SQDSolver** | ✅ IBM, BlueQubit | - | ✅ Statevector |
| **ExcitedStatesSolver** | ✅ QPE, VQE, SQD | - | ✅ CIS, TDDFT |
| **UV-Vis Spectroscopy** | ✅ quantum_sqd | - | ✅ TDA, TDDFT |
| **DOS Calculator** | ✅ compute_quantum_dos | - | ✅ Classical DOS |
| **NMR Calculator** | - | ✅ Quantum density | ✅ HF density |
| **Raman/IR** | ✅ Quantum polarizability | - | ✅ Classical |
| **Thermochemistry** | ✅ Quantum frequencies | - | ✅ Classical |
| **Property Calculator** | - | ✅ Quantum density | ✅ HF density |
| **Drug Discovery** | - | ✅ VQE/SQD binding | ✅ Force field |
| **Hartree-Fock** | - | - | ✅ Only classical |

**Legend:**
- ✅ = Supported
- - = Not applicable

---

## 9. FUTURE QUANTUM CAPABILITIES (Roadmap)

From `QUANTUM_ROADMAP_NEXT_STEPS.md`:

1. **Quantum Molecular Dynamics**: Time evolution on quantum hardware
2. **Quantum Machine Learning**: Feature extraction from quantum states
3. **Quantum Error Correction**: Logical qubits for longer circuits
4. **Quantum Chemistry ML**: Train on quantum data for classical predictions

---

## Conclusion

The Kanad framework implements a sophisticated 3-tier quantum architecture:
1. **Quantum solvers** (VQE, SQD, Hi-VQE) execute directly on IBM Quantum and BlueQubit hardware
2. **Analysis modules** seamlessly integrate quantum results when available
3. **Classical fallbacks** ensure functionality without quantum access

The `backend` parameter controls execution mode throughout the framework, with automatic optimizations for cloud backends (SPSA optimizer, error mitigation, job tracking).

This architecture enables **quantum advantage today** while maintaining **classical compatibility** for development and testing.
