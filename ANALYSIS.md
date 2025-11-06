# Gap Analysis: Current State vs Hi-VQE Requirements

**Date:** November 4, 2025

---

## Current VQE Flow (What We Have)

```
User Request (H2O, 20 iterations)
    ↓
1. Build Molecule → 10 electrons, 7 orbitals
    ↓
2. Build Hamiltonian → Jordan-Wigner → 1086 Pauli terms on 14 qubits
    ↓
3. Build Ansatz → Governance → 84 parameters, deep circuit
    ↓
4. VQE Loop (20 iterations):
    For each iteration:
        - Adjust 84 parameters (optimizer: SLSQP)
        - Build quantum circuit with new parameters
        - Measure ALL 1086 Pauli terms: <ψ|P_i|ψ> for each term
        - Sum up: E = Σ coeff_i * <ψ|P_i|ψ>
        - Optimizer uses E to update parameters (needs gradient!)
        - Gradient = N+1 more function evals (85 total per iteration!)
    ↓
5. Result: Energy = -74.82 Ha (WORSE than HF -74.96 Ha!)
   Total measurements: 20 iter × 85 evals × 1086 Pauli = 1.8 million measurements!
   Status: FAILED - got stuck
```

**Problems:**
- ❌ Must measure 1086 Pauli terms per function eval
- ❌ 85 function evals per iteration (gradient computation)
- ❌ Gets stuck at wrong energy
- ❌ 1.8 MILLION total measurements for H2O

---

## Hi-VQE Flow (What We Need)

```
User Request (H2O, 10 iterations)
    ↓
1. Build Molecule → 10 electrons, 7 orbitals
    ↓
2. Active Space Selection → 6 active orbitals (freeze O 1s)
    ↓
3. Build Hamiltonian → Jordan-Wigner → ~500 Pauli terms on 12 qubits
    ↓
4. Initialize Subspace → Start with HF configuration: |1100000000⟩
    ↓
5. Hi-VQE Loop (10 iterations):
    For each iteration:
        a) Prepare quantum state with current parameters
        b) Sample configurations (ONE measurement in Z basis):
           Measure: |ψ⟩ → bitstrings (e.g., |1100000000⟩, |1010000001⟩, ...)
           Get ~100 configurations from 1000 shots
        c) Filter valid configurations (correct electron count)
        d) Add to subspace: S = S ∪ {new configs}
        e) Classical diagonalization:
           - Project H into subspace S: H_sub = ⟨config_i|H|config_j⟩
           - Exact solve: E, ψ = diagonalize(H_sub)
           - NO quantum measurements needed!
        f) Identify important configs (high amplitude in ψ)
        g) Generate excitations from important configs
        h) Prune low-amplitude configs (amplitude < 10^-6)
        i) Optimize parameters to sample MORE important configs
    ↓
6. Result: Energy = -75.06 Ha (matches FCI!)
   Total measurements: 10 iterations × 1 measurement = 10 measurements!
   Status: SUCCESS - chemical accuracy
```

**Advantages:**
- ✅ Only 1 measurement per iteration (Z basis sampling)
- ✅ Exact energy from classical diagonalization
- ✅ 10 measurements total (vs 1.8 million!)
- ✅ 180,000x more efficient!

---

## Gap Analysis: What's Missing

### ❌ MISSING Component #1: Active Space Selection

**What we have:**
```python
# kanad/core/hamiltonians/openfermion_jw.py
def build_h2_hamiltonian(molecule):
    # Uses ALL orbitals
    n_orbitals = molecule.n_orbitals  # 7 for H2O
    n_qubits = 2 * n_orbitals  # 14 qubits
```

**What we need:**
```python
def build_h2_hamiltonian(molecule, active_space='auto'):
    # Auto-select active space
    if active_space == 'auto':
        frozen_orbitals, active_orbitals = select_active_space(molecule)
        # H2O: freeze [0] (O 1s), active [1,2,3,4,5,6]
        # Result: 6 active → 12 qubits (not 14)
```

**Where to add:**
- File: `kanad/core/hamiltonians/openfermion_jw.py`
- Add function: `select_active_space(molecule)`
- Modify: `build_h2_hamiltonian()` to accept `active_space` parameter

---

### ❌ MISSING Component #2: Configuration Sampling

**What we have:**
```python
# vqe_solver.py _objective_function()
def _objective_function(self, params):
    circuit = self.ansatz.build_circuit(params)
    # Measure energy by measuring ALL Pauli terms
    energy = sum(coeff * measure_pauli(circuit, pauli)
                 for pauli, coeff in hamiltonian)
    return energy
```

**What we need:**
```python
def sample_configurations(self, params, shots=1000):
    circuit = self.ansatz.build_circuit(params)
    # Measure in computational basis ONLY (Z basis)
    result = execute(circuit, backend, shots=shots)
    counts = result.get_counts()

    # Convert to configurations
    configs = []
    for bitstring, count in counts.items():
        configs.extend([bitstring] * count)

    # Filter: keep only configs with correct electron count
    n_electrons = self.molecule.n_electrons
    valid = [c for c in configs if c.count('1') == n_electrons]

    return valid
```

**Where to add:**
- New file: `kanad/core/configuration.py`
- Or add to: `kanad/utils/vqe_solver.py` as new method

---

### ❌ MISSING Component #3: Subspace Management

**What we have:**
- Nothing! We don't track configurations at all.

**What we need:**
```python
class ConfigurationSubspace:
    def __init__(self, n_qubits, n_electrons):
        # Start with HF configuration
        hf_config = '1' * n_electrons + '0' * (n_qubits - n_electrons)
        self.configs = [hf_config]
        self.config_to_index = {hf_config: 0}

    def add_configs(self, new_configs):
        for config in new_configs:
            if config not in self.config_to_index:
                idx = len(self.configs)
                self.configs.append(config)
                self.config_to_index[config] = idx

    def generate_single_excitations(self, config):
        """Generate all single excitations from config."""
        excitations = []
        for i in range(len(config)):
            if config[i] == '1':
                for j in range(len(config)):
                    if config[j] == '0':
                        new = list(config)
                        new[i], new[j] = '0', '1'
                        excitations.append(''.join(new))
        return excitations

    def prune(self, amplitudes, threshold=1e-6):
        """Remove low-amplitude configurations."""
        keep_indices = [i for i, amp in enumerate(amplitudes) if abs(amp) > threshold]
        self.configs = [self.configs[i] for i in keep_indices]
        self.config_to_index = {c: i for i, c in enumerate(self.configs)}
```

**Where to add:**
- New file: `kanad/core/subspace.py`

---

### ❌ MISSING Component #4: Classical Diagonalization

**What we have:**
- Nothing! We always use quantum measurement.

**What we need:**
```python
def project_hamiltonian(hamiltonian, subspace):
    """
    Project Hamiltonian into configuration subspace.

    Returns:
        H_sub: numpy array (n_configs × n_configs)
    """
    n = len(subspace.configs)
    H_sub = np.zeros((n, n), dtype=complex)

    for i, config_i in enumerate(subspace.configs):
        for j, config_j in enumerate(subspace.configs):
            # Compute <config_i|H|config_j>
            H_sub[i, j] = compute_matrix_element(hamiltonian, config_i, config_j)

    return H_sub

def compute_matrix_element(hamiltonian, bra, ket):
    """
    Compute <bra|H|ket> for two configurations.

    This is the KEY function - converts Pauli operators to configuration basis.
    """
    if bra == ket:
        # Diagonal: <config|H|config>
        # Sum all Z terms that match config
        energy = 0.0
        for pauli_string, coeff in hamiltonian:
            if applies_diagonal(pauli_string, bra):
                sign = get_sign(pauli_string, bra)
                energy += coeff * sign
        return energy
    else:
        # Off-diagonal: <config_i|H|config_j>
        # Only single/double excitation operators contribute
        energy = 0.0
        for pauli_string, coeff in hamiltonian:
            if connects_configs(pauli_string, bra, ket):
                matrix_elem = compute_pauli_element(pauli_string, bra, ket)
                energy += coeff * matrix_elem
        return energy

def diagonalize_subspace(H_sub):
    """
    Exact diagonalization of subspace Hamiltonian.

    Returns:
        energy: Ground state energy (exact!)
        amplitudes: Ground state wavefunction coefficients
    """
    eigenvalues, eigenvectors = np.linalg.eigh(H_sub)
    ground_energy = eigenvalues[0]
    ground_state = eigenvectors[:, 0]
    return ground_energy, ground_state
```

**Where to add:**
- New file: `kanad/core/classical_solver.py`
- This is CRITICAL - need to implement Pauli → configuration matrix elements correctly

---

### 🔧 Component to MODIFY: VQE Solver

**Current solve() method:**
```python
def solve(self, initial_parameters=None):
    # Standard VQE optimization
    result = minimize(
        self._objective_function,  # Measures all Pauli terms
        initial_parameters,
        method=self.optimizer_method
    )
    return result
```

**Need to add Hi-VQE mode:**
```python
def solve(self, initial_parameters=None, mode='standard'):
    if mode == 'hivqe':
        return self._solve_hivqe(initial_parameters)
    else:
        return self._solve_standard(initial_parameters)

def _solve_hivqe(self, initial_parameters):
    """Hi-VQE: Configuration sampling + classical diagonalization."""

    # Initialize subspace with HF configuration
    subspace = ConfigurationSubspace(self.n_qubits, self.n_electrons)

    # Get HF config
    hf_config = get_hf_configuration(self.molecule)

    params = initial_parameters if initial_parameters is not None else np.random.randn(self.n_parameters) * 0.1

    for iteration in range(self.max_iterations):
        # 1. Sample configurations (1 measurement!)
        configs = self.sample_configurations(params, shots=1000)

        # 2. Add to subspace
        subspace.add_configs(configs)

        # 3. Classical diagonalization (NO quantum measurements!)
        H_sub = project_hamiltonian(self.hamiltonian, subspace)
        energy, amplitudes = diagonalize_subspace(H_sub)

        print(f"Iteration {iteration+1}: E = {energy:.8f} Ha, subspace size = {len(subspace.configs)}")

        # 4. Generate excitations from high-amplitude configs
        important_configs = [subspace.configs[i] for i, a in enumerate(amplitudes)
                           if abs(a) > 0.1]
        for config in important_configs:
            excitations = subspace.generate_single_excitations(config)
            subspace.add_configs(excitations)

        # 5. Prune low-amplitude configs
        subspace.prune(amplitudes, threshold=1e-6)

        # 6. Optimize parameters to sample important configs better
        # This is where we'd update params, but for now keep simple

        # Check convergence
        if iteration > 0 and abs(energy - prev_energy) < self.conv_threshold:
            break
        prev_energy = energy

    return {
        'energy': energy,
        'iterations': iteration + 1,
        'subspace_size': len(subspace.configs),
        'amplitudes': amplitudes
    }
```

**Where to modify:**
- File: `kanad/utils/vqe_solver.py`
- Add: `_solve_hivqe()` method
- Add: `sample_configurations()` method
- Modify: `solve()` to accept `mode` parameter

---

### ❌ MISSING Component #5: Qubit Reduction Utilities

**What we need:**
```python
def select_active_space(molecule, method='valence'):
    """
    Automatically select active orbitals.

    For H2O:
        - Total: 7 orbitals (O: 1s, 2s, 2px, 2py, 2pz; H: 1s, 1s)
        - Frozen: O 1s (core, doesn't participate in bonding)
        - Active: 6 orbitals → 12 qubits
    """
    if method == 'valence':
        # Freeze core orbitals
        core_orbitals = identify_core_orbitals(molecule)
        all_orbitals = list(range(molecule.n_orbitals))
        active_orbitals = [i for i in all_orbitals if i not in core_orbitals]
        return core_orbitals, active_orbitals

def identify_core_orbitals(molecule):
    """Identify core orbitals to freeze."""
    core = []
    for atom in molecule.atoms:
        if atom.symbol == 'O':
            # Oxygen: freeze 1s
            core.append(0)  # Assuming first orbital is O 1s
        # H has no core electrons
    return core

def apply_frozen_core(hamiltonian, frozen_orbitals):
    """
    Remove frozen orbital terms from Hamiltonian.

    Returns:
        reduced_hamiltonian: Hamiltonian in active space only
        frozen_core_energy: Classical contribution from frozen electrons
    """
    # This requires understanding Hamiltonian structure
    pass
```

**Where to add:**
- New file: `kanad/core/active_space.py`

---

## What to Remove

### 🗑️ Remove #1: Specialized Hamiltonians (Not Used)
```bash
kanad/core/hamiltonians/ionic_hamiltonian.py      # Not used, adds complexity
kanad/core/hamiltonians/metallic_hamiltonian.py   # Not used
kanad/core/hamiltonians/periodic_hamiltonian.py   # Not used
```

### 🗑️ Remove #2: Broken Ansatz
```bash
kanad/ansatze/governance_optimized.py  # Broken MP2 initialization
```

### 🗑️ Remove #3: Test Files Already Archived
- Already moved to `archive/old_tests/` ✓

### 🗑️ Remove #4: Documentation Bloat
- Already moved to `archive/old_docs/` ✓

---

## What to Change

### 🔧 Change #1: Hamiltonian Builder
**File:** `kanad/core/hamiltonians/openfermion_jw.py`

**Current:**
```python
def build_h2_hamiltonian(molecule):
    n_orbitals = molecule.n_orbitals  # Uses ALL orbitals
```

**Change to:**
```python
def build_h2_hamiltonian(molecule, active_space=None, frozen_core=True):
    if frozen_core and active_space is None:
        frozen, active = select_active_space(molecule)
        # Use only active orbitals
        n_orbitals = len(active)  # Reduced!
    else:
        n_orbitals = molecule.n_orbitals
```

### 🔧 Change #2: API to Accept Active Space
**File:** `api/services/experiment_service.py`

**Add to config:**
```python
config = {
    'molecule': 'H2O',
    'basis': 'sto-3g',
    'ansatz': 'governance',
    'optimizer': 'SLSQP',
    'active_space': 'auto',  # NEW: auto, valence, or manual
    'frozen_core': True,      # NEW: freeze core electrons
}
```

### 🔧 Change #3: VQE Solver to Support Hi-VQE Mode
**File:** `kanad/utils/vqe_solver.py`

**Add parameter:**
```python
def solve(self, mode='standard'):  # or mode='hivqe'
```

---

## Implementation Priority

### Phase 1: Qubit Reduction (Immediate Win)
1. Add `select_active_space()` function
2. Modify `build_h2_hamiltonian()` to use active space
3. Test: H2O should use 12 qubits, not 14

**Impact:** Reduces circuit depth, fewer Pauli terms

### Phase 2: Configuration Sampling
1. Add `ConfigurationSubspace` class
2. Add `sample_configurations()` to VQE solver
3. Test: Can we sample valid configurations?

**Impact:** Foundation for Hi-VQE

### Phase 3: Classical Diagonalization
1. Add `project_hamiltonian()` function
2. Implement `compute_matrix_element()` (Pauli → config basis)
3. Add `diagonalize_subspace()`
4. Test: Can we get exact energy in small subspace?

**Impact:** This is THE key to Hi-VQE efficiency

### Phase 4: Hi-VQE Integration
1. Add `_solve_hivqe()` method to VQE solver
2. Iterative subspace growth
3. Test: H2O with 10 iterations, chemical accuracy

**Impact:** 180,000x efficiency improvement!

---

## Key Technical Challenges

### Challenge #1: Pauli → Configuration Matrix Elements
**Problem:** Need to compute `<config_i|Pauli_string|config_j>`

**Example:** For Pauli string `X₀Y₁Z₂`:
```
X₀ flips qubit 0
Y₁ flips qubit 1 (with i phase)
Z₂ measures qubit 2 (±1 depending on state)

For |110⟩ → |config⟩:
X₀|110⟩ = |010⟩
Y₁|010⟩ = i|000⟩
Z₂|000⟩ = |000⟩

So: <110|X₀Y₁Z₂|110⟩ connects |110⟩ to |000⟩
```

**Need to implement:** Proper Pauli operator application in configuration basis

### Challenge #2: Efficient Subspace Growth
**Problem:** Subspace can grow exponentially if not pruned

**Solution:**
- Prune configs with amplitude < 10⁻⁶
- Only generate excitations from high-amplitude configs
- Typical size: 100-1000 configs for H2O

### Challenge #3: Parameter Optimization in Hi-VQE
**Problem:** What should circuit parameters optimize for?

**Answer:** Optimize to sample MORE of the important configurations
- Use amplitudes from current iteration
- Adjust params so circuit has higher probability for high-amplitude configs
- This is different from minimizing energy!

---

## Testing Strategy

### Test 1: Active Space
```python
# H2O should use 12 qubits, not 14
molecule = Molecule.from_smiles("O")
H = build_h2_hamiltonian(molecule, frozen_core=True)
assert H.n_qubits == 12  # Not 14!
```

### Test 2: Configuration Sampling
```python
# Can we sample valid configurations?
configs = solver.sample_configurations(params, shots=1000)
assert all(c.count('1') == molecule.n_electrons for c in configs)
```

### Test 3: Subspace Diagonalization
```python
# Small test: H2 with 2 configurations
subspace = ConfigurationSubspace(4, 2)
subspace.add_configs(['1100', '1001'])  # HF and single excitation
H_sub = project_hamiltonian(H, subspace)
energy, amps = diagonalize_subspace(H_sub)
# Should get exact energy in this 2-config subspace
```

### Test 4: Full Hi-VQE
```python
# H2O with Hi-VQE mode
solver = VQESolver(molecule=h2o, mode='hivqe', max_iterations=10)
result = solver.solve()
assert abs(result['energy'] - fci_energy) < 0.0016  # Chemical accuracy
assert result['iterations'] <= 10
```

---

## Summary

### What We Have:
- ✅ Basic VQE framework
- ✅ Ansatz library
- ✅ Hamiltonian builders
- ✅ Standard optimizers

### What We're Missing:
- ❌ Active space selection → 14 qubits instead of 12 for H2O
- ❌ Configuration sampling → Can't do Hi-VQE without this
- ❌ Subspace management → No tracking of configurations
- ❌ Classical diagonalization → Can't get exact energy in subspace
- ❌ Hi-VQE mode in solver → Just has standard VQE

### Implementation Order:
1. Active space (reduce qubits) - **EASY, BIG IMPACT**
2. Configuration sampling - **MODERATE DIFFICULTY**
3. Classical diagonalization - **HARD (need Pauli matrix elements)**
4. Hi-VQE integration - **MODERATE (glue it together)**

### Expected Results After Implementation:
- H2O: 12 qubits (not 14)
- H2O: 10 iterations (not 20+)
- H2O: 10 measurements (not 1.8 million!)
- H2O: -75.06 Ha ± 0.001 Ha (chemical accuracy)
- Larger molecules: Actually feasible on cloud quantum hardware

---

**End of Analysis**
