# Hi-VQE Implementation Plan - Production Ready for Large Molecules

**Goal:** Achieve Hi-VQE efficiency for large molecules on cloud quantum hardware

---

## Current State (H2O):
- ❌ 14 qubits (should be 12)
- ❌ 1086 Pauli terms measured per eval
- ❌ 84 parameters
- ❌ ~1800 function evals
- ❌ Energy stuck at -74.82 Ha (worse than HF -74.96 Ha!)
- ❌ Timeouts after 3 minutes

## Target State (H2O):
- ✅ 12 qubits (active space + symmetries)
- ✅ 1 computational basis measurement per iteration
- ✅ 10-20 parameters
- ✅ 5-10 iterations
- ✅ Chemical accuracy (< 1 kcal/mol error)
- ✅ Complete in < 1 minute

---

## Phase 1: IMMEDIATE WINS (Today - 2 hours)

### 1.1 Add SPSA Optimizer (30 min)
**File:** `kanad/utils/vqe_solver.py`
**Impact:** 42x reduction in function evals for cloud

```python
def _spsa_minimize(self, initial_params):
    """SPSA: Only 2 evals per iteration regardless of parameter count!"""
    params = initial_params.copy()

    for iteration in range(self.max_iterations):
        # Random perturbation
        delta = 2 * np.random.randint(0, 2, size=len(params)) - 1

        # Only 2 function evaluations
        c = 0.1 / (iteration + 1)**0.602  # Decay step size
        a = 0.16 / (iteration + 1)**0.101  # Decay perturbation

        params_plus = params + a * delta
        params_minus = params - a * delta

        f_plus = self._objective_function(params_plus)
        f_minus = self._objective_function(params_minus)

        # Gradient estimate
        grad = (f_plus - f_minus) / (2 * a) * delta

        # Update
        params = params - c * grad

        if abs(f_plus - f_minus) < self.conv_threshold:
            break

    return params
```

**Test:** Run H2 with SPSA, verify 2 evals/iteration

### 1.2 Fix Active Space for H2O (30 min)
**File:** `kanad/bonds/covalent_bond.py` or molecule creation

**H2O Active Space:**
- Total: 7 spatial orbitals (14 spin-orbitals)
- Core (freeze): O 1s (2 electrons) → Classical treatment
- Active: O 2s, 2px, 2py, 2pz, H 1s (2 orbitals) = 6 spatial → **12 qubits**

```python
def get_active_space(molecule):
    """Select active orbitals for molecule."""
    if molecule.name == "H2O":
        # Freeze O 1s core electrons
        frozen_orbitals = [0]  # O 1s orbital
        active_orbitals = [1, 2, 3, 4, 5, 6]  # Valence orbitals
        return frozen_orbitals, active_orbitals
```

**Test:** Build H2O Hamiltonian, verify 12 qubits (not 14)

### 1.3 Backend-Specific Optimizer Selection (15 min)
**File:** `api/services/experiment_service.py`

```python
# Smart optimizer selection
if backend in ['ibm', 'bluequbit']:
    if 'optimizer' not in config:
        config['optimizer'] = 'SPSA'  # 2 evals/iter for cloud
        print(f"📊 Cloud backend: Using SPSA optimizer (2 evals/iteration)")
else:
    if 'optimizer' not in config:
        config['optimizer'] = 'SLSQP'  # Gradient-based for fast classical
```

**Test:** API request with IBM backend uses SPSA automatically

### 1.4 Test Phase 1 (45 min)
```bash
# Test 1: H2 with SPSA
python test_h2_spsa.py  # Should use ~20 function evals for 10 iterations

# Test 2: H2O with 12 qubits
python test_h2o_active_space.py  # Should show 12 qubits, not 14

# Test 3: Cloud backend auto-selects SPSA
curl -X POST /api/experiments -d '{"backend": "ibm", "molecule": "H2"}'
```

---

## Phase 2: CONFIGURATION SAMPLING (Tomorrow - 4 hours)

### 2.1 Add Configuration Sampler (1 hour)
**New file:** `kanad/core/configuration_sampling.py`

```python
class ConfigurationSampler:
    def sample_configs(self, circuit, shots=1000):
        """Sample computational basis states."""
        # Measure in Z basis only (no Pauli rotations!)
        result = execute(circuit, backend, shots=shots).result()
        counts = result.get_counts()

        configs = []
        for bitstring, count in counts.items():
            configs.extend([bitstring] * count)

        return configs

    def filter_valid(self, configs, n_electrons):
        """Keep only configs with correct electron count."""
        valid = []
        for config in configs:
            if config.count('1') == n_electrons:
                valid.append(config)
        return valid
```

### 2.2 Add Subspace Manager (1 hour)
**New file:** `kanad/core/subspace.py`

```python
class ConfigurationSubspace:
    def __init__(self, hf_config):
        self.configs = [hf_config]

    def add_configs(self, new_configs):
        """Add new configurations to subspace."""
        for config in new_configs:
            if config not in self.configs:
                self.configs.append(config)

    def generate_excitations(self, config):
        """Generate single/double excitations from config."""
        excitations = []
        # Single excitations
        for i in range(len(config)):
            if config[i] == '1':
                for j in range(len(config)):
                    if config[j] == '0':
                        new_config = list(config)
                        new_config[i], new_config[j] = '0', '1'
                        excitations.append(''.join(new_config))
        return excitations

    def prune(self, amplitudes, threshold=1e-6):
        """Remove low-amplitude configurations."""
        keep = []
        for i, amp in enumerate(amplitudes):
            if abs(amp) > threshold:
                keep.append(self.configs[i])
        self.configs = keep
```

### 2.3 Add Classical Diagonalization (1 hour)
**New file:** `kanad/core/classical_diag.py`

```python
def project_hamiltonian(hamiltonian, subspace):
    """Project Hamiltonian into configuration subspace."""
    n = len(subspace)
    H_sub = np.zeros((n, n))

    for i, config_i in enumerate(subspace):
        for j, config_j in enumerate(subspace):
            # Compute <config_i|H|config_j>
            H_sub[i, j] = hamiltonian_element(hamiltonian, config_i, config_j)

    return H_sub

def hamiltonian_element(hamiltonian, bra, ket):
    """Compute <bra|H|ket> matrix element."""
    if bra == ket:
        # Diagonal: Sum single-particle energies
        return sum_diagonal_terms(hamiltonian, bra)
    else:
        # Off-diagonal: Check if single/double excitation
        diff = hamming_distance(bra, ket)
        if diff == 2:  # Single excitation
            return single_excitation_element(hamiltonian, bra, ket)
        elif diff == 4:  # Double excitation
            return double_excitation_element(hamiltonian, bra, ket)
        else:
            return 0.0  # Not connected
```

### 2.4 Integrate into VQE Solver (1 hour)
**File:** `kanad/utils/vqe_solver.py`

Add mode selection:
```python
def solve(self, mode='standard'):
    """
    mode='standard': Traditional VQE (measure all Pauli terms)
    mode='qsci': Quantum Selected CI (configuration sampling)
    """
    if mode == 'qsci':
        return self._solve_qsci()
    else:
        return self._solve_standard()

def _solve_qsci(self):
    """Hi-VQE style: Configuration sampling + subspace diagonalization."""
    subspace = ConfigurationSubspace(hf_config)

    for iteration in range(self.max_iterations):
        # Sample configurations
        configs = self.sampler.sample_configs(self.circuit, shots=1000)
        valid_configs = self.sampler.filter_valid(configs, self.n_electrons)

        # Add to subspace
        subspace.add_configs(valid_configs)

        # Classical diagonalization
        H_sub = project_hamiltonian(self.hamiltonian, subspace.configs)
        energy, amplitudes = np.linalg.eigh(H_sub)
        energy = energy[0]
        amplitudes = amplitudes[:, 0]

        # Generate classical excitations
        important_configs = [subspace.configs[i] for i, a in enumerate(amplitudes) if abs(a) > 0.1]
        for config in important_configs:
            excitations = subspace.generate_excitations(config)
            subspace.add_configs(excitations)

        # Prune
        subspace.prune(amplitudes)

        # Optimize sampling
        self._optimize_sampling(important_configs)

        # Check convergence
        if abs(energy - prev_energy) < self.conv_threshold:
            break
        prev_energy = energy

    return {'energy': energy, 'iterations': iteration+1}
```

---

## Phase 3: PRODUCTION READY (Week 2)

### 3.1 Add Z2 Symmetry Tapering
- Reduce H2O from 12 → 10 qubits using symmetries

### 3.2 Add Noise Mitigation
- Zero-noise extrapolation
- Readout error mitigation
- For IBM/BlueQubit backends

### 3.3 Circuit Optimization
- Use hardware-native gates
- Minimize circuit depth
- Reduce gate count

### 3.4 Benchmarking
- Compare vs Hi-VQE paper results
- Test on: H2, LiH, H2O, NH3, CH4
- Verify chemical accuracy

---

## Success Metrics

### H2O Performance:
- **Qubits:** 12 (down from 14) ✅
- **Iterations:** 5-10 (down from 20+)
- **Function evals:** 10-20 (down from 1800)
- **Accuracy:** < 1 kcal/mol error
- **Time:** < 1 minute (down from timeout)

### Cloud Cost (IBM):
- **Before:** 1800 jobs × $0.01 = $18 per H2O experiment
- **After:** 20 jobs × $0.01 = $0.20 per H2O experiment
- **Savings:** 90x reduction

### Scalability:
- **H2O (10e):** 5-10 iterations, 12 qubits ✅
- **NH3 (10e):** 5-10 iterations, 12 qubits
- **CH4 (10e):** 5-10 iterations, 12 qubits
- **Benzene (30e):** 10-15 iterations, 18-20 qubits (with active space)

---

## Implementation Order:

**Today (2 hours):**
1. SPSA optimizer (30 min)
2. Active space for H2O (30 min)
3. Backend-specific optimizer selection (15 min)
4. Testing (45 min)

**Tomorrow (4 hours):**
1. Configuration sampling (1 hour)
2. Subspace management (1 hour)
3. Classical diagonalization (1 hour)
4. VQE integration (1 hour)

**This Week:**
- Validate on H2, LiH, H2O
- Fix any bugs
- Optimize performance

**Next Week:**
- Symmetry reduction
- Noise mitigation
- Production deployment

---

**End of Plan**
