# Governance + Hi-VQE Integration Plan

**Date:** November 4, 2025
**Purpose:** Rebuild KANAD to combine governance protocols with Hi-VQE efficiency

---

## Governance Module Architecture (What We Have)

### Core Concept:
**Governance = Physics-aware constraints on quantum circuits**

Different bond types have different quantum mechanical behaviors:
- **Covalent:** Orbital hybridization (sp³), paired electron sharing, molecular orbitals
- **Ionic:** Charge transfer, localized electrons, minimal entanglement
- **Metallic:** Delocalized electrons, collective states, band structure

### Key Components:

1. **Governance Protocols** (`kanad/governance/protocols/`)
   - `covalent_protocol.py` - Enforces hybridization, MO formation, paired entanglement
   - `ionic_protocol.py` - Enforces charge localization, transfer operators
   - `metallic_protocol.py` - Enforces delocalization, collective states

2. **Governance-Aware Hamiltonians** (`kanad/core/hamiltonians/`)
   - `covalent_hamiltonian.py` - Builds H with hybridization physics
   - `ionic_hamiltonian.py` - Builds H with charge transfer terms
   - `metallic_hamiltonian.py` - Builds H with band structure

3. **Governance-Aware Ansatze** (`kanad/ansatze/`)
   - `governance_aware_ansatz.py` - Circuit respects bonding physics
   - Uses protocols to validate/construct physically meaningful gates

### How It Works:

```python
# Example: H2O (covalent bonding)
molecule = Molecule.from_smiles("O")
protocol = CovalentGovernanceProtocol()  # Knows about sp³ hybridization

# Protocol enforces:
# 1. Hybridization transformations first (sp³ for O)
# 2. Bonding/antibonding MO pairs
# 3. Electron pair entanglement (Bell states)
# 4. NO long-range entanglement (only bonding pairs)

hamiltonian = CovalentHamiltonian(molecule, use_governance=True)
ansatz = CovalentGovernanceAnsatz(n_qubits=12, n_electrons=10, protocol=protocol)
```

**Key Insight:**
Governance protocols KNOW THE PHYSICS - they can guide:
- Which orbitals to include in active space
- How to construct physically meaningful configurations
- Which excitations are important for bonding

---

## Hi-VQE Architecture (What We Need)

### Core Concept:
**Don't measure energy - sample configurations and solve exactly in subspace**

```
Standard VQE:
  Adjust params → Measure <ψ|H|ψ> (1000s of Pauli terms) → Update → Repeat
  Problem: Expensive measurements, gets stuck

Hi-VQE:
  Adjust params → Sample bitstrings (1 Z measurement) → Classical solve → Update → Repeat
  Advantage: 1 measurement, exact energy, faster convergence
```

### Key Components Needed:

1. **Configuration Sampler**
   - Measure quantum state in computational basis
   - Filter configs (correct electron count)
   - Track sampled configurations

2. **Subspace Manager**
   - Store configuration subspace
   - Generate excitations from important configs
   - Prune low-amplitude configs

3. **Classical Diagonalizer**
   - Project Hamiltonian into config subspace
   - Exact eigensolve (no quantum measurements!)
   - Extract amplitudes

---

## THE KEY INSIGHT: Governance + Hi-VQE Synergy

**Governance protocols can GUIDE Hi-VQE!**

### Problem Hi-VQE Faces:
- Which configurations should we sample?
- Which excitations should we generate?
- How to avoid exponential growth of subspace?

### How Governance Helps:

1. **Physics-Guided Configuration Selection**
   ```python
   # Covalent protocol KNOWS:
   # - Bonding configs have paired electrons
   # - Important configs respect hybridization
   # - Single/double excitations between bonding pairs matter

   protocol = CovalentGovernanceProtocol()

   # Generate only PHYSICALLY MEANINGFUL excitations
   def generate_covalent_excitations(config, protocol):
       excitations = []
       # Only excitations that preserve:
       # - Hybridization character
       # - Bonding/antibonding structure
       # - Spin pairing
       bonding_pairs = protocol.get_bonding_pairs()
       for (i, j) in bonding_pairs:
           # Single excitation within bond
           excitations.append(single_excitation(config, i, j))
       return excitations
   ```

2. **Smart Active Space from Governance**
   ```python
   # Covalent protocol knows H2O has:
   # - O: 1s (core, freeze!), 2s, 2px, 2py, 2pz (valence, active!)
   # - H: 1s (valence, active!)

   def get_governance_active_space(molecule, protocol):
       if protocol.bond_type == 'covalent':
           # Freeze core, keep valence orbitals involved in bonding
           core = protocol.identify_core_orbitals(molecule)
           bonding = protocol.identify_bonding_orbitals(molecule)
           return core, bonding  # H2O: freeze [0], active [1,2,3,4,5,6]
   ```

3. **Efficient Subspace Growth**
   ```python
   # Instead of generating ALL single/double excitations (exponential!),
   # generate only excitations that respect bonding physics

   for config in subspace:
       if protocol.is_bonding_config(config):
           # Generate excitations respecting MO character
           excitations = protocol.generate_physical_excitations(config)
           # Much smaller set than all possible excitations!
           subspace.add(excitations)
   ```

---

## Integration Architecture

### New Component: `GovernanceHiVQESolver`

```python
class GovernanceHiVQESolver(VQESolver):
    """
    Hi-VQE solver with governance protocol integration.

    Combines:
    - Hi-VQE efficiency (config sampling + classical diag)
    - Governance physics (smart active space, guided excitations)
    """

    def __init__(self, molecule, bond_type='auto', ...):
        # Auto-detect bond type or use specified
        if bond_type == 'auto':
            self.protocol = detect_bond_protocol(molecule)
        else:
            self.protocol = get_protocol(bond_type)

        # Get physics-guided active space
        frozen, active = self.protocol.get_active_space(molecule)

        # Build Hamiltonian with governance
        self.hamiltonian = build_governed_hamiltonian(
            molecule, protocol=self.protocol,
            frozen_core=frozen, active_space=active
        )

        # Build governance-aware ansatz
        self.ansatz = build_governed_ansatz(
            n_qubits=len(active)*2,
            n_electrons=molecule.n_electrons - len(frozen)*2,
            protocol=self.protocol
        )

    def solve(self):
        """Hi-VQE with governance guidance."""
        # Initialize with HF configuration respecting governance
        hf_config = self.protocol.get_hf_configuration()
        subspace = ConfigurationSubspace(hf_config)

        for iteration in range(self.max_iterations):
            # 1. Sample configurations
            configs = self.sample_configurations(self.params, shots=1000)

            # 2. Filter using governance (e.g., spin symmetry, charge conservation)
            valid_configs = self.protocol.filter_valid_configs(configs)
            subspace.add_configs(valid_configs)

            # 3. Classical diagonalization
            H_sub = project_hamiltonian(self.hamiltonian, subspace)
            energy, amplitudes = diagonalize(H_sub)

            # 4. GOVERNANCE-GUIDED excitation generation
            important_configs = get_high_amplitude_configs(amplitudes)
            for config in important_configs:
                # Use protocol to generate PHYSICAL excitations only!
                excitations = self.protocol.generate_excitations(config)
                subspace.add_configs(excitations)

            # 5. Prune (can use governance to identify unphysical configs)
            subspace = self.protocol.prune_subspace(subspace, amplitudes)

            # 6. Optimize parameters
            # ... update params to sample important configs better

            if converged:
                break

        return energy, subspace
```

---

## Implementation Plan

### Phase 1: Governance-Aware Active Space

**File:** `kanad/core/active_space.py` (NEW)

```python
def get_governance_active_space(molecule, protocol):
    """
    Use governance protocol to select active orbitals.

    For H2O with CovalentGovernanceProtocol:
        - Identifies O 1s as core (freeze)
        - Identifies O 2s, 2p and H 1s as valence (active)
        - Returns: frozen=[0], active=[1,2,3,4,5,6]
        - Result: 12 qubits instead of 14
    """
    if protocol.bond_type == 'covalent':
        # Covalent: freeze core, active = valence orbitals
        core_orbitals = []
        for atom in molecule.atoms:
            if atom.symbol == 'O':
                core_orbitals.append(0)  # O 1s

        all_orbitals = list(range(molecule.n_orbitals))
        active_orbitals = [i for i in all_orbitals if i not in core_orbitals]

        return core_orbitals, active_orbitals

    elif protocol.bond_type == 'ionic':
        # Ionic: different active space strategy
        # Keep orbitals involved in charge transfer
        pass

    elif protocol.bond_type == 'metallic':
        # Metallic: need conduction band orbitals
        pass
```

**Test:**
```python
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol

protocol = CovalentGovernanceProtocol()
frozen, active = get_governance_active_space(h2o_molecule, protocol)

assert len(active) == 6  # H2O: 6 active orbitals
assert len(frozen) == 1  # O 1s frozen
# Result: 12 qubits instead of 14!
```

### Phase 2: Physics-Guided Configuration Sampling

**File:** `kanad/core/configuration.py` (NEW)

```python
class GovernanceConfigurationManager:
    """Manage configurations with governance protocol guidance."""

    def __init__(self, protocol, n_qubits, n_electrons):
        self.protocol = protocol
        self.subspace = ConfigurationSubspace(n_qubits, n_electrons)

    def filter_configs(self, sampled_configs):
        """Filter using governance rules."""
        valid = []
        for config in sampled_configs:
            # Basic: correct electron count
            if config.count('1') != self.n_electrons:
                continue

            # Governance: spin symmetry, charge conservation, etc.
            if self.protocol.is_valid_configuration(config):
                valid.append(config)

        return valid

    def generate_excitations(self, config):
        """Generate excitations respecting physics."""
        # Use protocol to generate only meaningful excitations
        if self.protocol.bond_type == 'covalent':
            # Only excitations within/between bonding pairs
            return self.protocol.generate_covalent_excitations(config)
        else:
            # Other bond types have different rules
            return self.protocol.generate_excitations(config)
```

### Phase 3: Classical Diagonalization with Governance

**File:** `kanad/core/classical_solver.py` (NEW)

```python
def project_hamiltonian_with_governance(hamiltonian, subspace, protocol):
    """
    Project Hamiltonian with governance awareness.

    Can use protocol to:
    - Skip obviously zero matrix elements (symmetry)
    - Use physical approximations (e.g., σ/π separation)
    - Cache important matrix elements
    """
    n = len(subspace.configs)
    H_sub = np.zeros((n, n), dtype=complex)

    for i, config_i in enumerate(subspace.configs):
        for j, config_j in enumerate(subspace.configs):
            # Use governance to check if element is zero by symmetry
            if protocol.matrix_element_is_zero(config_i, config_j):
                H_sub[i, j] = 0.0
                continue

            # Compute matrix element
            H_sub[i, j] = compute_matrix_element(hamiltonian, config_i, config_j)

    return H_sub
```

### Phase 4: Full Integration

**Modify:** `kanad/utils/vqe_solver.py`

```python
class VQESolver:
    def solve(self, mode='standard'):
        if mode == 'governance_hivqe':
            return self._solve_governance_hivqe()
        else:
            return self._solve_standard()

    def _solve_governance_hivqe(self):
        """Governance-guided Hi-VQE."""
        # Get protocol
        protocol = self.get_governance_protocol()

        # Initialize with governance-aware HF config
        hf_config = protocol.get_hf_configuration(self.n_qubits, self.n_electrons)
        subspace = GovernanceConfigurationManager(protocol, self.n_qubits, self.n_electrons)
        subspace.add_configs([hf_config])

        for iteration in range(self.max_iterations):
            # Sample
            configs = self.sample_configurations(self.params, shots=1000)
            valid = subspace.filter_configs(configs)
            subspace.add_configs(valid)

            # Classical solve
            H_sub = project_hamiltonian_with_governance(
                self.hamiltonian, subspace, protocol
            )
            energy, amplitudes = np.linalg.eigh(H_sub)
            energy = energy[0]
            amplitudes = amplitudes[:, 0]

            # Governance-guided expansion
            important = get_high_amplitude_configs(subspace.configs, amplitudes)
            for config in important:
                excitations = subspace.generate_excitations(config)
                subspace.add_configs(excitations)

            # Prune with governance
            subspace = protocol.prune_subspace(subspace, amplitudes)

            if converged:
                break

        return energy
```

---

## Expected Benefits

### Compared to Standard VQE:
- ✅ 180,000x fewer measurements (10 vs 1.8M for H2O)
- ✅ Exact energy in subspace (no measurement noise)
- ✅ Faster convergence (5-10 iterations)

### Compared to Generic Hi-VQE:
- ✅ Physics-guided active space (12 vs 14 qubits for H2O)
- ✅ Smarter excitation generation (100s vs 1000s of configs)
- ✅ Faster subspace growth (use bonding physics)
- ✅ Better initial guess (governance-aware HF)

### Unique KANAD Advantage:
- ✅ **Governance protocols encode domain knowledge**
- ✅ **Automatic adaptation to bond type**
- ✅ **Physics-constrained = more efficient**
- ✅ **Publishable novel approach**

---

## Testing Strategy

### Test 1: Governance Active Space
```python
# H2O should use 12 qubits with covalent protocol
protocol = CovalentGovernanceProtocol()
frozen, active = get_governance_active_space(h2o, protocol)
assert len(active) * 2 == 12  # 6 active orbitals = 12 qubits
```

### Test 2: Configuration Filtering
```python
# Invalid configs should be rejected
config_manager = GovernanceConfigurationManager(protocol, 12, 10)
sampled = ['110000000000', '010000000001', ...]  # Some invalid
valid = config_manager.filter_configs(sampled)
assert all(c.count('1') == 10 for c in valid)  # Correct electron count
assert all(protocol.is_valid_configuration(c) for c in valid)  # Physics valid
```

### Test 3: Excitation Generation
```python
# Should generate fewer, more meaningful excitations
config = '110000000000'  # HF config
all_excitations = generate_all_single_excitations(config)  # ~60 excitations
gov_excitations = protocol.generate_covalent_excitations(config)  # ~10-20 excitations
assert len(gov_excitations) < len(all_excitations)
assert all(protocol.is_physical_excitation(e) for e in gov_excitations)
```

### Test 4: Full Governance Hi-VQE
```python
# H2O with governance Hi-VQE
solver = VQESolver(molecule=h2o, bond_type='covalent')
result = solver.solve(mode='governance_hivqe', max_iterations=10)

# Should achieve chemical accuracy
assert abs(result['energy'] - h2o.fci_energy) < 0.0016  # < 1 kcal/mol
assert result['iterations'] <= 10
assert result['n_qubits'] == 12  # Not 14!
assert result['subspace_size'] < 1000  # Manageable subspace
```

---

## Summary

**Governance Module:** Physics-aware constraints for different bond types
**Hi-VQE:** Efficient VQE via configuration sampling + classical solve
**Integration:** Use governance to GUIDE Hi-VQE configuration selection

**Result:** Best of both worlds:
- Hi-VQE efficiency (10 measurements vs millions)
- Governance physics (smart active space, guided excitations)
- KANAD uniqueness (automated bond-type adaptation)

**This is the path to production-ready quantum chemistry for large molecules!**

---

**End of Integration Plan**
