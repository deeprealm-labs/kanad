# Governance Implementation Plan - Making Governance ACTIVE

**Date:** October 7, 2025
**Objective:** Make governance protocols the ACTIVE physics engine in Kanad (the core innovation)

---

## Executive Summary

**Current State:** Governance is metadata/validation layer
**Target State:** Governance actively controls Hamiltonian, ansatz, and solver behavior
**Philosophy:** "Physics of bonding dictates quantum representation"

---

## Part 1: Architecture Analysis

### What Governance SHOULD Control:

```
┌─────────────────────────────────────────────────────────────┐
│                    GOVERNANCE LAYER                          │
│  ┌────────────┐  ┌────────────┐  ┌────────────┐            │
│  │  Ionic     │  │ Covalent   │  │ Metallic   │            │
│  │  Protocol  │  │ Protocol   │  │ Protocol   │            │
│  └────────────┘  └────────────┘  └────────────┘            │
│         ↓                ↓                ↓                 │
└─────────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────────┐
│              REPRESENTATION SELECTION                        │
│  Ionic: Localized AOs  │ Covalent: MOs  │ Metallic: Bloch  │
└─────────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────────┐
│            HAMILTONIAN CONSTRUCTION                          │
│  Protocol.construct_hamiltonian(atoms, basis)                │
│  - Selects appropriate operators                             │
│  - Builds basis functions per governance                     │
│  - Applies constraints                                       │
└─────────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────────┐
│               ANSATZ CONSTRUCTION                            │
│  Protocol.construct_ansatz(n_qubits, n_electrons)            │
│  - Enforces entanglement topology                            │
│  - Applies gate ordering rules                               │
│  - Preserves required symmetries                             │
└─────────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────────┐
│                 SOLVER BEHAVIOR                              │
│  - Initial guess from protocol                               │
│  - Convergence criteria                                      │
│  - Observable measurement                                    │
└─────────────────────────────────────────────────────────────┘
```

### Current Implementation Status:

| Component | Current | Target | Priority |
|-----------|---------|--------|----------|
| Protocol Attachment | ✅ Works | ✅ Keep | - |
| Rule Definitions | ✅ Works | ✅ Keep | - |
| Ansatz Construction | ⚠️ Partial | ✅ Fix | **HIGH** |
| Hamiltonian Construction | ❌ Not Used | ✅ Implement | **HIGH** |
| Representation Selection | ❌ Not Used | ✅ Implement | **MEDIUM** |
| SCF Integration | ❌ Not Used | ⚠️ Optional | LOW |

---

## Part 2: Implementation Plan

### Phase 1: Fix Governance Ansatze (HIGH PRIORITY) ⚡

**Goal:** Make governance ansatze work with VQE solver

#### Step 1.1: Fix CovalentGovernanceAnsatz

**File:** `kanad/ansatze/governance_aware_ansatz.py`

**Issues to Fix:**
1. Add `n_parameters` property
2. Make hamiltonian parameter optional
3. Ensure `build_circuit()` returns proper circuit object

**Implementation:**

```python
class CovalentGovernanceAnsatz(BaseAnsatz):
    def __init__(
        self,
        n_qubits: int,
        n_electrons: int,
        n_layers: int = 2,
        hybridization: str = 'sp3',
        protocol: Optional[CovalentGovernanceProtocol] = None,
        hamiltonian: Optional[Any] = None  # Make optional
    ):
        super().__init__(n_qubits, n_electrons)

        # Initialize protocol if not provided
        if protocol is None:
            from kanad.governance.protocols import CovalentGovernanceProtocol
            protocol = CovalentGovernanceProtocol()

        self.protocol = protocol
        self.n_layers = n_layers
        self.hybridization = hybridization
        self.hamiltonian = hamiltonian  # Store if provided

        # Build circuit using governance
        self._circuit_state = None
        self._parameter_count = 0

    @property
    def n_parameters(self) -> int:
        """Return number of variational parameters."""
        if self._parameter_count == 0 and self._circuit_state is None:
            # Build circuit to count parameters
            self.build_circuit()
        return self._parameter_count

    def build_circuit(self):
        """Build governance-aware circuit."""
        from kanad.governance.protocols.base_protocol import QuantumCircuitState

        # Create circuit state
        circuit_state = QuantumCircuitState(self.n_qubits)

        # Apply governance rules
        context = {
            'n_electrons': self.n_electrons,
            'n_layers': self.n_layers,
            'hybridization': self.hybridization
        }

        # Use protocol to construct ansatz
        governed_state = self.protocol.construct_ansatz(
            representation={'n_qubits': self.n_qubits, 'n_electrons': self.n_electrons}
        )

        # Count parameters from circuit
        self._parameter_count = sum(
            len(gate['params']) for gate in governed_state.gates
        )

        # Convert to actual quantum circuit
        self.circuit = self._convert_to_circuit(governed_state)
        self._circuit_state = governed_state

        return self.circuit
```

#### Step 1.2: Test Governance Ansatz

**Create test:** `tests/validation/test_governance_ansatz_active.py`

```python
def test_covalent_ansatz_uses_governance():
    """Verify covalent ansatz applies governance rules."""
    from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz
    from kanad.governance.protocols import CovalentGovernanceProtocol

    # Create ansatz
    ansatz = CovalentGovernanceAnsatz(n_qubits=4, n_electrons=2)

    # Build circuit
    circuit = ansatz.build_circuit()

    # Check governance was applied
    assert ansatz._circuit_state is not None, "Circuit state not created"
    assert ansatz._circuit_state.is_hybridized, "Hybridization not applied"
    assert ansatz._circuit_state.has_mo_pairs, "MO pairs not formed"
    assert ansatz._circuit_state.is_paired, "Electron pairing not applied"

    # Check circuit structure follows covalent rules
    gates = ansatz._circuit_state.gates

    # Rule 1: Hybridization comes first
    first_gates = gates[:4]  # First n_qubits gates
    assert all(g['type'] in ['ry', 'rx', 'rz'] for g in first_gates), \
        "First gates should be single-qubit rotations (hybridization)"

    # Rule 2: Has 2-qubit gates for pairing
    two_qubit_gates = [g for g in gates if len(g['qubits']) == 2]
    assert len(two_qubit_gates) > 0, "No 2-qubit gates (pairing missing)"

    # Rule 3: Sparse entanglement (covalent is local)
    assert ansatz._circuit_state.is_sparse(), \
        "Entanglement not sparse (violates covalent locality)"

    print("✅ Covalent governance ACTIVELY applied to ansatz!")
```

---

### Phase 2: Implement Governance in Hamiltonian Construction (HIGH PRIORITY) 🔥

**Goal:** Make Hamiltonian construction use governance protocols

#### Step 2.1: Add Protocol-Driven Hamiltonian Construction

**File:** `kanad/core/hamiltonians/covalent_hamiltonian.py`

**New Method:**

```python
def _construct_with_governance(self):
    """Construct Hamiltonian using governance protocol guidance."""
    if not self.use_governance or not self.governance_protocol:
        # Fall back to standard construction
        return self._construct_standard()

    logger.info("🔥 Constructing Hamiltonian with ACTIVE governance")

    # Step 1: Let protocol choose representation
    representation_type = self.governance_protocol.get_representation_type()
    logger.info(f"  Governance selected: {representation_type} representation")

    # Step 2: Let protocol guide basis function construction
    if representation_type == 'molecular_orbital':
        # Covalent bonding uses MO basis
        self._construct_mo_basis()
    elif representation_type == 'atomic_orbital':
        # Ionic bonding uses AO basis
        self._construct_ao_basis()
    elif representation_type == 'bloch':
        # Metallic bonding uses Bloch basis
        self._construct_bloch_basis()

    # Step 3: Let protocol filter allowed operators
    allowed_ops = self.governance_protocol.get_allowed_operators()
    forbidden_ops = self.governance_protocol.get_forbidden_operators()

    logger.info(f"  Allowed operators: {allowed_ops}")
    logger.info(f"  Forbidden operators: {forbidden_ops}")

    # Step 4: Construct matrices with governance constraints
    self._build_h_core_governed(allowed_ops)
    self._build_eri_governed(allowed_ops)

    return self
```

**Add to CovalentGovernanceProtocol:**

```python
def get_representation_type(self) -> str:
    """Get recommended representation for covalent bonding."""
    return 'molecular_orbital'  # Covalent uses MO basis

def get_allowed_operators(self) -> List[str]:
    """Operators allowed in covalent Hamiltonian."""
    return [
        'kinetic',           # Kinetic energy
        'nuclear_attraction', # Electron-nuclear attraction
        'coulomb',           # Electron-electron repulsion
        'exchange',          # Exchange interaction
        'hybridization',     # Orbital mixing
    ]

def get_forbidden_operators(self) -> List[str]:
    """Operators forbidden in covalent Hamiltonian."""
    return [
        'long_range_coulomb',  # No long-range in localized covalent
        'collective_exchange', # No collective effects
    ]
```

#### Step 2.2: Test Governance-Driven Hamiltonian

**Test:**

```python
def test_hamiltonian_uses_governance():
    """Verify Hamiltonian construction respects governance."""
    from kanad.bonds import BondFactory

    # Create H2 bond
    h2 = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
    ham = h2.hamiltonian

    # Verify governance is active
    assert ham.use_governance, "Governance not enabled"
    assert ham.governance_protocol is not None, "No protocol attached"

    # Check that representation was selected by protocol
    protocol = ham.governance_protocol
    rep_type = protocol.get_representation_type()
    assert rep_type == 'molecular_orbital', \
        f"Expected MO representation for covalent, got {rep_type}"

    # Check matrix construction followed governance
    allowed_ops = protocol.get_allowed_operators()
    forbidden_ops = protocol.get_forbidden_operators()

    # Verify h_core only contains allowed terms
    # (Would need to track which operators contributed to h_core)

    print("✅ Hamiltonian construction ACTIVELY used governance!")
```

---

### Phase 3: Implement Representation Selection (MEDIUM PRIORITY)

**Goal:** Let governance choose quantum representation

#### Step 3.1: Create Governance-Aware Representations

**New File:** `kanad/core/representations/governance_representations.py`

```python
class GovernanceAwareRepresentation(BaseRepresentation):
    """Base class for governance-driven representations."""

    def __init__(self, protocol: BaseGovernanceProtocol):
        self.protocol = protocol
        super().__init__()

    @abstractmethod
    def construct_from_governance(self, atoms, basis):
        """Construct representation using governance rules."""
        pass

class MolecularOrbitalRepresentation(GovernanceAwareRepresentation):
    """MO representation for covalent bonding."""

    def construct_from_governance(self, atoms, basis):
        """Build MO basis using covalent governance."""
        # 1. Form atomic orbitals
        aos = self._build_atomic_orbitals(atoms, basis)

        # 2. Apply hybridization (from governance)
        hybrid_aos = self.protocol.apply_hybridization(aos)

        # 3. Form molecular orbitals (bonding/antibonding pairs)
        mos = self.protocol.form_molecular_orbitals(hybrid_aos)

        return mos

class AtomicOrbitalRepresentation(GovernanceAwareRepresentation):
    """AO representation for ionic bonding."""

    def construct_from_governance(self, atoms, basis):
        """Build localized AO basis using ionic governance."""
        # Ionic bonding keeps orbitals localized
        aos = self._build_atomic_orbitals(atoms, basis)

        # No hybridization (ionic keeps s, p separate)
        # No delocalization (electrons stay on atoms)

        return aos
```

---

### Phase 4: Full Integration Test (CRITICAL) ⭐

**Goal:** Test complete governance workflow end-to-end

#### Step 4.1: Comprehensive Governance Workflow Test

**Test File:** `tests/validation/test_full_governance_workflow.py`

```python
def test_h2_with_active_governance():
    """
    Test H2 with ACTIVE governance at every step.

    Workflow:
    1. Create bond → Covalent protocol selected
    2. Hamiltonian → MO representation chosen
    3. Ansatz → Covalent rules applied
    4. VQE → Governance-constrained optimization
    """
    from kanad.bonds import BondFactory
    from kanad.solvers.vqe_solver import VQESolver

    print("\n" + "="*80)
    print("FULL GOVERNANCE WORKFLOW TEST")
    print("="*80)

    # Step 1: Create bond with governance
    print("\n1. Creating H2 bond...")
    h2 = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')

    assert h2.hamiltonian.use_governance, "Governance not enabled!"
    protocol = h2.hamiltonian.governance_protocol
    print(f"   ✅ Protocol: {type(protocol).__name__}")
    print(f"   ✅ Rules: {len(protocol.rules)}")

    # Step 2: Check Hamiltonian used governance
    print("\n2. Checking Hamiltonian construction...")
    rep_type = protocol.get_representation_type()
    print(f"   ✅ Representation: {rep_type}")

    # Step 3: Create governance ansatz
    print("\n3. Creating governance ansatz...")
    solver = VQESolver(
        bond=h2,
        ansatz_type='governance',  # Use governance ansatz
        mapper_type='jordan_wigner',
        max_iterations=50,
        backend='statevector'
    )

    # Verify ansatz used governance
    ansatz = solver.ansatz
    print(f"   ✅ Ansatz type: {type(ansatz).__name__}")

    if hasattr(ansatz, '_circuit_state'):
        state = ansatz._circuit_state
        print(f"   ✅ Hybridized: {state.is_hybridized}")
        print(f"   ✅ MO pairs: {state.has_mo_pairs}")
        print(f"   ✅ Paired: {state.is_paired}")
        print(f"   ✅ Sparse: {state.is_sparse()}")

    # Step 4: Run VQE with governance
    print("\n4. Running VQE with governance...")
    result = solver.solve()

    energy_gov = result['energy']
    print(f"   ✅ Energy (governance): {energy_gov:.6f} Ha")

    # Step 5: Compare with non-governance
    print("\n5. Comparing with standard VQE (no governance)...")
    solver_std = VQESolver(
        bond=h2,
        ansatz_type='hardware_efficient',  # Standard ansatz
        mapper_type='jordan_wigner',
        max_iterations=50,
        backend='statevector'
    )
    result_std = solver_std.solve()
    energy_std = result_std['energy']
    print(f"   Standard energy: {energy_std:.6f} Ha")
    print(f"   Governance energy: {energy_gov:.6f} Ha")
    print(f"   Difference: {abs(energy_gov - energy_std):.6f} Ha")

    # Step 6: Verify governance improved or maintained accuracy
    _, scf_energy = h2.hamiltonian.solve_scf()
    print(f"\n6. Comparing with SCF...")
    print(f"   SCF energy: {scf_energy:.6f} Ha")
    print(f"   Governance error: {abs(energy_gov - scf_energy):.6f} Ha")
    print(f"   Standard error: {abs(energy_std - scf_energy):.6f} Ha")

    # Governance should be at least as good as standard
    assert abs(energy_gov - scf_energy) <= abs(energy_std - scf_energy) * 1.1, \
        "Governance ansatz worse than standard!"

    print("\n" + "="*80)
    print("✅ FULL GOVERNANCE WORKFLOW SUCCESSFUL!")
    print("="*80)
```

---

## Part 3: Expected Outcomes

### What Governance WILL Change:

1. **Ansatz Structure:**
   - **Covalent:** Hybridization gates → Pairing gates → Correlation
   - **Ionic:** Local gates only, minimal entanglement
   - **Metallic:** QFT-like collective gates, maximal entanglement

2. **Circuit Topology:**
   - **Covalent:** Paired connectivity (qubits i, i+1 for bonding orbitals)
   - **Ionic:** Nearest-neighbor only
   - **Metallic:** All-to-all connectivity

3. **Parameter Count:**
   - **Covalent:** Moderate (paired excitations)
   - **Ionic:** Few (local rotations only)
   - **Metallic:** Many (collective transformations)

4. **Convergence:**
   - Governance-constrained parameter space is smaller
   - May converge faster (fewer bad local minima)
   - May be more physically meaningful

---

## Part 4: Implementation Timeline

### Week 1: Core Fixes
- [ ] Day 1-2: Fix CovalentGovernanceAnsatz integration
- [ ] Day 3-4: Fix IonicGovernanceAnsatz integration
- [ ] Day 5: Test both ansatze with VQE

### Week 2: Hamiltonian Integration
- [ ] Day 1-2: Implement `_construct_with_governance()`
- [ ] Day 3: Add representation selection to protocols
- [ ] Day 4-5: Test governance-driven Hamiltonian construction

### Week 3: Full Integration
- [ ] Day 1-2: Create governance-aware representations
- [ ] Day 3-4: Full workflow tests
- [ ] Day 5: Performance comparison

### Week 4: Validation & Documentation
- [ ] Day 1-2: Test multiple molecules (H2, H2O, LiH, etc.)
- [ ] Day 3: Benchmark governance vs standard
- [ ] Day 4-5: Documentation and examples

---

## Part 5: Success Metrics

### Quantitative Metrics:

1. **Ansatz Verification:**
   - ✅ `ansatz._circuit_state.is_hybridized == True` for covalent
   - ✅ `ansatz._circuit_state.is_sparse() == True` for ionic
   - ✅ `ansatz._circuit_state.is_collectively_entangled == True` for metallic

2. **Accuracy:**
   - Governance VQE error ≤ Standard VQE error
   - Both within 1% of SCF energy

3. **Efficiency:**
   - Governance uses ≤ parameters than standard (more constrained)
   - May require fewer iterations (smaller search space)

### Qualitative Metrics:

1. **Physical Meaning:**
   - Circuit structure reflects bonding physics
   - Gates correspond to physical processes (hybridization, pairing, etc.)

2. **Interpretability:**
   - Can explain why certain gates are present
   - Can validate circuit against chemistry knowledge

3. **Differentiation:**
   - Clearly different from Qiskit Nature (no governance there)
   - Novel approach backed by physics

---

## Part 6: Risk Mitigation

### Risk 1: Governance Constraints Too Restrictive

**Symptom:** Governance VQE significantly worse than standard VQE

**Mitigation:**
- Add "relaxation" mode to governance
- Allow optional rules vs required rules
- Provide override mechanism for specific cases

### Risk 2: Implementation Complexity

**Symptom:** Hard to integrate governance everywhere

**Mitigation:**
- Start with ansatz only (highest impact)
- Add Hamiltonian governance incrementally
- Keep standard path as fallback

### Risk 3: Performance Overhead

**Symptom:** Governance makes code slower

**Mitigation:**
- Governance rules run once during construction, not per iteration
- Cache governance decisions
- Profile and optimize hot paths

---

## Conclusion

This plan makes governance the **ACTIVE CORE** of Kanad, transforming it from metadata to the physics engine that drives:
- Representation selection
- Hamiltonian construction
- Ansatz structure
- Circuit topology

**This is THE differentiation from Qiskit Nature** - physics-driven quantum chemistry, not just generic quantum algorithms applied to molecules.

**Start with Phase 1** (fix ansatze) - highest impact, lowest risk, fastest path to demonstrating active governance.
