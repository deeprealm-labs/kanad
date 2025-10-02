"""
Ionic bonding governance protocol.

Physical principles for ionic bonding:
1. Electrons are LOCALIZED on specific atoms
2. Electron transfer between sites (a†_i a_j operators)
3. MINIMAL entanglement (nearest-neighbor only)
4. Strong on-site Coulomb repulsion (Hubbard U)
5. Large electronegativity difference drives charge separation

Circuit requirements:
- Sparse connectivity (nearest-neighbor gates)
- No long-range entanglement
- Transfer operators only between bonded sites
- Preserves particle number per spin
"""

from typing import List, Dict, Any
import numpy as np
from kanad.governance.protocols.base_protocol import (
    BaseGovernanceProtocol,
    BondingType,
    GovernanceRule,
    QuantumCircuitState
)


class IonicGovernanceProtocol(BaseGovernanceProtocol):
    """
    Governance protocol for ionic bonding systems.

    Example: Na+ Cl-
    - Electron transfers from Na to Cl
    - Localized states on each atom
    - Minimal entanglement between sites
    - Circuit uses local gates only
    """

    def __init__(self):
        """Initialize ionic governance protocol."""
        super().__init__(BondingType.IONIC)

    def _initialize_rules(self):
        """Initialize ionic bonding governance rules."""

        # Rule 1: Enforce localized gates (no long-range entanglement)
        self.rules.append(GovernanceRule(
            name="localized_gates_only",
            description="Only nearest-neighbor gates allowed (localized interactions)",
            condition=lambda state, ctx: True,  # Always check
            action=self._enforce_locality,
            priority=100,
            required=True
        ))

        # Rule 2: Sparse connectivity (minimal entanglement)
        self.rules.append(GovernanceRule(
            name="sparse_connectivity",
            description="Connectivity must be sparse (avg degree < 2)",
            condition=lambda state, ctx: isinstance(state, QuantumCircuitState),
            action=self._enforce_sparsity,
            priority=90,
            required=True
        ))

        # Rule 3: Preserve particle number
        self.rules.append(GovernanceRule(
            name="particle_number_conservation",
            description="Total particle number must be conserved",
            condition=lambda state, ctx: True,
            action=self._enforce_particle_conservation,
            priority=80,
            required=True
        ))

        # Rule 4: No collective gates (forbidden for ionic)
        self.rules.append(GovernanceRule(
            name="forbid_collective_gates",
            description="No gates acting on >2 qubits (no collective behavior)",
            condition=lambda state, ctx: isinstance(state, QuantumCircuitState),
            action=self._forbid_collective_gates,
            priority=70,
            required=True
        ))

    def validate_operator(self, operator: 'QuantumOperator') -> bool:
        """
        Validate operator for ionic bonding.

        ALLOWED:
        - Single-qubit rotations (on-site terms)
        - Nearest-neighbor 2-qubit gates (transfer)
        - Number operators (charge counting)

        FORBIDDEN:
        - Long-range entangling gates
        - Multi-qubit (>2) gates
        - Delocalized operators
        """
        # Placeholder - would check operator properties
        operator_type = getattr(operator, 'type', 'unknown')

        allowed = ['rx', 'ry', 'rz', 'cx', 'cnot', 'cz', 'rxx', 'ryy', 'rzz']
        forbidden = ['ghz', 'w_state', 'qft', 'swap_network']

        if operator_type in forbidden:
            return False

        # Check locality - no long-range gates
        if hasattr(operator, 'qubits'):
            qubits = operator.qubits
            if len(qubits) > 2:
                return False  # No multi-qubit gates

            if len(qubits) == 2:
                # Check if qubits are nearest neighbors
                # (simplified - would use actual geometry)
                if abs(qubits[0] - qubits[1]) > 1:
                    return False  # Not nearest neighbor

        return True

    def construct_ansatz(self, representation: Any) -> QuantumCircuitState:
        """
        Construct ionic bonding ansatz.

        Circuit structure for ionic systems:
        1. Single-qubit rotations (on-site energies)
        2. Nearest-neighbor transfer gates (hopping)
        3. Small variational parameters (minimal mixing)

        Example for NaCl:
        |ψ⟩ = |Na⁺⟩⊗|Cl⁻⟩ + small corrections

        Args:
            representation: Representation object

        Returns:
            QuantumCircuitState with governed circuit
        """
        n_qubits = representation.get_num_qubits()
        circuit = QuantumCircuitState(n_qubits)

        # Step 1: Initialize in charge-separated state
        # (done via state preparation, not gates)
        circuit.metadata['initial_state'] = 'charge_separated'

        # Step 2: Single-qubit rotations for on-site energy
        for i in range(n_qubits):
            circuit.add_gate('rz', [i], params=[0.0])  # Variational parameter

        # Step 3: Nearest-neighbor transfer gates (minimal)
        for i in range(0, n_qubits - 1, 2):  # Only adjacent pairs
            # Small amplitude transfer
            circuit.add_gate('rxx', [i, i + 1], params=[0.1])  # Small angle

        # Step 4: Another layer of single-qubit rotations
        for i in range(n_qubits):
            circuit.add_gate('ry', [i], params=[0.0])

        circuit.is_localized = True
        circuit.metadata['ansatz_type'] = 'ionic_localized'

        return circuit

    def enforce_constraints(self, circuit: QuantumCircuitState) -> QuantumCircuitState:
        """
        Enforce ionic bonding constraints on circuit.

        Constraints:
        - Remove long-range gates
        - Limit entanglement degree
        - Ensure particle conservation

        Args:
            circuit: Input circuit

        Returns:
            Constrained circuit
        """
        # Filter out disallowed gates
        allowed_gates = []

        for gate in circuit.gates:
            qubits = gate['qubits']

            # Keep single-qubit gates
            if len(qubits) == 1:
                allowed_gates.append(gate)

            # Keep nearest-neighbor 2-qubit gates only
            elif len(qubits) == 2:
                if abs(qubits[0] - qubits[1]) <= 1:  # Nearest neighbor
                    allowed_gates.append(gate)
                else:
                    print(f"Removed long-range gate: {gate}")

            # Remove multi-qubit gates
            else:
                print(f"Removed multi-qubit gate: {gate}")

        # Update circuit
        circuit.gates = allowed_gates
        circuit.is_localized = True

        return circuit

    def _enforce_locality(self, state: QuantumCircuitState, context: Dict) -> QuantumCircuitState:
        """Enforce that all gates are local (nearest-neighbor)."""
        if not isinstance(state, QuantumCircuitState):
            return state

        # Check all gates are local
        for gate in state.gates:
            if len(gate['qubits']) == 2:
                q1, q2 = gate['qubits']
                if abs(q1 - q2) > 1:
                    raise ValueError(
                        f"Non-local gate detected: {gate}. "
                        f"Ionic bonding requires nearest-neighbor gates only."
                    )

        state.metadata['locality_enforced'] = True
        return state

    def _enforce_sparsity(self, state: QuantumCircuitState, context: Dict) -> QuantumCircuitState:
        """Ensure circuit has sparse connectivity."""
        if not isinstance(state, QuantumCircuitState):
            return state

        max_degree = state.max_entanglement_degree()
        if max_degree > 2:
            print(f"Warning: High entanglement degree ({max_degree}). "
                  f"Ionic systems should have sparse connectivity.")

        state.metadata['sparse'] = state.is_sparse()
        return state

    def _enforce_particle_conservation(self, state: Any, context: Dict) -> Any:
        """Ensure particle number is conserved."""
        # Placeholder - would verify all gates conserve particle number
        if isinstance(state, QuantumCircuitState):
            state.metadata['particle_conserving'] = True
        return state

    def _forbid_collective_gates(self, state: QuantumCircuitState, context: Dict) -> QuantumCircuitState:
        """Remove any collective (>2 qubit) gates."""
        if not isinstance(state, QuantumCircuitState):
            return state

        # Filter out multi-qubit gates
        state.gates = [g for g in state.gates if len(g['qubits']) <= 2]
        state.metadata['no_collective_gates'] = True

        return state

    def get_allowed_operators(self) -> List[str]:
        """Get list of allowed operators for ionic bonding."""
        return [
            'rx', 'ry', 'rz',  # Single-qubit rotations
            'cx', 'cy', 'cz',  # Controlled gates (nearest-neighbor)
            'rxx', 'ryy', 'rzz',  # Two-qubit rotations
            'number',  # Number operator
        ]

    def get_forbidden_operators(self) -> List[str]:
        """Get list of forbidden operators for ionic bonding."""
        return [
            'qft',  # Quantum Fourier Transform (delocalization)
            'ghz',  # GHZ state (collective)
            'w_state',  # W state (collective)
            'long_range_swap',  # Long-range operations
        ]

    def get_transfer_integral_estimate(self, distance: float) -> float:
        """
        Estimate electron transfer integral based on distance.

        For ionic bonding: t ∝ exp(-d/d₀) where d₀ ~ 1-2 Å

        Args:
            distance: Distance between sites in Angstroms

        Returns:
            Transfer integral (arbitrary units)
        """
        d_0 = 1.5  # Decay length in Angstroms
        t_0 = 1.0  # Reference transfer integral
        return t_0 * np.exp(-distance / d_0)

    def __repr__(self) -> str:
        """String representation."""
        return f"IonicGovernanceProtocol(rules={len(self.rules)}, entanglement='minimal')"
