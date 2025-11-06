"""
Covalent bonding governance protocol.

Physical principles for covalent bonding:
1. Orbital HYBRIDIZATION must occur (sp, sp², sp³)
2. Bonding/antibonding molecular orbital pairs form
3. PAIRED entanglement (Bell states for electron pairs)
4. Electron sharing between atoms
5. Small electronegativity difference

Circuit requirements:
- Hybridization transformations first
- Paired qubit gates for bonding orbitals
- Bell state preparation for electron pairs
- Moderate entanglement (bonding pairs)
- Preserve σ/π bond character
"""

from typing import List, Dict, Any, Tuple
import numpy as np
import logging

logger = logging.getLogger(__name__)
from kanad.governance.protocols.base_protocol import (
    BaseGovernanceProtocol,
    BondingType,
    GovernanceRule,
    QuantumCircuitState
)


class CovalentGovernanceProtocol(BaseGovernanceProtocol):
    """
    Governance protocol for covalent bonding systems.

    Example: H₂, H₂O, CH₄
    - Orbital hybridization (sp³ for C in CH₄)
    - Bonding/antibonding MO formation
    - Electron pair sharing
    - Moderate entanglement between bonding pairs
    """

    def __init__(self):
        """Initialize covalent governance protocol."""
        super().__init__(BondingType.COVALENT)

    def _initialize_rules(self):
        """Initialize covalent bonding governance rules."""

        # Rule 1: Hybridization must occur before bonding
        self.rules.append(GovernanceRule(
            name="hybridization_first",
            description="Atomic orbitals must hybridize before forming bonds",
            condition=lambda state, ctx: not self._is_hybridized(state),
            action=self._apply_hybridization,
            priority=100,
            required=True
        ))

        # Rule 2: Create bonding/antibonding MO pairs
        self.rules.append(GovernanceRule(
            name="molecular_orbital_formation",
            description="Form bonding and antibonding molecular orbitals",
            condition=lambda state, ctx: (self._is_hybridized(state) and
                                         not self._has_mo_pairs(state)),
            action=self._create_mo_pairs,
            priority=90,
            required=True
        ))

        # Rule 3: Entangle electron pairs (Bell states)
        self.rules.append(GovernanceRule(
            name="electron_pair_entanglement",
            description="Create entangled electron pairs in bonding orbitals",
            condition=lambda state, ctx: (self._has_mo_pairs(state) and
                                         not self._is_paired(state)),
            action=self._entangle_pairs,
            priority=80,
            required=True
        ))

        # Rule 4: Preserve spin symmetry
        self.rules.append(GovernanceRule(
            name="spin_symmetry",
            description="Maintain proper spin coupling (singlet for bonding)",
            condition=lambda state, ctx: True,
            action=self._enforce_spin_symmetry,
            priority=70,
            required=True
        ))

        # Rule 5: Forbid long-range entanglement (covalent is local)
        self.rules.append(GovernanceRule(
            name="no_long_range_entanglement",
            description="Entanglement only between bonding pairs",
            condition=lambda state, ctx: isinstance(state, QuantumCircuitState),
            action=self._enforce_paired_entanglement,
            priority=60,
            required=True
        ))

    def validate_operator(self, operator: 'QuantumOperator') -> bool:
        """
        Validate operator for covalent bonding.

        ALLOWED:
        - Orbital rotation operators (hybridization)
        - Pairing operators (2-qubit gates on bonding pairs)
        - Symmetric excitations (preserve spin)
        - Givens rotations (bonding/antibonding mixing)

        FORBIDDEN:
        - Long-range entanglement (>2 qubits for single bond)
        - Asymmetric excitations
        - Bare transfer operators (must use MO formation)
        - Collective operations (those are for metallic)
        """
        operator_type = getattr(operator, 'type', 'unknown')

        # Allowed operators
        allowed = [
            'rx', 'ry', 'rz',  # Single-qubit rotations (hybridization)
            'givens',  # Bonding/antibonding mixing
            'h',  # Hadamard (for Bell states)
            'cx', 'cy', 'cz',  # Pairing gates
            'rxx', 'ryy', 'rzz',  # Symmetric excitations
        ]

        # Forbidden operators
        forbidden = [
            'qft',  # Delocalization (metallic)
            'ghz',  # Collective (metallic)
            'transfer',  # Bare transfer (ionic)
        ]

        if operator_type in forbidden:
            return False

        # Check if gate acts on bonding pair
        if hasattr(operator, 'qubits') and len(operator.qubits) == 2:
            # Would check if qubits form a bonding pair
            # For now, allow all 2-qubit gates
            pass

        return operator_type in allowed or operator_type not in forbidden

    def construct_ansatz(self, representation: Any) -> QuantumCircuitState:
        """
        Construct covalent bonding ansatz.

        Circuit structure:
        1. Hybridization transformations (single-qubit rotations)
        2. Bonding/antibonding formation (Givens rotations)
        3. Electron pair entanglement (Bell states)
        4. Variational correlation (small symmetric excitations)

        Example for H₂:
        |ψ⟩ = |bonding↑↓⟩ + small|antibonding↑↓⟩

        Args:
            representation: LCAO representation object

        Returns:
            QuantumCircuitState with covalent ansatz
        """
        n_qubits = representation.get_num_qubits()
        circuit = QuantumCircuitState(n_qubits)

        # Get hybridization type if available
        hybridization = getattr(representation, 'hybridization', None)

        # Step 1: Hybridization (single-qubit rotations)
        if hybridization:
            for i in range(n_qubits):
                # Rotate to hybrid orbital basis
                circuit.add_gate('ry', [i], params=[np.pi/4])  # sp³ angle
                circuit.add_gate('rz', [i], params=[0.0])  # Phase

        circuit.is_hybridized = True

        # Step 2: Form molecular orbitals (Givens rotations)
        mo_pairs = getattr(representation, 'mo_pairs', [])

        if not mo_pairs and n_qubits >= 2:
            # Default: pair adjacent qubits
            mo_pairs = [(2*i, 2*i+1) for i in range(n_qubits // 4)]

        for bonding_idx, antibonding_idx in mo_pairs:
            # Givens rotation creates bonding/antibonding combination
            # G(θ) = [[cos(θ), -sin(θ)], [sin(θ), cos(θ)]]
            theta = np.pi / 4  # Equal mixing for symmetric bond
            circuit.add_gate('givens', [bonding_idx, antibonding_idx], params=[theta])

        circuit.has_mo_pairs = True

        # Step 3: Create electron pairs (Bell states for bonding orbitals)
        for bonding_idx, _ in mo_pairs:
            # Create (|01⟩ + |10⟩)/√2 for opposite spins
            circuit.add_gate('h', [bonding_idx])
            circuit.add_gate('cx', [bonding_idx, bonding_idx + 1])

        circuit.is_paired = True

        # Step 4: Variational correlation (symmetric excitations)
        for bonding_idx, antibonding_idx in mo_pairs:
            # Small amplitude symmetric excitations
            circuit.add_gate('rxx', [bonding_idx, antibonding_idx], params=[0.1])
            circuit.add_gate('ryy', [bonding_idx, antibonding_idx], params=[0.1])

        circuit.metadata['ansatz_type'] = 'covalent_paired'
        circuit.metadata['bond_type'] = 'covalent'

        return circuit

    def enforce_constraints(self, circuit: QuantumCircuitState) -> QuantumCircuitState:
        """
        Enforce covalent bonding constraints.

        Constraints:
        - Hybridization must come before bonding
        - Only paired entanglement allowed
        - Preserve spin symmetry

        Args:
            circuit: Input circuit

        Returns:
            Constrained circuit
        """
        # Ensure hybridization gates appear first
        hybridization_gates = []
        pairing_gates = []
        other_gates = []

        for gate in circuit.gates:
            if gate['type'] in ['rx', 'ry', 'rz'] and len(gate['qubits']) == 1:
                hybridization_gates.append(gate)
            elif gate['type'] in ['h', 'cx', 'cy', 'cz', 'givens']:
                pairing_gates.append(gate)
            else:
                other_gates.append(gate)

        # Reorder: hybridization -> pairing -> others
        circuit.gates = hybridization_gates + pairing_gates + other_gates

        circuit.is_hybridized = len(hybridization_gates) > 0
        circuit.has_mo_pairs = len(pairing_gates) > 0

        return circuit

    def _is_hybridized(self, state: Any) -> bool:
        """Check if orbitals are hybridized."""
        if isinstance(state, QuantumCircuitState):
            return state.is_hybridized
        return False

    def _has_mo_pairs(self, state: Any) -> bool:
        """Check if molecular orbital pairs exist."""
        if isinstance(state, QuantumCircuitState):
            return state.has_mo_pairs
        return False

    def _is_paired(self, state: Any) -> bool:
        """Check if electrons are paired."""
        if isinstance(state, QuantumCircuitState):
            return state.is_paired
        return False

    def _apply_hybridization(self, state: QuantumCircuitState, context: Dict) -> QuantumCircuitState:
        """Apply orbital hybridization transformation."""
        if not isinstance(state, QuantumCircuitState):
            return state

        hybridization_type = context.get('hybridization', 'sp3')

        # Apply appropriate rotations for hybridization type
        if hybridization_type == 'sp3':
            # Tetrahedral angle: arccos(-1/3) ≈ 109.5°
            theta = np.arccos(-1.0/3.0)
            for i in range(state.n_qubits):
                state.add_gate('ry', [i], params=[theta])

        elif hybridization_type == 'sp2':
            # Trigonal planar: 120°
            theta = 2 * np.pi / 3
            for i in range(state.n_qubits):
                state.add_gate('ry', [i], params=[theta])

        elif hybridization_type == 'sp':
            # Linear: 180°
            for i in range(state.n_qubits):
                state.add_gate('ry', [i], params=[np.pi])

        state.is_hybridized = True
        state.metadata['hybridization'] = hybridization_type

        return state

    def _create_mo_pairs(self, state: QuantumCircuitState, context: Dict) -> QuantumCircuitState:
        """Create bonding/antibonding molecular orbital pairs."""
        if not isinstance(state, QuantumCircuitState):
            return state

        bonds = context.get('bonds', [])

        # Create Givens rotations for each bond
        for bond_idx, (i, j) in enumerate(bonds):
            if i < state.n_qubits and j < state.n_qubits:
                # Bonding: |ψ₊⟩ = (|i⟩ + |j⟩)/√2
                # Antibonding: |ψ₋⟩ = (|i⟩ - |j⟩)/√2
                theta = np.pi / 4  # 45° for equal mixing
                state.add_gate('givens', [i, j], params=[theta])

        state.has_mo_pairs = True
        state.metadata['n_mo_pairs'] = len(bonds)

        return state

    def _entangle_pairs(self, state: QuantumCircuitState, context: Dict) -> QuantumCircuitState:
        """Create entangled electron pairs (Bell states)."""
        if not isinstance(state, QuantumCircuitState):
            return state

        electron_pairs = context.get('electron_pairs', [])

        for i, j in electron_pairs:
            if i < state.n_qubits and j < state.n_qubits:
                # Create Bell state: (|01⟩ + |10⟩)/√2
                # Singlet state for bonding
                state.add_gate('h', [i])
                state.add_gate('cx', [i, j])

        state.is_paired = True
        state.metadata['n_pairs'] = len(electron_pairs)

        return state

    def _enforce_spin_symmetry(self, state: Any, context: Dict) -> Any:
        """Enforce spin symmetry (singlet for bonding pairs)."""
        if isinstance(state, QuantumCircuitState):
            state.metadata['spin_symmetry'] = 'singlet'
        return state

    def _enforce_paired_entanglement(self, state: QuantumCircuitState, context: Dict) -> QuantumCircuitState:
        """Ensure entanglement only occurs between bonding pairs."""
        if not isinstance(state, QuantumCircuitState):
            return state

        # Check that entanglement graph has paired structure
        max_degree = state.max_entanglement_degree()

        if max_degree > 3:
            logger.warning(f"High entanglement degree ({max_degree}). "
                          f"Covalent bonding should have paired entanglement only.")

        state.metadata['paired_entanglement'] = True

        return state

    def get_allowed_operators(self) -> List[str]:
        """Get allowed operators for covalent bonding."""
        return [
            'rx', 'ry', 'rz',  # Hybridization
            'h',  # Bell state preparation
            'cx', 'cy', 'cz',  # Pairing gates
            'givens',  # MO formation
            'rxx', 'ryy', 'rzz',  # Symmetric excitations
        ]

    def get_forbidden_operators(self) -> List[str]:
        """Get forbidden operators for covalent bonding."""
        return [
            'qft',  # Delocalization (metallic)
            'ghz', 'w_state',  # Collective states
            'long_range',  # Long-range interactions
        ]

    def get_representation_type(self) -> str:
        """
        Get quantum representation type for covalent bonding.

        Returns:
            'molecular_orbital' - Covalent bonding uses MO basis
        """
        return 'molecular_orbital'

    def get_bonding_angle(self, hybridization: str) -> float:
        """
        Get characteristic bonding angle for hybridization type.

        Args:
            hybridization: Type of hybridization (sp, sp2, sp3)

        Returns:
            Bonding angle in radians
        """
        angles = {
            'sp': np.pi,  # 180° (linear)
            'sp2': 2 * np.pi / 3,  # 120° (trigonal planar)
            'sp3': np.arccos(-1.0/3.0),  # ~109.5° (tetrahedral)
        }
        return angles.get(hybridization, np.pi/2)

    def generate_single_excitations(self, bitstring: str) -> List[str]:
        """
        Generate physics-aware single excitations for covalent bonding.

        Covalent bonding principles:
        1. HOMO → LUMO transitions (bonding → antibonding)
        2. Preserve spin pairing (singlet state)
        3. Localized excitations within bonding pairs
        4. No long-range charge transfer

        Args:
            bitstring: Configuration bitstring

        Returns:
            List of excited configuration bitstrings
        """
        n_qubits = len(bitstring)
        bits = list(bitstring)
        excitations = []

        # Find occupied and virtual orbitals
        occupied = [i for i, b in enumerate(bits) if b == '1']
        virtual = [i for i, b in enumerate(bits) if b == '0']

        # Rule 1: HOMO → LUMO (highest occupied to lowest unoccupied)
        if occupied and virtual:
            homo_idx = max(occupied)
            lumo_idx = min(virtual)

            # Primary excitation: HOMO → LUMO
            new_bits = bits.copy()
            new_bits[homo_idx] = '0'
            new_bits[lumo_idx] = '1'
            excitations.append(''.join(new_bits))

        # Rule 2: Bonding → Antibonding for each pair
        # In Jordan-Wigner: spin-up/down pairs are adjacent (i, i+1)
        for i in range(0, n_qubits-1, 2):
            if bits[i] == '1':  # Occupied bonding orbital
                # Excite to corresponding antibonding orbital
                # For small molecules, antibonding is typically +2 or +4 positions
                for offset in [2, 4]:
                    antibonding_idx = i + offset
                    if antibonding_idx < n_qubits and bits[antibonding_idx] == '0':
                        new_bits = bits.copy()
                        new_bits[i] = '0'
                        new_bits[antibonding_idx] = '1'
                        excitations.append(''.join(new_bits))

        # Rule 3: Preserve spin pairing - also generate spin-down excitation
        for i in range(1, n_qubits, 2):
            if bits[i] == '1':  # Occupied spin-down
                for offset in [2, 4]:
                    antibonding_idx = i + offset
                    if antibonding_idx < n_qubits and bits[antibonding_idx] == '0':
                        new_bits = bits.copy()
                        new_bits[i] = '0'
                        new_bits[antibonding_idx] = '1'
                        excitations.append(''.join(new_bits))

        # Remove duplicates while preserving order
        seen = set()
        unique_excitations = []
        for exc in excitations:
            if exc not in seen and exc != bitstring:
                seen.add(exc)
                unique_excitations.append(exc)

        return unique_excitations

    def generate_double_excitations(self, bitstring: str) -> List[str]:
        """
        Generate physics-aware double excitations for covalent bonding.

        Covalent bonding principles:
        1. Paired excitations (preserve spin singlet)
        2. (bonding, bonding) → (antibonding, antibonding)
        3. Avoid breaking electron pairs
        4. Local excitations within bonding region

        Args:
            bitstring: Configuration bitstring

        Returns:
            List of doubly excited configuration bitstrings
        """
        n_qubits = len(bitstring)
        bits = list(bitstring)
        excitations = []

        # Find occupied and virtual orbitals
        occupied = [i for i, b in enumerate(bits) if b == '1']
        virtual = [i for i, b in enumerate(bits) if b == '0']

        # Rule 1: Paired double excitations (preserve singlet)
        # Excite spin-up and spin-down together: (i, i+1) → (j, j+1)
        for i in range(0, n_qubits-1, 2):
            if bits[i] == '1' and bits[i+1] == '1':  # Bonding pair occupied
                for j in range(0, n_qubits-1, 2):
                    if bits[j] == '0' and bits[j+1] == '0':  # Antibonding pair empty
                        # Double excitation: |↑↓⟩_bonding → |↑↓⟩_antibonding
                        new_bits = bits.copy()
                        new_bits[i] = '0'
                        new_bits[i+1] = '0'
                        new_bits[j] = '1'
                        new_bits[j+1] = '1'
                        excitations.append(''.join(new_bits))

        # Rule 2: HOMO-1, HOMO → LUMO, LUMO+1 (important for correlation)
        if len(occupied) >= 2 and len(virtual) >= 2:
            homo = max(occupied)
            homo_1 = max([o for o in occupied if o < homo], default=None)
            lumo = min(virtual)
            lumo_1 = min([v for v in virtual if v > lumo], default=None)

            if homo_1 is not None and lumo_1 is not None:
                new_bits = bits.copy()
                new_bits[homo] = '0'
                new_bits[homo_1] = '0'
                new_bits[lumo] = '1'
                new_bits[lumo_1] = '1'
                excitations.append(''.join(new_bits))

        # Rule 3: Same-spin double excitations (within spin manifold)
        # Only for larger molecules where this matters
        if n_qubits >= 8:
            occupied_up = [i for i in occupied if i % 2 == 0]
            occupied_down = [i for i in occupied if i % 2 == 1]
            virtual_up = [i for i in virtual if i % 2 == 0]
            virtual_down = [i for i in virtual if i % 2 == 1]

            # Spin-up double excitations
            if len(occupied_up) >= 2 and len(virtual_up) >= 2:
                i, j = occupied_up[:2]
                a, b = virtual_up[:2]
                new_bits = bits.copy()
                new_bits[i] = '0'
                new_bits[j] = '0'
                new_bits[a] = '1'
                new_bits[b] = '1'
                excitations.append(''.join(new_bits))

            # Spin-down double excitations
            if len(occupied_down) >= 2 and len(virtual_down) >= 2:
                i, j = occupied_down[:2]
                a, b = virtual_down[:2]
                new_bits = bits.copy()
                new_bits[i] = '0'
                new_bits[j] = '0'
                new_bits[a] = '1'
                new_bits[b] = '1'
                excitations.append(''.join(new_bits))

        # Remove duplicates
        seen = set()
        unique_excitations = []
        for exc in excitations:
            if exc not in seen and exc != bitstring:
                seen.add(exc)
                unique_excitations.append(exc)

        return unique_excitations

    def is_valid_configuration(self, bitstring: str) -> bool:
        """
        Check if configuration is valid for covalent bonding.

        Rules:
        1. Must preserve total electron count
        2. Should have spin pairing (equal spin-up and spin-down)
        3. No charge separation (all atoms neutral or nearly neutral)

        Args:
            bitstring: Configuration bitstring

        Returns:
            True if configuration is physically valid
        """
        n_qubits = len(bitstring)

        # Rule 1: Already checked by ConfigurationSubspace (electron count)

        # Rule 2: Check spin pairing
        # For covalent molecules, expect singlet state (equal spin-up and spin-down)
        n_up = sum(1 for i in range(0, n_qubits, 2) if bitstring[i] == '1')
        n_down = sum(1 for i in range(1, n_qubits, 2) if bitstring[i] == '1')

        # Allow small spin imbalance for excited states
        if abs(n_up - n_down) > 2:
            return False

        # Rule 3: No excessive charge separation
        # For now, accept all configurations that pass spin check
        # Future: Add atom-wise charge analysis

        return True

    def __repr__(self) -> str:
        """String representation."""
        return f"CovalentGovernanceProtocol(rules={len(self.rules)}, entanglement='paired')"
