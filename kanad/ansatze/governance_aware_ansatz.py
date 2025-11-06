"""
Governance-aware ansätze using bonding-specific protocols.

These ansätze integrate with governance protocols to enforce
bonding-specific physical constraints in circuit construction.
"""

from typing import List, Optional, Dict
import numpy as np
from kanad.ansatze.base_ansatz import BaseAnsatz, QuantumCircuit, Parameter
from kanad.governance.protocols.ionic_protocol import IonicGovernanceProtocol
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol
from kanad.ansatze.hardware_efficient_ansatz import get_hf_state_qubits


class IonicGovernanceAnsatz(BaseAnsatz):
    """
    Ansatz governed by ionic bonding protocol.

    Physical Model:
        - Localized atomic orbitals
        - Minimal entanglement (charge transfer only)
        - On-site rotations dominant
        - Long-range interactions via protocol

    ADVANTAGES:
        - Physically motivated for ionic systems
        - Enforces localization
        - Efficient for charge-separated states
        - Natural for NaCl, LiF, etc.

    DISADVANTAGES:
        - Not suitable for covalent systems
        - Limited expressivity
    """

    def __init__(
        self,
        n_qubits: int,
        n_electrons: int,
        n_layers: int = 2,
        protocol: Optional[IonicGovernanceProtocol] = None,
        mapper: str = 'jordan_wigner'
    ):
        """
        Initialize ionic governance ansatz.

        Args:
            n_qubits: Number of qubits
            n_electrons: Number of electrons
            n_layers: Number of ansatz layers
            protocol: IonicGovernanceProtocol (optional, will create default)
            mapper: Fermion-to-qubit mapper type ('jordan_wigner', 'bravyi_kitaev', 'parity')
        """
        super().__init__(n_qubits, n_electrons)
        self.n_layers = n_layers
        self.protocol = protocol
        self.mapper = mapper.lower()
        self._built = False

    @property
    def n_parameters(self) -> int:
        """Return number of variational parameters."""
        if not self._built:
            # Build circuit to count parameters
            self.build_circuit()
        # Ionic: 2 parameters per qubit per layer (Rz, Ry)
        return 2 * self.n_qubits * self.n_layers

    def build_circuit(
        self,
        initial_state: Optional[List[int]] = None,
        ionization_threshold: float = 1.5
    ) -> QuantumCircuit:
        """
        Build ionic-governed circuit.

        Structure:
        1. Initial ionic state (charge-separated)
        2. On-site rotations (local excitations)
        3. Minimal charge-transfer entanglement (governed by protocol)

        Args:
            initial_state: Initial charge distribution
            ionization_threshold: Electronegativity difference for ionic bonding

        Returns:
            Ionic ansatz circuit
        """
        circuit = QuantumCircuit(self.n_qubits)

        # 1. Prepare initial ionic state
        if initial_state is not None:
            for qubit, occupation in enumerate(initial_state):
                if occupation == 1:
                    circuit.x(qubit)
            circuit.barrier()
        else:
            # Default: Hartree-Fock state with MAPPER-AWARE preparation
            # Different fermion-to-qubit mappings use different HF state encodings:
            #   - Jordan-Wigner: |1100⟩ for H2 (qubits [2, 3])
            #   - Bravyi-Kitaev: |1000⟩ for H2 (qubit [3] only)
            #   - Parity: Similar to JW

            # CRITICAL: Use mapper-aware HF state preparation
            hf_qubits = get_hf_state_qubits(self.n_qubits, self.n_electrons, self.mapper)
            for qubit in hf_qubits:
                circuit.x(qubit)

            circuit.barrier()

        # 2. Apply ionic-governed layers
        for layer_idx in range(self.n_layers):
            self._apply_ionic_layer(circuit, layer_idx, ionization_threshold)

        self.circuit = circuit
        self._built = True
        return circuit

    def _apply_ionic_layer(
        self,
        circuit: QuantumCircuit,
        layer_idx: int,
        ionization_threshold: float
    ):
        """
        Apply one layer with ionic governance.

        Layer structure:
        1. On-site rotations (RZ dominant for charge state)
        2. Charge-transfer entanglement (governed)
        3. Local RY rotations

        Args:
            circuit: Circuit to modify
            layer_idx: Layer index
            ionization_threshold: Ionization threshold
        """
        # On-site RZ rotations (charge state adjustments)
        for qubit in range(self.n_qubits):
            param_z = Parameter(f'θ_z_{layer_idx}_{qubit}')
            circuit.rz(param_z, qubit)

        # Charge-transfer entanglement (governed by protocol)
        # Only connect qubits that satisfy ionic bonding criteria
        if self.protocol is not None:
            allowed_transfers = self._get_allowed_charge_transfers(ionization_threshold)
        else:
            # Default: linear chain
            allowed_transfers = [(i, i+1) for i in range(self.n_qubits - 1)]

        for (donor, acceptor) in allowed_transfers:
            # Charge transfer = CNOT (donor controls acceptor)
            circuit.cx(donor, acceptor)

        # Local RY rotations (superposition within atomic levels)
        for qubit in range(self.n_qubits):
            param_y = Parameter(f'θ_y_{layer_idx}_{qubit}')
            circuit.ry(param_y, qubit)

    def _get_allowed_charge_transfers(
        self,
        ionization_threshold: float
    ) -> List[tuple]:
        """
        Get allowed charge transfer pairs based on governance protocol.

        Args:
            ionization_threshold: Electronegativity difference threshold

        Returns:
            List of (donor, acceptor) qubit pairs
        """
        # Simplified: assume adjacent pairs with governance check
        # In practice, would use protocol.validate_ionic_bonding()
        allowed = []

        for i in range(self.n_qubits - 1):
            # Assume governance allows adjacent transfers
            # Real implementation would check electronegativity differences
            allowed.append((i, i + 1))

        return allowed

    def __repr__(self) -> str:
        return (f"IonicGovernanceAnsatz(n_qubits={self.n_qubits}, "
                f"layers={self.n_layers}, governed=True)")


class CovalentGovernanceAnsatz(BaseAnsatz):
    """
    Ansatz governed by covalent bonding protocol.

    Physical Model:
        - Hybridized orbital pairs
        - Bonding/antibonding structure
        - Paired entanglement
        - Symmetric electron sharing

    ADVANTAGES:
        - Physically motivated for covalent systems
        - Enforces hybridization
        - Natural for H2, CH4, benzene
        - Captures bond order

    DISADVANTAGES:
        - Not suitable for ionic systems
        - Requires paired qubit structure
    """

    def __init__(
        self,
        n_qubits: int,
        n_electrons: int,
        n_layers: int = 2,
        hybridization: str = 'sp3',
        protocol: Optional[CovalentGovernanceProtocol] = None,
        mapper: str = 'jordan_wigner'
    ):
        """
        Initialize covalent governance ansatz.

        Args:
            n_qubits: Number of qubits (must be even for pairing)
            n_electrons: Number of electrons
            n_layers: Number of ansatz layers
            hybridization: Hybridization type ('sp', 'sp2', 'sp3')
            protocol: CovalentGovernanceProtocol (optional)
            mapper: Fermion-to-qubit mapper type ('jordan_wigner', 'bravyi_kitaev', 'parity')
        """
        super().__init__(n_qubits, n_electrons)
        self.n_layers = n_layers
        self.hybridization = hybridization
        self.protocol = protocol
        self.mapper = mapper.lower()
        self._built = False

        if n_qubits % 2 != 0:
            raise ValueError("Covalent ansatz requires even number of qubits for MO pairing")

    @property
    def n_parameters(self) -> int:
        """Return number of variational parameters."""
        if not self._built:
            # Build circuit to count parameters
            self.build_circuit()
        # Covalent: 2 parameters per qubit (hybridization) + 2 per pair (bonding)
        n_pairs = self.n_qubits // 2
        return (2 * self.n_qubits + 2 * n_pairs) * self.n_layers

    def build_circuit(
        self,
        initial_state: Optional[List[int]] = None,
        overlap_threshold: float = 0.3
    ) -> QuantumCircuit:
        """
        Build covalent-governed circuit.

        Structure:
        1. Initial Hartree-Fock state
        2. Hybridization rotations (create hybrid orbitals)
        3. Bonding/antibonding pair entanglement
        4. Correlation corrections

        Args:
            initial_state: Initial orbital occupation
            overlap_threshold: Minimum orbital overlap for bonding

        Returns:
            Covalent ansatz circuit
        """
        circuit = QuantumCircuit(self.n_qubits)

        # 1. Prepare initial state (Hartree-Fock)
        if initial_state is not None:
            for qubit, occupation in enumerate(initial_state):
                if occupation == 1:
                    circuit.x(qubit)
            circuit.barrier()
        else:
            # Default: Hartree-Fock state with MAPPER-AWARE preparation
            # Different fermion-to-qubit mappings use different HF state encodings:
            #   - Jordan-Wigner: |1100⟩ for H2 (qubits [2, 3])
            #   - Bravyi-Kitaev: |1000⟩ for H2 (qubit [3] only)
            #   - Parity: Similar to JW

            # CRITICAL: Use mapper-aware HF state preparation
            hf_qubits = get_hf_state_qubits(self.n_qubits, self.n_electrons, self.mapper)
            for qubit in hf_qubits:
                circuit.x(qubit)

            circuit.barrier()

        # 2. Apply covalent-governed layers
        for layer_idx in range(self.n_layers):
            self._apply_covalent_layer(circuit, layer_idx, overlap_threshold)

        self.circuit = circuit
        self._built = True
        return circuit

    def _apply_covalent_layer(
        self,
        circuit: QuantumCircuit,
        layer_idx: int,
        overlap_threshold: float
    ):
        """
        Apply one layer with covalent governance.

        Layer structure:
        1. Hybridization rotations (RY + RZ)
        2. Bonding/antibonding pair entanglement (RXX + RYY)
        3. Intra-pair correlations

        Args:
            circuit: Circuit to modify
            layer_idx: Layer index
            overlap_threshold: Overlap threshold for bonding
        """
        # Hybridization rotations (sp3 uses both RY and RZ)
        for qubit in range(self.n_qubits):
            param_y = Parameter(f'θ_hyb_y_{layer_idx}_{qubit}')
            param_z = Parameter(f'θ_hyb_z_{layer_idx}_{qubit}')
            circuit.ry(param_y, qubit)
            circuit.rz(param_z, qubit)

        # Bonding/antibonding pair entanglement
        n_pairs = self.n_qubits // 2

        for pair_idx in range(n_pairs):
            bonding_qubit = 2 * pair_idx
            antibonding_qubit = 2 * pair_idx + 1

            # Check governance (would use protocol.validate_orbital_hybridization)
            if self._is_bonding_allowed(pair_idx, overlap_threshold):
                # Bonding entanglement: RXX + RYY (real amplitudes)
                param_xx = Parameter(f'θ_bond_xx_{layer_idx}_{pair_idx}')
                param_yy = Parameter(f'θ_bond_yy_{layer_idx}_{pair_idx}')

                circuit.rxx(param_xx, bonding_qubit, antibonding_qubit)
                circuit.ryy(param_yy, bonding_qubit, antibonding_qubit)

        # Inter-pair correlations (between different bonds)
        if self.protocol is not None:
            allowed_correlations = self._get_allowed_correlations()
        else:
            # Default: adjacent pairs
            allowed_correlations = [(2*i, 2*(i+1)) for i in range(n_pairs - 1)]

        for (qubit_i, qubit_j) in allowed_correlations:
            circuit.cz(qubit_i, qubit_j)

    def _is_bonding_allowed(
        self,
        pair_idx: int,
        overlap_threshold: float
    ) -> bool:
        """
        Check if bonding is allowed for this orbital pair.

        Args:
            pair_idx: MO pair index
            overlap_threshold: Minimum overlap required

        Returns:
            True if bonding allowed
        """
        # Simplified: allow all pairs
        # Real implementation would use protocol.validate_orbital_hybridization()
        # with actual overlap integrals
        return True

    def _get_allowed_correlations(self) -> List[tuple]:
        """
        Get allowed inter-bond correlations.

        Returns:
            List of (qubit_i, qubit_j) pairs for correlations
        """
        # Simplified: connect bonding orbitals of adjacent pairs
        n_pairs = self.n_qubits // 2
        correlations = []

        for i in range(n_pairs - 1):
            # Connect bonding orbital of pair i to bonding orbital of pair i+1
            correlations.append((2 * i, 2 * (i + 1)))

        return correlations

    def __repr__(self) -> str:
        return (f"CovalentGovernanceAnsatz(n_qubits={self.n_qubits}, "
                f"layers={self.n_layers}, hybridization='{self.hybridization}', governed=True)")


class AdaptiveGovernanceAnsatz(BaseAnsatz):
    """
    Adaptive ansatz that switches governance based on bonding type.

    Automatically selects ionic or covalent governance based on
    molecular properties.
    """

    def __init__(
        self,
        n_qubits: int,
        n_electrons: int,
        n_layers: int = 2,
        bonding_type: str = 'auto'
    ):
        """
        Initialize adaptive governance ansatz.

        Args:
            n_qubits: Number of qubits
            n_electrons: Number of electrons
            n_layers: Number of ansatz layers
            bonding_type: 'ionic', 'covalent', or 'auto'
        """
        super().__init__(n_qubits, n_electrons)
        self.n_layers = n_layers
        self.bonding_type = bonding_type
        self._delegate_ansatz = None

    @property
    def n_parameters(self) -> int:
        """
        Return number of variational parameters.

        Since this is an adaptive ansatz that delegates, we need to build
        the circuit first to determine the parameter count.
        """
        if self._delegate_ansatz is not None:
            return self._delegate_ansatz.n_parameters
        else:
            # If not yet built, estimate based on covalent (larger parameter count)
            # Covalent: 3 * n_qubits * n_layers
            # This will be updated when build_circuit() is called
            return 3 * self.n_qubits * self.n_layers

    def build_circuit(
        self,
        initial_state: Optional[List[int]] = None,
        electronegativity_diff: Optional[float] = None,
        overlap_matrix: Optional[np.ndarray] = None
    ) -> QuantumCircuit:
        """
        Build adaptive circuit.

        Selects ionic or covalent governance based on chemical properties.

        Args:
            initial_state: Initial state
            electronegativity_diff: Electronegativity difference between atoms
            overlap_matrix: Orbital overlap matrix

        Returns:
            Adaptive ansatz circuit
        """
        # Determine bonding type
        if self.bonding_type == 'auto':
            bonding_type = self._determine_bonding_type(
                electronegativity_diff,
                overlap_matrix
            )
        else:
            bonding_type = self.bonding_type

        # Delegate to appropriate ansatz
        if bonding_type == 'ionic':
            self._delegate_ansatz = IonicGovernanceAnsatz(
                self.n_qubits,
                self.n_electrons,
                self.n_layers
            )
            circuit = self._delegate_ansatz.build_circuit(initial_state)
        else:  # covalent
            # Ensure even qubits for covalent
            if self.n_qubits % 2 != 0:
                raise ValueError("Covalent bonding requires even number of qubits")

            self._delegate_ansatz = CovalentGovernanceAnsatz(
                self.n_qubits,
                self.n_electrons,
                self.n_layers
            )
            circuit = self._delegate_ansatz.build_circuit(initial_state)

        self.circuit = circuit
        return circuit

    def _determine_bonding_type(
        self,
        electronegativity_diff: Optional[float],
        overlap_matrix: Optional[np.ndarray]
    ) -> str:
        """
        Automatically determine bonding type.

        Criteria:
        - Ionic: ΔEN > 1.7, low overlap
        - Covalent: ΔEN < 1.7, high overlap

        Args:
            electronegativity_diff: Electronegativity difference
            overlap_matrix: Overlap matrix

        Returns:
            'ionic' or 'covalent'
        """
        # Use electronegativity if available
        if electronegativity_diff is not None:
            if electronegativity_diff > 1.7:
                return 'ionic'
            else:
                return 'covalent'

        # Use overlap matrix if available
        if overlap_matrix is not None:
            # High off-diagonal overlap → covalent
            off_diag = overlap_matrix - np.diag(np.diag(overlap_matrix))
            max_overlap = np.max(np.abs(off_diag))

            if max_overlap > 0.3:
                return 'covalent'
            else:
                return 'ionic'

        # Default: covalent
        return 'covalent'

    def get_num_parameters(self) -> int:
        """
        Get number of parameters.

        Returns number from delegate ansatz if built, otherwise estimate.
        """
        if self._delegate_ansatz is not None:
            return self._delegate_ansatz.get_num_parameters()
        else:
            # Estimate before circuit is built
            # Assume covalent (more parameters)
            return self.n_layers * self.n_qubits * 3

    def __repr__(self) -> str:
        return (f"AdaptiveGovernanceAnsatz(n_qubits={self.n_qubits}, "
                f"bonding_type='{self.bonding_type}', adaptive=True)")
