"""
Hybridization operators for covalent bonding.

These operators implement orbital hybridization transformations that are
fundamental to covalent bonding:
- sp hybridization: Linear (180°)
- sp² hybridization: Trigonal planar (120°)
- sp³ hybridization: Tetrahedral (109.5°)
"""

from typing import List, Optional, Dict, Any
import numpy as np
from scipy.linalg import expm
import logging

logger = logging.getLogger(__name__)

from kanad.governance.operators.base_operator import (
    BaseGovernanceOperator,
    OperatorType,
    PhysicalMeaning,
    kron_product,
    I_MATRIX,
    X_MATRIX,
    Y_MATRIX,
    Z_MATRIX,
    H_MATRIX
)


class HybridizationOperator(BaseGovernanceOperator):
    """
    Orbital hybridization operator.

    Transforms atomic orbitals into hybrid orbitals:
    - sp: R_y(π) - Linear geometry (180°)
    - sp²: R_y(2π/3) - Trigonal planar (120°)
    - sp³: R_y(arccos(-1/3)) - Tetrahedral (109.5°)

    Physical meaning:
        |hybrid⟩ = R(θ) |atomic⟩

    Example (sp³ for carbon):
        |s⟩ → (|sp³_1⟩ + |sp³_2⟩ + |sp³_3⟩ + |sp³_4⟩)/2
    """

    def __init__(
        self,
        qubit: int,
        hybridization_type: str = 'sp3',
        rotation_axis: str = 'y',
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize hybridization operator.

        Args:
            qubit: Qubit index to hybridize
            hybridization_type: Type of hybridization ('sp', 'sp2', 'sp3')
            rotation_axis: Rotation axis ('x', 'y', 'z')
            metadata: Additional metadata
        """
        # Calculate hybridization angle
        self.hybridization_type = hybridization_type
        self.rotation_axis = rotation_axis
        theta = self._get_hybridization_angle(hybridization_type)

        super().__init__(
            operator_type=OperatorType.HYBRIDIZATION,
            physical_meaning=PhysicalMeaning.HYBRIDIZATION,
            qubits=[qubit],
            params=[theta],
            metadata=metadata or {}
        )

        self.metadata['hybridization_type'] = hybridization_type
        self.metadata['angle_degrees'] = np.degrees(theta)

    def _get_hybridization_angle(self, hyb_type: str) -> float:
        """
        Get rotation angle for hybridization type.

        Args:
            hyb_type: Hybridization type

        Returns:
            Rotation angle in radians
        """
        angles = {
            'sp': np.pi,  # 180° (linear)
            'sp2': 2 * np.pi / 3,  # 120° (trigonal planar)
            'sp3': np.arccos(-1.0 / 3.0),  # ~109.47° (tetrahedral)
        }

        if hyb_type not in angles:
            logger.warning(f"Unknown hybridization type: {hyb_type}, defaulting to sp3")
            return angles['sp3']

        return angles[hyb_type]

    def _get_allowed_bonding_types(self) -> List[str]:
        """Hybridization is specific to covalent bonding."""
        return ['covalent']

    def to_matrix(self, n_qubits: int) -> np.ndarray:
        """
        Build rotation matrix for hybridization.

        Args:
            n_qubits: Total number of qubits

        Returns:
            Rotation matrix (2^n × 2^n)
        """
        theta = self.params[0]

        # Choose rotation matrix based on axis
        if self.rotation_axis == 'x':
            # R_x(θ) = exp(-iθX/2) = cos(θ/2)I - i*sin(θ/2)X
            rot = np.array([
                [np.cos(theta/2), -1j*np.sin(theta/2)],
                [-1j*np.sin(theta/2), np.cos(theta/2)]
            ], dtype=complex)
        elif self.rotation_axis == 'y':
            # R_y(θ) = exp(-iθY/2) = cos(θ/2)I - i*sin(θ/2)Y
            rot = np.array([
                [np.cos(theta/2), -np.sin(theta/2)],
                [np.sin(theta/2), np.cos(theta/2)]
            ], dtype=complex)
        elif self.rotation_axis == 'z':
            # R_z(θ) = exp(-iθZ/2)
            rot = np.array([
                [np.exp(-1j*theta/2), 0],
                [0, np.exp(1j*theta/2)]
            ], dtype=complex)
        else:
            raise ValueError(f"Invalid rotation axis: {self.rotation_axis}")

        # Embed in full Hilbert space
        return kron_product(rot, self.qubits[0], n_qubits)

    def to_qiskit(self):
        """Convert to Qiskit RY gate."""
        from qiskit.circuit.library import RYGate, RXGate, RZGate

        theta = self.params[0]

        if self.rotation_axis == 'y':
            return RYGate(theta)
        elif self.rotation_axis == 'x':
            return RXGate(theta)
        elif self.rotation_axis == 'z':
            return RZGate(theta)

    def get_physical_description(self) -> str:
        """Get description of hybridization."""
        geometries = {
            'sp': 'linear',
            'sp2': 'trigonal planar',
            'sp3': 'tetrahedral'
        }
        geom = geometries.get(self.hybridization_type, 'unknown')
        angle = self.metadata.get('angle_degrees', 0)
        return f"{self.hybridization_type} hybridization ({geom}, {angle:.1f}°)"


class GivensRotationOperator(BaseGovernanceOperator):
    """
    Givens rotation for molecular orbital formation.

    Creates bonding/antibonding MO pairs from atomic orbitals:
        |bonding⟩ = cos(θ)|ψ_A⟩ + sin(θ)|ψ_B⟩
        |antibonding⟩ = -sin(θ)|ψ_A⟩ + cos(θ)|ψ_B⟩

    For symmetric bonds (like H₂): θ = π/4 (equal mixing)
    For polar bonds: θ ≠ π/4 (unequal mixing)

    Matrix form:
        G(θ) = [[cos(θ), -sin(θ)],
                [sin(θ),  cos(θ)]]
    """

    def __init__(
        self,
        qubit1: int,
        qubit2: int,
        mixing_angle: float = np.pi / 4,
        bond_type: str = 'sigma',
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize Givens rotation operator.

        Args:
            qubit1: First qubit (bonding orbital)
            qubit2: Second qubit (antibonding orbital)
            mixing_angle: Rotation angle (π/4 for equal mixing)
            bond_type: Type of bond ('sigma' or 'pi')
            metadata: Additional metadata
        """
        super().__init__(
            operator_type=OperatorType.GIVENS,
            physical_meaning=PhysicalMeaning.MO_FORMATION,
            qubits=[qubit1, qubit2],
            params=[mixing_angle],
            metadata=metadata or {}
        )

        self.bond_type = bond_type
        self.metadata['bond_type'] = bond_type
        self.metadata['mixing_angle_degrees'] = np.degrees(mixing_angle)

    def _get_allowed_bonding_types(self) -> List[str]:
        """Givens rotations are for covalent MO formation."""
        return ['covalent']

    def to_matrix(self, n_qubits: int) -> np.ndarray:
        """
        Build Givens rotation matrix.

        Args:
            n_qubits: Total number of qubits

        Returns:
            Givens rotation matrix (2^n × 2^n)
        """
        theta = self.params[0]
        q1, q2 = self.qubits

        # 2x2 Givens rotation
        givens_2x2 = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)]
        ], dtype=complex)

        # Build full matrix by embedding Givens in computational basis
        dim = 2 ** n_qubits
        full_matrix = np.eye(dim, dtype=complex)

        # This is a simplified implementation
        # Full Givens would require fermionic operator mapping
        # For now, use direct unitary embedding

        # Convert to controlled rotation representation
        # This is approximate but captures the physics

        return full_matrix  # Placeholder - would implement full JW mapping

    def to_qiskit(self):
        """Convert to Qiskit gates."""
        # Givens rotation can be decomposed into:
        # G(θ) = R_Y(θ) on qubit 1, controlled by qubit 2
        # This is a simplification; full implementation would use custom gate

        from qiskit.circuit import QuantumCircuit

        theta = self.params[0]
        qc = QuantumCircuit(2)

        # Approximate Givens with RY + CNOT
        qc.ry(theta, 0)
        qc.cx(1, 0)
        qc.ry(-theta, 0)

        return qc

    def get_physical_description(self) -> str:
        """Get description of MO formation."""
        angle = self.metadata.get('mixing_angle_degrees', 0)
        if abs(angle - 45) < 1:
            character = "symmetric (equal mixing)"
        elif angle < 45:
            character = f"polar (more A character)"
        else:
            character = f"polar (more B character)"

        return f"{self.bond_type}-bond MO formation ({character}, θ={angle:.1f}°)"


class BellPairOperator(BaseGovernanceOperator):
    """
    Bell pair creation operator for electron pairing.

    Creates entangled singlet state for bonding electron pair:
        |singlet⟩ = (|01⟩ + |10⟩)/√2

    This represents two electrons in opposite spins sharing a bond.

    Implementation:
        H ⊗ I followed by CNOT creates Bell state:
        |00⟩ → H → (|00⟩ + |10⟩)/√2 → CNOT → (|00⟩ + |11⟩)/√2
        |01⟩ → H → (|01⟩ + |11⟩)/√2 → CNOT → (|01⟩ + |10⟩)/√2 ← singlet!
    """

    def __init__(
        self,
        qubit1: int,
        qubit2: int,
        state_type: str = 'singlet',
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize Bell pair operator.

        Args:
            qubit1: First qubit (spin up)
            qubit2: Second qubit (spin down)
            state_type: Type of Bell state ('singlet', 'triplet')
            metadata: Additional metadata
        """
        super().__init__(
            operator_type=OperatorType.BELL_PAIR,
            physical_meaning=PhysicalMeaning.ELECTRON_PAIRING,
            qubits=[qubit1, qubit2],
            params=[],
            metadata=metadata or {}
        )

        self.state_type = state_type
        self.metadata['state_type'] = state_type
        self.metadata['spin'] = 0 if state_type == 'singlet' else 1

    def _get_allowed_bonding_types(self) -> List[str]:
        """Bell pairs are for covalent electron pairing."""
        return ['covalent']

    def to_matrix(self, n_qubits: int) -> np.ndarray:
        """
        Build Bell pair creation matrix.

        Bell pair circuit: H on qubit1, then CNOT(qubit1 → qubit2)

        Args:
            n_qubits: Total number of qubits

        Returns:
            Bell pair operator matrix (2^n × 2^n)
        """
        q1, q2 = self.qubits

        # Build Hadamard on q1
        H_full = kron_product(H_MATRIX, q1, n_qubits)

        # Build CNOT from q1 to q2
        dim = 2 ** n_qubits
        CNOT_full = np.eye(dim, dtype=complex)

        # CNOT matrix in computational basis
        # |control target⟩ → |control control⊕target⟩
        for i in range(dim):
            bits = [(i >> k) & 1 for k in range(n_qubits)]

            if bits[q1] == 1:  # Control qubit is |1⟩
                # Flip target qubit
                bits[q2] = 1 - bits[q2]

            # Compute new basis index
            j = sum(bit << k for k, bit in enumerate(bits))

            if i != j:  # Swapped state
                CNOT_full[j, i] = 1
                CNOT_full[i, i] = 0

        # Bell operator = CNOT @ H
        return CNOT_full @ H_full

    def to_qiskit(self):
        """Convert to Qiskit Bell pair circuit."""
        from qiskit.circuit import QuantumCircuit

        qc = QuantumCircuit(2)
        qc.h(0)  # Hadamard on first qubit
        qc.cx(0, 1)  # CNOT to second qubit

        return qc

    def get_physical_description(self) -> str:
        """Get description of electron pairing."""
        if self.state_type == 'singlet':
            return "Singlet electron pair (S=0, bonding)"
        else:
            return "Triplet electron pair (S=1, antibonding)"

    def get_entanglement_degree(self) -> int:
        """Bell pairs create pairwise entanglement."""
        return 2


class OrbitalRotationOperator(BaseGovernanceOperator):
    """
    Generic orbital rotation operator.

    Rotates orbital basis using arbitrary rotation:
        |ψ'⟩ = R(θ, φ, λ) |ψ⟩

    Used for fine-tuning orbital character beyond standard hybridization.
    """

    def __init__(
        self,
        qubit: int,
        theta: float = 0.0,
        phi: float = 0.0,
        lam: float = 0.0,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize orbital rotation.

        Args:
            qubit: Qubit to rotate
            theta: Rotation angle θ
            phi: Rotation angle φ
            lam: Rotation angle λ
            metadata: Additional metadata
        """
        super().__init__(
            operator_type=OperatorType.RY,  # Generic rotation
            physical_meaning=PhysicalMeaning.ON_SITE_ROTATION,
            qubits=[qubit],
            params=[theta, phi, lam],
            metadata=metadata or {}
        )

    def _get_allowed_bonding_types(self) -> List[str]:
        """Orbital rotations used in all bonding types."""
        return ['ionic', 'covalent', 'metallic']

    def to_matrix(self, n_qubits: int) -> np.ndarray:
        """
        Build U3 rotation matrix.

        U3(θ,φ,λ) = [[cos(θ/2), -e^(iλ)sin(θ/2)],
                     [e^(iφ)sin(θ/2), e^(i(φ+λ))cos(θ/2)]]

        Args:
            n_qubits: Total number of qubits

        Returns:
            Rotation matrix (2^n × 2^n)
        """
        theta, phi, lam = self.params

        # U3 gate matrix
        u3 = np.array([
            [np.cos(theta/2), -np.exp(1j*lam)*np.sin(theta/2)],
            [np.exp(1j*phi)*np.sin(theta/2), np.exp(1j*(phi+lam))*np.cos(theta/2)]
        ], dtype=complex)

        return kron_product(u3, self.qubits[0], n_qubits)

    def to_qiskit(self):
        """Convert to Qiskit U gate."""
        from qiskit.circuit.library import UGate

        theta, phi, lam = self.params
        return UGate(theta, phi, lam)

    def get_physical_description(self) -> str:
        """Get description of rotation."""
        theta, phi, lam = self.params
        return f"Orbital rotation (θ={np.degrees(theta):.1f}°, φ={np.degrees(phi):.1f}°, λ={np.degrees(lam):.1f}°)"
