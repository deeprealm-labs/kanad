"""
Transfer operators for ionic bonding.

These operators implement electron transfer between localized atomic sites,
which is the key process in ionic bonding:
- Electron hopping: a†_i a_j
- On-site repulsion: U n_i↑ n_i↓
- Charge localization: minimal entanglement
"""

from typing import List, Optional, Dict, Any
import numpy as np
import logging

logger = logging.getLogger(__name__)

from kanad.governance.operators.base_operator import (
    BaseGovernanceOperator,
    OperatorType,
    PhysicalMeaning,
    kron_product,
    two_qubit_operator_matrix,
    I_MATRIX,
    X_MATRIX,
    Y_MATRIX,
    Z_MATRIX
)


class ElectronTransferOperator(BaseGovernanceOperator):
    """
    Electron transfer operator for ionic bonding.

    Implements: t_ij (a†_i a_j + a†_j a_i)

    Physical meaning:
        Electron hops from site j to site i with amplitude t_ij

    In Jordan-Wigner representation:
        a†_i a_j = 1/4 [(X_i - iY_i) Z_{i+1}...Z_{j-1} (X_j + iY_j)]

    Simplification for nearest neighbors (|i-j| = 1):
        ≈ RXX(t) or parametrized X⊗X rotation
    """

    def __init__(
        self,
        site_i: int,
        site_j: int,
        transfer_amplitude: float,
        distance: float,
        enforce_locality: bool = True,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize electron transfer operator.

        Args:
            site_i: Source site (qubit index)
            site_j: Target site (qubit index)
            transfer_amplitude: Transfer integral t_ij (in Hartree)
            distance: Physical distance between sites (Angstrom)
            enforce_locality: If True, reject long-range transfers
            metadata: Additional metadata
        """
        self.distance = distance
        self.enforce_locality = enforce_locality

        # Check locality constraint for ionic bonding
        if enforce_locality and abs(site_i - site_j) > 1:
            raise ValueError(
                f"Ionic bonding requires nearest-neighbor transfer only. "
                f"Distance between qubits {site_i} and {site_j} is {abs(site_i - site_j)}"
            )

        super().__init__(
            operator_type=OperatorType.TRANSFER,
            physical_meaning=PhysicalMeaning.ELECTRON_TRANSFER,
            qubits=[site_i, site_j],
            params=[transfer_amplitude],
            metadata=metadata or {}
        )

        self.metadata['distance_angstrom'] = distance
        self.metadata['transfer_amplitude_hartree'] = transfer_amplitude
        self.metadata['transfer_amplitude_ev'] = transfer_amplitude * 27.211
        self.metadata['locality_enforced'] = enforce_locality

    def _get_allowed_bonding_types(self) -> List[str]:
        """Transfer operators are for ionic bonding."""
        return ['ionic']

    def to_matrix(self, n_qubits: int) -> np.ndarray:
        """
        Build transfer operator matrix using Jordan-Wigner transformation.

        For nearest neighbors, use RXX approximation:
            exp(-i*t*X⊗X) ≈ transfer operator

        Args:
            n_qubits: Total number of qubits

        Returns:
            Transfer operator matrix (2^n × 2^n)
        """
        t = self.params[0]  # Transfer amplitude
        q1, q2 = self.qubits

        # For nearest neighbors, use RXX(2*t) as approximation
        # RXX(θ) = exp(-i*θ/2 * X⊗X)
        #        = cos(θ/2)*I⊗I - i*sin(θ/2)*X⊗X

        theta = 2 * t  # Factor of 2 from definition

        # Build X⊗X operator
        XX = two_qubit_operator_matrix(X_MATRIX, X_MATRIX, q1, q2, n_qubits)

        # Build identity
        dim = 2 ** n_qubits
        II = np.eye(dim, dtype=complex)

        # RXX(θ) = cos(θ/2)*I - i*sin(θ/2)*XX
        rxx = np.cos(theta/2) * II - 1j * np.sin(theta/2) * XX

        return rxx

    def to_qiskit(self):
        """Convert to Qiskit RXX gate."""
        from qiskit.circuit.library import RXXGate

        t = self.params[0]
        theta = 2 * t

        return RXXGate(theta)

    def get_physical_description(self) -> str:
        """Get description of electron transfer."""
        t_ev = self.metadata.get('transfer_amplitude_ev', 0)
        d = self.metadata.get('distance_angstrom', 0)

        return f"Electron transfer (t={t_ev:.2f} eV, d={d:.2f} Å)"

    def get_entanglement_degree(self) -> int:
        """Transfer creates pairwise entanglement."""
        return 2

    def is_local(self, max_distance: int = 1) -> bool:
        """Check if transfer is local (nearest-neighbor)."""
        return abs(self.qubits[0] - self.qubits[1]) <= max_distance


class HubbardInteractionOperator(BaseGovernanceOperator):
    """
    Hubbard U on-site repulsion operator.

    Implements: U n_i↑ n_i↓

    Physical meaning:
        Energy penalty for double occupancy (two electrons on same site)

    In spin-blocked ordering [0↑, 1↑, ..., 0↓, 1↓, ...]:
        n_i↑ n_i↓ = ((I - Z_i)/2) ⊗ ((I - Z_{i+N})/2)
                   = (I - Z_i - Z_{i+N} + Z_i Z_{i+N}) / 4

    This operator is diagonal in computational basis.
    """

    def __init__(
        self,
        site: int,
        hubbard_u: float,
        n_orbitals: int,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize Hubbard interaction operator.

        Args:
            site: Orbital index (0 to n_orbitals-1)
            hubbard_u: Hubbard U parameter (in Hartree)
            n_orbitals: Total number of spatial orbitals
            metadata: Additional metadata
        """
        self.n_orbitals = n_orbitals
        self.site = site

        # In blocked spin ordering: alpha qubits [0, n_orb), beta qubits [n_orb, 2*n_orb)
        qubit_alpha = site
        qubit_beta = site + n_orbitals

        super().__init__(
            operator_type=OperatorType.RZZ,  # Diagonal operator
            physical_meaning=PhysicalMeaning.ON_SITE_ROTATION,
            qubits=[qubit_alpha, qubit_beta],
            params=[hubbard_u],
            metadata=metadata or {}
        )

        self.metadata['hubbard_u_hartree'] = hubbard_u
        self.metadata['hubbard_u_ev'] = hubbard_u * 27.211
        self.metadata['site_index'] = site

    def _get_allowed_bonding_types(self) -> List[str]:
        """Hubbard U is for ionic/correlated systems."""
        return ['ionic', 'covalent']  # Also relevant for covalent with correlation

    def to_matrix(self, n_qubits: int) -> np.ndarray:
        """
        Build Hubbard U operator matrix.

        U n↑ n↓ = U (I - Z↑)(I - Z↓) / 4
               = U (I - Z↑ - Z↓ + Z↑Z↓) / 4

        This is a diagonal operator.

        Args:
            n_qubits: Total number of qubits

        Returns:
            Hubbard operator matrix (2^n × 2^n)
        """
        U = self.params[0]
        q_alpha, q_beta = self.qubits

        # Build Z operators
        Z_alpha = kron_product(Z_MATRIX, q_alpha, n_qubits)
        Z_beta = kron_product(Z_MATRIX, q_beta, n_qubits)

        # Build ZZ operator
        ZZ = two_qubit_operator_matrix(Z_MATRIX, Z_MATRIX, q_alpha, q_beta, n_qubits)

        # Identity
        dim = 2 ** n_qubits
        II = np.eye(dim, dtype=complex)

        # U n↑ n↓ = U (I - Z↑ - Z↓ + Z↑Z↓) / 4
        hubbard_op = U * (II - Z_alpha - Z_beta + ZZ) / 4

        return hubbard_op

    def to_qiskit(self):
        """
        Convert to Qiskit diagonal gate.

        Hubbard U is diagonal, so use RZZ approximation.
        """
        from qiskit.circuit.library import RZZGate

        U = self.params[0]

        # RZZ implements exp(-i θ/2 Z⊗Z)
        # For diagonal Hubbard, use U as rotation angle
        return RZZGate(U)

    def get_physical_description(self) -> str:
        """Get description of Hubbard interaction."""
        U_ev = self.metadata.get('hubbard_u_ev', 0)
        site = self.metadata.get('site_index', 0)

        return f"On-site Coulomb repulsion (U={U_ev:.1f} eV) at site {site}"

    def conserves_particle_number(self) -> bool:
        """Hubbard U conserves particle number (diagonal operator)."""
        return True


class LocalDensityOperator(BaseGovernanceOperator):
    """
    Local electron density operator n_i = a†_i a_i.

    Physical meaning:
        Number of electrons on site i

    In Jordan-Wigner:
        n_i = (I - Z_i) / 2

    Used for measuring charge distribution in ionic systems.
    """

    def __init__(
        self,
        site: int,
        spin: str = 'up',
        n_orbitals: Optional[int] = None,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize local density operator.

        Args:
            site: Orbital index
            spin: Spin ('up' or 'down')
            n_orbitals: Total spatial orbitals (for spin mapping)
            metadata: Additional metadata
        """
        # Map to qubit index
        if n_orbitals is not None and spin == 'down':
            qubit = site + n_orbitals  # Beta spin in blocked ordering
        else:
            qubit = site  # Alpha spin

        super().__init__(
            operator_type=OperatorType.RZ,  # Diagonal
            physical_meaning=PhysicalMeaning.ON_SITE_ROTATION,
            qubits=[qubit],
            params=[],
            metadata=metadata or {}
        )

        self.metadata['site'] = site
        self.metadata['spin'] = spin

    def _get_allowed_bonding_types(self) -> List[str]:
        """Density operators universal."""
        return ['ionic', 'covalent', 'metallic']

    def to_matrix(self, n_qubits: int) -> np.ndarray:
        """
        Build density operator matrix.

        n_i = (I - Z_i) / 2

        Args:
            n_qubits: Total number of qubits

        Returns:
            Density operator (2^n × 2^n)
        """
        qubit = self.qubits[0]

        # Build Z operator
        Z = kron_product(Z_MATRIX, qubit, n_qubits)

        # Identity
        dim = 2 ** n_qubits
        I = np.eye(dim, dtype=complex)

        # n = (I - Z) / 2
        return (I - Z) / 2

    def to_qiskit(self):
        """Convert to measurement operator."""
        from qiskit.circuit.library import ZGate

        # Density is (I - Z)/2, use Z gate with scaling
        return ZGate()

    def get_physical_description(self) -> str:
        """Get description of density."""
        site = self.metadata.get('site', 0)
        spin = self.metadata.get('spin', 'up')

        return f"Electron density at site {site} (spin {spin})"

    def conserves_particle_number(self) -> bool:
        """Density measurement conserves N (diagonal)."""
        return True


class NearestNeighborHoppingOperator(BaseGovernanceOperator):
    """
    Simplified nearest-neighbor hopping for tight-binding models.

    Implements: t (|01⟩⟨10| + |10⟩⟨01|)

    This is the simplest form of electron transfer between adjacent sites,
    commonly used in Hubbard-like models for ionic systems.

    More efficient than full Jordan-Wigner transfer operator.
    """

    def __init__(
        self,
        site_i: int,
        site_j: int,
        hopping_amplitude: float,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize hopping operator.

        Args:
            site_i: First site
            site_j: Second site (must be adjacent to site_i)
            hopping_amplitude: Hopping parameter t
            metadata: Additional metadata
        """
        if abs(site_i - site_j) != 1:
            raise ValueError(
                f"Nearest-neighbor hopping requires adjacent sites. "
                f"Got distance {abs(site_i - site_j)}"
            )

        super().__init__(
            operator_type=OperatorType.SWAP,  # Hopping is like partial SWAP
            physical_meaning=PhysicalMeaning.ELECTRON_TRANSFER,
            qubits=[site_i, site_j],
            params=[hopping_amplitude],
            metadata=metadata or {}
        )

        self.metadata['hopping_amplitude'] = hopping_amplitude

    def _get_allowed_bonding_types(self) -> List[str]:
        """Hopping for ionic tight-binding."""
        return ['ionic']

    def to_matrix(self, n_qubits: int) -> np.ndarray:
        """
        Build hopping operator.

        Hopping: t(|01⟩⟨10| + |10⟩⟨01|) = t(X⊗X + Y⊗Y)/2

        Args:
            n_qubits: Total number of qubits

        Returns:
            Hopping operator (2^n × 2^n)
        """
        t = self.params[0]
        q1, q2 = self.qubits

        # Build XX and YY
        XX = two_qubit_operator_matrix(X_MATRIX, X_MATRIX, q1, q2, n_qubits)
        YY = two_qubit_operator_matrix(Y_MATRIX, Y_MATRIX, q1, q2, n_qubits)

        # Hopping = t(XX + YY)/2
        return t * (XX + YY) / 2

    def to_qiskit(self):
        """Convert to Qiskit circuit."""
        from qiskit.circuit import QuantumCircuit

        t = self.params[0]
        qc = QuantumCircuit(2)

        # Decompose hopping into RXX + RYY
        qc.rxx(t, 0, 1)
        qc.ryy(t, 0, 1)

        return qc

    def get_physical_description(self) -> str:
        """Get description of hopping."""
        t = self.metadata.get('hopping_amplitude', 0)
        t_ev = t * 27.211

        return f"Nearest-neighbor hopping (t={t_ev:.2f} eV)"

    def is_local(self, max_distance: int = 1) -> bool:
        """Hopping is always local (nearest-neighbor by construction)."""
        return True
