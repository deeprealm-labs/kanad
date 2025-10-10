"""
Base classes for governance-aware quantum operators.

This module provides the foundation for physics-driven operator selection
in quantum circuits. Each operator carries metadata about its physical meaning
and can be validated against governance protocols.
"""

from abc import ABC, abstractmethod
from typing import List, Optional, Dict, Any, Tuple
from enum import Enum
import numpy as np
import logging

logger = logging.getLogger(__name__)


class OperatorType(Enum):
    """Types of quantum operators in chemistry."""
    # Single-qubit rotations
    RX = "rx"
    RY = "ry"
    RZ = "rz"
    H = "h"
    S = "s"
    T = "t"

    # Two-qubit gates
    CNOT = "cx"
    CZ = "cz"
    CY = "cy"
    SWAP = "swap"

    # Parameterized two-qubit
    RXX = "rxx"
    RYY = "ryy"
    RZZ = "rzz"

    # Chemistry-specific
    GIVENS = "givens"  # MO formation
    HYBRIDIZATION = "hybridization"  # Orbital hybridization
    TRANSFER = "transfer"  # Electron transfer
    BELL_PAIR = "bell_pair"  # Electron pairing

    # Collective (metallic)
    QFT = "qft"
    GHZ = "ghz"
    BLOCH_TRANSFORM = "bloch"


class PhysicalMeaning(Enum):
    """Physical interpretation of operators."""
    HYBRIDIZATION = "orbital_hybridization"  # sp, sp2, sp3
    MO_FORMATION = "molecular_orbital_formation"  # bonding/antibonding
    ELECTRON_PAIRING = "electron_pair_formation"  # singlet state
    ELECTRON_TRANSFER = "electron_transfer"  # hopping
    ON_SITE_ROTATION = "on_site_rotation"  # local energy
    SPIN_ROTATION = "spin_rotation"  # spin state
    ENTANGLEMENT = "entanglement_creation"  # correlation
    DELOCALIZATION = "delocalization"  # metallic

class BaseGovernanceOperator(ABC):
    """
    Abstract base class for governance-aware quantum operators.

    Each operator knows:
    1. Its quantum mechanical action (unitary matrix)
    2. Its physical meaning (what chemistry it represents)
    3. Which bonding types allow it
    4. How to validate against governance rules
    """

    def __init__(
        self,
        operator_type: OperatorType,
        physical_meaning: PhysicalMeaning,
        qubits: List[int],
        params: Optional[List[float]] = None,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize governance operator.

        Args:
            operator_type: Type of quantum operator
            physical_meaning: Physical interpretation
            qubits: List of qubit indices this operator acts on
            params: Optional parameters (angles, amplitudes)
            metadata: Additional metadata about physical context
        """
        self.operator_type = operator_type
        self.physical_meaning = physical_meaning
        self.qubits = qubits
        self.params = params or []
        self.metadata = metadata or {}

        # Bonding types that allow this operator
        self.allowed_bonding_types = self._get_allowed_bonding_types()

    @abstractmethod
    def _get_allowed_bonding_types(self) -> List[str]:
        """
        Get list of bonding types that allow this operator.

        Returns:
            List of bonding type names ('ionic', 'covalent', 'metallic')
        """
        pass

    @abstractmethod
    def to_matrix(self, n_qubits: int) -> np.ndarray:
        """
        Convert operator to matrix representation.

        Args:
            n_qubits: Total number of qubits in system

        Returns:
            Unitary matrix (2^n × 2^n)
        """
        pass

    @abstractmethod
    def to_qiskit(self):
        """
        Convert to Qiskit gate representation.

        Returns:
            Qiskit gate object
        """
        pass

    def validate_bonding_type(self, bonding_type: str) -> bool:
        """
        Check if operator is allowed for given bonding type.

        Args:
            bonding_type: Bonding type ('ionic', 'covalent', 'metallic')

        Returns:
            True if operator is valid for this bonding type
        """
        return bonding_type in self.allowed_bonding_types

    def get_physical_description(self) -> str:
        """
        Get human-readable description of operator's physical meaning.

        Returns:
            Description string
        """
        descriptions = {
            PhysicalMeaning.HYBRIDIZATION: "Orbital hybridization (sp/sp²/sp³)",
            PhysicalMeaning.MO_FORMATION: "Molecular orbital formation (bonding/antibonding)",
            PhysicalMeaning.ELECTRON_PAIRING: "Electron pair formation (singlet state)",
            PhysicalMeaning.ELECTRON_TRANSFER: "Electron transfer between sites",
            PhysicalMeaning.ON_SITE_ROTATION: "On-site energy rotation",
            PhysicalMeaning.SPIN_ROTATION: "Spin state rotation",
            PhysicalMeaning.ENTANGLEMENT: "Quantum correlation creation",
            PhysicalMeaning.DELOCALIZATION: "Electron delocalization (metallic)"
        }
        return descriptions.get(self.physical_meaning, "Unknown")

    def get_entanglement_degree(self) -> int:
        """
        Get number of qubits this operator entangles.

        Returns:
            0: No entanglement (single-qubit)
            2: Pairwise entanglement
            >2: Multi-qubit entanglement
        """
        if len(self.qubits) == 1:
            return 0
        elif len(self.qubits) == 2:
            # Check if gate actually entangles (e.g., CNOT does, but ZZ⊗I doesn't)
            entangling_types = [
                OperatorType.CNOT, OperatorType.CZ, OperatorType.CY,
                OperatorType.GIVENS, OperatorType.BELL_PAIR,
                OperatorType.RXX, OperatorType.RYY, OperatorType.RZZ
            ]
            return 2 if self.operator_type in entangling_types else 0
        else:
            return len(self.qubits)

    def is_local(self, max_distance: int = 1) -> bool:
        """
        Check if operator acts on neighboring qubits only.

        Args:
            max_distance: Maximum allowed qubit distance

        Returns:
            True if operator is local
        """
        if len(self.qubits) <= 1:
            return True

        # Check pairwise distances
        for i in range(len(self.qubits)):
            for j in range(i + 1, len(self.qubits)):
                if abs(self.qubits[i] - self.qubits[j]) > max_distance:
                    return False

        return True

    def conserves_particle_number(self) -> bool:
        """
        Check if operator conserves total particle number.

        All unitary operators in quantum chemistry should conserve N,
        but this method can be overridden for specific checks.

        Returns:
            True if operator conserves particle number
        """
        # By default, assume properly constructed operators conserve N
        # Subclasses can override for specific validation
        return True

    def __repr__(self) -> str:
        """String representation."""
        return (f"{self.__class__.__name__}("
                f"type={self.operator_type.value}, "
                f"qubits={self.qubits}, "
                f"meaning={self.physical_meaning.value})")


class QuantumOperator:
    """
    Simple quantum operator without governance (for backwards compatibility).

    Used in validation: protocol.validate_operator(QuantumOperator(...))
    """

    def __init__(self, operator_type: str, qubits: List[int], params: Optional[List[float]] = None):
        """
        Initialize basic quantum operator.

        Args:
            operator_type: Type of operator (string)
            qubits: Qubit indices
            params: Optional parameters
        """
        self.type = operator_type
        self.qubits = qubits
        self.params = params or []

    def __repr__(self) -> str:
        return f"QuantumOperator(type={self.type}, qubits={self.qubits})"


def kron_product(operator: np.ndarray, position: int, n_qubits: int) -> np.ndarray:
    """
    Create tensor product with single-qubit operator at specified position.

    Helper function for building full operator matrices.

    Args:
        operator: 2x2 operator matrix
        position: Qubit position (0-indexed)
        n_qubits: Total number of qubits

    Returns:
        Full operator matrix (2^n × 2^n)
    """
    I = np.eye(2, dtype=complex)

    # Build tensor product: I ⊗ I ⊗ ... ⊗ operator ⊗ ... ⊗ I
    # IMPORTANT: Build from left to right (MSB to LSB)
    result = np.array([[1.0]], dtype=complex)

    for i in range(n_qubits - 1, -1, -1):  # Reverse order
        if i == position:
            result = np.kron(result, operator)
        else:
            result = np.kron(result, I)

    return result


def two_qubit_operator_matrix(
    op1: np.ndarray,
    op2: np.ndarray,
    qubit1: int,
    qubit2: int,
    n_qubits: int
) -> np.ndarray:
    """
    Create matrix for two-qubit operator op1 ⊗ op2.

    Args:
        op1: Operator on first qubit (2x2)
        op2: Operator on second qubit (2x2)
        qubit1: First qubit index
        qubit2: Second qubit index
        n_qubits: Total number of qubits

    Returns:
        Full operator matrix (2^n × 2^n)
    """
    I = np.eye(2, dtype=complex)

    # Build tensor product
    result = np.array([[1.0]], dtype=complex)

    for i in range(n_qubits - 1, -1, -1):
        if i == qubit1:
            result = np.kron(result, op1)
        elif i == qubit2:
            result = np.kron(result, op2)
        else:
            result = np.kron(result, I)

    return result


# Pauli matrices for common use
I_MATRIX = np.array([[1, 0], [0, 1]], dtype=complex)
X_MATRIX = np.array([[0, 1], [1, 0]], dtype=complex)
Y_MATRIX = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z_MATRIX = np.array([[1, 0], [0, -1]], dtype=complex)
H_MATRIX = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S_MATRIX = np.array([[1, 0], [0, 1j]], dtype=complex)
T_MATRIX = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)
