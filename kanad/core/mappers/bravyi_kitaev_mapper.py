"""
Bravyi-Kitaev transformation for fermionic-to-qubit mapping.

More efficient than Jordan-Wigner for large systems with non-local interactions.
Uses OpenFermion's validated implementation for correctness.
"""

from typing import Dict, List
import numpy as np
from kanad.core.mappers.base_mapper import BaseMapper
import logging

logger = logging.getLogger(__name__)


class BravyiKitaevMapper(BaseMapper):
    """
    Bravyi-Kitaev transformation using OpenFermion.

    Uses a tree-based encoding that balances locality and operator weight:
    - Better than JW for long-range interactions
    - Requires O(log n) operator weight vs O(n) for JW
    - More complex encoding but more efficient gates

    Encoding idea:
    - Instead of sequential Z string (JW), use binary tree structure
    - Each qubit encodes partial occupation and parity information
    - Reduces Pauli weight of excitation operators

    APPLICATIONS:
    - General-purpose mapper
    - Good for molecules with many orbitals
    - Useful when inter-orbital distances vary

    NOTE: This implementation uses OpenFermion's validated BK transformation
    to ensure correctness. The BK algorithm is complex and error-prone when
    implemented from scratch.
    """

    def __init__(self):
        """Initialize Bravyi-Kitaev mapper."""
        self._n_qubits_cache = None
        self._check_openfermion()

    def _check_openfermion(self):
        """Check if OpenFermion is available."""
        try:
            import openfermion
            self._openfermion_available = True
        except ImportError:
            logger.warning("OpenFermion not available. BK mapper may not work correctly.")
            self._openfermion_available = False

    def n_qubits(self, n_spin_orbitals: int) -> int:
        """
        One qubit per spin orbital (same as JW).

        Args:
            n_spin_orbitals: Number of spin orbitals

        Returns:
            Number of qubits
        """
        return n_spin_orbitals

    def map_number_operator(self, orbital: int, n_orbitals: int) -> Dict[str, complex]:
        """
        Map number operator n_i = a†_i a_i.

        Uses OpenFermion's BK transformation for correctness.

        Args:
            orbital: Orbital index
            n_orbitals: Total orbitals

        Returns:
            Pauli operator dictionary
        """
        if not self._openfermion_available:
            raise ImportError("OpenFermion required for BK mapper. Install with: pip install openfermion")

        from openfermion import FermionOperator, bravyi_kitaev

        # Create number operator: a†_i a_i
        fermion_op = FermionOperator(f'{orbital}^ {orbital}')

        # Apply BK transformation
        qubit_op = bravyi_kitaev(fermion_op, n_qubits=n_orbitals)

        # Convert to our dictionary format
        return self._qubit_op_to_dict(qubit_op, n_orbitals)

    def map_excitation_operator(
        self,
        orbital_from: int,
        orbital_to: int,
        n_orbitals: int
    ) -> Dict[str, complex]:
        """
        Map excitation operator a†_i a_j.

        Uses OpenFermion's BK transformation for correctness.

        Args:
            orbital_from: Initial orbital (j)
            orbital_to: Final orbital (i)
            n_orbitals: Total orbitals

        Returns:
            Pauli operator dictionary
        """
        if not self._openfermion_available:
            raise ImportError("OpenFermion required for BK mapper. Install with: pip install openfermion")

        from openfermion import FermionOperator, bravyi_kitaev

        i, j = orbital_to, orbital_from

        if i == j:
            return self.map_number_operator(i, n_orbitals)

        # Create excitation operator: a†_i a_j
        fermion_op = FermionOperator(f'{i}^ {j}')

        # Apply BK transformation
        qubit_op = bravyi_kitaev(fermion_op, n_qubits=n_orbitals)

        # Convert to our dictionary format
        return self._qubit_op_to_dict(qubit_op, n_orbitals)

    def map_creation_operator(self, orbital: int, n_orbitals: int) -> Dict[str, complex]:
        """
        Map creation operator a†_i.

        Args:
            orbital: Orbital index
            n_orbitals: Total orbitals

        Returns:
            Pauli operator dictionary
        """
        if not self._openfermion_available:
            raise ImportError("OpenFermion required for BK mapper. Install with: pip install openfermion")

        from openfermion import FermionOperator, bravyi_kitaev

        # Create creation operator: a†_i
        fermion_op = FermionOperator(f'{orbital}^')

        # Apply BK transformation
        qubit_op = bravyi_kitaev(fermion_op, n_qubits=n_orbitals)

        # Convert to our dictionary format
        return self._qubit_op_to_dict(qubit_op, n_orbitals)

    def map_annihilation_operator(self, orbital: int, n_orbitals: int) -> Dict[str, complex]:
        """
        Map annihilation operator a_i.

        Args:
            orbital: Orbital index
            n_orbitals: Total orbitals

        Returns:
            Pauli operator dictionary
        """
        if not self._openfermion_available:
            raise ImportError("OpenFermion required for BK mapper. Install with: pip install openfermion")

        from openfermion import FermionOperator, bravyi_kitaev

        # Create annihilation operator: a_i
        fermion_op = FermionOperator(f'{orbital}')

        # Apply BK transformation
        qubit_op = bravyi_kitaev(fermion_op, n_qubits=n_orbitals)

        # Convert to our dictionary format
        return self._qubit_op_to_dict(qubit_op, n_orbitals)

    def map_double_excitation(
        self,
        orb_from_1: int,
        orb_from_2: int,
        orb_to_1: int,
        orb_to_2: int,
        n_orbitals: int
    ) -> Dict[str, complex]:
        """
        Map double excitation operator a†_i a†_k a_l a_j.

        Args:
            orb_from_1: j (first annihilation)
            orb_from_2: l (second annihilation)
            orb_to_1: i (first creation)
            orb_to_2: k (second creation)
            n_orbitals: Total orbitals

        Returns:
            Pauli operator dictionary
        """
        if not self._openfermion_available:
            raise ImportError("OpenFermion required for BK mapper. Install with: pip install openfermion")

        from openfermion import FermionOperator, bravyi_kitaev

        j, l, i, k = orb_from_1, orb_from_2, orb_to_1, orb_to_2

        # Create double excitation: a†_i a†_k a_l a_j
        fermion_op = FermionOperator(f'{i}^ {k}^ {l} {j}')

        # Apply BK transformation
        qubit_op = bravyi_kitaev(fermion_op, n_qubits=n_orbitals)

        # Convert to our dictionary format
        return self._qubit_op_to_dict(qubit_op, n_orbitals)

    def _qubit_op_to_dict(self, qubit_op, n_qubits: int) -> Dict[str, complex]:
        """
        Convert OpenFermion QubitOperator to our dictionary format.

        Args:
            qubit_op: OpenFermion QubitOperator
            n_qubits: Number of qubits

        Returns:
            Dictionary mapping Pauli strings to coefficients
        """
        result = {}

        for term, coeff in qubit_op.terms.items():
            # term is a tuple of (qubit_index, pauli_char) pairs
            # e.g., ((0, 'X'), (1, 'Y'))

            # Build Pauli string
            pauli_str = ['I'] * n_qubits
            for qubit_idx, pauli_char in term:
                pauli_str[qubit_idx] = pauli_char

            pauli_string = ''.join(pauli_str)

            if pauli_string in result:
                result[pauli_string] += coeff
            else:
                result[pauli_string] = coeff

        # Filter negligible terms
        return {k: v for k, v in result.items() if abs(v) > 1e-12}

    def get_operator_weight(self, orbital_from: int, orbital_to: int, n_orbitals: int) -> int:
        """
        Get the weight (number of non-I Paulis) for an excitation operator.

        BK: O(log n) weight (theoretically)
        JW: O(n) weight

        Args:
            orbital_from: Initial orbital
            orbital_to: Final orbital
            n_orbitals: Total orbitals

        Returns:
            Operator weight (actual, not theoretical)
        """
        pauli_dict = self.map_excitation_operator(orbital_from, orbital_to, n_orbitals)

        # Find the maximum weight across all Pauli strings
        max_weight = 0
        for pauli_str in pauli_dict.keys():
            weight = sum(1 for p in pauli_str if p != 'I')
            max_weight = max(max_weight, weight)

        return max_weight

    def __repr__(self) -> str:
        """String representation."""
        status = "available" if self._openfermion_available else "unavailable (OpenFermion not installed)"
        return f"BravyiKitaevMapper(weight='O(log n)', status='{status}')"
