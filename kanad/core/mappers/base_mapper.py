"""
Base class for fermionic-to-qubit mappers.

Different bonding types use different mapping strategies:
- Ionic: Jordan-Wigner (local, sequential)
- Covalent: Paired mapping (orbital-centric)
- Metallic: Momentum-space mapping (collective)
"""

from abc import ABC, abstractmethod
from typing import Dict, Tuple, List
import numpy as np


class BaseMapper(ABC):
    """
    Abstract base class for fermionic-to-qubit mappers.

    Transforms fermionic operators (creation/annihilation) to
    qubit operators (Pauli matrices).
    """

    @abstractmethod
    def n_qubits(self, n_spin_orbitals: int) -> int:
        """
        Compute number of qubits needed for n spin orbitals.

        Args:
            n_spin_orbitals: Number of spin orbitals

        Returns:
            Number of qubits required
        """
        pass

    @abstractmethod
    def map_number_operator(self, orbital: int, n_orbitals: int) -> Dict[str, complex]:
        """
        Map fermionic number operator n_i = a†_i a_i to Pauli operators.

        Args:
            orbital: Orbital index
            n_orbitals: Total number of orbitals

        Returns:
            Dictionary mapping Pauli strings to coefficients
        """
        pass

    @abstractmethod
    def map_excitation_operator(
        self,
        orbital_from: int,
        orbital_to: int,
        n_orbitals: int
    ) -> Dict[str, complex]:
        """
        Map fermionic excitation operator a†_i a_j to Pauli operators.

        Args:
            orbital_from: Initial orbital
            orbital_to: Final orbital
            n_orbitals: Total number of orbitals

        Returns:
            Dictionary mapping Pauli strings to coefficients
        """
        pass

    def map_hamiltonian_term(
        self,
        indices: Tuple[int, ...],
        coeff: complex,
        n_orbitals: int
    ) -> Dict[str, complex]:
        """
        Map a general Hamiltonian term to Pauli operators.

        Args:
            indices: Tuple of orbital indices
            coeff: Coefficient
            n_orbitals: Total number of orbitals

        Returns:
            Dictionary mapping Pauli strings to coefficients
        """
        if len(indices) == 2:
            # One-body term: h_ij a†_i a_j
            i, j = indices
            if i == j:
                # Number operator
                return {k: coeff * v for k, v in self.map_number_operator(i, n_orbitals).items()}
            else:
                # Excitation operator
                return {k: coeff * v for k, v in self.map_excitation_operator(j, i, n_orbitals).items()}
        elif len(indices) == 4:
            # Two-body term: g_ijkl a†_i a†_j a_l a_k
            # This requires combining multiple Pauli strings
            return self._map_two_body_term(indices, coeff, n_orbitals)
        else:
            raise ValueError(f"Unsupported number of indices: {len(indices)}")

    def _map_two_body_term(
        self,
        indices: Tuple[int, int, int, int],
        coeff: complex,
        n_orbitals: int
    ) -> Dict[str, complex]:
        """
        Map two-body term (default implementation).

        Subclasses can override for more efficient implementations.
        """
        # This is complex and depends on the specific mapping
        # For now, return empty dict (will be implemented by subclasses)
        return {}

    def pauli_string_multiply(
        self,
        pauli1: Dict[str, complex],
        pauli2: Dict[str, complex]
    ) -> Dict[str, complex]:
        """
        Multiply two Pauli operator dictionaries.

        Args:
            pauli1: First Pauli operator
            pauli2: Second Pauli operator

        Returns:
            Product of operators
        """
        result = {}

        for p1, c1 in pauli1.items():
            for p2, c2 in pauli2.items():
                # Multiply Pauli strings
                p_prod, phase = self._multiply_pauli_strings(p1, p2)
                coeff = c1 * c2 * phase

                if p_prod in result:
                    result[p_prod] += coeff
                else:
                    result[p_prod] = coeff

        # Remove near-zero terms
        return {k: v for k, v in result.items() if abs(v) > 1e-12}

    def _multiply_pauli_strings(self, p1: str, p2: str) -> Tuple[str, complex]:
        """
        Multiply two Pauli strings.

        Args:
            p1: First Pauli string (e.g., 'XIYZ')
            p2: Second Pauli string

        Returns:
            (product_string, phase)
        """
        if len(p1) != len(p2):
            raise ValueError("Pauli strings must have same length")

        result = []
        phase = 1.0

        for a, b in zip(p1, p2):
            prod, ph = self._multiply_single_pauli(a, b)
            result.append(prod)
            phase *= ph

        return ''.join(result), phase

    def _multiply_single_pauli(self, a: str, b: str) -> Tuple[str, complex]:
        """
        Multiply two single-qubit Pauli operators.

        Multiplication table:
            I * X = X,  X * I = X,  X * X = I
            I * Y = Y,  Y * I = Y,  Y * Y = I
            I * Z = Z,  Z * I = Z,  Z * Z = I
            X * Y = iZ, Y * X = -iZ
            Y * Z = iX, Z * Y = -iX
            Z * X = iY, X * Z = -iY

        Returns:
            (product, phase)
        """
        if a == 'I':
            return b, 1.0
        if b == 'I':
            return a, 1.0
        if a == b:
            return 'I', 1.0

        # X * Y = iZ
        if (a, b) == ('X', 'Y'):
            return 'Z', 1j
        if (a, b) == ('Y', 'X'):
            return 'Z', -1j

        # Y * Z = iX
        if (a, b) == ('Y', 'Z'):
            return 'X', 1j
        if (a, b) == ('Z', 'Y'):
            return 'X', -1j

        # Z * X = iY
        if (a, b) == ('Z', 'X'):
            return 'Y', 1j
        if (a, b) == ('X', 'Z'):
            return 'Y', -1j

        raise ValueError(f"Invalid Pauli operators: {a}, {b}")

    def __repr__(self) -> str:
        """String representation."""
        return f"{self.__class__.__name__}()"
