"""
Jordan-Wigner transformation for fermionic-to-qubit mapping.

Best suited for IONIC bonding where electrons are localized.
"""

from typing import Dict, Tuple, List
import numpy as np
from kanad.core.mappers.base_mapper import BaseMapper


class JordanWignerMapper(BaseMapper):
    """
    Jordan-Wigner transformation.

    Maps fermionic operators to qubits in a LOCAL, SEQUENTIAL manner:
        |0⟩ → empty orbital
        |1⟩ → occupied orbital

    Fermionic operators:
        a†_j = (X_j - iY_j)/2 ⊗ Z_0 ⊗ Z_1 ⊗ ... ⊗ Z_{j-1}
        a_j  = (X_j + iY_j)/2 ⊗ Z_0 ⊗ Z_1 ⊗ ... ⊗ Z_{j-1}

    The Z string enforces fermionic anticommutation relations.

    BEST FOR IONIC BONDING:
    - Electrons are localized on atoms
    - Site-to-site ordering is natural
    - No need for collective transformations
    """

    def n_qubits(self, n_spin_orbitals: int) -> int:
        """
        One qubit per spin orbital.

        Args:
            n_spin_orbitals: Number of spin orbitals

        Returns:
            Number of qubits = n_spin_orbitals
        """
        return n_spin_orbitals

    def map_number_operator(self, orbital: int, n_orbitals: int) -> Dict[str, complex]:
        """
        Map number operator n_i = a†_i a_i to Pauli operators.

        n_i → (I - Z_i) / 2

        Args:
            orbital: Orbital index
            n_orbitals: Total number of orbitals

        Returns:
            Pauli operator dictionary
        """
        # Build Pauli string
        pauli_I = 'I' * n_orbitals
        pauli_Z = list('I' * n_orbitals)
        pauli_Z[orbital] = 'Z'
        pauli_Z = ''.join(pauli_Z)

        return {
            pauli_I: 0.5,
            pauli_Z: -0.5
        }

    def map_excitation_operator(
        self,
        orbital_from: int,
        orbital_to: int,
        n_orbitals: int
    ) -> Dict[str, complex]:
        """
        Map excitation operator a†_i a_j to Pauli operators.

        a†_i a_j → ½[(X_i X_j + Y_i Y_j) ⊗ Z_string] + i/2[(X_i Y_j - Y_i X_j) ⊗ Z_string]

        where Z_string = Z_{j+1} ⊗ Z_{j+2} ⊗ ... ⊗ Z_{i-1} if j < i

        Args:
            orbital_from: Initial orbital (j)
            orbital_to: Final orbital (i)
            n_orbitals: Total number of orbitals

        Returns:
            Pauli operator dictionary
        """
        i, j = orbital_to, orbital_from

        if i == j:
            # Number operator
            return self.map_number_operator(i, n_orbitals)

        # Build Z string for anticommutation
        z_string = self._build_z_string(j, i, n_orbitals)

        # Build Pauli strings for excitation
        xx_string = self._build_pauli_string('X', 'X', j, i, z_string, n_orbitals)
        yy_string = self._build_pauli_string('Y', 'Y', j, i, z_string, n_orbitals)
        xy_string = self._build_pauli_string('X', 'Y', j, i, z_string, n_orbitals)
        yx_string = self._build_pauli_string('Y', 'X', j, i, z_string, n_orbitals)

        return {
            xx_string: 0.25,
            yy_string: 0.25,
            xy_string: 0.25j,
            yx_string: -0.25j
        }

    def _build_z_string(self, j: int, i: int, n_orbitals: int) -> List[int]:
        """
        Build list of qubit indices for Z string.

        For j < i: apply Z on qubits (j+1, j+2, ..., i-1)
        For i < j: apply Z on qubits (i+1, i+2, ..., j-1)

        Args:
            j: First orbital
            i: Second orbital
            n_orbitals: Total orbitals

        Returns:
            List of indices where Z should be applied
        """
        if i > j:
            return list(range(j + 1, i))
        elif j > i:
            return list(range(i + 1, j))
        else:
            return []

    def _build_pauli_string(
        self,
        pauli_j: str,
        pauli_i: str,
        j: int,
        i: int,
        z_indices: List[int],
        n_orbitals: int
    ) -> str:
        """
        Build Pauli string for two-site operator with Z string.

        Args:
            pauli_j: Pauli operator at site j ('X', 'Y', or 'Z')
            pauli_i: Pauli operator at site i
            j: First site
            i: Second site
            z_indices: Indices for Z operators
            n_orbitals: Total number of orbitals

        Returns:
            Pauli string
        """
        pauli = ['I'] * n_orbitals
        pauli[j] = pauli_j
        pauli[i] = pauli_i

        for idx in z_indices:
            pauli[idx] = 'Z'

        return ''.join(pauli)

    def map_double_excitation(
        self,
        orb_from_1: int,
        orb_from_2: int,
        orb_to_1: int,
        orb_to_2: int,
        n_orbitals: int
    ) -> Dict[str, complex]:
        """
        Map double excitation operator a†_i a†_j a_l a_k.

        This is used for two-electron integrals in Hamiltonians.

        Args:
            orb_from_1: First annihilation orbital (k)
            orb_from_2: Second annihilation orbital (l)
            orb_to_1: First creation orbital (i)
            orb_to_2: Second creation orbital (j)
            n_orbitals: Total orbitals

        Returns:
            Pauli operator dictionary
        """
        # This is complex - for now use product of single excitations
        # More efficient implementations exist
        exc1 = self.map_excitation_operator(orb_from_1, orb_to_1, n_orbitals)
        exc2 = self.map_excitation_operator(orb_from_2, orb_to_2, n_orbitals)

        return self.pauli_string_multiply(exc1, exc2)

    def __repr__(self) -> str:
        """String representation."""
        return "JordanWignerMapper(locality='sequential', best_for='ionic_bonding')"
