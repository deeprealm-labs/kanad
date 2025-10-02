"""
Hybrid Orbital Mapper for covalent bonding.

Maps bonding/antibonding MO PAIRS to qubits instead of individual orbitals.
"""

from typing import Dict, Tuple, List
import numpy as np
from kanad.core.mappers.base_mapper import BaseMapper


class HybridOrbitalMapper(BaseMapper):
    """
    Custom mapper for covalent bonding with hybrid orbitals.

    KEY INNOVATION: Maps bonding/antibonding PAIRS to qubits

    Instead of Jordan-Wigner's one-to-one mapping:
        |orbital_i⟩ → |qubit_i⟩

    We use molecular orbital pairs:
        (|bonding⟩, |antibonding⟩) → (|qubit_2i⟩, |qubit_{2i+1}⟩)

    Encoding for each MO pair:
        |00⟩ → both empty
        |01⟩ → antibonding occupied
        |10⟩ → bonding occupied (ground state for filled bond)
        |11⟩ → both occupied (excited state)

    BEST FOR COVALENT BONDING:
    - Natural encoding of bonding/antibonding character
    - Excitations are local within MO pairs
    - Hybridization structure is explicit
    """

    def __init__(self, mo_pairs: List[Tuple[int, int]]):
        """
        Initialize with molecular orbital pairs.

        Args:
            mo_pairs: List of (bonding_orbital_idx, antibonding_orbital_idx) tuples
        """
        self.mo_pairs = mo_pairs
        self.n_pairs = len(mo_pairs)

        # Build orbital → pair mapping
        self.orbital_to_pair = {}
        for pair_idx, (bonding, antibonding) in enumerate(mo_pairs):
            self.orbital_to_pair[bonding] = (pair_idx, 'bonding')
            self.orbital_to_pair[antibonding] = (pair_idx, 'antibonding')

    def n_qubits(self, n_spin_orbitals: int) -> int:
        """
        Two qubits per MO pair.

        For N hybrid orbitals forming N/2 bonds → N qubits total

        Args:
            n_spin_orbitals: Number of spin orbitals

        Returns:
            Number of qubits
        """
        # Each pair needs 2 qubits
        return 2 * self.n_pairs

    def map_number_operator(self, orbital: int, n_orbitals: int) -> Dict[str, complex]:
        """
        Map number operator n_i = a†_i a_i.

        For bonding orbital:    n_b → (I - Z_{2p}) / 2
        For antibonding orbital: n_a → (I - Z_{2p+1}) / 2

        where p is the pair index.

        Args:
            orbital: Orbital index
            n_orbitals: Total number of orbitals

        Returns:
            Pauli operator dictionary
        """
        if orbital not in self.orbital_to_pair:
            # Not part of a pair - fall back to JW
            return self._jordan_wigner_number(orbital, n_orbitals)

        pair_idx, mo_type = self.orbital_to_pair[orbital]

        # Qubit index for this orbital
        if mo_type == 'bonding':
            qubit_idx = 2 * pair_idx
        else:  # antibonding
            qubit_idx = 2 * pair_idx + 1

        # Build Pauli strings
        n_qubits = self.n_qubits(n_orbitals)
        pauli_I = 'I' * n_qubits
        pauli_Z = list('I' * n_qubits)
        pauli_Z[qubit_idx] = 'Z'
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
        Map excitation operator a†_i a_j.

        SPECIAL CASE: Bonding → Antibonding within same pair
            a†_antibonding a_bonding → (X_{2p} X_{2p+1} + Y_{2p} Y_{2p+1}) / 2

        This is the key covalent excitation!

        General case: Use modified Jordan-Wigner

        Args:
            orbital_from: Initial orbital
            orbital_to: Final orbital
            n_orbitals: Total orbitals

        Returns:
            Pauli operator dictionary
        """
        if orbital_from == orbital_to:
            return self.map_number_operator(orbital_from, n_orbitals)

        # Check if both orbitals are in the same pair
        if self._are_mo_pair(orbital_from, orbital_to):
            return self._map_mo_pair_excitation(orbital_from, orbital_to, n_orbitals)
        else:
            # Inter-pair excitation - use standard approach
            return self._map_inter_pair_excitation(orbital_from, orbital_to, n_orbitals)

    def _are_mo_pair(self, orb1: int, orb2: int) -> bool:
        """Check if two orbitals form a bonding/antibonding pair."""
        for bonding, antibonding in self.mo_pairs:
            if (orb1, orb2) == (bonding, antibonding) or (orb2, orb1) == (bonding, antibonding):
                return True
        return False

    def _map_mo_pair_excitation(
        self,
        orbital_from: int,
        orbital_to: int,
        n_orbitals: int
    ) -> Dict[str, complex]:
        """
        Map excitation within an MO pair.

        Bonding → Antibonding: Key covalent excitation
            (X_i X_j + Y_i Y_j) / 2

        where i, j are the qubit indices for the MO pair.
        """
        # Get pair index
        pair_idx_from = self.orbital_to_pair[orbital_from][0]
        pair_idx_to = self.orbital_to_pair[orbital_to][0]

        assert pair_idx_from == pair_idx_to, "Should be same pair"

        # Qubit indices
        qubit_i = 2 * pair_idx_from
        qubit_j = 2 * pair_idx_from + 1

        n_qubits = self.n_qubits(n_orbitals)

        # Build XX string
        xx_string = list('I' * n_qubits)
        xx_string[qubit_i] = 'X'
        xx_string[qubit_j] = 'X'
        xx_string = ''.join(xx_string)

        # Build YY string
        yy_string = list('I' * n_qubits)
        yy_string[qubit_i] = 'Y'
        yy_string[qubit_j] = 'Y'
        yy_string = ''.join(yy_string)

        return {
            xx_string: 0.25,
            yy_string: 0.25
        }

    def _map_inter_pair_excitation(
        self,
        orbital_from: int,
        orbital_to: int,
        n_orbitals: int
    ) -> Dict[str, complex]:
        """
        Map excitation between different MO pairs.

        Use Jordan-Wigner-like approach with Z string.
        """
        # For simplicity, use the same JW transformation
        # but on the qubit indices corresponding to the orbitals

        if orbital_from not in self.orbital_to_pair:
            from_qubit = orbital_from
        else:
            pair_idx, mo_type = self.orbital_to_pair[orbital_from]
            from_qubit = 2 * pair_idx + (0 if mo_type == 'bonding' else 1)

        if orbital_to not in self.orbital_to_pair:
            to_qubit = orbital_to
        else:
            pair_idx, mo_type = self.orbital_to_pair[orbital_to]
            to_qubit = 2 * pair_idx + (0 if mo_type == 'bonding' else 1)

        n_qubits = self.n_qubits(n_orbitals)

        # Build Z string
        if to_qubit > from_qubit:
            z_indices = list(range(from_qubit + 1, to_qubit))
        else:
            z_indices = list(range(to_qubit + 1, from_qubit))

        # Build Pauli strings
        xx_string = self._build_pauli_string('X', 'X', from_qubit, to_qubit, z_indices, n_qubits)
        yy_string = self._build_pauli_string('Y', 'Y', from_qubit, to_qubit, z_indices, n_qubits)
        xy_string = self._build_pauli_string('X', 'Y', from_qubit, to_qubit, z_indices, n_qubits)
        yx_string = self._build_pauli_string('Y', 'X', from_qubit, to_qubit, z_indices, n_qubits)

        return {
            xx_string: 0.25,
            yy_string: 0.25,
            xy_string: 0.25j,
            yx_string: -0.25j
        }

    def _build_pauli_string(
        self,
        pauli_j: str,
        pauli_i: str,
        j: int,
        i: int,
        z_indices: List[int],
        n_qubits: int
    ) -> str:
        """Build Pauli string with Z operators."""
        pauli = ['I'] * n_qubits
        pauli[j] = pauli_j
        pauli[i] = pauli_i

        for idx in z_indices:
            pauli[idx] = 'Z'

        return ''.join(pauli)

    def _jordan_wigner_number(self, orbital: int, n_orbitals: int) -> Dict[str, complex]:
        """Fallback to JW for orbitals not in pairs."""
        n_qubits = self.n_qubits(n_orbitals)
        pauli_I = 'I' * n_qubits
        pauli_Z = list('I' * n_qubits)
        pauli_Z[orbital] = 'Z'
        pauli_Z = ''.join(pauli_Z)

        return {
            pauli_I: 0.5,
            pauli_Z: -0.5
        }

    def create_bond_excitation_operator(self, pair_idx: int) -> Dict[str, complex]:
        """
        Create bonding → antibonding excitation for a specific MO pair.

        This is the primary excitation in covalent molecules.

        Args:
            pair_idx: Index of the MO pair

        Returns:
            Pauli operator dictionary
        """
        qubit_bonding = 2 * pair_idx
        qubit_antibonding = 2 * pair_idx + 1

        n_qubits = 2 * self.n_pairs

        # Build XX and YY strings
        xx_string = list('I' * n_qubits)
        xx_string[qubit_bonding] = 'X'
        xx_string[qubit_antibonding] = 'X'
        xx_string = ''.join(xx_string)

        yy_string = list('I' * n_qubits)
        yy_string[qubit_bonding] = 'Y'
        yy_string[qubit_antibonding] = 'Y'
        yy_string = ''.join(yy_string)

        return {
            xx_string: 0.5,
            yy_string: 0.5
        }

    def __repr__(self) -> str:
        """String representation."""
        return f"HybridOrbitalMapper(n_pairs={self.n_pairs}, locality='paired', best_for='covalent_bonding')"
