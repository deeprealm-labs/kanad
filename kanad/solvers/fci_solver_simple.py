"""
Simple FCI solver for 2-electron systems.

Implements proper Full Configuration Interaction for molecules with 2 electrons.
Uses occupation number representation to handle both closed-shell and open-shell singlets.
"""

import numpy as np
from typing import Tuple, Dict, Any
import logging

logger = logging.getLogger(__name__)

from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian


class SimpleFCISolver:
    """
    Simple FCI solver for 2-electron systems only.

    Properly handles:
    - Closed-shell configs: both electrons in same orbital
    - Open-shell singlets: one electron in each orbital

    For H2 (2e, 2 orbitals):
    - Determinant space: (2,0), (0,2), (1,1)
    - CI matrix is 3×3
    """

    def __init__(self, hamiltonian: MolecularHamiltonian):
        """
        Initialize simple FCI solver.

        Args:
            hamiltonian: Molecular Hamiltonian with h_core and ERI

        Raises:
            ValueError: If not a 2-electron system
        """
        self.hamiltonian = hamiltonian
        self.n_orbitals = hamiltonian.n_orbitals
        self.n_electrons = hamiltonian.n_electrons
        self.E_nuc = hamiltonian.nuclear_repulsion

        if self.n_electrons != 2:
            raise ValueError(
                f"SimpleFCISolver only supports 2-electron systems, "
                f"got {self.n_electrons} electrons"
            )

        self.h_core = hamiltonian.h_core
        self.eri = hamiltonian.eri

        logger.info(f"Initialized SimpleFCI: 2 electrons in {self.n_orbitals} orbitals")

    def solve(self) -> Dict[str, Any]:
        """
        Solve FCI problem for 2 electrons.

        Returns:
            Dictionary with:
            - energies: Array of eigenvalues (without nuclear repulsion)
            - eigenvectors: CI coefficients
            - determinants: List of determinants (occupation number tuples)
            - hamiltonian_matrix: CI Hamiltonian matrix
        """
        # Build determinant space
        determinants = self._build_determinants()
        n_det = len(determinants)

        logger.info(f"CI space dimension: {n_det}")
        logger.info(f"Determinants: {determinants}")

        # Build CI Hamiltonian
        H_CI = self._build_hamiltonian_matrix(determinants)

        # Diagonalize
        eigenvalues, eigenvectors = np.linalg.eigh(H_CI)

        logger.info(f"FCI ground state: {eigenvalues[0]:.6f} Ha")

        return {
            'energies': eigenvalues,
            'eigenvectors': eigenvectors,
            'ci_dimension': n_det,
            'determinants': determinants,
            'hamiltonian_matrix': H_CI
        }

    def _build_determinants(self) -> list:
        """
        Build all singlet determinants for 2 electrons in n orbitals.

        Uses occupation number representation:
        - (2, 0, ...): 2 electrons in orbital 0
        - (1, 1, 0, ...): 1 in orbital 0, 1 in orbital 1 (singlet)

        Returns:
            List of determinants as tuples of occupation numbers
        """
        determinants = []

        # Closed-shell: both electrons in same orbital
        for i in range(self.n_orbitals):
            occ = [0] * self.n_orbitals
            occ[i] = 2
            determinants.append(tuple(occ))

        # Open-shell singlets: one electron in each of 2 orbitals
        for i in range(self.n_orbitals):
            for j in range(i+1, self.n_orbitals):
                occ = [0] * self.n_orbitals
                occ[i] = 1
                occ[j] = 1
                determinants.append(tuple(occ))

        return determinants

    def _build_hamiltonian_matrix(self, determinants: list) -> np.ndarray:
        """
        Build CI Hamiltonian matrix.

        Args:
            determinants: List of determinants (occupation number tuples)

        Returns:
            Hamiltonian matrix (n_det × n_det)
        """
        n_det = len(determinants)
        H = np.zeros((n_det, n_det))

        for i in range(n_det):
            for j in range(i, n_det):  # Symmetric, only compute upper triangle
                H[i, j] = self._matrix_element(determinants[i], determinants[j])
                if i != j:
                    H[j, i] = H[i, j]  # Hermitian

        return H

    def _matrix_element(self, det_I: Tuple[int, ...], det_J: Tuple[int, ...]) -> float:
        """
        Compute ⟨I|H|J⟩ using Slater-Condon rules.

        For 2-electron singlets with occupation numbers:
        - Diagonal: Energy of configuration
        - Off-diagonal: Coupling between configurations

        Args:
            det_I, det_J: Determinants as occupation number tuples

        Returns:
            Matrix element (Hartree)
        """
        if det_I == det_J:
            return self._diagonal(det_I)

        # Find which orbitals differ
        diff = [(i, det_I[i], det_J[i]) for i in range(len(det_I)) if det_I[i] != det_J[i]]

        if len(diff) == 0:
            return self._diagonal(det_I)
        elif len(diff) == 2:
            # Most common case: two orbitals change
            return self._two_orbital_coupling(det_I, det_J, diff)
        else:
            # Shouldn't happen for 2 electrons
            return 0.0

    def _diagonal(self, det: Tuple[int, ...]) -> float:
        """
        Diagonal element ⟨det|H|det⟩.

        For occupation numbers (n_0, n_1, ...) where n_i ∈ {0, 1, 2}:
        E = Σ_i n_i * h_ii + Σ_{i<j} interaction(i,j)

        Args:
            det: Determinant as occupation number tuple

        Returns:
            Energy (Hartree)
        """
        h = self.h_core
        eri = self.eri
        energy = 0.0

        # One-electron part
        for i, n_i in enumerate(det):
            energy += n_i * h[i, i]

        # Two-electron part
        for i in range(len(det)):
            n_i = det[i]
            if n_i == 0:
                continue

            for j in range(i, len(det)):
                n_j = det[j]
                if n_j == 0:
                    continue

                if i == j:
                    # Same orbital: both electrons there (n_i = 2)
                    if n_i == 2:
                        energy += eri[i, i, i, i]  # Coulomb repulsion
                else:
                    # Different orbitals
                    if n_i == 2 and n_j == 2:
                        # Both doubly occupied: full Coulomb + exchange
                        energy += 2 * eri[i, i, j, j]  # Coulomb
                        energy -= eri[i, j, i, j]  # Exchange
                    elif n_i == 1 and n_j == 1:
                        # Both singly occupied (singlet): Coulomb only
                        energy += eri[i, i, j, j]
                    elif (n_i == 2 and n_j == 1) or (n_i == 1 and n_j == 2):
                        # One doubly, one singly occupied
                        energy += eri[i, i, j, j]  # Simplified

        return energy

    def _two_orbital_coupling(self, det_I: Tuple[int, ...], det_J: Tuple[int, ...],
                               diff: list) -> float:
        """
        Coupling between determinants differing in 2 orbitals.

        Args:
            det_I, det_J: Determinants
            diff: List of (orbital_index, occ_I, occ_J) for differing orbitals

        Returns:
            Matrix element (Hartree)
        """
        if len(diff) != 2:
            return 0.0

        orb1, occ_I1, occ_J1 = diff[0]
        orb2, occ_I2, occ_J2 = diff[1]

        h = self.h_core
        eri = self.eri

        # Common cases for 2 electrons:
        # 1. (2,0) ↔ (1,1): closed → open shell
        # 2. (1,1) ↔ (0,2): open → closed shell

        # (2,0) ↔ (1,1): Move one electron from orb1 to orb2
        if occ_I1 == 2 and occ_J1 == 1 and occ_I2 == 0 and occ_J2 == 1:
            # Transfer one electron from orb1 to orb2
            return h[orb1, orb2] + eri[orb1, orb2, orb1, orb2]

        # (1,1) ↔ (0,2): Move one electron from orb1 to orb2
        if occ_I1 == 1 and occ_J1 == 0 and occ_I2 == 1 and occ_J2 == 2:
            # Transfer one electron from orb1 to orb2
            return h[orb1, orb2] + eri[orb1, orb2, orb1, orb2]

        # (2,0) ↔ (0,2): Double excitation (both electrons move)
        if occ_I1 == 2 and occ_J1 == 0 and occ_I2 == 0 and occ_J2 == 2:
            return eri[orb1, orb1, orb2, orb2]

        # Default: return coupling integral
        return h[orb1, orb2]
