"""
Full Configuration Interaction (FCI) Solver.

Exact solution of electronic structure within a given basis set.
Builds and diagonalizes the CI Hamiltonian in the space of all possible
Slater determinants.
"""

from typing import Tuple, List, Dict, Any
import numpy as np
import logging
from itertools import combinations

logger = logging.getLogger(__name__)

from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian


class FCISolver:
    """
    Full Configuration Interaction solver.

    Exact diagonalization in the space of all Slater determinants
    for a given number of electrons and spatial orbitals.

    For spin-restricted case (singlets), uses spatial orbital formulation.
    """

    def __init__(
        self,
        hamiltonian: MolecularHamiltonian,
        spin_multiplicity: int = 1
    ):
        """
        Initialize FCI solver.

        Args:
            hamiltonian: Molecular Hamiltonian with h_core and ERI
            spin_multiplicity: 2S+1 (1=singlet, 2=doublet, 3=triplet, etc.)
        """
        self.hamiltonian = hamiltonian
        self.n_orbitals = hamiltonian.n_orbitals
        self.n_electrons = hamiltonian.n_electrons
        self.spin_multiplicity = spin_multiplicity

        # Nuclear repulsion
        self.E_nuc = hamiltonian.nuclear_repulsion

        # One- and two-electron integrals
        self.h_core = hamiltonian.h_core
        self.eri = hamiltonian.eri  # In physicist's notation: (ik|jl)

        logger.info(f"Initialized FCI: {self.n_electrons} electrons in {self.n_orbitals} orbitals")

    def solve(self, n_states: int = 1) -> Dict[str, Any]:
        """
        Solve FCI problem.

        Args:
            n_states: Number of states to compute (ground + excited)

        Returns:
            Dictionary with:
                - energies: List of energies (Hartree)
                - eigenvectors: CI coefficient matrices
                - ci_dimension: Size of CI space
        """
        logger.info("Building FCI Hamiltonian...")

        # Build CI space (all possible determinants)
        determinants = self._build_determinant_space()
        n_det = len(determinants)

        logger.info(f"CI space dimension: {n_det}")

        # Build Hamiltonian matrix
        H_CI = self._build_ci_hamiltonian(determinants)

        # Diagonalize
        logger.info("Diagonalizing CI Hamiltonian...")
        eigenvalues, eigenvectors = np.linalg.eigh(H_CI)

        # Add nuclear repulsion
        energies = eigenvalues[:n_states] + self.E_nuc

        logger.info(f"FCI ground state: {energies[0]:.6f} Ha")

        return {
            'energies': energies,
            'eigenvectors': eigenvectors[:, :n_states],
            'ci_dimension': n_det,
            'determinants': determinants,
            'hamiltonian_matrix': H_CI
        }

    def _build_determinant_space(self) -> List[Tuple[int, ...]]:
        """
        Build complete FCI space for singlets using occupation number representation.

        Each determinant is a tuple of length n_orbitals containing occupation numbers:
        - 0: orbital is empty
        - 1: orbital is singly occupied (only in singlet open-shell configs)
        - 2: orbital is doubly occupied (closed-shell)

        For H2 (2e, 2 orb):
        - (2, 0): Both electrons in orbital 0
        - (0, 2): Both electrons in orbital 1
        - (1, 1): One electron in each orbital (singlet coupled)

        Returns:
            List of determinants as tuples of occupation numbers
        """
        if self.spin_multiplicity != 1:
            raise NotImplementedError("Only singlet FCI implemented for now")

        if self.n_electrons % 2 != 0:
            raise ValueError("Odd number of electrons not supported in RHF-based FCI")

        determinants = []

        # Generate all ways to distribute n_electrons among n_orbitals
        # with allowed occupations: 0, 1 (open-shell), or 2 (closed-shell)
        # Subject to: sum of occupations = n_electrons

        # For 2 electrons specifically (general case is more complex):
        if self.n_electrons == 2:
            # Closed-shell configs: both electrons in same orbital
            for i in range(self.n_orbitals):
                occ = [0] * self.n_orbitals
                occ[i] = 2
                determinants.append(tuple(occ))

            # Open-shell singlets: one electron in each of two orbitals
            for i in range(self.n_orbitals):
                for j in range(i+1, self.n_orbitals):
                    occ = [0] * self.n_orbitals
                    occ[i] = 1
                    occ[j] = 1
                    determinants.append(tuple(occ))

        else:
            # For more electrons, use recursive generation
            # For now, raise error
            raise NotImplementedError(
                f"FCI for {self.n_electrons} electrons not yet implemented. "
                "Only 2-electron systems supported."
            )

        determinants.sort()
        return determinants

    def _build_ci_hamiltonian(self, determinants: List[Tuple[int, ...]]) -> np.ndarray:
        """
        Build full CI Hamiltonian matrix.

        H_IJ = ⟨Φ_I|H|Φ_J⟩

        where Φ_I are Slater determinants.

        Args:
            determinants: List of determinants

        Returns:
            CI Hamiltonian matrix (n_det × n_det)
        """
        n_det = len(determinants)
        H_CI = np.zeros((n_det, n_det), dtype=float)

        for I in range(n_det):
            for J in range(I, n_det):  # Symmetric matrix
                # Compute matrix element
                H_IJ = self._ci_matrix_element(determinants[I], determinants[J])
                H_CI[I, J] = H_IJ
                H_CI[J, I] = H_IJ  # Hermitian

        return H_CI

    def _ci_matrix_element(self, det_I: Tuple[int, ...], det_J: Tuple[int, ...]) -> float:
        """
        Compute Hamiltonian matrix element between two determinants.

        Uses Slater-Condon rules:
        - Identical determinants: ⟨Φ|H|Φ⟩ = Σ_i h_ii + Σ_{i<j} (2J_ij - K_ij)
        - Single excitation: ⟨Φ_i^a|H|Φ⟩ = h_ia + Σ_j (2⟨ij|aj⟩ - ⟨ij|ja⟩)
        - Double excitation: ⟨Φ_ij^ab|H|Φ⟩ = ⟨ij|ab⟩ - ⟨ij|ba⟩
        - Higher excitations: 0

        For RHF-based FCI (closed-shell singlets):
        Each spatial orbital is doubly occupied.

        Args:
            det_I, det_J: Determinants as tuples of occupied orbital indices

        Returns:
            Matrix element H_IJ (Hartree)
        """
        set_I = set(det_I)
        set_J = set(det_J)

        # Find excitation level
        diff_I = set_I - set_J  # Orbitals in I but not J
        diff_J = set_J - set_I  # Orbitals in J but not I

        n_excitation = len(diff_I)

        if n_excitation == 0:
            # Identical determinants
            return self._diagonal_element(det_I)
        elif n_excitation == 1:
            # Single excitation
            i = list(diff_I)[0]
            a = list(diff_J)[0]
            return self._single_excitation_element(det_I, det_J, i, a)
        elif n_excitation == 2:
            # Double excitation
            orbs_I = list(diff_I)
            orbs_J = list(diff_J)
            return self._double_excitation_element(orbs_I[0], orbs_I[1], orbs_J[0], orbs_J[1])
        else:
            # Higher excitations
            return 0.0

    def _diagonal_element(self, det: Tuple[int, ...]) -> float:
        """
        Diagonal element: ⟨Φ|H|Φ⟩

        Handles both:
        - Closed-shell: len(det) = n_electrons/2, each orbital doubly occupied
        - Open-shell: len(det) = n_electrons, each orbital singly occupied (singlet)

        For H2:
        - det = (0,): Both electrons in orbital 0 → E = 2*h[0,0] + eri[0,0,0,0]
        - det = (0,1): One in 0, one in 1 (singlet) → E = h[0,0] + h[1,1] + eri[0,1,0,1]
        """
        energy = 0.0
        n_orb_in_det = len(det)

        # Determine if closed-shell (doubly occupied) or open-shell (singly occupied)
        n_expected_pairs = self.n_electrons // 2
        is_closed_shell = (n_orb_in_det == n_expected_pairs)

        if is_closed_shell:
            # Each orbital in det is doubly occupied
            occ = list(det)

            # One-electron part (factor of 2 for both spins)
            for i in occ:
                energy += 2.0 * self.h_core[i, i]

            # Two-electron part
            for i in occ:
                for j in occ:
                    # Coulomb: both electrons in same spatial orbital
                    J_ij = self.eri[i, i, j, j]
                    energy += 2.0 * J_ij  # Factor of 2 for spin

                    # Exchange
                    K_ij = self.eri[i, j, i, j]
                    energy -= K_ij
        else:
            # Open-shell: each orbital singly occupied
            # For singlet: one α, one β (or symmetric combination)
            occ = list(det)

            # One-electron part (one electron per orbital)
            for i in occ:
                energy += self.h_core[i, i]

            # Two-electron part: interaction between singly-occupied orbitals
            for i in range(len(occ)):
                for j in range(i+1, len(occ)):
                    orb_i = occ[i]
                    orb_j = occ[j]
                    # Coulomb: electrons in different spatial orbitals
                    J = self.eri[orb_i, orb_i, orb_j, orb_j]
                    # Exchange: for singlet, only α-β interaction (no exchange)
                    # Actually for singlet coupled: E = J_ij (Coulomb only, no exchange)
                    energy += J

        return energy

    def _single_excitation_element(
        self,
        det_I: Tuple[int, ...],
        det_J: Tuple[int, ...],
        i: int,
        a: int
    ) -> float:
        """
        Single excitation in spatial orbital representation: i → a

        Handles different cases:
        1. Closed-shell → Closed-shell: Both electrons move from i to a
        2. Open-shell → Open-shell: One electron moves from i to a

        For H2 example:
        - (0,) → (0,1): excite one electron from 0 to 1 (0 becomes singly occupied)
        - (0,1) → (1,): excite from 0 to 1 (1 becomes doubly occupied)
        """
        n_expected_pairs = self.n_electrons // 2
        is_I_closed = (len(det_I) == n_expected_pairs)
        is_J_closed = (len(det_J) == n_expected_pairs)

        energy = 0.0

        if is_I_closed and is_J_closed:
            # Both closed-shell: both electrons move i→a
            energy = 2.0 * self.h_core[i, a]
            energy += self.eri[i, i, a, a]  # Coulomb
            energy -= self.eri[i, a, i, a]  # Exchange

        elif not is_I_closed and not is_J_closed:
            # Both open-shell: one electron moves i→a
            energy = self.h_core[i, a]

            # Interaction with other singly-occupied orbitals
            common_occ = set(det_I) & set(det_J)
            for j in common_occ:
                # Coulomb interaction with electron in orbital j
                energy += self.eri[i, a, j, j]
                # No exchange for singlet-coupled electrons in different orbitals

        else:
            # Mixed: closed→open or open→closed
            # This is more complex - for now use simplified formula
            energy = self.h_core[i, a]

        # Phase factor from reordering
        phase = self._excitation_phase(det_I, det_J)

        return phase * energy

    def _double_excitation_element(
        self,
        i: int,
        j: int,
        a: int,
        b: int
    ) -> float:
        """
        Double excitation: ⟨Φ_ij^ab|H|Φ⟩

        For doubly-occupied orbitals:
        H_ijab = 4 ⟨ij|ab⟩ - 2 ⟨ij|ba⟩
        """
        # Two-electron integral: (ij|ab) - (ij|ba)
        coulomb = self.eri[i, j, a, b]
        exchange = self.eri[i, j, b, a]

        energy = 4.0 * coulomb - 2.0 * exchange

        return energy

    def _excitation_phase(self, det_I: Tuple[int, ...], det_J: Tuple[int, ...]) -> float:
        """
        Compute phase factor from orbital reordering.

        For single excitation i→a, phase = (-1)^m where m is the number
        of orbitals between i and a.

        For now, return +1 (proper phase calculation requires more logic).
        """
        # Simplified: assuming phase is handled implicitly
        # Proper implementation would track orbital ordering
        return 1.0

    def get_fci_hamiltonian_matrix(self) -> np.ndarray:
        """
        Get the full FCI Hamiltonian matrix (without nuclear repulsion).

        This is what VQE should use for energy calculations.

        Returns:
            FCI Hamiltonian matrix (electronic part only)
        """
        determinants = self._build_determinant_space()
        H_CI = self._build_ci_hamiltonian(determinants)
        return H_CI
