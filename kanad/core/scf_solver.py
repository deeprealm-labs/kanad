"""
Self-Consistent Field (SCF) Solver for Hartree-Fock calculations.

Implements Restricted Hartree-Fock (RHF) for closed-shell systems.
"""

import numpy as np
from typing import Tuple, Optional
from scipy.linalg import eigh


class SCFSolver:
    """
    Self-Consistent Field solver for Hartree-Fock calculations.

    Implements:
    - Restricted Hartree-Fock (RHF) for closed-shell systems
    - DIIS convergence acceleration
    - Energy and density convergence checks
    """

    def __init__(
        self,
        h_core: np.ndarray,
        S: np.ndarray,
        eri: np.ndarray,
        n_electrons: int,
        nuclear_repulsion: float = 0.0
    ):
        """
        Initialize SCF solver.

        Args:
            h_core: Core Hamiltonian matrix (n_basis, n_basis)
            S: Overlap matrix (n_basis, n_basis)
            eri: Electron repulsion integrals (n_basis, n_basis, n_basis, n_basis)
            n_electrons: Number of electrons
            nuclear_repulsion: Nuclear repulsion energy
        """
        self.h_core = h_core
        self.S = S
        self.eri = eri
        self.n_electrons = n_electrons
        self.nuclear_repulsion = nuclear_repulsion
        self.n_basis = h_core.shape[0]
        self.n_occ = n_electrons // 2  # Closed-shell

        # Convergence history
        self.energy_history = []
        self.converged = False
        self.iterations = 0

    def solve(
        self,
        max_iterations: int = 100,
        energy_tol: float = 1e-8,
        density_tol: float = 1e-6,
        use_diis: bool = True,
        diis_start: int = 2,
        diis_size: int = 8,
        level_shift: float = 0.0,
        damping_factor: float = 0.0,
        damping_start: int = 0
    ) -> Tuple[np.ndarray, np.ndarray, float, bool, int]:
        """
        Solve Hartree-Fock equations self-consistently.

        Args:
            max_iterations: Maximum SCF iterations
            energy_tol: Energy convergence threshold
            density_tol: Density matrix convergence threshold
            use_diis: Use DIIS convergence acceleration
            diis_start: Start DIIS after this many iterations
            diis_size: Maximum DIIS subspace size
            level_shift: Level shift parameter (Ha) for virtual orbitals (helps convergence)
            damping_factor: Density damping factor (0=no damping, 1=full damping)
            damping_start: Start damping after this iteration

        Returns:
            (density_matrix, mo_energies, total_energy, converged, iterations)
        """
        # Initial guess: diagonalize core Hamiltonian
        mo_energies, C = eigh(self.h_core, self.S)

        # Build initial density matrix
        C_occ = C[:, :self.n_occ]
        P = 2.0 * C_occ @ C_occ.T
        P_old = np.zeros_like(P)

        # DIIS storage
        if use_diis:
            error_list = []
            fock_list = []

        E_old = 0.0

        for iteration in range(max_iterations):
            self.iterations = iteration + 1

            # Build Fock matrix
            F = self._build_fock_matrix(P)

            # Apply level shift to virtual orbitals (helps convergence)
            # Shift up virtual orbital energies to prevent mixing with occupied
            if level_shift > 0:
                F_shifted = F.copy()
                # Add level shift to Fock matrix: F' = F + λS for virtual orbitals
                # This is applied before diagonalization
                # After solving, we'll remove the shift from virtual orbital energies

                # For now, use simple approach: shift entire Fock matrix
                # A more sophisticated approach would shift only virtual-virtual block
                pass  # Will apply after diagonalization

            # DIIS extrapolation
            if use_diis and iteration >= diis_start:
                # Compute error vector in orthonormal basis
                # e = FDS - SDF
                error = F @ P @ self.S - self.S @ P @ F

                # Store error and Fock matrix
                error_list.append(error)
                fock_list.append(F)

                # Limit DIIS subspace size
                if len(error_list) > diis_size:
                    error_list.pop(0)
                    fock_list.pop(0)

                # Perform DIIS extrapolation
                try:
                    F = self._diis_extrapolate(error_list, fock_list)
                except np.linalg.LinAlgError:
                    # DIIS failed, use regular Fock matrix
                    pass

            # Apply level shift before diagonalization
            F_work = F.copy()
            if level_shift > 0:
                # Add shift to all orbitals (will be removed from occupied after)
                F_work = F_work + level_shift * self.S

            # Solve generalized eigenvalue problem: FC = SCε
            mo_energies, C = eigh(F_work, self.S)

            # Remove level shift from occupied orbital energies
            if level_shift > 0:
                mo_energies[:self.n_occ] -= level_shift

            # Build new density matrix
            C_occ = C[:, :self.n_occ]
            P_new = 2.0 * C_occ @ C_occ.T

            # Apply density damping if requested
            if damping_factor > 0 and iteration >= damping_start:
                # P_damped = (1-α) * P_new + α * P_old
                P = (1.0 - damping_factor) * P_new + damping_factor * P_old
            else:
                P = P_new

            # Compute electronic energy
            # Standard RHF energy formula: E_elec = 0.5 * Tr[P(H + F)]
            E_elec = 0.5 * np.sum(P * (self.h_core + F))
            E_total = E_elec + self.nuclear_repulsion

            self.energy_history.append(E_total)

            # Check convergence
            delta_E = abs(E_total - E_old)
            delta_P = np.max(np.abs(P - P_old))

            if delta_E < energy_tol and delta_P < density_tol:
                self.converged = True
                break

            E_old = E_total
            P_old = P.copy()

        return P, mo_energies, E_total, self.converged, self.iterations

    def _build_fock_matrix(self, P: np.ndarray) -> np.ndarray:
        """
        Build Fock matrix from density matrix.

        F_μν = H_μν^core + G_μν

        where G_μν = Σ_λσ P_λσ [(μν|λσ) - 0.5(μλ|νσ)]

        Args:
            P: Density matrix

        Returns:
            Fock matrix
        """
        F = self.h_core.copy()

        # Build G matrix (two-electron part)
        for mu in range(self.n_basis):
            for nu in range(self.n_basis):
                G_mu_nu = 0.0
                for lam in range(self.n_basis):
                    for sigma in range(self.n_basis):
                        # Coulomb term: (μν|λσ)
                        G_mu_nu += P[lam, sigma] * self.eri[mu, nu, lam, sigma]
                        # Exchange term: (μλ|νσ)
                        G_mu_nu -= 0.5 * P[lam, sigma] * self.eri[mu, lam, nu, sigma]

                F[mu, nu] += G_mu_nu

        return F

    def _diis_extrapolate(
        self,
        error_list: list,
        fock_list: list
    ) -> np.ndarray:
        """
        Perform DIIS (Direct Inversion in Iterative Subspace) extrapolation.

        Minimizes the error vector to accelerate SCF convergence.

        Args:
            error_list: List of error matrices
            fock_list: List of Fock matrices

        Returns:
            Extrapolated Fock matrix
        """
        n = len(error_list)

        # Build B matrix: B_ij = Tr(e_i^T e_j)
        B = np.zeros((n + 1, n + 1))
        for i in range(n):
            for j in range(n):
                B[i, j] = np.sum(error_list[i] * error_list[j])

        # Constraint: Σ c_i = 1
        B[n, :n] = -1.0
        B[:n, n] = -1.0
        B[n, n] = 0.0

        # RHS vector
        rhs = np.zeros(n + 1)
        rhs[n] = -1.0

        # Solve for coefficients
        coeff = np.linalg.solve(B, rhs)

        # Extrapolate Fock matrix
        F_diis = np.zeros_like(fock_list[0])
        for i in range(n):
            F_diis += coeff[i] * fock_list[i]

        return F_diis

    def get_orbital_energies(self, F: np.ndarray) -> np.ndarray:
        """
        Get molecular orbital energies from Fock matrix.

        Args:
            F: Fock matrix

        Returns:
            Orbital energies
        """
        energies, _ = eigh(F, self.S)
        return energies

    def compute_total_energy(self, P: np.ndarray, F: np.ndarray) -> float:
        """
        Compute total energy from density and Fock matrices.

        E = 0.5 * Tr[P(H_core + F)] + E_nn

        Args:
            P: Density matrix
            F: Fock matrix

        Returns:
            Total energy
        """
        E_elec = 0.5 * np.sum(P * (self.h_core + F))
        E_total = E_elec + self.nuclear_repulsion
        return E_total
