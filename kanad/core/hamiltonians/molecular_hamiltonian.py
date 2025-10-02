"""
Base class for molecular Hamiltonians.

Provides common functionality for all bonding-type-specific Hamiltonians.
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Tuple, Optional
import numpy as np


class MolecularHamiltonian(ABC):
    """
    Base class for molecular Hamiltonians.

    General molecular electronic Hamiltonian:
        H = H_core + H_ee + H_nn

    where:
        H_core = T + V_ne (kinetic + nuclear-electron attraction)
        H_ee = electron-electron repulsion
        H_nn = nuclear-nuclear repulsion (constant energy shift)

    Each bonding type (ionic, covalent, metallic) emphasizes different terms.
    """

    def __init__(
        self,
        n_orbitals: int,
        n_electrons: int,
        nuclear_repulsion: float = 0.0
    ):
        """
        Initialize molecular Hamiltonian.

        Args:
            n_orbitals: Number of spatial orbitals
            n_electrons: Number of electrons
            nuclear_repulsion: Nuclear repulsion energy (constant)
        """
        self.n_orbitals = n_orbitals
        self.n_electrons = n_electrons
        self.nuclear_repulsion = nuclear_repulsion

        # Integral matrices (filled by subclasses)
        self.h_core: Optional[np.ndarray] = None  # (n_orbitals, n_orbitals)
        self.eri: Optional[np.ndarray] = None     # (n_orbitals, n_orbitals, n_orbitals, n_orbitals)

    @abstractmethod
    def to_matrix(self) -> np.ndarray:
        """
        Convert Hamiltonian to matrix representation.

        Returns:
            Hamiltonian matrix in the chosen basis
        """
        pass

    @abstractmethod
    def compute_energy(self, density_matrix: np.ndarray) -> float:
        """
        Compute total electronic energy from density matrix.

        E = Tr[P * H_core] + ½ Tr[P * G] + E_nn

        where G = 2J - K (Coulomb - Exchange)

        Args:
            density_matrix: One-particle density matrix

        Returns:
            Total electronic energy
        """
        pass

    def get_hf_energy(self, max_iter: int = 50, tol: float = 1e-6) -> Tuple[float, np.ndarray, np.ndarray]:
        """
        Compute Hartree-Fock energy via SCF procedure.

        Args:
            max_iter: Maximum SCF iterations
            tol: Energy convergence tolerance

        Returns:
            (energy, density_matrix, mo_coefficients)
        """
        if self.h_core is None or self.eri is None:
            raise ValueError("Hamiltonian not built - h_core or ERI is None")

        from kanad.core.integrals.two_electron import TwoElectronIntegrals
        from kanad.core.integrals.basis_sets import ContractedGaussian

        # Initial guess: core Hamiltonian eigenvalues
        eigenvalues, eigenvectors = np.linalg.eigh(self.h_core)

        # Build initial density matrix (fill lowest n_electrons/2 orbitals)
        n_occ = self.n_electrons // 2
        C_occ = eigenvectors[:, :n_occ]
        P = 2 * C_occ @ C_occ.T  # Factor of 2 for closed-shell

        energy_old = 0.0

        for iteration in range(max_iter):
            # Build Fock matrix
            # Create dummy basis functions for ERI calc
            basis_funcs = [ContractedGaussian([], np.zeros(3)) for _ in range(self.n_orbitals)]
            eri_calc = TwoElectronIntegrals(basis_funcs)
            eri_calc.basis_functions = basis_funcs
            eri_calc.n_basis = self.n_orbitals

            # Use pre-computed ERI
            J = np.zeros_like(self.h_core)
            K = np.zeros_like(self.h_core)

            for i in range(self.n_orbitals):
                for j in range(self.n_orbitals):
                    for k in range(self.n_orbitals):
                        for l in range(self.n_orbitals):
                            J[i, j] += P[k, l] * self.eri[i, j, k, l]
                            K[i, j] += P[k, l] * self.eri[i, k, j, l]

            G = 2 * J - K
            F = self.h_core + G

            # Solve Fock eigenvalue problem
            eigenvalues, eigenvectors = np.linalg.eigh(F)

            # Update density
            C_occ = eigenvectors[:, :n_occ]
            P = 2 * C_occ @ C_occ.T

            # Compute energy
            energy = self.compute_energy(P)

            # Check convergence
            if abs(energy - energy_old) < tol:
                return energy, P, eigenvectors

            energy_old = energy

        raise RuntimeError(f"SCF did not converge in {max_iter} iterations")

    def get_one_body_tensor(self) -> np.ndarray:
        """
        Get one-body (core) Hamiltonian tensor.

        Returns:
            h_core matrix
        """
        if self.h_core is None:
            raise ValueError("Core Hamiltonian not initialized")
        return self.h_core.copy()

    def get_two_body_tensor(self) -> np.ndarray:
        """
        Get two-body (electron repulsion) tensor.

        Returns:
            ERI tensor (n, n, n, n)
        """
        if self.eri is None:
            raise ValueError("ERI tensor not initialized")
        return self.eri.copy()

    def get_nuclear_repulsion(self) -> float:
        """Get nuclear repulsion energy."""
        return self.nuclear_repulsion

    def __repr__(self) -> str:
        """String representation."""
        return (f"{self.__class__.__name__}("
                f"n_orbitals={self.n_orbitals}, "
                f"n_electrons={self.n_electrons}, "
                f"E_nn={self.nuclear_repulsion:.6f})")
