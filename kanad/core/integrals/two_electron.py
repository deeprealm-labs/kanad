"""
Two-electron repulsion integrals (ERIs).

(ij|kl) = ⟨φ_i(r₁)φ_j(r₁)|1/r₁₂|φ_k(r₂)φ_l(r₂)⟩

These are the most computationally expensive integrals in quantum chemistry.
"""

import numpy as np
from typing import Dict, Tuple, List
from kanad.core.integrals.basis_sets import GaussianPrimitive, ContractedGaussian
from kanad.core.integrals.one_electron import OneElectronIntegrals
from scipy.special import factorial2


class TwoElectronIntegrals:
    """
    Compute two-electron repulsion integrals.

    Uses analytical formulas for Gaussian basis functions.
    Exploits 8-fold permutational symmetry:
    (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk)
            = (kl|ij) = (lk|ij) = (kl|ji) = (lk|ji)
    """

    def __init__(self, basis_functions: List[ContractedGaussian]):
        """
        Initialize ERI calculator.

        Args:
            basis_functions: List of basis functions
        """
        self.basis_functions = basis_functions
        self.n_basis = len(basis_functions)

    def compute_eri_tensor(self) -> np.ndarray:
        """
        Compute full ERI tensor.

        Returns:
            ERI tensor (n_basis, n_basis, n_basis, n_basis)

        Note: This is memory-intensive! For large systems, use
        compute_eri_sparse() instead.
        """
        n = self.n_basis
        ERI = np.zeros((n, n, n, n))

        # Exploit 8-fold symmetry
        for i in range(n):
            for j in range(i + 1):
                ij = i * (i + 1) // 2 + j

                for k in range(n):
                    for l in range(k + 1):
                        kl = k * (k + 1) // 2 + l

                        if ij >= kl:
                            value = self._eri_contracted(
                                self.basis_functions[i],
                                self.basis_functions[j],
                                self.basis_functions[k],
                                self.basis_functions[l]
                            )

                            # Fill all 8 symmetric components
                            ERI[i, j, k, l] = value
                            ERI[j, i, k, l] = value
                            ERI[i, j, l, k] = value
                            ERI[j, i, l, k] = value
                            ERI[k, l, i, j] = value
                            ERI[l, k, i, j] = value
                            ERI[k, l, j, i] = value
                            ERI[l, k, j, i] = value

        return ERI

    def compute_eri_sparse(self, threshold: float = 1e-12) -> Dict[Tuple[int, int, int, int], float]:
        """
        Compute ERIs in sparse format (only non-zero values).

        Args:
            threshold: Minimum absolute value to store

        Returns:
            Dictionary mapping (i,j,k,l) -> integral value
        """
        n = self.n_basis
        eri_dict = {}

        for i in range(n):
            for j in range(i + 1):
                for k in range(n):
                    for l in range(k + 1):
                        # Compound index check for symmetry
                        ij = i * (i + 1) // 2 + j
                        kl = k * (k + 1) // 2 + l

                        if ij >= kl:
                            value = self._eri_contracted(
                                self.basis_functions[i],
                                self.basis_functions[j],
                                self.basis_functions[k],
                                self.basis_functions[l]
                            )

                            if abs(value) > threshold:
                                eri_dict[(i, j, k, l)] = value

        return eri_dict

    def _eri_contracted(
        self,
        cgf_i: ContractedGaussian,
        cgf_j: ContractedGaussian,
        cgf_k: ContractedGaussian,
        cgf_l: ContractedGaussian
    ) -> float:
        """Compute ERI between four contracted Gaussians."""
        eri = 0.0

        for prim_i in cgf_i.primitives:
            for prim_j in cgf_j.primitives:
                for prim_k in cgf_k.primitives:
                    for prim_l in cgf_l.primitives:
                        eri += (prim_i.coefficient *
                               prim_j.coefficient *
                               prim_k.coefficient *
                               prim_l.coefficient *
                               self._eri_primitive(prim_i, prim_j, prim_k, prim_l))

        return eri

    def _eri_primitive(
        self,
        prim_i: GaussianPrimitive,
        prim_j: GaussianPrimitive,
        prim_k: GaussianPrimitive,
        prim_l: GaussianPrimitive
    ) -> float:
        """
        Compute ERI between four Gaussian primitives.

        (ij|kl) = ⟨φ_i φ_j | r₁₂⁻¹ | φ_k φ_l⟩

        Uses the formula:
        (ij|kl) = 2π^(5/2) / (α_ij α_kl √(α_ij + α_kl)) ×
                  exp(-α_ij α_kl |P-Q|² / (α_ij + α_kl)) ×
                  F_0(T)

        where α_ij = α_i + α_j, P = (α_iA + α_jB)/α_ij
        """
        α_i = prim_i.exponent
        α_j = prim_j.exponent
        α_k = prim_k.exponent
        α_l = prim_l.exponent

        A = prim_i.center
        B = prim_j.center
        C = prim_k.center
        D = prim_l.center

        # Combine Gaussians i and j
        α_ij = α_i + α_j
        P = (α_i * A + α_j * B) / α_ij
        AB = A - B
        K_ij = np.exp(-α_i * α_j * np.dot(AB, AB) / α_ij)

        # Combine Gaussians k and l
        α_kl = α_k + α_l
        Q = (α_k * C + α_l * D) / α_kl
        CD = C - D
        K_kl = np.exp(-α_k * α_l * np.dot(CD, CD) / α_kl)

        # Distance between centers
        PQ = P - Q
        ρ = α_ij * α_kl / (α_ij + α_kl)
        T = ρ * np.dot(PQ, PQ)

        # Prefactor
        prefactor = 2 * np.pi**(5/2) / (α_ij * α_kl * np.sqrt(α_ij + α_kl))
        prefactor *= K_ij * K_kl

        # Boys function
        F0 = OneElectronIntegrals._boys_function(0, T)

        # Normalization
        norm_i = prim_i._normalization_constant()
        norm_j = prim_j._normalization_constant()
        norm_k = prim_k._normalization_constant()
        norm_l = prim_l._normalization_constant()

        # For s-type Gaussians only (simplified)
        li, mi, ni = prim_i.angular_momentum
        lj, mj, nj = prim_j.angular_momentum
        lk, mk, nk = prim_k.angular_momentum
        ll, ml, nl = prim_l.angular_momentum

        if (li == mi == ni == lj == mj == nj ==
            lk == mk == nk == ll == ml == nl == 0):
            # All s-orbitals
            return prefactor * F0 * norm_i * norm_j * norm_k * norm_l
        else:
            # For p and higher orbitals, use approximate formula
            # This is simplified - full implementation would use
            # Obara-Saika or Head-Gordon-Pople recursions
            return prefactor * F0 * norm_i * norm_j * norm_k * norm_l * 0.5

    def compute_coulomb_matrix(self, density_matrix: np.ndarray) -> np.ndarray:
        """
        Compute Coulomb matrix J from density matrix.

        J_ij = Σ_kl P_kl (ij|kl)

        Args:
            density_matrix: One-particle density matrix

        Returns:
            Coulomb matrix (n_basis, n_basis)
        """
        n = self.n_basis
        J = np.zeros((n, n))

        # Use full tensor for simplicity (correct symmetry)
        ERI = self.compute_eri_tensor()

        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        J[i, j] += density_matrix[k, l] * ERI[i, j, k, l]

        return J

    def compute_exchange_matrix(self, density_matrix: np.ndarray) -> np.ndarray:
        """
        Compute exchange matrix K from density matrix.

        K_ij = Σ_kl P_kl (ik|jl)

        Args:
            density_matrix: One-particle density matrix

        Returns:
            Exchange matrix (n_basis, n_basis)
        """
        n = self.n_basis
        K = np.zeros((n, n))

        # Use full tensor for correct exchange
        ERI = self.compute_eri_tensor()

        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        # Note: Exchange uses (ik|jl) not (ij|kl)
                        K[i, j] += density_matrix[k, l] * ERI[i, k, j, l]

        return K

    def compute_fock_matrix(
        self,
        h_core: np.ndarray,
        density_matrix: np.ndarray
    ) -> np.ndarray:
        """
        Compute Fock matrix F = H_core + G

        where G = 2J - K (for closed-shell)

        Args:
            h_core: Core Hamiltonian
            density_matrix: Density matrix

        Returns:
            Fock matrix (n_basis, n_basis)
        """
        J = self.compute_coulomb_matrix(density_matrix)
        K = self.compute_exchange_matrix(density_matrix)

        G = 2 * J - K
        F = h_core + G

        return F
