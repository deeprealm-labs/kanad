"""
Overlap integral computation for Gaussian basis functions.

Overlap integrals: S_ij = ⟨φ_i|φ_j⟩

These are fundamental for normalizing molecular orbitals and are
used in the generalized eigenvalue problem HC = SCE.
"""

import numpy as np
from typing import Tuple
from scipy.special import factorial2, erf
from kanad.core.integrals.basis_sets import GaussianPrimitive, ContractedGaussian


class OverlapIntegrals:
    """
    Compute overlap integrals between Gaussian basis functions.

    Uses analytical formulas for Gaussian products.
    """

    @staticmethod
    def overlap_primitive(
        prim_a: GaussianPrimitive,
        prim_b: GaussianPrimitive
    ) -> float:
        """
        Compute overlap between two Gaussian primitives.

        S = ⟨φ_a|φ_b⟩ = ∫ φ_a(r) φ_b(r) dr

        Uses the Gaussian product theorem:
        exp(-α|r-A|²) exp(-β|r-B|²) = K exp(-γ|r-P|²)

        where:
            γ = α + β
            P = (αA + βB)/(α + β)
            K = exp(-αβ|A-B|²/γ)

        Args:
            prim_a: First Gaussian primitive
            prim_b: Second Gaussian primitive

        Returns:
            Overlap integral value
        """
        α = prim_a.exponent
        β = prim_b.exponent
        A = prim_a.center
        B = prim_b.center
        la, ma, na = prim_a.angular_momentum
        lb, mb, nb = prim_b.angular_momentum

        # Gaussian product parameters
        γ = α + β
        P = (α * A + β * B) / γ
        AB = A - B

        # Gaussian product prefactor
        K = np.exp(-α * β * np.dot(AB, AB) / γ)

        # Compute overlap for each Cartesian direction
        Sx = OverlapIntegrals._overlap_1d(la, lb, A[0] - P[0], B[0] - P[0], γ)
        Sy = OverlapIntegrals._overlap_1d(ma, mb, A[1] - P[1], B[1] - P[1], γ)
        Sz = OverlapIntegrals._overlap_1d(na, nb, A[2] - P[2], B[2] - P[2], γ)

        # Coefficients already include normalization in STO-3G
        return K * Sx * Sy * Sz

    @staticmethod
    def _overlap_1d(l1: int, l2: int, PA: float, PB: float, γ: float) -> float:
        """
        1D overlap integral for Cartesian Gaussians.

        S_l1,l2 = ∫ x^l1 exp(-α(x-A)²) x^l2 exp(-β(x-B)²) dx

        Simplified for s and p orbitals.

        Args:
            l1: Angular momentum quantum number of first Gaussian
            l2: Angular momentum quantum number of second Gaussian
            PA: P - A (center difference)
            PB: P - B (center difference)
            γ: Combined exponent α + β

        Returns:
            1D overlap integral
        """
        # For s-s overlap (l1=l2=0)
        if l1 == 0 and l2 == 0:
            return np.sqrt(np.pi / γ)

        # For s-p or p-s overlap (l1=0,l2=1 or l1=1,l2=0)
        elif l1 == 0 and l2 == 1:
            return PB * np.sqrt(np.pi / γ)
        elif l1 == 1 and l2 == 0:
            return PA * np.sqrt(np.pi / γ)

        # For p-p overlap (l1=l2=1)
        elif l1 == 1 and l2 == 1:
            return (PA * PB + 1/(2*γ)) * np.sqrt(np.pi / γ)

        # For higher angular momentum (d, f, etc.) - use approximation
        else:
            return np.sqrt(np.pi / γ) * (1.0 / (l1 + l2 + 1))

    @staticmethod
    def _binomial_prefactor(l1: int, l2: int, i: int, j: int, PA: float, PB: float) -> float:
        """
        Compute binomial prefactor for overlap integrals.

        Uses binomial expansion of (x-A)^l1 and (x-B)^l2.
        """
        from scipy.special import comb

        result = 0.0
        for k in range(max(0, i + j - l2), min(i, l1 - 2*i + 1) + 1):
            result += (comb(l1, k, exact=True) *
                      comb(l2, i + j - k, exact=True) *
                      PA**(l1 - 2*i - k) *
                      PB**(l2 - 2*j - i - j + k))

        return result / (4**i * 4**j) if i + j > 0 else 1.0

    @staticmethod
    def overlap_contracted(
        cgf_a: ContractedGaussian,
        cgf_b: ContractedGaussian
    ) -> float:
        """
        Compute overlap between two contracted Gaussian functions.

        S = Σᵢⱼ cᵢ cⱼ Nᵢ Nⱼ ⟨gᵢ|gⱼ⟩

        where:
        - cᵢ, cⱼ are contraction coefficients
        - Nᵢ, Nⱼ are normalization constants for primitives
        - ⟨gᵢ|gⱼ⟩ is the unnormalized Gaussian overlap integral

        Args:
            cgf_a: First contracted Gaussian
            cgf_b: Second contracted Gaussian

        Returns:
            Overlap integral value
        """
        overlap = 0.0

        for prim_a in cgf_a.primitives:
            for prim_b in cgf_b.primitives:
                # Overlap of normalized primitives:
                # ⟨φ_a|φ_b⟩ = N_a × N_b × ⟨g_a|g_b⟩
                overlap += (prim_a.coefficient *
                           prim_b.coefficient *
                           prim_a._normalization_constant() *
                           prim_b._normalization_constant() *
                           OverlapIntegrals.overlap_primitive(prim_a, prim_b))

        return overlap

    @staticmethod
    def build_overlap_matrix(basis_functions: list) -> np.ndarray:
        """
        Build the full overlap matrix S for a basis set.

        S_ij = ⟨φ_i|φ_j⟩

        Args:
            basis_functions: List of ContractedGaussian objects

        Returns:
            Overlap matrix (n_basis, n_basis)
        """
        n = len(basis_functions)
        S = np.zeros((n, n))

        for i in range(n):
            for j in range(i, n):  # Use symmetry
                S[i, j] = OverlapIntegrals.overlap_contracted(
                    basis_functions[i],
                    basis_functions[j]
                )
                S[j, i] = S[i, j]  # Symmetric

        return S
