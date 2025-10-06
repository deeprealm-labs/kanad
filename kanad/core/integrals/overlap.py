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
        # Use XPA = P - A and XPB = P - B (standard convention from Szabo-Ostlund)
        Sx = OverlapIntegrals._overlap_1d(la, lb, P[0] - A[0], P[0] - B[0], γ)
        Sy = OverlapIntegrals._overlap_1d(ma, mb, P[1] - A[1], P[1] - B[1], γ)
        Sz = OverlapIntegrals._overlap_1d(na, nb, P[2] - A[2], P[2] - B[2], γ)

        # Coefficients already include normalization in STO-3G
        return K * Sx * Sy * Sz

    @staticmethod
    def _overlap_1d(l1: int, l2: int, XPA: float, XPB: float, γ: float) -> float:
        """
        1D overlap integral for Cartesian Gaussians using recursion relations.

        S_l1,l2 = ∫ x^l1 exp(-α(x-A)²) x^l2 exp(-β(x-B)²) dx

        Uses Obara-Saika recursion relations for arbitrary angular momentum.

        Args:
            l1: Angular momentum quantum number of first Gaussian
            l2: Angular momentum quantum number of second Gaussian
            XPA: P - A (where P is the Gaussian product center)
            XPB: P - B
            γ: Combined exponent α + β

        Returns:
            1D overlap integral
        """
        # Base case: S(0,0)
        if l1 == 0 and l2 == 0:
            return np.sqrt(np.pi / γ)

        # Use recursion to build up from s-s
        # Recursion relation (increasing l2):
        # S(l1, l2+1) = XPB·S(l1,l2) + l1/(2γ)·S(l1-1,l2) + l2/(2γ)·S(l1,l2-1)

        # For efficiency, handle common cases explicitly
        if l1 == 0 and l2 == 1:
            return XPB * np.sqrt(np.pi / γ)
        elif l1 == 1 and l2 == 0:
            return XPA * np.sqrt(np.pi / γ)
        elif l1 == 1 and l2 == 1:
            return (XPA * XPB + 1/(2*γ)) * np.sqrt(np.pi / γ)

        # For l2 = 2, 3, etc., use recursion
        # First build S(l1, 0) through S(l1, l2) using recursion

        # Build table S[i,j] for i=0..l1, j=0..l2
        S = {}

        # Base: S(0,0)
        S[(0,0)] = np.sqrt(np.pi / γ)

        # Build column 0: S(i, 0) for i = 1..l1
        for i in range(1, l1 + 1):
            # S(i,0) = XPA·S(i-1,0) + (i-1)/(2γ)·S(i-2,0)
            S[(i,0)] = XPA * S[(i-1,0)]
            if i >= 2:
                S[(i,0)] += (i-1)/(2*γ) * S[(i-2,0)]

        # Build row 0: S(0, j) for j = 1..l2
        for j in range(1, l2 + 1):
            # S(0,j) = XPB·S(0,j-1) + (j-1)/(2γ)·S(0,j-2)
            S[(0,j)] = XPB * S[(0,j-1)]
            if j >= 2:
                S[(0,j)] += (j-1)/(2*γ) * S[(0,j-2)]

        # Build remaining entries: S(i,j) for i=1..l1, j=1..l2
        for i in range(1, l1 + 1):
            for j in range(1, l2 + 1):
                # S(i,j) = XPB·S(i,j-1) + i/(2γ)·S(i-1,j-1) + (j-1)/(2γ)·S(i,j-2)
                S[(i,j)] = XPB * S[(i,j-1)]
                S[(i,j)] += i/(2*γ) * S[(i-1,j-1)]
                if j >= 2:
                    S[(i,j)] += (j-1)/(2*γ) * S[(i,j-2)]

        return S[(l1, l2)]

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
