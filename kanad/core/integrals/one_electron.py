"""
One-electron integrals for molecular Hamiltonians.

Includes:
- Kinetic energy: T_ij = ⟨φ_i|-½∇²|φ_j⟩
- Nuclear attraction: V_ij = ⟨φ_i|Σ_A -Z_A/r_iA|φ_j⟩
- Core Hamiltonian: H_core = T + V
"""

import numpy as np
from typing import List
from kanad.core.integrals.basis_sets import GaussianPrimitive, ContractedGaussian
from kanad.core.integrals.overlap import OverlapIntegrals
from kanad.core.atom import Atom
from scipy.special import erf


class OneElectronIntegrals:
    """
    Compute one-electron integrals for molecular systems.
    """

    def __init__(self, atoms: List[Atom], basis_functions: List[ContractedGaussian]):
        """
        Initialize one-electron integral calculator.

        Args:
            atoms: List of Atom objects
            basis_functions: List of basis functions
        """
        self.atoms = atoms
        self.basis_functions = basis_functions
        self.n_basis = len(basis_functions)

    def compute_kinetic(self) -> np.ndarray:
        """
        Compute kinetic energy matrix T_ij = ⟨φ_i|-½∇²|φ_j⟩

        For Gaussians, this can be computed analytically using:
        ∇² exp(-αr²) = (2α)² exp(-αr²) - 6α exp(-αr²)

        Returns:
            Kinetic energy matrix (n_basis, n_basis)
        """
        T = np.zeros((self.n_basis, self.n_basis))

        for i in range(self.n_basis):
            for j in range(i, self.n_basis):
                T[i, j] = self._kinetic_contracted(
                    self.basis_functions[i],
                    self.basis_functions[j]
                )
                T[j, i] = T[i, j]  # Symmetric

        return T

    def _kinetic_contracted(
        self,
        cgf_a: ContractedGaussian,
        cgf_b: ContractedGaussian
    ) -> float:
        """Compute kinetic integral between contracted Gaussians."""
        kinetic = 0.0

        for prim_a in cgf_a.primitives:
            for prim_b in cgf_b.primitives:
                kinetic += (prim_a.coefficient *
                           prim_b.coefficient *
                           self._kinetic_primitive(prim_a, prim_b))

        return kinetic

    def _kinetic_primitive(
        self,
        prim_a: GaussianPrimitive,
        prim_b: GaussianPrimitive
    ) -> float:
        """
        Compute kinetic integral between two primitives.

        Simplified formula for s-orbitals.
        """
        α = prim_a.exponent
        β = prim_b.exponent
        la, ma, na = prim_a.angular_momentum
        lb, mb, nb = prim_b.angular_momentum

        # Get basic overlap
        S_ab = OverlapIntegrals.overlap_primitive(prim_a, prim_b)

        # For s-s overlap (simplified kinetic formula)
        if la == ma == na == lb == mb == nb == 0:
            A = prim_a.center
            B = prim_b.center
            AB = np.linalg.norm(A - B)

            # Kinetic energy: T = (3αβ/(α+β) - 2αβ²R²/(α+β)) × S
            γ = α + β
            T = (3 * α * β / γ - 2 * α * β**2 * AB**2 / γ) * S_ab
            return T
        else:
            # Simplified for p-orbitals
            return β * 3 * S_ab

    def compute_nuclear_attraction(self) -> np.ndarray:
        """
        Compute nuclear-electron attraction matrix.

        V_ij = ⟨φ_i|Σ_A -Z_A/|r-R_A||φ_j⟩

        Returns:
            Nuclear attraction matrix (n_basis, n_basis)
        """
        V = np.zeros((self.n_basis, self.n_basis))

        for i in range(self.n_basis):
            for j in range(i, self.n_basis):
                V[i, j] = self._nuclear_contracted(
                    self.basis_functions[i],
                    self.basis_functions[j]
                )
                V[j, i] = V[i, j]  # Symmetric

        return V

    def _nuclear_contracted(
        self,
        cgf_a: ContractedGaussian,
        cgf_b: ContractedGaussian
    ) -> float:
        """Compute nuclear attraction between contracted Gaussians."""
        nuclear = 0.0

        for prim_a in cgf_a.primitives:
            for prim_b in cgf_b.primitives:
                for atom in self.atoms:
                    nuclear += (prim_a.coefficient *
                               prim_b.coefficient *
                               self._nuclear_primitive(prim_a, prim_b, atom))

        return nuclear

    def _nuclear_primitive(
        self,
        prim_a: GaussianPrimitive,
        prim_b: GaussianPrimitive,
        atom: Atom
    ) -> float:
        """
        Compute nuclear attraction integral for one nucleus.

        V_A = -Z_A ⟨φ_a|1/|r-R_A||φ_b⟩

        Uses Boys function for analytical evaluation.
        """
        α = prim_a.exponent
        β = prim_b.exponent
        A = prim_a.center
        B = prim_b.center
        C = atom.position
        Z = atom.atomic_number

        # Gaussian product
        γ = α + β
        P = (α * A + β * B) / γ
        AB = A - B
        PC = P - C

        # Prefactor
        K = np.exp(-α * β * np.dot(AB, AB) / γ)
        prefactor = -Z * 2 * np.pi / γ * K

        # Boys function evaluation
        T = γ * np.dot(PC, PC)
        F0 = self._boys_function(0, T)

        # For s-type Gaussians (simplified)
        la, ma, na = prim_a.angular_momentum
        lb, mb, nb = prim_b.angular_momentum

        # Simplified - assume coefficients include normalization
        if la == ma == na == lb == mb == nb == 0:
            # Both are s-orbitals
            return prefactor * F0
        else:
            # For p and higher orbitals, use approximate formula
            return prefactor * F0 * 0.5

    @staticmethod
    def _boys_function(n: int, T: float) -> float:
        """
        Boys function F_n(T) used in nuclear attraction integrals.

        F_n(T) = ∫₀¹ t^(2n) exp(-Tt²) dt

        For n=0: F_0(T) = √(π/4T) erf(√T) for T > 0
                         = 1 for T = 0
        """
        if T < 1e-10:
            return 1.0 / (2 * n + 1)

        if n == 0:
            return 0.5 * np.sqrt(np.pi / T) * erf(np.sqrt(T))

        # Recursion for higher n
        # F_n(T) = [(2n-1)F_{n-1}(T) - exp(-T)] / (2T)
        F_prev = OneElectronIntegrals._boys_function(0, T)
        for i in range(1, n + 1):
            F_curr = ((2 * i - 1) * F_prev - np.exp(-T)) / (2 * T)
            F_prev = F_curr

        return F_prev

    def compute_core_hamiltonian(self) -> np.ndarray:
        """
        Compute core Hamiltonian H_core = T + V_ne

        Returns:
            Core Hamiltonian matrix (n_basis, n_basis)
        """
        T = self.compute_kinetic()
        V = self.compute_nuclear_attraction()
        return T + V

    def compute_nuclear_repulsion(self) -> float:
        """
        Compute nuclear-nuclear repulsion energy.

        E_nn = Σ_{A<B} Z_A Z_B / |R_A - R_B|

        Returns:
            Nuclear repulsion energy in Hartree
        """
        E_nn = 0.0

        for i, atom_i in enumerate(self.atoms):
            for j, atom_j in enumerate(self.atoms):
                if i < j:
                    R_ij = atom_i.distance_to(atom_j)
                    if R_ij > 1e-10:  # Avoid division by zero
                        # Convert to Bohr for atomic units
                        from kanad.core.constants.conversion_factors import ConversionFactors
                        R_ij_bohr = ConversionFactors.length_to_bohr(R_ij, 'angstrom')
                        E_nn += atom_i.atomic_number * atom_j.atomic_number / R_ij_bohr

        return E_nn
