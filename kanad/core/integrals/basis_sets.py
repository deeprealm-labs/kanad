"""
Gaussian basis sets for molecular orbital calculations.

Implements commonly used basis sets:
- STO-3G: Minimal basis (3 Gaussians per STO)
- 3-21G: Split valence
- 6-31G: Split valence with polarization options
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Dict
from kanad.core.constants.atomic_data import PeriodicTable


@dataclass
class GaussianPrimitive:
    """
    A single Gaussian primitive function.

    ψ(r) = N * x^l * y^m * z^n * exp(-α * r²)

    where N is the normalization constant.
    """

    exponent: float  # α (zeta)
    coefficient: float  # Contraction coefficient
    angular_momentum: Tuple[int, int, int]  # (l, m, n)
    center: np.ndarray  # Atomic position (x, y, z)

    def __post_init__(self):
        """Ensure center is numpy array."""
        if not isinstance(self.center, np.ndarray):
            self.center = np.ndarray(self.center)

    @property
    def l(self) -> int:
        """Total angular momentum quantum number."""
        return sum(self.angular_momentum)

    def evaluate(self, point: np.ndarray) -> float:
        """
        Evaluate Gaussian at a point in space.

        Args:
            point: 3D coordinate (x, y, z)

        Returns:
            Value of Gaussian primitive at point
        """
        r = point - self.center
        x, y, z = r
        lx, ly, lz = self.angular_momentum

        # Gaussian part
        r_squared = np.dot(r, r)
        gaussian = np.exp(-self.exponent * r_squared)

        # Angular part
        angular = (x ** lx) * (y ** ly) * (z ** lz)

        # Normalization
        norm = self._normalization_constant()

        return norm * self.coefficient * angular * gaussian

    def _normalization_constant(self) -> float:
        """
        Compute normalization constant for Gaussian primitive.

        For a 3D Cartesian Gaussian:
        φ(r) = N × x^lx × y^ly × z^lz × exp(-α|r-R|²)

        Normalization:
        N = (2α/π)^(3/4) × [(8α)^lx / (2lx-1)!!]^(1/2)
                         × [(8α)^ly / (2ly-1)!!]^(1/2)
                         × [(8α)^lz / (2lz-1)!!]^(1/2)

        Note: The (2α/π)^(3/4) factor appears ONCE for the 3D Gaussian,
        not once per dimension.
        """
        α = self.exponent
        lx, ly, lz = self.angular_momentum
        from scipy.special import factorial2

        # 3D Gaussian base normalization (appears once)
        base_norm = (2 * α / np.pi) ** 0.75

        # Angular momentum normalization for each component
        # For each direction: [(8α)^l / (2l-1)!!]^(1/2)
        norm_x = np.sqrt((8 * α) ** lx / (factorial2(2 * lx - 1, exact=True) if lx > 0 else 1))
        norm_y = np.sqrt((8 * α) ** ly / (factorial2(2 * ly - 1, exact=True) if ly > 0 else 1))
        norm_z = np.sqrt((8 * α) ** lz / (factorial2(2 * lz - 1, exact=True) if lz > 0 else 1))

        return base_norm * norm_x * norm_y * norm_z

    @staticmethod
    def _1d_normalization(α: float, l: int) -> float:
        """
        DEPRECATED: This method had incorrect normalization.
        Use _normalization_constant() instead.
        """
        # Kept for backwards compatibility but should not be used
        from scipy.special import factorial2
        numerator = (2 * α / np.pi) ** 0.75 * (4 * α) ** (l / 2)
        denominator = np.sqrt(factorial2(2 * l - 1, exact=True) if l > 0 else 1)
        return numerator / denominator


@dataclass
class ContractedGaussian:
    """
    A contracted Gaussian function (linear combination of primitives).

    ψ_contracted = Σᵢ cᵢ * ψᵢ(primitive)
    """

    primitives: List[GaussianPrimitive]
    shell_type: str  # 's', 'p', 'd', etc.

    def evaluate(self, point: np.ndarray) -> float:
        """Evaluate contracted Gaussian at point."""
        return sum(prim.evaluate(point) for prim in self.primitives)

    @property
    def center(self) -> np.ndarray:
        """Get center of contracted function."""
        return self.primitives[0].center if self.primitives else np.zeros(3)


class BasisSet:
    """
    Complete basis set for a molecule.

    Manages all basis functions for all atoms.
    """

    # STO-3G basis set data (exponents and coefficients)
    STO3G_DATA = {
        'H': {
            's': [
                (3.42525091, 0.15432897),
                (0.62391373, 0.53532814),
                (0.16885540, 0.44463454),
            ]
        },
        'C': {
            's': [
                (71.6168370, 0.15432897),
                (13.0450960, 0.53532814),
                (3.53051220, 0.44463454),
            ],
            's_valence': [
                (2.94124940, -0.09996723),
                (0.68348310, 0.39951283),
                (0.22228990, 0.70011547),
            ],
            'p': [
                (2.94124940, 0.15591627),
                (0.68348310, 0.60768372),
                (0.22228990, 0.39195739),
            ],
        },
        'N': {
            's': [
                (99.1061690, 0.15432897),
                (18.0523120, 0.53532814),
                (4.88566020, 0.44463454),
            ],
            's_valence': [
                (3.78045590, -0.09996723),
                (0.87849660, 0.39951283),
                (0.28571440, 0.70011547),
            ],
            'p': [
                (3.78045590, 0.15591627),
                (0.87849660, 0.60768372),
                (0.28571440, 0.39195739),
            ],
        },
        'O': {
            's': [
                (130.7093200, 0.15432897),
                (23.8088610, 0.53532814),
                (6.44360830, 0.44463454),
            ],
            's_valence': [
                (5.03315130, -0.09996723),
                (1.16959610, 0.39951283),
                (0.38038900, 0.70011547),
            ],
            'p': [
                (5.03315130, 0.15591627),
                (1.16959610, 0.60768372),
                (0.38038900, 0.39195739),
            ],
        },
        'F': {
            's': [
                (166.6791300, 0.15432897),
                (30.3608120, 0.53532814),
                (8.21682070, 0.44463454),
            ],
            's_valence': [
                (6.46480320, -0.09996723),
                (1.50228120, 0.39951283),
                (0.48858850, 0.70011547),
            ],
            'p': [
                (6.46480320, 0.15591627),
                (1.50228120, 0.60768372),
                (0.48858850, 0.39195739),
            ],
        },
    }

    def __init__(self, basis_name: str = 'sto-3g'):
        """
        Initialize basis set.

        Args:
            basis_name: Name of basis set ('sto-3g', '3-21g', '6-31g')
        """
        self.basis_name = basis_name.lower()
        self.basis_functions: List[ContractedGaussian] = []

    def build_basis(self, atoms: List['Atom']) -> None:
        """
        Build basis functions for a list of atoms.

        Args:
            atoms: List of Atom objects
        """
        self.basis_functions = []

        for atom in atoms:
            if self.basis_name == 'sto-3g':
                self._add_sto3g_functions(atom)
            else:
                raise NotImplementedError(f"Basis set '{self.basis_name}' not implemented yet")

    def _add_sto3g_functions(self, atom: 'Atom') -> None:
        """
        Add STO-3G basis functions for an atom.

        IMPORTANT: STO-3G exponents are in atomic units (bohr^-2).
        Atom positions are in Angstroms, so we convert to Bohr here.
        """
        from kanad.core.constants.conversion_factors import ConversionFactors

        symbol = atom.symbol
        # Convert position from Angstroms to Bohr (atomic units)
        position_angstrom = atom.position
        position = position_angstrom * ConversionFactors.ANGSTROM_TO_BOHR

        if symbol not in self.STO3G_DATA:
            raise ValueError(f"STO-3G basis not available for element '{symbol}'")

        basis_data = self.STO3G_DATA[symbol]

        # Add s orbital (core)
        if 's' in basis_data:
            s_primitives = []
            for exp, coeff in basis_data['s']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff,
                    angular_momentum=(0, 0, 0),  # s orbital
                    center=position
                )
                s_primitives.append(prim)

            self.basis_functions.append(
                ContractedGaussian(s_primitives, shell_type='s')
            )

        # Add valence s orbital
        if 's_valence' in basis_data:
            s_val_primitives = []
            for exp, coeff in basis_data['s_valence']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff,
                    angular_momentum=(0, 0, 0),
                    center=position
                )
                s_val_primitives.append(prim)

            self.basis_functions.append(
                ContractedGaussian(s_val_primitives, shell_type='s')
            )

        # Add p orbitals (px, py, pz)
        if 'p' in basis_data:
            # Normalization fix for Cartesian p orbitals
            # STO-3G coefficients need 1/√2 factor for proper normalization
            p_norm_factor = 1.0 / np.sqrt(2.0)

            # px orbital
            px_primitives = []
            for exp, coeff in basis_data['p']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff * p_norm_factor,
                    angular_momentum=(1, 0, 0),  # px
                    center=position
                )
                px_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(px_primitives, shell_type='px')
            )

            # py orbital
            py_primitives = []
            for exp, coeff in basis_data['p']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff * p_norm_factor,
                    angular_momentum=(0, 1, 0),  # py
                    center=position
                )
                py_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(py_primitives, shell_type='py')
            )

            # pz orbital
            pz_primitives = []
            for exp, coeff in basis_data['p']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff * p_norm_factor,
                    angular_momentum=(0, 0, 1),  # pz
                    center=position
                )
                pz_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(pz_primitives, shell_type='pz')
            )

    @property
    def n_basis_functions(self) -> int:
        """Get total number of basis functions."""
        return len(self.basis_functions)

    def get_function(self, idx: int) -> ContractedGaussian:
        """Get basis function by index."""
        return self.basis_functions[idx]
