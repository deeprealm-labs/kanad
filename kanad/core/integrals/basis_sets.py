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
        # For each direction: [(4α)^l / (2l-1)!!]^(1/2)
        norm_x = np.sqrt((4 * α) ** lx / (factorial2(2 * lx - 1, exact=True) if lx > 0 else 1))
        norm_y = np.sqrt((4 * α) ** ly / (factorial2(2 * ly - 1, exact=True) if ly > 0 else 1))
        norm_z = np.sqrt((4 * α) ** lz / (factorial2(2 * lz - 1, exact=True) if lz > 0 else 1))

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
        'He': {
            's': [
                (6.36242139, 0.15432897),
                (1.15892300, 0.53532814),
                (0.31364979, 0.44463454),
            ]
        },
        'Li': {
            's': [
                (16.1195750, 0.15432897),
                (2.93620070, 0.53532814),
                (0.79465050, 0.44463454),
            ],
            's_valence': [
                (0.63628970, -0.09996723),
                (0.14786010, 0.39951283),
                (0.04808870, 0.70011547),
            ],
            'p': [
                (0.63628970, 0.15591627),
                (0.14786010, 0.60768372),
                (0.04808870, 0.39195739),
            ],
        },
        'Be': {
            's': [
                (30.1678710, 0.15432897),
                (5.49515020, 0.53532814),
                (1.48731180, 0.44463454),
            ],
            's_valence': [
                (1.31405750, -0.09996723),
                (0.30538190, 0.39951283),
                (0.09937410, 0.70011547),
            ],
            'p': [
                (1.31405750, 0.15591627),
                (0.30538190, 0.60768372),
                (0.09937410, 0.39195739),
            ],
        },
        'B': {
            's': [
                (48.7914840, 0.15432897),
                (8.88676660, 0.53532814),
                (2.40525900, 0.44463454),
            ],
            's_valence': [
                (2.23627170, -0.09996723),
                (0.51982050, 0.39951283),
                (0.16906180, 0.70011547),
            ],
            'p': [
                (2.23627170, 0.15591627),
                (0.51982050, 0.60768372),
                (0.16906180, 0.39195739),
            ],
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
        'Ne': {
            's': [
                (207.0156540, 0.15432897),
                (37.7084590, 0.53532814),
                (10.2055870, 0.44463454),
            ],
            's_valence': [
                (8.02461490, -0.09996723),
                (1.86606560, 0.39951283),
                (0.60698070, 0.70011547),
            ],
            'p': [
                (8.02461490, 0.15591627),
                (1.86606560, 0.60768372),
                (0.60698070, 0.39195739),
            ],
        },
        'Na': {
            's': [
                (251.2357930, 0.15432897),
                (45.7956490, 0.53532814),
                (12.3906680, 0.44463454),
            ],
            's_valence': [
                (9.80765210, -0.09996723),
                (2.28124780, 0.39951283),
                (0.74239960, 0.70011547),
            ],
            's_valence2': [  # 3s valence
                (0.97616300, -0.22995122),
                (0.22684400, 0.28640156),
                (0.07388800, 1.16998209),
            ],
            'p': [
                (9.80765210, 0.15591627),
                (2.28124780, 0.60768372),
                (0.74239960, 0.39195739),
            ],
            'p_valence': [  # 3p
                (0.97616300, 0.07087427),
                (0.22684400, 0.33975129),
                (0.07388800, 0.72715858),
            ],
        },
        'Mg': {
            's': [
                (299.4392410, 0.15432897),
                (54.5252290, 0.53532814),
                (14.7566330, 0.44463454),
            ],
            's_valence': [
                (11.7864540, -0.09996723),
                (2.74349760, 0.39951283),
                (0.89225340, 0.70011547),
            ],
            's_valence2': [  # 3s
                (1.38970720, -0.22995122),
                (0.32321650, 0.28640156),
                (0.10520050, 1.16998209),
            ],
            'p': [
                (11.7864540, 0.15591627),
                (2.74349760, 0.60768372),
                (0.89225340, 0.39195739),
            ],
            'p_valence': [  # 3p
                (1.38970720, 0.07087427),
                (0.32321650, 0.33975129),
                (0.10520050, 0.72715858),
            ],
        },
        'Al': {
            's': [
                (352.7451490, 0.15432897),
                (64.2529460, 0.53532814),
                (17.3831440, 0.44463454),
            ],
            's_valence': [
                (13.9511730, -0.09996723),
                (3.24749210, 0.39951283),
                (1.05672200, 0.70011547),
            ],
            's_valence2': [  # 3s
                (1.87572640, -0.22995122),
                (0.43653090, 0.28640156),
                (0.14200540, 1.16998209),
            ],
            'p': [
                (13.9511730, 0.15591627),
                (3.24749210, 0.60768372),
                (1.05672200, 0.39195739),
            ],
            'p_valence': [  # 3p
                (1.87572640, 0.07087427),
                (0.43653090, 0.33975129),
                (0.14200540, 0.72715858),
            ],
        },
        'Si': {
            's': [
                (411.1966060, 0.15432897),
                (74.9337050, 0.53532814),
                (20.2726580, 0.44463454),
            ],
            's_valence': [
                (16.3092900, -0.09996723),
                (3.79787040, 0.39951283),
                (1.23567320, 0.70011547),
            ],
            's_valence2': [  # 3s
                (2.43226860, -0.22995122),
                (0.56585030, 0.28640156),
                (0.18417620, 1.16998209),
            ],
            'p': [
                (16.3092900, 0.15591627),
                (3.79787040, 0.60768372),
                (1.23567320, 0.39195739),
            ],
            'p_valence': [  # 3p
                (2.43226860, 0.07087427),
                (0.56585030, 0.33975129),
                (0.18417620, 0.72715858),
            ],
        },
        'P': {
            's': [
                (475.0453630, 0.15432897),
                (86.5890100, 0.53532814),
                (23.4244420, 0.44463454),
            ],
            's_valence': [
                (18.8632830, -0.09996723),
                (4.39205100, 0.39951283),
                (1.42869200, 0.70011547),
            ],
            's_valence2': [  # 3s
                (3.06697910, -0.22995122),
                (0.71391200, 0.28640156),
                (0.23236980, 1.16998209),
            ],
            'p': [
                (18.8632830, 0.15591627),
                (4.39205100, 0.60768372),
                (1.42869200, 0.39195739),
            ],
            'p_valence': [  # 3p
                (3.06697910, 0.07087427),
                (0.71391200, 0.33975129),
                (0.23236980, 0.72715858),
            ],
        },
        'S': {
            's': [
                (544.4067680, 0.15432897),
                (99.2198140, 0.53532814),
                (26.8355560, 0.44463454),
            ],
            's_valence': [
                (21.6144680, -0.09996723),
                (5.03456150, 0.39951283),
                (1.63831540, 0.70011547),
            ],
            's_valence2': [  # 3s
                (3.78654150, -0.22995122),
                (0.88204660, 0.28640156),
                (0.28727710, 1.16998209),
            ],
            'p': [
                (21.6144680, 0.15591627),
                (5.03456150, 0.60768372),
                (1.63831540, 0.39195739),
            ],
            'p_valence': [  # 3p
                (3.78654150, 0.07087427),
                (0.88204660, 0.33975129),
                (0.28727710, 0.72715858),
            ],
        },
        'Cl': {
            's': [
                (619.4118600, 0.15432897),
                (112.8203000, 0.53532814),
                (30.5231900, 0.44463454),
            ],
            's_valence': [
                (24.5549300, -0.09996723),
                (5.7172140, 0.39951283),
                (1.8604370, 0.70011547),
            ],
            's_valence2': [  # 3s
                (4.5846780, -0.22995122),
                (1.0681850, 0.28640156),
                (0.3478640, 1.16998209),
            ],
            'p': [
                (24.5549300, 0.15591627),
                (5.7172140, 0.60768372),
                (1.8604370, 0.39195739),
            ],
            'p_valence': [  # 3p
                (4.5846780, 0.07087427),
                (1.0681850, 0.33975129),
                (0.3478640, 0.72715858),
            ],
        },
        'Ar': {
            's': [
                (700.2896050, 0.15432897),
                (127.6109300, 0.53532814),
                (34.5287110, 0.44463454),
            ],
            's_valence': [
                (27.7034990, -0.09996723),
                (6.45044490, 0.39951283),
                (2.09824100, 0.70011547),
            ],
            's_valence2': [  # 3s
                (5.47067450, -0.22995122),
                (1.27438980, 0.28640156),
                (0.41487480, 1.16998209),
            ],
            'p': [
                (27.7034990, 0.15591627),
                (6.45044490, 0.60768372),
                (2.09824100, 0.39195739),
            ],
            'p_valence': [  # 3p
                (5.47067450, 0.07087427),
                (1.27438980, 0.33975129),
                (0.41487480, 0.72715858),
            ],
        },
    }

    # 6-31G basis set data (split-valence)
    # Core: 6 Gaussians contracted, Valence: 3 + 1 Gaussians
    G31G_DATA = {
        'H': {
            's': [  # Core (contracted from 3 primitives)
                (18.7311370, 0.03349460),
                (2.8253937, 0.23472695),
                (0.6401217, 0.81375733),
            ],
            's_valence': [  # Valence (single primitive)
                (0.1612778, 1.0),
            ],
        },
        'C': {
            's': [  # Core 1s (6 primitives)
                (3047.5249000, 0.00183470),
                (457.3695100, 0.01403730),
                (103.9486900, 0.06884260),
                (29.2101550, 0.23218440),
                (9.2866630, 0.46794130),
                (3.1639270, 0.36231200),
            ],
            's_valence': [  # Valence 2s (3 primitives)
                (7.8682724, -0.11933240),
                (1.8812885, -0.16085420),
                (0.5442493, 1.14345640),
            ],
            's_valence2': [  # Valence 2s' (single primitive)
                (0.1687144, 1.0),
            ],
            'p': [  # Valence 2p (3 primitives)
                (7.8682724, 0.06899910),
                (1.8812885, 0.31642400),
                (0.5442493, 0.74430830),
            ],
            'p_valence': [  # Valence 2p' (single primitive)
                (0.1687144, 1.0),
            ],
        },
        'N': {
            's': [  # Core 1s
                (4173.5110000, 0.00183480),
                (627.4579000, 0.01403730),
                (142.9021000, 0.06878660),
                (40.2343300, 0.23218440),
                (12.8202100, 0.46794130),
                (4.3906440, 0.36231200),
            ],
            's_valence': [  # Valence 2s
                (11.6263580, -0.11496220),
                (2.7162800, -0.16923060),
                (0.7722180, 1.14585610),
            ],
            's_valence2': [  # Valence 2s'
                (0.2120313, 1.0),
            ],
            'p': [  # Valence 2p
                (11.6263580, 0.06758920),
                (2.7162800, 0.32390260),
                (0.7722180, 0.74089990),
            ],
            'p_valence': [  # Valence 2p'
                (0.2120313, 1.0),
            ],
        },
        'O': {
            's': [  # Core 1s
                (5484.6717000, 0.00183110),
                (825.2349500, 0.01395010),
                (188.0469600, 0.06844510),
                (52.9645000, 0.23271430),
                (16.8975700, 0.47019300),
                (5.7996353, 0.35852090),
            ],
            's_valence': [  # Valence 2s
                (15.5396160, -0.11077750),
                (3.5999336, -0.14802630),
                (1.0137618, 1.13076700),
            ],
            's_valence2': [  # Valence 2s'
                (0.2700058, 1.0),
            ],
            'p': [  # Valence 2p
                (15.5396160, 0.07087427),
                (3.5999336, 0.33975129),
                (1.0137618, 0.72715858),
            ],
            'p_valence': [  # Valence 2p'
                (0.2700058, 1.0),
            ],
        },
        'F': {
            's': [  # Core 1s
                (7001.7130900, 0.00182938),
                (1051.3660900, 0.01395017),
                (239.2856900, 0.06844508),
                (67.3974453, 0.23271434),
                (21.5199573, 0.47019290),
                (7.40310130, 0.35852085),
            ],
            's_valence': [  # Valence 2s
                (20.8479528, -0.10852694),
                (4.80830834, -0.14606829),
                (1.34406986, 1.12887134),
            ],
            's_valence2': [  # Valence 2s'
                (0.35815139, 1.0),
            ],
            'p': [  # Valence 2p
                (20.8479528, 0.07160729),
                (4.80830834, 0.34591193),
                (1.34406986, 0.72216143),
            ],
            'p_valence': [  # Valence 2p'
                (0.35815139, 1.0),
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
            elif self.basis_name == '6-31g':
                self._add_6_31g_functions(atom)
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
            # px orbital
            px_primitives = []
            for exp, coeff in basis_data['p']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff,  # Use coefficient as-is (normalization handled by _normalization_constant)
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
                    coefficient=coeff,  # Use coefficient as-is
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
                    coefficient=coeff,  # Use coefficient as-is
                    angular_momentum=(0, 0, 1),  # pz
                    center=position
                )
                pz_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(pz_primitives, shell_type='pz')
            )

        # Add second valence s orbital (for period 3 elements: Na, Mg, Al, etc.)
        if 's_valence2' in basis_data:
            s_val2_primitives = []
            for exp, coeff in basis_data['s_valence2']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff,
                    angular_momentum=(0, 0, 0),
                    center=position
                )
                s_val2_primitives.append(prim)

            self.basis_functions.append(
                ContractedGaussian(s_val2_primitives, shell_type='s')
            )

        # Add valence p orbitals (for period 3 elements)
        if 'p_valence' in basis_data:
            # px valence
            px_val_primitives = []
            for exp, coeff in basis_data['p_valence']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff,  # Use coefficient as-is
                    angular_momentum=(1, 0, 0),
                    center=position
                )
                px_val_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(px_val_primitives, shell_type='px')
            )

            # py valence
            py_val_primitives = []
            for exp, coeff in basis_data['p_valence']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff,  # Use coefficient as-is
                    angular_momentum=(0, 1, 0),
                    center=position
                )
                py_val_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(py_val_primitives, shell_type='py')
            )

            # pz valence
            pz_val_primitives = []
            for exp, coeff in basis_data['p_valence']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff,  # Use coefficient as-is
                    angular_momentum=(0, 0, 1),
                    center=position
                )
                pz_val_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(pz_val_primitives, shell_type='pz')
            )

    def _add_6_31g_functions(self, atom: 'Atom') -> None:
        """
        Add 6-31G basis functions for an atom.

        6-31G is a split-valence basis:
        - Core: 6 Gaussians contracted to 1 function
        - Valence: 3 Gaussians + 1 Gaussian (split into 2 functions)
        """
        from kanad.core.constants.conversion_factors import ConversionFactors

        symbol = atom.symbol
        position_angstrom = atom.position
        position = position_angstrom * ConversionFactors.ANGSTROM_TO_BOHR

        if symbol not in self.G31G_DATA:
            raise ValueError(f"6-31G basis not available for element '{symbol}'")

        basis_data = self.G31G_DATA[symbol]

        # Add core s orbital (contracted)
        if 's' in basis_data:
            s_core_primitives = []
            for exp, coeff in basis_data['s']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff,
                    angular_momentum=(0, 0, 0),
                    center=position
                )
                s_core_primitives.append(prim)

            self.basis_functions.append(
                ContractedGaussian(s_core_primitives, shell_type='s')
            )

        # Add valence s orbital (contracted, 3 Gaussians)
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

        # Add valence s' orbital (single Gaussian)
        if 's_valence2' in basis_data:
            s_val2_primitives = []
            for exp, coeff in basis_data['s_valence2']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff,
                    angular_momentum=(0, 0, 0),
                    center=position
                )
                s_val2_primitives.append(prim)

            self.basis_functions.append(
                ContractedGaussian(s_val2_primitives, shell_type='s')
            )

        # Add p orbitals (contracted, 3 Gaussians)
        if 'p' in basis_data:
            p_norm_factor = 1.0 / np.sqrt(2.0)

            # px
            px_primitives = []
            for exp, coeff in basis_data['p']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff * p_norm_factor,
                    angular_momentum=(1, 0, 0),
                    center=position
                )
                px_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(px_primitives, shell_type='px')
            )

            # py
            py_primitives = []
            for exp, coeff in basis_data['p']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff * p_norm_factor,
                    angular_momentum=(0, 1, 0),
                    center=position
                )
                py_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(py_primitives, shell_type='py')
            )

            # pz
            pz_primitives = []
            for exp, coeff in basis_data['p']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff * p_norm_factor,
                    angular_momentum=(0, 0, 1),
                    center=position
                )
                pz_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(pz_primitives, shell_type='pz')
            )

        # Add p' orbitals (single Gaussian)
        if 'p_valence' in basis_data:
            p_norm_factor = 1.0 / np.sqrt(2.0)

            # px'
            px_val_primitives = []
            for exp, coeff in basis_data['p_valence']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff * p_norm_factor,
                    angular_momentum=(1, 0, 0),
                    center=position
                )
                px_val_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(px_val_primitives, shell_type='px')
            )

            # py'
            py_val_primitives = []
            for exp, coeff in basis_data['p_valence']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff * p_norm_factor,
                    angular_momentum=(0, 1, 0),
                    center=position
                )
                py_val_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(py_val_primitives, shell_type='py')
            )

            # pz'
            pz_val_primitives = []
            for exp, coeff in basis_data['p_valence']:
                prim = GaussianPrimitive(
                    exponent=exp,
                    coefficient=coeff * p_norm_factor,
                    angular_momentum=(0, 0, 1),
                    center=position
                )
                pz_val_primitives.append(prim)
            self.basis_functions.append(
                ContractedGaussian(pz_val_primitives, shell_type='pz')
            )

    @property
    def n_basis_functions(self) -> int:
        """Get total number of basis functions."""
        return len(self.basis_functions)

    def get_function(self, idx: int) -> ContractedGaussian:
        """Get basis function by index."""
        return self.basis_functions[idx]
