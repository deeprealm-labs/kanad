"""
Atomic element properties and periodic table data.

Data sources:
- Covalent radii: Cordero et al. (2008)
- Electronegativity: Pauling scale
- Ionization energies: NIST database
"""

from dataclasses import dataclass
from typing import Dict, Optional


@dataclass
class AtomicProperties:
    """Properties of an atomic element."""

    symbol: str
    atomic_number: int
    atomic_mass: float  # amu (atomic mass units)
    covalent_radius: float  # Angstrom
    van_der_waals_radius: float  # Angstrom
    electronegativity: float  # Pauling scale
    valence_electrons: int
    ionization_energy: float  # eV (first ionization energy)
    electron_affinity: float  # eV
    oxidation_states: tuple  # Common oxidation states
    is_metal: bool
    group: int  # Periodic table group
    period: int  # Periodic table period


class PeriodicTable:
    """
    Complete periodic table with bonding-relevant properties.

    Provides atomic data for quantum chemistry calculations with focus
    on bonding analysis.
    """

    ELEMENTS: Dict[str, AtomicProperties] = {
        # Period 1
        'H': AtomicProperties('H', 1, 1.008, 0.31, 1.20, 2.20, 1, 13.598, 0.754, (-1, 1), False, 1, 1),
        'He': AtomicProperties('He', 2, 4.003, 0.28, 1.40, 0.00, 2, 24.587, 0.00, (0,), False, 18, 1),

        # Period 2
        'Li': AtomicProperties('Li', 3, 6.941, 1.28, 1.82, 0.98, 1, 5.392, 0.618, (1,), True, 1, 2),
        'Be': AtomicProperties('Be', 4, 9.012, 0.96, 1.53, 1.57, 2, 9.323, 0.00, (2,), True, 2, 2),
        'B': AtomicProperties('B', 5, 10.81, 0.84, 1.92, 2.04, 3, 8.298, 0.277, (3,), False, 13, 2),
        'C': AtomicProperties('C', 6, 12.011, 0.76, 1.70, 2.55, 4, 11.260, 1.263, (-4, 2, 4), False, 14, 2),
        'N': AtomicProperties('N', 7, 14.007, 0.71, 1.55, 3.04, 5, 14.534, 0.07, (-3, 3, 5), False, 15, 2),
        'O': AtomicProperties('O', 8, 15.999, 0.66, 1.52, 3.44, 6, 13.618, 1.461, (-2,), False, 16, 2),
        'F': AtomicProperties('F', 9, 18.998, 0.57, 1.47, 3.98, 7, 17.423, 3.401, (-1,), False, 17, 2),
        'Ne': AtomicProperties('Ne', 10, 20.180, 0.58, 1.54, 0.00, 8, 21.565, 0.00, (0,), False, 18, 2),

        # Period 3
        'Na': AtomicProperties('Na', 11, 22.990, 1.66, 2.27, 0.93, 1, 5.139, 0.548, (1,), True, 1, 3),
        'Mg': AtomicProperties('Mg', 12, 24.305, 1.41, 1.73, 1.31, 2, 7.646, 0.00, (2,), True, 2, 3),
        'Al': AtomicProperties('Al', 13, 26.982, 1.21, 1.84, 1.61, 3, 5.986, 0.433, (3,), True, 13, 3),
        'Si': AtomicProperties('Si', 14, 28.085, 1.11, 2.10, 1.90, 4, 8.152, 1.385, (-4, 4), False, 14, 3),
        'P': AtomicProperties('P', 15, 30.974, 1.07, 1.80, 2.19, 5, 10.487, 0.747, (-3, 3, 5), False, 15, 3),
        'S': AtomicProperties('S', 16, 32.06, 1.05, 1.80, 2.58, 6, 10.360, 2.077, (-2, 4, 6), False, 16, 3),
        'Cl': AtomicProperties('Cl', 17, 35.45, 1.02, 1.75, 3.16, 7, 12.968, 3.617, (-1, 1, 3, 5, 7), False, 17, 3),
        'Ar': AtomicProperties('Ar', 18, 39.948, 1.06, 1.88, 0.00, 8, 15.760, 0.00, (0,), False, 18, 3),

        # Period 4 - Transition metals and selected elements
        'K': AtomicProperties('K', 19, 39.098, 2.03, 2.75, 0.82, 1, 4.341, 0.501, (1,), True, 1, 4),
        'Ca': AtomicProperties('Ca', 20, 40.078, 1.76, 2.31, 1.00, 2, 6.113, 0.025, (2,), True, 2, 4),
        'Sc': AtomicProperties('Sc', 21, 44.956, 1.70, 2.11, 1.36, 3, 6.561, 0.188, (3,), True, 3, 4),
        'Ti': AtomicProperties('Ti', 22, 47.867, 1.60, 2.07, 1.54, 4, 6.828, 0.084, (2, 3, 4), True, 4, 4),
        'V': AtomicProperties('V', 23, 50.942, 1.53, 2.05, 1.63, 5, 6.746, 0.525, (2, 3, 4, 5), True, 5, 4),
        'Cr': AtomicProperties('Cr', 24, 51.996, 1.39, 2.05, 1.66, 6, 6.767, 0.666, (2, 3, 6), True, 6, 4),
        'Mn': AtomicProperties('Mn', 25, 54.938, 1.39, 2.05, 1.55, 7, 7.434, 0.00, (2, 3, 4, 7), True, 7, 4),
        'Fe': AtomicProperties('Fe', 26, 55.845, 1.32, 2.04, 1.83, 8, 7.902, 0.163, (2, 3), True, 8, 4),
        'Co': AtomicProperties('Co', 27, 58.933, 1.26, 2.00, 1.88, 9, 7.881, 0.661, (2, 3), True, 9, 4),
        'Ni': AtomicProperties('Ni', 28, 58.693, 1.24, 1.63, 1.91, 10, 7.640, 1.156, (2, 3), True, 10, 4),
        'Cu': AtomicProperties('Cu', 29, 63.546, 1.32, 1.40, 1.90, 11, 7.726, 1.228, (1, 2), True, 11, 4),
        'Zn': AtomicProperties('Zn', 30, 65.38, 1.22, 1.39, 1.65, 12, 9.394, 0.00, (2,), True, 12, 4),
        'Ga': AtomicProperties('Ga', 31, 69.723, 1.22, 1.87, 1.81, 3, 5.999, 0.300, (3,), True, 13, 4),
        'As': AtomicProperties('As', 33, 74.922, 1.19, 1.85, 2.18, 5, 9.815, 0.810, (-3, 3, 5), False, 15, 4),
        'Br': AtomicProperties('Br', 35, 79.904, 1.20, 1.85, 2.96, 7, 11.814, 3.365, (-1, 1, 3, 5), False, 17, 4),

        # Period 5 - Selected elements
        'Rb': AtomicProperties('Rb', 37, 85.468, 2.20, 3.03, 0.82, 1, 4.177, 0.486, (1,), True, 1, 5),
        'Sr': AtomicProperties('Sr', 38, 87.62, 1.95, 2.49, 0.95, 2, 5.695, 0.048, (2,), True, 2, 5),
        'Ag': AtomicProperties('Ag', 47, 107.87, 1.45, 1.72, 1.93, 11, 7.576, 1.302, (1,), True, 11, 5),
        'Cd': AtomicProperties('Cd', 48, 112.41, 1.44, 1.58, 1.69, 12, 8.994, 0.00, (2,), True, 12, 5),
        'I': AtomicProperties('I', 53, 126.90, 1.39, 1.98, 2.66, 7, 10.451, 3.059, (-1, 1, 3, 5, 7), False, 17, 5),

        # Period 6 - Selected elements
        'Cs': AtomicProperties('Cs', 55, 132.91, 2.44, 3.43, 0.79, 1, 3.894, 0.472, (1,), True, 1, 6),
        'Ba': AtomicProperties('Ba', 56, 137.33, 2.15, 2.68, 0.89, 2, 5.212, 0.145, (2,), True, 2, 6),
        'Pt': AtomicProperties('Pt', 78, 195.08, 1.36, 1.75, 2.28, 10, 9.000, 2.128, (2, 4), True, 10, 6),
        'Au': AtomicProperties('Au', 79, 196.97, 1.36, 1.66, 2.54, 11, 9.226, 2.309, (1, 3), True, 11, 6),
        'Hg': AtomicProperties('Hg', 80, 200.59, 1.32, 1.55, 2.00, 12, 10.438, 0.00, (1, 2), True, 12, 6),
    }

    @classmethod
    def get_element(cls, symbol: str) -> AtomicProperties:
        """
        Get atomic properties for an element.

        Args:
            symbol: Chemical symbol (e.g., 'H', 'C', 'Na')

        Returns:
            AtomicProperties for the element

        Raises:
            KeyError: If element not found
        """
        if symbol not in cls.ELEMENTS:
            raise KeyError(f"Element '{symbol}' not found in periodic table")
        return cls.ELEMENTS[symbol]

    @classmethod
    def electronegativity_difference(cls, sym1: str, sym2: str) -> float:
        """
        Calculate electronegativity difference for bonding classification.

        Args:
            sym1: First element symbol
            sym2: Second element symbol

        Returns:
            Absolute electronegativity difference (Pauling scale)
        """
        en1 = cls.ELEMENTS[sym1].electronegativity
        en2 = cls.ELEMENTS[sym2].electronegativity
        return abs(en1 - en2)

    @classmethod
    def is_metal(cls, symbol: str) -> bool:
        """Check if element is a metal."""
        return cls.ELEMENTS[symbol].is_metal

    @classmethod
    def get_valence_electrons(cls, symbol: str) -> int:
        """Get number of valence electrons."""
        return cls.ELEMENTS[symbol].valence_electrons

    @classmethod
    def get_covalent_radius(cls, symbol: str) -> float:
        """Get covalent radius in Angstroms."""
        return cls.ELEMENTS[symbol].covalent_radius

    @classmethod
    def estimate_bond_length(cls, sym1: str, sym2: str) -> float:
        """
        Estimate equilibrium bond length from covalent radii.

        Args:
            sym1: First element symbol
            sym2: Second element symbol

        Returns:
            Estimated bond length in Angstroms
        """
        r1 = cls.get_covalent_radius(sym1)
        r2 = cls.get_covalent_radius(sym2)
        return r1 + r2

    @classmethod
    def list_elements(cls) -> list:
        """Get list of all available element symbols."""
        return sorted(cls.ELEMENTS.keys())

    @classmethod
    def get_by_atomic_number(cls, z: int) -> Optional[AtomicProperties]:
        """
        Get element by atomic number.

        Args:
            z: Atomic number

        Returns:
            AtomicProperties or None if not found
        """
        for elem in cls.ELEMENTS.values():
            if elem.atomic_number == z:
                return elem
        return None
