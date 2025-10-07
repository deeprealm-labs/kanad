"""
Atom class representing individual atoms in molecules.
"""

import numpy as np
from typing import Optional
from kanad.core.constants.atomic_data import PeriodicTable, AtomicProperties


class Atom:
    """
    Represents an atom in a molecular system.

    Attributes:
        symbol: Chemical symbol (e.g., 'H', 'C', 'Na')
        position: 3D coordinates in Angstroms
        properties: Atomic properties from periodic table
    """

    def __init__(
        self,
        symbol: str,
        position: Optional[np.ndarray] = None,
        charge: int = 0
    ):
        """
        Initialize an atom.

        Args:
            symbol: Chemical symbol
            position: Position in 3D space (Angstroms), defaults to origin
            charge: Formal charge on the atom
        """
        self.symbol = symbol
        # Always store position as numpy array
        if position is not None:
            self.position = np.array(position) if not isinstance(position, np.ndarray) else position
        else:
            self.position = np.zeros(3)
        self.charge = charge

        # Get atomic properties from periodic table
        self.properties = PeriodicTable.get_element(symbol)

    @property
    def atomic_number(self) -> int:
        """Get atomic number (Z)."""
        return self.properties.atomic_number

    @property
    def atomic_mass(self) -> float:
        """Get atomic mass in amu."""
        return self.properties.atomic_mass

    @property
    def n_electrons(self) -> int:
        """Get number of electrons (considering charge)."""
        return self.atomic_number - self.charge

    @property
    def n_valence(self) -> int:
        """Get number of valence electrons."""
        return self.properties.valence_electrons

    @property
    def electronegativity(self) -> float:
        """Get Pauling electronegativity."""
        return self.properties.electronegativity

    @property
    def covalent_radius(self) -> float:
        """Get covalent radius in Angstroms."""
        return self.properties.covalent_radius

    @property
    def is_metal(self) -> bool:
        """Check if atom is a metal."""
        return self.properties.is_metal

    def distance_to(self, other: 'Atom') -> float:
        """
        Calculate distance to another atom in Angstroms.

        Args:
            other: Another Atom object

        Returns:
            Distance in Angstroms
        """
        return np.linalg.norm(self.position - other.position)

    def __repr__(self) -> str:
        """String representation."""
        x, y, z = self.position
        return f"Atom({self.symbol}, pos=[{x:.3f}, {y:.3f}, {z:.3f}])"

    def __str__(self) -> str:
        """Human-readable string."""
        return f"{self.symbol} at ({self.position[0]:.3f}, {self.position[1]:.3f}, {self.position[2]:.3f}) Å"
