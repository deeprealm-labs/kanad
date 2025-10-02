"""
Base class for all bond types.

Provides common interface for ionic, covalent, and metallic bonds.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List
import numpy as np

from kanad.core.atom import Atom


class BaseBond(ABC):
    """
    Abstract base class for chemical bonds.

    All bond types (ionic, covalent, metallic) inherit from this class
    and implement the compute_energy and analyze methods.
    """

    def __init__(
        self,
        atoms: List[Atom],
        bond_type: str,
        distance: Optional[float] = None
    ):
        """
        Initialize base bond.

        Args:
            atoms: List of atoms in the bond
            bond_type: Type of bond ('ionic', 'covalent', 'metallic')
            distance: Bond distance in Angstroms (optional)
        """
        self.atoms = atoms
        self.bond_type = bond_type
        self.distance = distance
        self.hamiltonian = None
        self.governance = None
        self.mapper = None
        self.representation = None

    @abstractmethod
    def compute_energy(
        self,
        method: str = 'VQE',
        **kwargs
    ) -> Dict[str, Any]:
        """
        Compute bond energy.

        Args:
            method: Computational method ('VQE', 'exact', 'HF')
            **kwargs: Method-specific parameters

        Returns:
            Dictionary with energy and analysis results
        """
        pass

    @abstractmethod
    def analyze(self) -> Dict[str, Any]:
        """
        Analyze bond properties.

        Returns:
            Dictionary with bond analysis (bond order, charge transfer, etc.)
        """
        pass

    def get_bond_length(self) -> float:
        """
        Get bond length.

        Returns:
            Bond length in Angstroms
        """
        if self.distance is not None:
            return self.distance
        elif len(self.atoms) >= 2:
            # Calculate from atom positions
            return np.linalg.norm(
                self.atoms[0].position - self.atoms[1].position
            ) * 0.529177  # Bohr to Angstrom
        else:
            return 0.0

    def get_atoms(self) -> List[Atom]:
        """Get atoms in the bond."""
        return self.atoms

    def get_bond_type(self) -> str:
        """Get bond type."""
        return self.bond_type

    def __repr__(self) -> str:
        """String representation."""
        atom_symbols = '-'.join([atom.symbol for atom in self.atoms])
        return f"{self.__class__.__name__}({atom_symbols}, {self.bond_type})"
