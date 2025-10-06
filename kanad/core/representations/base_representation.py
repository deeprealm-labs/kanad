"""
Base representation for quantum chemical systems.

Each bonding type (ionic, covalent, metallic) requires a different
quantum representation strategy.
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
import numpy as np



class BaseRepresentation(ABC):
    """
    Abstract base class for quantum representations of bonding.

    Different bonding types require fundamentally different representations:
    - Ionic: Localized, minimal entanglement
    - Covalent: Hybrid orbitals, paired entanglement
    - Metallic: Delocalized k-space, collective entanglement
    """

    def __init__(self, molecule: 'Molecule'):
        """
        Initialize representation.

        Args:
            molecule: Molecule object with atoms and bonds
        """
        self.molecule = molecule
        self.n_qubits: Optional[int] = None
        self.n_electrons = sum(atom.n_electrons for atom in molecule.atoms)
        self.n_orbitals: Optional[int] = None

    @abstractmethod
    def build_hamiltonian(self) -> 'Hamiltonian':
        """
        Build the Hamiltonian in this representation.

        Returns:
            Hamiltonian object specific to this representation type
        """
        pass

    @abstractmethod
    def get_reference_state(self) -> np.ndarray:
        """
        Get the reference state (e.g., Hartree-Fock state).

        Returns:
            Reference state vector
        """
        pass

    @abstractmethod
    def compute_observables(self, state: np.ndarray) -> Dict[str, float]:
        """
        Compute physical observables from quantum state.

        Args:
            state: Quantum state vector

        Returns:
            Dictionary of observable names to values
        """
        pass

    @abstractmethod
    def to_qubit_operator(self) -> 'QubitOperator':
        """
        Map representation to qubit operators.

        Returns:
            QubitOperator representing the Hamiltonian
        """
        pass

    @abstractmethod
    def get_num_qubits(self) -> int:
        """
        Get number of qubits required for this representation.

        Returns:
            Number of qubits
        """
        pass

    def get_num_electrons(self) -> int:
        """Get total number of electrons."""
        return self.n_electrons

    def get_num_orbitals(self) -> int:
        """Get number of spatial orbitals."""
        if self.n_orbitals is None:
            raise ValueError("Number of orbitals not set")
        return self.n_orbitals


class Molecule:
    """
    Simple molecule class for holding atoms and geometry.

    This is a minimal implementation for Phase 3.
    """

    def __init__(self, atoms: List['Atom'], bonds: Optional[List] = None, spin: int = 0):
        """
        Initialize molecule.

        Args:
            atoms: List of Atom objects
            bonds: Optional list of Bond objects
            spin: Spin multiplicity (2S, where S is total spin)
        """
        self.atoms = atoms
        self.bonds = bonds or []
        self.n_atoms = len(atoms)
        self.spin = spin

    @property
    def n_electrons(self) -> int:
        """Total number of electrons."""
        return sum(atom.n_electrons for atom in self.atoms)

    @property
    def symbols(self) -> List[str]:
        """List of atomic symbols."""
        return [atom.symbol for atom in self.atoms]

    @property
    def positions(self) -> np.ndarray:
        """Atomic positions as (n_atoms, 3) array."""
        return np.array([atom.position for atom in self.atoms])

    def distance_matrix(self) -> np.ndarray:
        """
        Compute distance matrix between all atoms.

        Returns:
            (n_atoms, n_atoms) distance matrix in Angstroms
        """
        n = self.n_atoms
        D = np.zeros((n, n))

        for i in range(n):
            for j in range(i + 1, n):
                d = self.atoms[i].distance_to(self.atoms[j])
                D[i, j] = d
                D[j, i] = d

        return D

    def __repr__(self) -> str:
        """String representation."""
        formula = ''.join(self.symbols)
        return f"Molecule({formula}, n_atoms={self.n_atoms}, n_electrons={self.n_electrons})"
