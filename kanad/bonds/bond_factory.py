"""
Bond Factory - User-facing API for creating bonds.

Provides simple interface for researchers to create bonds without
worrying about low-level details. Automatically determines bond type
from electronegativity and applies appropriate governance.
"""

from typing import Union, List, Dict, Any, Optional
from enum import Enum
import numpy as np

from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule


class BondType(Enum):
    """Types of chemical bonds."""
    IONIC = "ionic"
    COVALENT = "covalent"
    METALLIC = "metallic"
    AUTO = "auto"  # Automatically determine from properties


class BondFactory:
    """
    Factory for creating bonds with automatic governance.

    USER-FACING API - This is what researchers interact with!

    Features:
    - Automatic bond type detection from electronegativity
    - Governance protocol selection
    - Simplified molecule creation
    - Preset geometries for common molecules
    """

    # Electronegativity thresholds
    EN_IONIC_THRESHOLD = 1.7  # ΔEN > 1.7 → ionic
    EN_POLAR_THRESHOLD = 0.4  # ΔEN > 0.4 → polar covalent

    @staticmethod
    def create_bond(
        atom_1: Union[str, Atom],
        atom_2: Union[str, Atom],
        bond_type: Union[BondType, str] = BondType.AUTO,
        distance: Optional[float] = None,
        **kwargs
    ) -> 'BaseBond':
        """
        Create a bond between two atoms.

        Args:
            atom_1: First atom (symbol like 'H' or Atom object)
            atom_2: Second atom (symbol like 'Cl' or Atom object)
            bond_type: Type of bond ('ionic', 'covalent', 'metallic', or 'auto')
            distance: Bond distance in Angstroms (optional, will use defaults)
            **kwargs: Additional parameters for specific bond types

        Returns:
            Bond object (IonicBond, CovalentBond, or MetallicBond)

        Examples:
            >>> # Auto-detect bond type
            >>> bond = BondFactory.create_bond('Na', 'Cl')
            >>> print(bond.bond_type)  # 'ionic'

            >>> # Explicit bond type
            >>> bond = BondFactory.create_bond('H', 'H', bond_type='covalent')

            >>> # With specific distance
            >>> bond = BondFactory.create_bond('C', 'C', distance=1.54)
        """
        # Import here to avoid circular imports
        from kanad.bonds.ionic_bond import IonicBond
        from kanad.bonds.covalent_bond import CovalentBond
        from kanad.bonds.metallic_bond import MetallicBond

        # Convert strings to Atom objects if needed
        if isinstance(atom_1, str):
            atom_1 = Atom(atom_1, position=np.array([0.0, 0.0, 0.0]))
        if isinstance(atom_2, str):
            # Create temporary atom for distance estimation
            atom_2_temp = Atom(atom_2, position=np.zeros(3))
            # Default: place second atom along z-axis
            if distance is None:
                distance = BondFactory._estimate_bond_length(atom_1, atom_2_temp)
            atom_2 = Atom(atom_2, position=np.array([0.0, 0.0, distance]))

        # Convert string bond type to enum
        if isinstance(bond_type, str):
            bond_type = BondType(bond_type.lower())

        # Auto-determine bond type if requested
        if bond_type == BondType.AUTO:
            bond_type = BondFactory._determine_bond_type(atom_1, atom_2)

        # Create appropriate bond with governance
        if bond_type == BondType.IONIC:
            return IonicBond(atom_1, atom_2, distance=distance, **kwargs)
        elif bond_type == BondType.COVALENT:
            return CovalentBond(atom_1, atom_2, distance=distance, **kwargs)
        elif bond_type == BondType.METALLIC:
            return MetallicBond([atom_1, atom_2], **kwargs)
        else:
            raise ValueError(f"Unknown bond type: {bond_type}")

    @staticmethod
    def _determine_bond_type(atom_1: Atom, atom_2: Atom) -> BondType:
        """
        Automatically determine bond type from atomic properties.

        Rules (Pauling electronegativity scale):
            - ΔEN > 1.7: Ionic bonding
            - 0.4 < ΔEN ≤ 1.7: Polar/Nonpolar covalent
            - Both metals: Metallic bonding

        Args:
            atom_1: First atom
            atom_2: Second atom

        Returns:
            BondType enum (IONIC, COVALENT, or METALLIC)
        """
        # Get electronegativity values
        en1 = atom_1.electronegativity
        en2 = atom_2.electronegativity
        delta_en = abs(en1 - en2)

        # Check if both are metals
        is_metal_1 = atom_1.is_metal
        is_metal_2 = atom_2.is_metal

        if is_metal_1 and is_metal_2:
            # Both metals → metallic bonding
            return BondType.METALLIC
        elif delta_en > BondFactory.EN_IONIC_THRESHOLD:
            # Large EN difference → ionic bonding
            return BondType.IONIC
        else:
            # Small EN difference → covalent bonding
            return BondType.COVALENT

    @staticmethod
    def _estimate_bond_length(atom_1: Atom, atom_2: Atom) -> float:
        """
        Estimate bond length from covalent radii.

        Args:
            atom_1: First atom
            atom_2: Second atom

        Returns:
            Estimated bond length in Angstroms
        """
        # Sum of covalent radii (already in Angstroms)
        return atom_1.covalent_radius + atom_2.covalent_radius

    @staticmethod
    def create_molecule(
        atoms: List[Union[str, Atom]],
        geometry: Union[str, np.ndarray] = 'auto',
        bond_types: Optional[List[Union[BondType, str]]] = None,
        basis: str = 'sto-3g'
    ) -> Molecule:
        """
        Create a molecule with multiple atoms.

        Args:
            atoms: List of atoms (symbols or Atom objects)
            geometry: Geometry specification:
                - 'auto': Optimize geometry
                - 'linear': Linear arrangement
                - 'water', 'methane', etc.: Preset geometries
                - np.ndarray: Custom positions (n_atoms, 3)
            bond_types: Optional list of bond types for each bond
            basis: Basis set name (default: 'sto-3g')

        Returns:
            Molecule object

        Examples:
            >>> # Simple diatomic
            >>> mol = BondFactory.create_molecule(['H', 'H'])

            >>> # Water molecule with preset geometry
            >>> mol = BondFactory.create_molecule(['H', 'O', 'H'], geometry='water')

            >>> # Custom positions
            >>> positions = np.array([[0,0,0], [0,0,1.4], [0,0,2.8]])
            >>> mol = BondFactory.create_molecule(['H', 'H', 'H'], geometry=positions)
        """
        # Convert strings to Atom objects
        atom_objects = []
        for i, a in enumerate(atoms):
            if isinstance(a, str):
                # Temporary position (will be updated by geometry)
                atom_objects.append(Atom(a, position=np.zeros(3)))
            else:
                atom_objects.append(a)

        # Determine geometry
        if isinstance(geometry, np.ndarray):
            # Custom positions provided
            positions = geometry
        elif geometry == 'linear':
            # Linear arrangement along z-axis
            positions = BondFactory._generate_linear_geometry(atom_objects)
        elif geometry in ['water', 'h2o']:
            # Water molecule (H-O-H)
            positions = BondFactory._generate_water_geometry()
        elif geometry in ['methane', 'ch4']:
            # Methane (tetrahedral)
            positions = BondFactory._generate_methane_geometry()
        elif geometry == 'auto':
            # Simple linear for now
            positions = BondFactory._generate_linear_geometry(atom_objects)
        else:
            raise ValueError(f"Unknown geometry: {geometry}")

        # Assign positions
        if len(positions) != len(atom_objects):
            raise ValueError(
                f"Geometry has {len(positions)} positions but {len(atom_objects)} atoms"
            )

        for atom, pos in zip(atom_objects, positions):
            atom.position = pos

        # Create Molecule object
        return Molecule(atom_objects)

    @staticmethod
    def _generate_linear_geometry(atoms: List[Atom]) -> np.ndarray:
        """Generate linear geometry along z-axis."""
        n_atoms = len(atoms)
        positions = np.zeros((n_atoms, 3))

        # Space atoms along z-axis
        cumulative_distance = 0.0
        for i in range(1, n_atoms):
            # Estimate bond length between consecutive atoms
            bond_length = BondFactory._estimate_bond_length(atoms[i-1], atoms[i])
            cumulative_distance += bond_length
            positions[i, 2] = cumulative_distance

        return positions

    @staticmethod
    def _generate_water_geometry() -> np.ndarray:
        """
        Generate water molecule geometry.

        H-O-H angle: 104.5 degrees
        O-H bond length: 0.96 Å
        """
        angle_rad = np.radians(104.5 / 2)  # Half angle from z-axis
        oh_distance = 0.96  # Angstroms

        positions = np.array([
            [oh_distance * np.sin(angle_rad), 0.0, oh_distance * np.cos(angle_rad)],  # H1
            [0.0, 0.0, 0.0],  # O at origin
            [-oh_distance * np.sin(angle_rad), 0.0, oh_distance * np.cos(angle_rad)]  # H2
        ])

        return positions

    @staticmethod
    def _generate_methane_geometry() -> np.ndarray:
        """
        Generate methane geometry (tetrahedral).

        Tetrahedral angle: 109.47 degrees
        C-H bond length: 1.09 Å
        """
        ch_distance = 1.09  # Angstroms

        # Tetrahedral vertices (normalized)
        # Place C at origin, H at tetrahedral vertices
        positions = np.array([
            [0.0, 0.0, 0.0],  # C at origin
            [1.0, 1.0, 1.0],   # H1
            [1.0, -1.0, -1.0], # H2
            [-1.0, 1.0, -1.0], # H3
            [-1.0, -1.0, 1.0]  # H4
        ])

        # Normalize and scale H positions
        for i in range(1, 5):
            positions[i] = positions[i] / np.linalg.norm(positions[i]) * ch_distance

        return positions

    @staticmethod
    def quick_bond_info(atom_1: str, atom_2: str) -> Dict[str, Any]:
        """
        Get quick information about a potential bond WITHOUT computing energy.

        Useful for exploring bonding possibilities before expensive calculations.

        Args:
            atom_1: First atom symbol
            atom_2: Second atom symbol

        Returns:
            Dictionary with bond information

        Example:
            >>> info = BondFactory.quick_bond_info('Na', 'Cl')
            >>> print(info['predicted_type'])  # 'ionic'
            >>> print(info['electronegativity_difference'])  # 2.23
        """
        # Create temporary atoms
        a1 = Atom(atom_1, position=np.zeros(3))
        a2 = Atom(atom_2, position=np.zeros(3))

        # Determine bond type
        bond_type = BondFactory._determine_bond_type(a1, a2)

        # Calculate properties
        en_diff = abs(a1.electronegativity - a2.electronegativity)
        estimated_length = BondFactory._estimate_bond_length(a1, a2)

        return {
            'atom_1': atom_1,
            'atom_2': atom_2,
            'predicted_type': bond_type.value,
            'electronegativity_difference': en_diff,
            'estimated_bond_length': estimated_length,
            'electronegativity_1': a1.electronegativity,
            'electronegativity_2': a2.electronegativity,
            'is_metal_1': a1.is_metal,
            'is_metal_2': a2.is_metal,
            'rationale': BondFactory._explain_bond_type(a1, a2, bond_type)
        }

    @staticmethod
    def _explain_bond_type(
        atom_1: Atom,
        atom_2: Atom,
        bond_type: BondType
    ) -> str:
        """Generate human-readable explanation of bond type selection."""
        en_diff = abs(atom_1.electronegativity - atom_2.electronegativity)

        if bond_type == BondType.METALLIC:
            return (f"Both {atom_1.symbol} and {atom_2.symbol} are metals, "
                   f"forming metallic bonding with delocalized electrons.")
        elif bond_type == BondType.IONIC:
            return (f"Electronegativity difference ({en_diff:.2f}) exceeds {BondFactory.EN_IONIC_THRESHOLD}, "
                   f"indicating significant electron transfer (ionic bonding).")
        else:  # COVALENT
            return (f"Electronegativity difference ({en_diff:.2f}) is moderate, "
                   f"indicating electron sharing (covalent bonding).")
