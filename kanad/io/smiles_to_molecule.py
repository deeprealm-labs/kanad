"""
SMILES to Molecule Converter

Converts SMILES strings to Kanad Molecule objects using RDKit.
"""

import logging
from typing import Dict, List, Tuple, Optional
import numpy as np

logger = logging.getLogger(__name__)


class SMILESConverter:
    """
    Convert SMILES strings to Kanad molecules.

    Uses RDKit for:
    - SMILES parsing
    - 3D structure generation
    - Geometry optimization

    Usage:
        converter = SMILESConverter()
        molecule = converter.smiles_to_molecule('O', charge=0)
    """

    def __init__(self, optimize_geometry: bool = True):
        """
        Initialize SMILES converter.

        Args:
            optimize_geometry: Whether to optimize 3D geometry
        """
        self.optimize_geometry = optimize_geometry

        # Try to import RDKit
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            self.Chem = Chem
            self.AllChem = AllChem
            logger.info("RDKit initialized for SMILES conversion")
        except ImportError:
            raise ImportError(
                "RDKit required for SMILES input. Install with: pip install rdkit"
            )

    def smiles_to_molecule(
        self,
        smiles: str,
        charge: int = 0,
        multiplicity: int = 1,
        name: Optional[str] = None
    ) -> 'Molecule':
        """
        Convert SMILES string to Kanad Molecule.

        Args:
            smiles: SMILES string (e.g., 'O' for water, 'CCO' for ethanol)
            charge: Molecular charge
            multiplicity: Spin multiplicity (1=singlet, 2=doublet, etc.)
            name: Optional molecule name

        Returns:
            Kanad Molecule object

        Examples:
            >>> converter = SMILESConverter()
            >>> h2o = converter.smiles_to_molecule('O', name='water')
            >>> nacl = converter.smiles_to_molecule('[Na+].[Cl-]', name='NaCl')
            >>> ethanol = converter.smiles_to_molecule('CCO', name='ethanol')
        """
        logger.info(f"Converting SMILES: {smiles}")

        # Parse SMILES
        mol = self.Chem.MolFromSmiles(smiles)

        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        # Add hydrogens
        mol = self.Chem.AddHs(mol)

        # Generate 3D coordinates
        self.AllChem.EmbedMolecule(mol, randomSeed=42)

        # Optimize geometry if requested
        if self.optimize_geometry:
            self.AllChem.MMFFOptimizeMolecule(mol)
            logger.info("Geometry optimized with MMFF force field")

        # Extract atomic data
        atoms = self._extract_atoms(mol)

        # Create Kanad Molecule
        from kanad.core.molecule import Molecule

        # Convert multiplicity to spin (2S)
        # multiplicity = 2S + 1, so spin = multiplicity - 1
        spin = multiplicity - 1

        molecule = Molecule(
            atoms=atoms,
            charge=charge,
            spin=spin
        )

        # Store name as attribute if provided
        if name:
            molecule.name = name

        logger.info(f"Molecule created: {len(atoms)} atoms, charge={charge}, spin={spin}")

        return molecule

    def _extract_atoms(self, rdkit_mol) -> List['Atom']:
        """Extract atom data from RDKit molecule."""
        from kanad.core.atom import Atom

        atoms = []

        conf = rdkit_mol.GetConformer()

        for atom in rdkit_mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)

            atoms.append(Atom(
                symbol=atom.GetSymbol(),
                position=[pos.x, pos.y, pos.z]
            ))

        return atoms

    def common_molecules(self) -> Dict[str, str]:
        """Return dictionary of common molecules and their SMILES."""
        return {
            # Simple molecules
            'H2': '[H][H]',
            'O2': 'O=O',
            'N2': 'N#N',
            'H2O': 'O',
            'H2O2': 'OO',
            'CO': '[C-]#[O+]',
            'CO2': 'O=C=O',
            'NH3': 'N',
            'CH4': 'C',

            # Inorganic
            'HCl': 'Cl',
            'HF': 'F',
            'NaCl': '[Na+].[Cl-]',

            # Organic
            'methanol': 'CO',
            'ethanol': 'CCO',
            'acetone': 'CC(=O)C',
            'benzene': 'c1ccccc1',
            'toluene': 'Cc1ccccc1',

            # Biomolecules (small)
            'glycine': 'C(C(=O)O)N',  # Simplest amino acid
            'alanine': 'C[C@@H](C(=O)O)N',
            'formaldehyde': 'C=O',
            'formic_acid': 'C(=O)O',
            'acetic_acid': 'CC(=O)O'
        }

    def batch_convert(
        self,
        smiles_list: List[str],
        **kwargs
    ) -> List['Molecule']:
        """
        Convert multiple SMILES strings to molecules.

        Args:
            smiles_list: List of SMILES strings
            **kwargs: Arguments passed to smiles_to_molecule

        Returns:
            List of Molecule objects
        """
        molecules = []

        for smiles in smiles_list:
            try:
                mol = self.smiles_to_molecule(smiles, **kwargs)
                molecules.append(mol)
            except Exception as e:
                logger.warning(f"Failed to convert {smiles}: {e}")

        return molecules


def smiles_to_bond(
    smiles: str,
    charge: int = 0,
    bond_type: Optional[str] = None
) -> 'BaseBond':
    """
    Convenience function: SMILES directly to Bond.

    Args:
        smiles: SMILES string
        charge: Molecular charge
        bond_type: Force specific bond type (or auto-detect)

    Returns:
        Bond object

    Examples:
        >>> bond = smiles_to_bond('O')  # Water
        >>> bond = smiles_to_bond('[Na+].[Cl-]')  # NaCl
    """
    from kanad.bonds import BondFactory

    # Convert SMILES to molecule
    converter = SMILESConverter()
    molecule = converter.smiles_to_molecule(smiles, charge=charge)

    # Create bond
    # For multi-atom molecules, we create a generalized molecular bond
    if len(molecule.atoms) == 2:
        # Diatomic - use standard bond factory
        atom1 = molecule.atoms[0]['symbol']
        atom2 = molecule.atoms[1]['symbol']
        distance = np.linalg.norm(
            np.array(molecule.atoms[0]['position']) -
            np.array(molecule.atoms[1]['position'])
        )

        bond = BondFactory.create_bond(
            atom1, atom2,
            distance=distance,
            bond_type=bond_type
        )
    else:
        # Polyatomic - create molecular bond
        # This requires a MolecularBond class (to be implemented)
        logger.warning("Multi-atom molecules not fully supported yet")

        # For now, create a covalent bond using the first two heavy atoms
        heavy_atoms = [a for a in molecule.atoms if a['symbol'] != 'H']

        if len(heavy_atoms) >= 2:
            atom1 = heavy_atoms[0]
            atom2 = heavy_atoms[1]
            distance = np.linalg.norm(
                np.array(atom1['position']) - np.array(atom2['position'])
            )

            bond = BondFactory.create_bond(
                atom1['symbol'],
                atom2['symbol'],
                distance=distance
            )
        else:
            raise ValueError(f"Cannot create bond from SMILES: {smiles}")

    return bond
