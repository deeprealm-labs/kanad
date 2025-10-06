"""
SMILES String Parser for Kanad

Converts SMILES strings to Kanad Molecule objects using RDKit.
Includes 2D→3D structure generation and charge detection.

Examples:
    mol = from_smiles("CCO")  # Ethanol
    mol = from_smiles("c1ccccc1")  # Benzene
    mol = from_smiles("[NH4+]")  # Ammonium ion
"""

import numpy as np
from typing import List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


def from_smiles(
    smiles: str,
    name: Optional[str] = None,
    basis: str = 'sto-3g',
    optimize_geometry: bool = True
) -> 'Molecule':
    """
    Create Kanad Molecule from SMILES string.

    Args:
        smiles: SMILES string (e.g., "CCO", "c1ccccc1")
        name: Optional molecule name
        basis: Basis set for calculations (default: 'sto-3g')
        optimize_geometry: Whether to optimize 3D geometry with force field

    Returns:
        Molecule: Kanad Molecule object

    Raises:
        ValueError: If SMILES string is invalid
        ImportError: If RDKit is not installed

    Examples:
        >>> mol = from_smiles("CCO")  # Ethanol
        >>> print(mol)
        Molecule(C2H6O, 3 atoms, 26 electrons)

        >>> benzene = from_smiles("c1ccccc1", name="Benzene")
        >>> benzene.compute_energy(method='HF')
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        raise ImportError(
            "RDKit is required for SMILES parsing. "
            "Install with: pip install rdkit"
        )

    # Parse SMILES
    mol_rdkit = Chem.MolFromSmiles(smiles)
    if mol_rdkit is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    # Add hydrogens
    mol_rdkit = Chem.AddHs(mol_rdkit)

    # Generate 3D coordinates
    if optimize_geometry:
        # ETKDG: Extended Topological Kernel Distance Geometry
        AllChem.EmbedMolecule(mol_rdkit, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol_rdkit)  # Force field optimization
    else:
        AllChem.EmbedMolecule(mol_rdkit, randomSeed=42)

    # Extract atoms and coordinates
    atoms = smiles_to_atoms(mol_rdkit)

    # Detect molecular charge
    charge = Chem.GetFormalCharge(mol_rdkit)

    # Detect spin (singlet=0, doublet=1, triplet=2, etc.)
    # Count total electrons: sum from atoms (charges already applied by RDKit)
    num_electrons = sum(atom.n_electrons for atom in atoms)

    # Spin multiplicity (2S) where S is total spin
    # Even electrons → spin=0 (singlet, closed-shell)
    # Odd electrons → spin=1 (doublet, one unpaired electron)
    spin = 1 if num_electrons % 2 == 1 else 0

    # Create Kanad Molecule
    from kanad.core.molecule import Molecule
    molecule = Molecule(
        atoms=atoms,
        charge=charge,
        spin=spin,
        basis=basis
    )

    # Store metadata
    molecule._smiles = smiles
    molecule._name = name or smiles

    logger.info(f"Created molecule from SMILES: {smiles}")
    logger.info(f"  Formula: {molecule.formula}")
    logger.info(f"  Atoms: {len(atoms)}, Electrons: {num_electrons}, Charge: {charge}, Spin: {spin}")

    return molecule


def smiles_to_atoms(mol_rdkit) -> List['Atom']:
    """
    Convert RDKit Mol object to list of Kanad Atoms.

    Args:
        mol_rdkit: RDKit Mol object with 3D coordinates

    Returns:
        List[Atom]: List of Kanad Atom objects
    """
    from kanad.core.atom import Atom

    atoms = []
    conformer = mol_rdkit.GetConformer()

    for rdkit_atom in mol_rdkit.GetAtoms():
        # Get atomic symbol
        symbol = rdkit_atom.GetSymbol()

        # Get 3D position (in Angstroms)
        idx = rdkit_atom.GetIdx()
        pos = conformer.GetAtomPosition(idx)
        position = np.array([pos.x, pos.y, pos.z])

        # Get formal charge on this atom
        atom_charge = rdkit_atom.GetFormalCharge()

        # Create Kanad Atom
        atom = Atom(symbol=symbol, position=position, charge=atom_charge)
        atoms.append(atom)

    return atoms


def smiles_to_formula(smiles: str) -> str:
    """
    Get molecular formula from SMILES.

    Args:
        smiles: SMILES string

    Returns:
        str: Molecular formula (e.g., "C2H6O")

    Examples:
        >>> smiles_to_formula("CCO")
        'C2H6O'
    """
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError("RDKit required for SMILES parsing")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    return Chem.rdMolDescriptors.CalcMolFormula(mol)


def validate_smiles(smiles: str) -> Tuple[bool, str]:
    """
    Validate SMILES string.

    Args:
        smiles: SMILES string to validate

    Returns:
        Tuple[bool, str]: (is_valid, error_message)

    Examples:
        >>> validate_smiles("CCO")
        (True, "")
        >>> validate_smiles("invalid")
        (False, "Invalid SMILES syntax")
    """
    try:
        from rdkit import Chem
    except ImportError:
        return False, "RDKit not installed"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES syntax"

    return True, ""
