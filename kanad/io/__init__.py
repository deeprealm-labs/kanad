"""
Kanad Molecular I/O Module

Provides parsers and writers for standard chemistry file formats:
- SMILES strings (from_smiles)
- InChI strings (from_inchi)
- XYZ coordinate files (from_xyz, to_xyz)
- PDB protein structures (from_pdb, to_pdb)
- MOL/SDF files (from_mol, to_mol)

Usage:
    from kanad.io import from_smiles, to_xyz

    # Parse SMILES
    mol = from_smiles("CCO")  # Ethanol

    # Load from file
    mol = from_xyz("molecule.xyz")

    # Save to file
    to_xyz(mol, "output.xyz")
"""

from kanad.io.smiles_parser import from_smiles, smiles_to_atoms
from kanad.io.xyz_io import from_xyz, to_xyz
from kanad.io.crystal_builder import build_crystal, build_binary_crystal, get_kpath, get_lattice_info

__all__ = [
    'from_smiles',
    'smiles_to_atoms',
    'from_xyz',
    'to_xyz',
    'build_crystal',
    'build_binary_crystal',
    'get_kpath',
    'get_lattice_info',
]
