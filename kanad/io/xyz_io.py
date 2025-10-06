"""
XYZ File Format I/O for Kanad

Read and write XYZ coordinate files.

XYZ Format:
    Line 1: Number of atoms
    Line 2: Comment line (molecule name, energy, etc.)
    Line 3+: Symbol X Y Z (Angstroms)

Example XYZ file:
    3
    Water molecule
    O   0.000   0.000   0.119
    H   0.000   0.763  -0.477
    H   0.000  -0.763  -0.477
"""

import numpy as np
from typing import List, Optional
import logging

logger = logging.getLogger(__name__)


def from_xyz(
    filename: str,
    charge: int = 0,
    spin: int = 0,
    basis: str = 'sto-3g'
) -> 'Molecule':
    """
    Read molecule from XYZ file.

    Args:
        filename: Path to XYZ file
        charge: Molecular charge (default: 0)
        spin: Spin multiplicity 2S (default: 0 for singlet)
        basis: Basis set (default: 'sto-3g')

    Returns:
        Molecule: Kanad Molecule object

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is invalid

    Examples:
        >>> mol = from_xyz("water.xyz")
        >>> mol = from_xyz("h2.xyz", charge=0, spin=0)
    """
    from kanad.core.atom import Atom
    from kanad.core.molecule import Molecule

    try:
        with open(filename, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
    except FileNotFoundError:
        raise FileNotFoundError(f"XYZ file not found: {filename}")

    if len(lines) < 3:
        raise ValueError(f"Invalid XYZ file: too few lines ({len(lines)})")

    # Parse header
    try:
        n_atoms = int(lines[0])
    except ValueError:
        raise ValueError(f"Invalid XYZ file: first line must be number of atoms")

    comment = lines[1] if len(lines) > 1 else ""

    # Parse atoms
    atoms = []
    for i, line in enumerate(lines[2:2+n_atoms], start=1):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ line {i+2}: expected 'Symbol X Y Z', got '{line}'")

        symbol = parts[0]
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        except ValueError:
            raise ValueError(f"Invalid coordinates on line {i+2}: {line}")

        position = np.array([x, y, z])
        atom = Atom(symbol=symbol, position=position)
        atoms.append(atom)

    if len(atoms) != n_atoms:
        logger.warning(f"Expected {n_atoms} atoms, found {len(atoms)}")

    # Create molecule
    molecule = Molecule(atoms=atoms, charge=charge, spin=spin, basis=basis)
    molecule._name = comment or filename
    molecule._xyz_file = filename

    logger.info(f"Loaded molecule from {filename}")
    logger.info(f"  {n_atoms} atoms, formula: {molecule.formula}")

    return molecule


def to_xyz(
    molecule: 'Molecule',
    filename: str,
    comment: Optional[str] = None,
    include_energy: bool = False
):
    """
    Write molecule to XYZ file.

    Args:
        molecule: Kanad Molecule object
        filename: Output file path
        comment: Comment line (default: molecule name or formula)
        include_energy: Whether to include energy in comment line

    Examples:
        >>> to_xyz(mol, "output.xyz")
        >>> to_xyz(mol, "optimized.xyz", comment="Optimized geometry")
        >>> to_xyz(mol, "result.xyz", include_energy=True)
    """
    n_atoms = len(molecule.atoms)

    # Generate comment line
    if comment is None:
        comment = getattr(molecule, '_name', molecule.formula)

    if include_energy and hasattr(molecule, '_last_energy'):
        comment += f" Energy: {molecule._last_energy:.8f} Ha"

    # Write file
    with open(filename, 'w') as f:
        # Header
        f.write(f"{n_atoms}\n")
        f.write(f"{comment}\n")

        # Atom coordinates
        for atom in molecule.atoms:
            x, y, z = atom.position
            f.write(f"{atom.symbol:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")

    logger.info(f"Wrote molecule to {filename}")


def xyz_to_string(molecule: 'Molecule', comment: str = "") -> str:
    """
    Convert molecule to XYZ format string.

    Args:
        molecule: Kanad Molecule object
        comment: Comment line

    Returns:
        str: XYZ format string

    Examples:
        >>> xyz_str = xyz_to_string(mol, "Water")
        >>> print(xyz_str)
    """
    n_atoms = len(molecule.atoms)
    lines = [f"{n_atoms}", comment]

    for atom in molecule.atoms:
        x, y, z = atom.position
        lines.append(f"{atom.symbol:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}")

    return "\n".join(lines)


def parse_xyz_string(xyz_string: str, charge: int = 0, spin: int = 0, basis: str = 'sto-3g') -> 'Molecule':
    """
    Parse XYZ format string to Molecule.

    Args:
        xyz_string: XYZ format string
        charge: Molecular charge
        spin: Spin multiplicity
        basis: Basis set

    Returns:
        Molecule: Kanad Molecule object
    """
    from kanad.core.atom import Atom
    from kanad.core.molecule import Molecule

    lines = [line.strip() for line in xyz_string.split('\n') if line.strip()]

    n_atoms = int(lines[0])
    comment = lines[1] if len(lines) > 1 else ""

    atoms = []
    for line in lines[2:2+n_atoms]:
        parts = line.split()
        symbol = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        position = np.array([x, y, z])
        atom = Atom(symbol=symbol, position=position)
        atoms.append(atom)

    molecule = Molecule(atoms=atoms, charge=charge, spin=spin, basis=basis)
    molecule._name = comment
    return molecule
