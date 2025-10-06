"""
Unit tests for XYZ file I/O

Tests cover:
- Reading XYZ files
- Writing XYZ files
- Round-trip (write then read)
- Invalid file handling
- String conversion
"""

import pytest
import numpy as np
import tempfile
import os


def test_xyz_write_read_roundtrip():
    """Test writing and reading XYZ file (round-trip)."""
    from kanad.core.atom import Atom
    from kanad.core.molecule import Molecule
    from kanad.io import to_xyz, from_xyz

    # Create H2 molecule
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
    mol_original = Molecule(atoms=[h1, h2], charge=0, spin=0)

    # Write to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        temp_file = f.name

    try:
        to_xyz(mol_original, temp_file, comment="Hydrogen molecule")

        # Read back
        mol_loaded = from_xyz(temp_file)

        # Verify
        assert len(mol_loaded.atoms) == 2
        assert mol_loaded.atoms[0].symbol == 'H'
        assert mol_loaded.atoms[1].symbol == 'H'

        # Check coordinates match (within tolerance)
        for i in range(2):
            np.testing.assert_array_almost_equal(
                mol_loaded.atoms[i].position,
                mol_original.atoms[i].position,
                decimal=5
            )

    finally:
        os.unlink(temp_file)


def test_xyz_read_water():
    """Test reading water XYZ file."""
    from kanad.io import from_xyz
    import tempfile

    xyz_content = """3
Water molecule
O   0.000   0.000   0.119
H   0.000   0.763  -0.477
H   0.000  -0.763  -0.477
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(xyz_content)
        temp_file = f.name

    try:
        mol = from_xyz(temp_file)

        assert len(mol.atoms) == 3
        assert mol.atoms[0].symbol == 'O'
        assert mol.atoms[1].symbol == 'H'
        assert mol.atoms[2].symbol == 'H'

        # Check oxygen position
        np.testing.assert_array_almost_equal(
            mol.atoms[0].position,
            np.array([0.0, 0.0, 0.119]),
            decimal=3
        )

    finally:
        os.unlink(temp_file)


def test_xyz_write_with_energy():
    """Test writing XYZ with energy in comment."""
    from kanad.core.atom import Atom
    from kanad.core.molecule import Molecule
    from kanad.io import to_xyz
    import tempfile

    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
    mol = Molecule(atoms=[h1, h2], charge=0, spin=0)
    mol._last_energy = -1.133  # Store energy

    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        temp_file = f.name

    try:
        to_xyz(mol, temp_file, include_energy=True)

        # Read file and check comment line
        with open(temp_file, 'r') as f:
            lines = f.readlines()

        assert "Energy:" in lines[1]
        assert "-1.133" in lines[1]

    finally:
        os.unlink(temp_file)


def test_xyz_to_string():
    """Test converting molecule to XYZ string."""
    from kanad.core.atom import Atom
    from kanad.core.molecule import Molecule
    from kanad.io.xyz_io import xyz_to_string

    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
    mol = Molecule(atoms=[h1, h2], charge=0, spin=0)

    xyz_str = xyz_to_string(mol, comment="Test molecule")

    lines = xyz_str.split('\n')
    assert lines[0] == "2"
    assert lines[1] == "Test molecule"
    assert "H" in lines[2]
    assert "0.740000" in lines[3]


def test_parse_xyz_string():
    """Test parsing XYZ from string."""
    from kanad.io.xyz_io import parse_xyz_string

    xyz_str = """2
Hydrogen
H   0.000   0.000   0.000
H   0.740   0.000   0.000
"""

    mol = parse_xyz_string(xyz_str)

    assert len(mol.atoms) == 2
    assert mol.atoms[0].symbol == 'H'
    np.testing.assert_array_almost_equal(
        mol.atoms[1].position,
        np.array([0.74, 0.0, 0.0])
    )


def test_xyz_invalid_file():
    """Test error handling for invalid XYZ file."""
    from kanad.io import from_xyz

    with pytest.raises(FileNotFoundError):
        from_xyz("nonexistent_file.xyz")


def test_xyz_invalid_format():
    """Test error for invalid XYZ format."""
    from kanad.io import from_xyz
    import tempfile

    # Too few lines
    xyz_content = "2\n"

    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(xyz_content)
        temp_file = f.name

    try:
        with pytest.raises(ValueError, match="too few lines"):
            from_xyz(temp_file)
    finally:
        os.unlink(temp_file)


def test_xyz_invalid_coordinates():
    """Test error for invalid coordinate values."""
    from kanad.io import from_xyz
    import tempfile

    xyz_content = """1
Test
H   invalid   0.0   0.0
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(xyz_content)
        temp_file = f.name

    try:
        with pytest.raises(ValueError, match="Invalid coordinates"):
            from_xyz(temp_file)
    finally:
        os.unlink(temp_file)


def test_xyz_with_charge_spin():
    """Test XYZ reading with custom charge and spin."""
    from kanad.io import from_xyz
    import tempfile

    xyz_content = """1
H atom
H   0.0   0.0   0.0
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(xyz_content)
        temp_file = f.name

    try:
        mol = from_xyz(temp_file, charge=1, spin=0)
        assert mol.charge == 1
        assert mol.spin == 0
    finally:
        os.unlink(temp_file)
