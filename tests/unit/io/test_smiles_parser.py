"""
Unit tests for SMILES parser

Tests cover:
- Simple molecules (H2, H2O, CH4)
- Aromatic molecules (benzene)
- Charged species ([NH4+], [Cl-])
- Invalid SMILES handling
- Formula extraction
- Validation
"""

import pytest
import numpy as np


def test_import():
    """Test that SMILES parser can be imported."""
    from kanad.io import from_smiles
    assert from_smiles is not None


def test_simple_smiles_h2():
    """Test H2 parsing from SMILES."""
    pytest.importorskip("rdkit")
    from kanad.io import from_smiles

    mol = from_smiles("[H][H]", name="Hydrogen")

    assert mol is not None
    assert len(mol.atoms) == 2
    assert mol.atoms[0].symbol == 'H'
    assert mol.atoms[1].symbol == 'H'
    assert mol.n_electrons == 2
    assert mol.charge == 0
    assert mol.spin == 0  # Singlet


def test_simple_smiles_water():
    """Test H2O parsing from SMILES."""
    pytest.importorskip("rdkit")
    from kanad.io import from_smiles

    mol = from_smiles("O", name="Water")  # SMILES for water

    assert len(mol.atoms) == 3  # O + 2H (implicit H added)
    assert mol.formula == "H2O"
    assert mol.n_electrons == 10
    assert mol.charge == 0

    # Check that geometry is reasonable (not all atoms at origin)
    positions = np.array([atom.position for atom in mol.atoms])
    distances = np.linalg.norm(positions[1:] - positions[0], axis=1)
    assert all(d > 0.5 for d in distances)  # O-H bond ~0.96 Å


def test_simple_smiles_ethanol():
    """Test ethanol (CCO) parsing."""
    pytest.importorskip("rdkit")
    from kanad.io import from_smiles

    mol = from_smiles("CCO", name="Ethanol")

    assert mol.formula == "C2H6O"
    assert len(mol.atoms) == 9  # 2C + 6H + 1O
    assert mol.n_electrons == 26


def test_simple_smiles_methane():
    """Test CH4 parsing."""
    pytest.importorskip("rdkit")
    from kanad.io import from_smiles

    mol = from_smiles("C", name="Methane")

    assert mol.formula == "CH4"
    assert len(mol.atoms) == 5
    assert mol.n_electrons == 10


def test_aromatic_benzene():
    """Test benzene aromatic SMILES."""
    pytest.importorskip("rdkit")
    from kanad.io import from_smiles

    mol = from_smiles("c1ccccc1", name="Benzene")

    assert mol.formula == "C6H6"
    assert len(mol.atoms) == 12
    assert mol.n_electrons == 42


def test_charged_ammonium():
    """Test charged species: NH4+."""
    pytest.importorskip("rdkit")
    from kanad.io import from_smiles

    mol = from_smiles("[NH4+]", name="Ammonium")

    assert mol.formula == "H4N+"  # RDKit format
    assert mol.charge == 1
    assert mol.n_electrons == 10  # 11 (N) + 4 (H) - 1 (charge)


def test_charged_chloride():
    """Test Cl- anion."""
    pytest.importorskip("rdkit")
    from kanad.io import from_smiles

    mol = from_smiles("[Cl-]", name="Chloride")

    assert mol.charge == -1
    assert mol.n_electrons == 18  # 17 (Cl) + 1 (extra electron)


def test_invalid_smiles():
    """Test that invalid SMILES raises error."""
    pytest.importorskip("rdkit")
    from kanad.io import from_smiles

    with pytest.raises(ValueError, match="Invalid SMILES"):
        from_smiles("invalid_smiles_xyz123")


def test_smiles_to_formula():
    """Test formula extraction."""
    pytest.importorskip("rdkit")
    from kanad.io.smiles_parser import smiles_to_formula

    assert smiles_to_formula("CCO") == "C2H6O"
    assert smiles_to_formula("C") == "CH4"
    assert smiles_to_formula("c1ccccc1") == "C6H6"


def test_validate_smiles():
    """Test SMILES validation."""
    pytest.importorskip("rdkit")
    from kanad.io.smiles_parser import validate_smiles

    # Valid SMILES
    valid, msg = validate_smiles("CCO")
    assert valid is True
    assert msg == ""

    # Invalid SMILES
    valid, msg = validate_smiles("invalid")
    assert valid is False
    assert "Invalid" in msg


def test_smiles_no_optimization():
    """Test SMILES parsing without geometry optimization."""
    pytest.importorskip("rdkit")
    from kanad.io import from_smiles

    mol = from_smiles("C", optimize_geometry=False)
    assert mol is not None
    assert len(mol.atoms) == 5


def test_smiles_metadata():
    """Test that SMILES metadata is stored."""
    pytest.importorskip("rdkit")
    from kanad.io import from_smiles

    mol = from_smiles("CCO", name="Ethanol")
    assert mol._smiles == "CCO"
    assert mol._name == "Ethanol"


def test_rdkit_not_installed():
    """Test error when RDKit not installed."""
    import sys
    import importlib

    # Mock RDKit not being installed
    rdkit_module = sys.modules.get('rdkit')
    if rdkit_module:
        # Can't properly test this if RDKit is installed
        pytest.skip("RDKit is installed, cannot test import error")
