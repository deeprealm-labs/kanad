"""
Test that all kanad modules can be imported successfully.
"""
import pytest


def test_main_import():
    """Test main kanad package import."""
    import kanad
    assert kanad.__version__ == "0.1.0"
    # Updated after cleanup: ExcitedStatesSolver + QPESolver/FCISolver placeholders = 32
    assert len(kanad.__all__) == 32


def test_molecule_import():
    """Test Molecule class import."""
    from kanad import Molecule
    assert Molecule is not None


def test_bonds_import():
    """Test bond classes import."""
    from kanad import BondFactory, IonicBond, CovalentBond, MetallicBond
    assert all([BondFactory, IonicBond, CovalentBond, MetallicBond])


def test_hamiltonians_import():
    """Test Hamiltonian classes import."""
    from kanad import (
        MolecularHamiltonian,
        IonicHamiltonian,
        CovalentHamiltonian,
        MetallicHamiltonian
    )
    assert all([MolecularHamiltonian, IonicHamiltonian,
                CovalentHamiltonian, MetallicHamiltonian])


def test_ansatze_import():
    """Test ansatz classes import."""
    from kanad import (
        BaseAnsatz,
        UCCAnsatz,
        HardwareEfficientAnsatz,
        IonicGovernanceAnsatz,
        CovalentGovernanceAnsatz,
        AdaptiveGovernanceAnsatz
    )
    assert all([BaseAnsatz, UCCAnsatz, HardwareEfficientAnsatz,
                IonicGovernanceAnsatz, CovalentGovernanceAnsatz,
                AdaptiveGovernanceAnsatz])


def test_solvers_import():
    """Test solver classes import."""
    from kanad import VQESolver, SQDSolver
    # QPESolver and FCISolver are legacy (not yet reimplemented)
    assert VQESolver is not None
    assert SQDSolver is not None


def test_backends_import():
    """Test backend classes import."""
    from kanad import BlueQubitBackend
    # QiskitBackend and IBMRuntimeBackend require optional dependencies
    # Only test BlueQubitBackend which is always available
    assert BlueQubitBackend is not None


def test_ibm_solvers_import():
    """Test IBM-specific solver classes import."""
    # IBM solvers were removed during framework cleanup
    # Tests now use main VQESolver with backend='ibm'
    from kanad import VQESolver
    assert VQESolver is not None


def test_governance_import():
    """Test governance protocol classes import."""
    from kanad import (
        BaseGovernanceProtocol,
        BondingType,
        GovernanceRule,
        IonicGovernanceProtocol,
        CovalentGovernanceProtocol,
        MetallicGovernanceProtocol
    )
    assert all([BaseGovernanceProtocol, BondingType, GovernanceRule,
                IonicGovernanceProtocol, CovalentGovernanceProtocol,
                MetallicGovernanceProtocol])
