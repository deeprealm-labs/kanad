"""
Unit tests for representation layer (Phase 3).

Tests SecondQuantization, LCAO, and base representation classes.
"""

import pytest
import numpy as np
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule
from kanad.core.representations.second_quantization import SecondQuantizationRepresentation
from kanad.core.representations.lcao_representation import (
    LCAORepresentation,
    HybridizationType
)


class TestMolecule:
    """Test Molecule class."""

    def test_molecule_creation(self):
        """Test creating a molecule."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        mol = Molecule([h1, h2])

        assert mol.n_atoms == 2
        assert mol.n_electrons == 2
        assert mol.symbols == ['H', 'H']

    def test_distance_matrix(self):
        """Test distance matrix computation."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([1.0, 0.0, 0.0]))

        mol = Molecule([h1, h2])
        D = mol.distance_matrix()

        assert D.shape == (2, 2)
        assert D[0, 0] == 0.0
        assert abs(D[0, 1] - 1.0) < 1e-10
        assert abs(D[1, 0] - 1.0) < 1e-10

    def test_water_molecule(self):
        """Test H2O molecule."""
        o = Atom('O', position=np.zeros(3))
        h1 = Atom('H', position=np.array([0.757, 0.586, 0.0]))
        h2 = Atom('H', position=np.array([-0.757, 0.586, 0.0]))

        mol = Molecule([o, h1, h2])

        assert mol.n_atoms == 3
        assert mol.n_electrons == 10  # 8 + 1 + 1
        assert 'O' in mol.symbols
        assert mol.symbols.count('H') == 2


class TestSecondQuantizationRepresentation:
    """Test SecondQuantization representation for ionic bonding."""

    def test_h2_second_quantization(self):
        """Test H2 in second quantization."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        mol = Molecule([h1, h2])

        rep = SecondQuantizationRepresentation(mol)

        assert rep.n_orbitals == 2  # One per atom
        assert rep.n_spin_orbitals == 4  # With spin
        assert rep.n_qubits == 4
        assert rep.n_electrons == 2

    def test_nacl_representation(self):
        """Test NaCl ionic system."""
        na = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))
        mol = Molecule([na, cl])

        rep = SecondQuantizationRepresentation(mol)

        assert rep.n_orbitals == 2
        assert rep.n_electrons == 11 + 17  # Na + Cl

    def test_reference_state(self):
        """Test reference state generation."""
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        mol = Molecule([h1, h2])

        rep = SecondQuantizationRepresentation(mol)
        ref_state = rep.get_reference_state()

        # Should be normalized
        assert abs(np.linalg.norm(ref_state) - 1.0) < 1e-10

        # Should be a valid state vector
        assert ref_state.shape == (2 ** rep.n_qubits,)

    def test_on_site_energies(self):
        """Test on-site energy calculation."""
        na = Atom('Na', position=np.zeros(3))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))
        mol = Molecule([na, cl])

        rep = SecondQuantizationRepresentation(mol)
        energies = rep.get_on_site_energies()

        assert len(energies) == 2
        # Cl should have lower energy (more electronegative)
        assert energies[1] < energies[0]

    def test_hopping_matrix(self):
        """Test hopping matrix generation."""
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([1.0, 0.0, 0.0]))
        mol = Molecule([h1, h2])

        rep = SecondQuantizationRepresentation(mol)
        t = rep.get_hopping_matrix()

        assert t.shape == (2, 2)
        # Should be symmetric
        assert abs(t[0, 1] - t[1, 0]) < 1e-10
        # Diagonal should be zero
        assert abs(t[0, 0]) < 1e-10
        # Off-diagonal should be positive
        assert t[0, 1] > 0


class TestLCAORepresentation:
    """Test LCAO representation for covalent bonding."""

    def test_h2_lcao(self):
        """Test H2 in LCAO representation."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        mol = Molecule([h1, h2])

        rep = LCAORepresentation(mol)

        assert rep.n_electrons == 2
        # Each H has 1 orbital (1s)
        assert rep.n_orbitals == 2
        assert rep.n_qubits > 0

    def test_hybridization_detection(self):
        """Test automatic hybridization detection."""
        c = Atom('C', position=np.zeros(3))
        h = Atom('H', position=np.array([1.0, 0.0, 0.0]))
        mol = Molecule([c, h])

        rep = LCAORepresentation(mol)

        # Carbon should be sp³
        assert rep.hybridization[0] == HybridizationType.SP3
        # Hydrogen should have no hybridization
        assert rep.hybridization[1] == HybridizationType.NONE

    def test_sp3_hybrid_orbitals(self):
        """Test sp³ hybrid orbital construction."""
        c = Atom('C', position=np.zeros(3))
        mol = Molecule([c])

        rep = LCAORepresentation(mol)

        # sp³ carbon should have 4 hybrid orbitals
        c_hybrids = [h for h in rep.hybrid_orbitals if h['atom_index'] == 0]
        assert len(c_hybrids) == 4

        # All should be sp³ type
        for hybrid in c_hybrids:
            assert hybrid['type'] == 'sp3'

    def test_water_lcao(self):
        """Test H2O in LCAO representation."""
        o = Atom('O', position=np.zeros(3))
        h1 = Atom('H', position=np.array([0.757, 0.586, 0.0]))
        h2 = Atom('H', position=np.array([-0.757, 0.586, 0.0]))
        mol = Molecule([o, h1, h2])

        rep = LCAORepresentation(mol)

        assert rep.n_electrons == 10
        # O has 4 sp³ orbitals, each H has 1
        assert rep.n_orbitals == 6  # 4 + 1 + 1

    def test_reference_state_lcao(self):
        """Test LCAO reference state."""
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        mol = Molecule([h1, h2])

        rep = LCAORepresentation(mol)
        ref_state = rep.get_reference_state()

        # Should be normalized
        assert abs(np.linalg.norm(ref_state) - 1.0) < 1e-10

    def test_mo_pairs(self):
        """Test molecular orbital pair construction."""
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        mol = Molecule([h1, h2])

        rep = LCAORepresentation(mol)

        # H2 should have one MO pair (bonding + antibonding)
        assert len(rep.mo_pairs) == 1


class TestRepresentationComparison:
    """Test comparisons between different representations."""

    def test_h2_qubit_counts(self):
        """Compare qubit requirements for H2."""
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        mol = Molecule([h1, h2])

        sq_rep = SecondQuantizationRepresentation(mol)
        lcao_rep = LCAORepresentation(mol)

        # Both should need same number of qubits for H2
        assert sq_rep.get_num_qubits() == lcao_rep.get_num_qubits()

    def test_observable_computation(self):
        """Test observable computation for different representations."""
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        mol = Molecule([h1, h2])

        sq_rep = SecondQuantizationRepresentation(mol)
        ref_state_sq = sq_rep.get_reference_state()

        # Should compute without error
        obs_sq = sq_rep.compute_observables(ref_state_sq)
        assert 'site_occupations' in obs_sq


class TestHamiltonianConstruction:
    """Test Hamiltonian construction (placeholders)."""

    def test_ionic_hamiltonian_creation(self):
        """Test creating ionic Hamiltonian."""
        na = Atom('Na', position=np.zeros(3))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))
        mol = Molecule([na, cl])

        rep = SecondQuantizationRepresentation(mol)
        ham = rep.build_hamiltonian()

        assert ham is not None
        assert hasattr(ham, 'molecule')
        assert hasattr(ham, 'representation')

    def test_covalent_hamiltonian_creation(self):
        """Test creating covalent Hamiltonian."""
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        mol = Molecule([h1, h2])

        rep = LCAORepresentation(mol)
        ham = rep.build_hamiltonian()

        assert ham is not None
        assert hasattr(ham, 'molecule')


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
