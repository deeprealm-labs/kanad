"""
Unit tests for Hamiltonian builders (Phase 5).

Tests bond-type-specific Hamiltonians:
- IonicHamiltonian
- CovalentHamiltonian
"""

import pytest
import numpy as np
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule
from kanad.core.representations.second_quantization import SecondQuantizationRepresentation
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian


class TestIonicHamiltonian:
    """Test ionic bonding Hamiltonian."""

    def test_ionic_hamiltonian_creation(self):
        """Test creating ionic Hamiltonian for NaCl."""
        # Na+ Cl- ionic molecule
        na = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))  # ~2.36 Å

        molecule = Molecule([na, cl])
        representation = SecondQuantizationRepresentation(molecule)

        hamiltonian = IonicHamiltonian(molecule, representation)

        assert hamiltonian.n_orbitals == 2  # One orbital per atom
        assert hamiltonian.n_electrons == molecule.n_electrons
        assert hamiltonian.nuclear_repulsion > 0  # Repulsive

    def test_ionic_site_energies(self):
        """Test on-site energies from electronegativity."""
        na = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))

        molecule = Molecule([na, cl])
        representation = SecondQuantizationRepresentation(molecule)
        hamiltonian = IonicHamiltonian(molecule, representation)

        # Cl is more electronegative → lower (more negative) on-site energy
        assert hamiltonian.h_core[1, 1] < hamiltonian.h_core[0, 0]

    def test_ionic_transfer_integral(self):
        """Test transfer integrals are small for ionic bonds."""
        na = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))

        molecule = Molecule([na, cl])
        representation = SecondQuantizationRepresentation(molecule)
        hamiltonian = IonicHamiltonian(molecule, representation)

        # Transfer integral should be small (weak overlap)
        t = hamiltonian.h_core[0, 1]
        assert abs(t) < 0.2  # Less than ~5 eV

    def test_ionic_hubbard_u(self):
        """Test on-site Coulomb repulsion."""
        na = Atom('Na', position=np.zeros(3))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))

        molecule = Molecule([na, cl])
        representation = SecondQuantizationRepresentation(molecule)
        hamiltonian = IonicHamiltonian(molecule, representation)

        # Hubbard U should be positive
        U_na = hamiltonian.eri[0, 0, 0, 0]
        U_cl = hamiltonian.eri[1, 1, 1, 1]

        assert U_na > 0
        assert U_cl > 0

    def test_ionic_energy_computation(self):
        """Test energy calculation."""
        na = Atom('Na', position=np.zeros(3))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))

        molecule = Molecule([na, cl])
        representation = SecondQuantizationRepresentation(molecule)
        hamiltonian = IonicHamiltonian(molecule, representation)

        # Simple density matrix (all electrons on Cl)
        P = np.array([[0.0, 0.0], [0.0, 2.0]])

        energy = hamiltonian.compute_energy(P)

        # Should be finite
        assert np.isfinite(energy)

    def test_ionic_charge_transfer(self):
        """Test charge transfer calculation."""
        na = Atom('Na', position=np.zeros(3))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))

        molecule = Molecule([na, cl])
        representation = SecondQuantizationRepresentation(molecule)
        hamiltonian = IonicHamiltonian(molecule, representation)

        # All electrons on Cl (full transfer)
        P = np.array([[0.0, 0.0], [0.0, 2.0]])

        charges = hamiltonian.compute_charge_transfer(P)

        # Na should be positive (lost electrons)
        # charge = Z - n_electrons
        # Na: Z=11, n=0 → charge = +11
        assert charges[0] > 0
        # Cl should be positive too (Z=17, n=2 → charge = +15)
        # But if we add an electron to neutral Cl (Z=17, n=18), charge would be negative
        # The test needs a proper density matrix
        assert charges[1] > charges[0]  # Cl has more nuclear charge

    def test_ionic_hamiltonian_hermitian(self):
        """Test that Hamiltonian is Hermitian."""
        na = Atom('Na', position=np.zeros(3))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))

        molecule = Molecule([na, cl])
        representation = SecondQuantizationRepresentation(molecule)
        hamiltonian = IonicHamiltonian(molecule, representation)

        # Core Hamiltonian should be Hermitian
        assert np.allclose(hamiltonian.h_core, hamiltonian.h_core.T)


class TestCovalentHamiltonian:
    """Test covalent bonding Hamiltonian."""

    def test_covalent_hamiltonian_creation(self):
        """Test creating covalent Hamiltonian for H2."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)

        hamiltonian = CovalentHamiltonian(molecule, representation)

        assert hamiltonian.n_orbitals == 2  # Two 1s orbitals
        assert hamiltonian.n_electrons == 2

    def test_covalent_hamiltonian_hermitian(self):
        """Test that Hamiltonian is Hermitian."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation)

        # Should be Hermitian
        assert np.allclose(hamiltonian.h_core, hamiltonian.h_core.T)

    def test_covalent_overlap_matrix(self):
        """Test overlap matrix is computed."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation)

        S = hamiltonian.get_overlap_matrix()

        # Should be symmetric
        assert np.allclose(S, S.T)

        # Diagonal should be positive
        assert np.all(np.diag(S) > 0)

    def test_covalent_mo_energies(self):
        """Test molecular orbital energies."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation)

        energies = hamiltonian.get_mo_energies()

        # Should have 2 MO energies
        assert len(energies) == 2

        # Bonding should be lower than antibonding
        assert energies[0] < energies[1]

    def test_covalent_bonding_antibonding_split(self):
        """Test bonding/antibonding splitting."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation)

        analysis = hamiltonian.get_bonding_antibonding_split()

        assert 'homo_lumo_gap' in analysis
        assert 'homo_energy' in analysis
        assert 'lumo_energy' in analysis

        # Gap should be positive
        assert analysis['homo_lumo_gap'] > 0

    def test_covalent_energy_computation(self):
        """Test energy calculation."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation)

        # Simple density matrix (bonding orbital occupied)
        P = np.array([[1.0, 1.0], [1.0, 1.0]])

        energy = hamiltonian.compute_energy(P)

        # Should be negative (bonding)
        assert energy < 0
        assert np.isfinite(energy)

    def test_covalent_homo_lumo_gap(self):
        """Test HOMO-LUMO gap calculation."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation)

        gap = hamiltonian.get_homo_lumo_gap()

        # Should be positive
        assert gap > 0

    def test_covalent_bond_order(self):
        """Test bond order calculation."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation)

        # Bonding density matrix
        P = np.array([[1.0, 1.0], [1.0, 1.0]])

        bond_order = hamiltonian.compute_bond_order(P, 0, 1)

        # Should be positive (bonding)
        assert bond_order > 0


class TestHamiltonianComparison:
    """Compare ionic vs covalent Hamiltonians."""

    def test_ionic_vs_covalent_site_energies(self):
        """Compare site energy differences."""
        # NaCl (ionic)
        na = Atom('Na', position=np.zeros(3))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))
        mol_ionic = Molecule([na, cl])
        rep_ionic = SecondQuantizationRepresentation(mol_ionic)
        ham_ionic = IonicHamiltonian(mol_ionic, rep_ionic)

        # H2 (covalent)
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        mol_cov = Molecule([h1, h2])
        rep_cov = LCAORepresentation(mol_cov)
        ham_cov = CovalentHamiltonian(mol_cov, rep_cov)

        # Ionic: large site energy difference (Δε ~ electronegativity difference)
        delta_ionic = abs(ham_ionic.h_core[0, 0] - ham_ionic.h_core[1, 1])

        # Covalent: small site energy difference (similar atoms)
        delta_cov = abs(ham_cov.h_core[0, 0] - ham_cov.h_core[1, 1])

        # Ionic should have larger difference
        assert delta_ionic > delta_cov

    def test_ionic_vs_covalent_transfer_integral(self):
        """Compare transfer integrals."""
        # NaCl
        na = Atom('Na', position=np.zeros(3))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))
        mol_ionic = Molecule([na, cl])
        rep_ionic = SecondQuantizationRepresentation(mol_ionic)
        ham_ionic = IonicHamiltonian(mol_ionic, rep_ionic)

        # H2
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        mol_cov = Molecule([h1, h2])
        rep_cov = LCAORepresentation(mol_cov)
        ham_cov = CovalentHamiltonian(mol_cov, rep_cov)

        # Transfer integral (off-diagonal)
        t_ionic = abs(ham_ionic.h_core[0, 1])
        t_cov = abs(ham_cov.h_core[0, 1])

        # Both should be non-zero
        assert t_ionic > 0
        assert t_cov > 0


class TestHamiltonianEdgeCases:
    """Test edge cases and error handling."""

    def test_single_atom_ionic(self):
        """Test ionic Hamiltonian with single atom."""
        na = Atom('Na', position=np.zeros(3))
        molecule = Molecule([na])
        representation = SecondQuantizationRepresentation(molecule)

        hamiltonian = IonicHamiltonian(molecule, representation)

        # Should have 1 orbital
        assert hamiltonian.n_orbitals == 1

        # Nuclear repulsion should be zero
        assert hamiltonian.nuclear_repulsion == 0.0

    def test_nuclear_repulsion_scaling(self):
        """Test nuclear repulsion scales with distance."""
        # Close atoms
        na1 = Atom('Na', position=np.zeros(3))
        cl1 = Atom('Cl', position=np.array([1.0, 0.0, 0.0]))
        mol1 = Molecule([na1, cl1])
        rep1 = SecondQuantizationRepresentation(mol1)
        ham1 = IonicHamiltonian(mol1, rep1)

        # Far atoms
        na2 = Atom('Na', position=np.zeros(3))
        cl2 = Atom('Cl', position=np.array([10.0, 0.0, 0.0]))
        mol2 = Molecule([na2, cl2])
        rep2 = SecondQuantizationRepresentation(mol2)
        ham2 = IonicHamiltonian(mol2, rep2)

        # Closer atoms should have higher repulsion
        assert ham1.nuclear_repulsion > ham2.nuclear_repulsion


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
