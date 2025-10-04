#!/usr/bin/env python3
"""
Core functionality tests - Atoms, Bonds, Hamiltonians
"""

import pytest
import numpy as np
import sys
sys.path.insert(0, '/home/mk/deeprealm/kanad')

from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.bonds.ionic_bond import IonicBond
from kanad.bonds.metallic_bond import MetallicBond


class TestAtom:
    """Test Atom class."""

    def test_atom_creation(self):
        """Test creating atoms."""
        h = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        assert h.symbol == 'H'
        assert h.atomic_number == 1
        assert np.allclose(h.position, [0.0, 0.0, 0.0])

    def test_atom_properties(self):
        """Test atomic properties."""
        c = Atom('C')
        assert c.atomic_number == 6
        assert c.n_electrons == 6
        assert c.atomic_mass > 12.0

    def test_distance_calculation(self):
        """Test distance between atoms."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        distance = h1.distance_to(h2)
        assert np.isclose(distance, 0.74)


class TestCovalentBond:
    """Test Covalent Bond class."""

    def test_h2_bond(self):
        """Test H2 molecule."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        bond = CovalentBond(h1, h2)

        assert bond.hamiltonian is not None
        assert bond.hamiltonian.n_electrons == 2
        assert bond.hamiltonian.n_orbitals >= 2

    def test_covalent_energy(self):
        """Test energy calculation."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        bond = CovalentBond(h1, h2)
        result = bond.compute_energy(method='HF')

        assert 'energy' in result
        assert isinstance(result['energy'], (int, float))
        assert result['energy'] < 0  # Bound system


class TestIonicBond:
    """Test Ionic Bond class."""

    def test_lih_bond(self):
        """Test LiH molecule."""
        li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
        h = Atom('H', position=np.array([1.59, 0.0, 0.0]))

        bond = IonicBond(li, h)

        assert bond.hamiltonian is not None
        assert bond.hamiltonian.n_electrons == 4

    def test_ionic_energy(self):
        """Test ionic energy calculation."""
        li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
        h = Atom('H', position=np.array([1.59, 0.0, 0.0]))

        bond = IonicBond(li, h)
        result = bond.compute_energy(method='HF')

        assert 'energy' in result
        assert isinstance(result['energy'], (int, float))


class TestMetallicBond:
    """Test Metallic Bond class."""

    def test_na2_bond(self):
        """Test Na2 metallic bond."""
        na1 = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
        na2 = Atom('Na', position=np.array([3.66, 0.0, 0.0]))

        bond = MetallicBond([na1, na2])

        assert bond.hamiltonian is not None

    def test_metallic_energy(self):
        """Test metallic energy calculation."""
        na1 = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
        na2 = Atom('Na', position=np.array([3.66, 0.0, 0.0]))

        bond = MetallicBond([na1, na2])
        result = bond.compute_energy(method='tight_binding')

        assert 'energy' in result
        assert isinstance(result['energy'], (int, float))


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
