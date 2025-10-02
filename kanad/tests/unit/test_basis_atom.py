"""
Unit tests for Atom and BasisSet classes.
"""

import pytest
import numpy as np
from kanad.core.atom import Atom
from kanad.core.integrals.basis_sets import (
    BasisSet,
    GaussianPrimitive,
    ContractedGaussian,
)


class TestAtom:
    """Test Atom class."""

    def test_atom_creation(self):
        """Test basic atom creation."""
        h = Atom('H')

        assert h.symbol == 'H'
        assert h.atomic_number == 1
        assert h.n_electrons == 1
        assert h.n_valence == 1
        assert h.is_metal() is False

    def test_atom_with_position(self):
        """Test atom with specified position."""
        pos = np.array([1.0, 2.0, 3.0])
        c = Atom('C', position=pos)

        assert np.allclose(c.position, pos)
        assert c.atomic_number == 6
        assert c.n_valence == 4

    def test_atom_distance(self):
        """Test distance calculation between atoms."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        distance = h1.distance_to(h2)
        assert abs(distance - 0.74) < 1e-10

    def test_charged_atom(self):
        """Test atom with formal charge."""
        na_plus = Atom('Na', charge=1)

        assert na_plus.atomic_number == 11
        assert na_plus.n_electrons == 10  # Lost 1 electron

    def test_metal_detection(self):
        """Test metal vs non-metal classification."""
        na = Atom('Na')
        cl = Atom('Cl')

        assert na.is_metal() is True
        assert cl.is_metal() is False


class TestGaussianPrimitive:
    """Test Gaussian primitive functions."""

    def test_primitive_creation(self):
        """Test creating a Gaussian primitive."""
        prim = GaussianPrimitive(
            exponent=1.0,
            coefficient=0.5,
            angular_momentum=(0, 0, 0),
            center=np.zeros(3)
        )

        assert prim.exponent == 1.0
        assert prim.coefficient == 0.5
        assert prim.l == 0  # s orbital

    def test_s_orbital_angular_momentum(self):
        """Test s orbital (l=0)."""
        s_prim = GaussianPrimitive(
            exponent=1.0,
            coefficient=1.0,
            angular_momentum=(0, 0, 0),
            center=np.zeros(3)
        )

        assert s_prim.l == 0

    def test_p_orbital_angular_momentum(self):
        """Test p orbital (l=1)."""
        px_prim = GaussianPrimitive(
            exponent=1.0,
            coefficient=1.0,
            angular_momentum=(1, 0, 0),
            center=np.zeros(3)
        )

        assert px_prim.l == 1

    def test_gaussian_evaluation_at_center(self):
        """Test Gaussian value at its center."""
        prim = GaussianPrimitive(
            exponent=1.0,
            coefficient=1.0,
            angular_momentum=(0, 0, 0),
            center=np.zeros(3)
        )

        value = prim.evaluate(np.zeros(3))
        assert value > 0  # Should be maximum at center

    def test_gaussian_decay(self):
        """Test Gaussian decays with distance."""
        prim = GaussianPrimitive(
            exponent=1.0,
            coefficient=1.0,
            angular_momentum=(0, 0, 0),
            center=np.zeros(3)
        )

        value_at_center = prim.evaluate(np.zeros(3))
        value_far = prim.evaluate(np.array([5.0, 0.0, 0.0]))

        assert value_at_center > value_far
        assert value_far > 0  # But not zero (Gaussian never truly zero)


class TestContractedGaussian:
    """Test contracted Gaussian functions."""

    def test_contracted_creation(self):
        """Test creating contracted Gaussian."""
        primitives = [
            GaussianPrimitive(1.0, 0.5, (0, 0, 0), np.zeros(3)),
            GaussianPrimitive(0.5, 0.5, (0, 0, 0), np.zeros(3)),
        ]

        contracted = ContractedGaussian(primitives, shell_type='s')

        assert len(contracted.primitives) == 2
        assert contracted.shell_type == 's'

    def test_contracted_evaluation(self):
        """Test evaluation of contracted function."""
        primitives = [
            GaussianPrimitive(1.0, 0.5, (0, 0, 0), np.zeros(3)),
            GaussianPrimitive(0.5, 0.5, (0, 0, 0), np.zeros(3)),
        ]

        contracted = ContractedGaussian(primitives, shell_type='s')
        value = contracted.evaluate(np.zeros(3))

        assert value > 0


class TestBasisSet:
    """Test basis set construction."""

    def test_basis_set_creation(self):
        """Test creating basis set."""
        basis = BasisSet('sto-3g')

        assert basis.basis_name == 'sto-3g'
        assert basis.n_basis_functions == 0  # Not built yet

    def test_h2_basis_construction(self):
        """Test building STO-3G basis for H2."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        # Each H atom gets 1 s orbital
        assert basis.n_basis_functions == 2

    def test_ch4_basis_construction(self):
        """Test building STO-3G basis for CH4."""
        c = Atom('C', position=np.zeros(3))

        # Tetrahedral hydrogen positions
        h1 = Atom('H', position=np.array([0.63, 0.63, 0.63]))
        h2 = Atom('H', position=np.array([-0.63, -0.63, 0.63]))
        h3 = Atom('H', position=np.array([-0.63, 0.63, -0.63]))
        h4 = Atom('H', position=np.array([0.63, -0.63, -0.63]))

        basis = BasisSet('sto-3g')
        basis.build_basis([c, h1, h2, h3, h4])

        # C: 2 s + 3 p = 5 orbitals
        # 4 H: 4 s = 4 orbitals
        # Total: 9 orbitals
        assert basis.n_basis_functions == 9

    def test_water_basis_construction(self):
        """Test building STO-3G basis for H2O."""
        o = Atom('O', position=np.zeros(3))
        h1 = Atom('H', position=np.array([0.757, 0.586, 0.0]))
        h2 = Atom('H', position=np.array([-0.757, 0.586, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([o, h1, h2])

        # O: 2 s + 3 p = 5 orbitals
        # 2 H: 2 s = 2 orbitals
        # Total: 7 orbitals
        assert basis.n_basis_functions == 7

    def test_unknown_element_raises_error(self):
        """Test that unknown elements raise error."""
        fake_atom = Atom('H')  # Valid atom
        fake_atom.symbol = 'Xx'  # Change to invalid symbol

        basis = BasisSet('sto-3g')

        with pytest.raises(ValueError):
            basis.build_basis([fake_atom])

    def test_unsupported_basis_raises_error(self):
        """Test that unsupported basis sets raise error."""
        h = Atom('H')

        with pytest.raises(NotImplementedError):
            basis = BasisSet('cc-pvdz')  # Not implemented yet
            basis.build_basis([h])


class TestSTO3GAccuracy:
    """Test STO-3G basis accuracy against known values."""

    def test_sto3g_hydrogen_exponents(self):
        """Test STO-3G hydrogen exponents match reference."""
        expected_exponents = [3.42525091, 0.62391373, 0.16885540]

        sto3g_data = BasisSet.STO3G_DATA['H']['s']
        exponents = [exp for exp, coeff in sto3g_data]

        for i, (exp, exp_ref) in enumerate(zip(exponents, expected_exponents)):
            assert abs(exp - exp_ref) < 1e-6, f"Exponent {i} mismatch"

    def test_sto3g_carbon_core_exponents(self):
        """Test STO-3G carbon core s exponents."""
        expected = [71.6168370, 13.0450960, 3.53051220]

        sto3g_data = BasisSet.STO3G_DATA['C']['s']
        exponents = [exp for exp, coeff in sto3g_data]

        for i, (exp, exp_ref) in enumerate(zip(exponents, expected)):
            assert abs(exp - exp_ref) < 1e-5


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
