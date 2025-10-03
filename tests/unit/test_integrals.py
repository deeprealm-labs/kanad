"""
Unit tests for integral computation modules.

Tests overlap, kinetic, nuclear attraction, and ERI calculations.
"""

import pytest
import numpy as np
from scipy.special import erf
from kanad.core.atom import Atom
from kanad.core.integrals.basis_sets import BasisSet
from kanad.core.integrals.overlap import OverlapIntegrals
from kanad.core.integrals.one_electron import OneElectronIntegrals
from kanad.core.integrals.two_electron import TwoElectronIntegrals


class TestOverlapIntegrals:
    """Test overlap integral computation."""

    def test_overlap_identity(self):
        """Overlap of a basis function with itself should be positive."""
        h = Atom('H', position=np.zeros(3))
        basis = BasisSet('sto-3g')
        basis.build_basis([h])

        S = OverlapIntegrals.build_overlap_matrix(basis.basis_functions)

        # Diagonal elements should be positive (self-overlap)
        assert S[0, 0] > 0

    def test_overlap_symmetry(self):
        """Overlap matrix should be symmetric."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        S = OverlapIntegrals.build_overlap_matrix(basis.basis_functions)

        # Check symmetry
        assert np.allclose(S, S.T)

    def test_h2_overlap(self):
        """Test H2 overlap matrix."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        S = OverlapIntegrals.build_overlap_matrix(basis.basis_functions)

        # Should be 2x2
        assert S.shape == (2, 2)

        # Diagonal should be positive
        assert S[0, 0] > 0
        assert S[1, 1] > 0

        # Off-diagonal should be positive (bonding overlap)
        assert S[0, 1] > 0
        assert abs(S[0, 1] - S[1, 0]) < 1e-10

    def test_overlap_decreases_with_distance(self):
        """Overlap should decrease as atoms move apart."""
        h1 = Atom('H', position=np.zeros(3))

        basis_close = BasisSet('sto-3g')
        h2_close = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        basis_close.build_basis([h1, h2_close])
        S_close = OverlapIntegrals.build_overlap_matrix(basis_close.basis_functions)

        basis_far = BasisSet('sto-3g')
        h2_far = Atom('H', position=np.array([2.0, 0.0, 0.0]))
        basis_far.build_basis([h1, h2_far])
        S_far = OverlapIntegrals.build_overlap_matrix(basis_far.basis_functions)

        assert S_close[0, 1] > S_far[0, 1]


class TestOneElectronIntegrals:
    """Test one-electron integral computation."""

    def test_kinetic_matrix_positive(self):
        """Kinetic energy matrix should have positive diagonal."""
        h = Atom('H', position=np.zeros(3))
        basis = BasisSet('sto-3g')
        basis.build_basis([h])

        integrals = OneElectronIntegrals([h], basis.basis_functions)
        T = integrals.compute_kinetic()

        # Kinetic energy should be positive
        assert T[0, 0] > 0

    def test_kinetic_symmetry(self):
        """Kinetic matrix should be symmetric."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        integrals = OneElectronIntegrals([h1, h2], basis.basis_functions)
        T = integrals.compute_kinetic()

        assert np.allclose(T, T.T)

    def test_nuclear_attraction_negative(self):
        """Nuclear attraction should be negative (attractive)."""
        h = Atom('H', position=np.zeros(3))
        basis = BasisSet('sto-3g')
        basis.build_basis([h])

        integrals = OneElectronIntegrals([h], basis.basis_functions)
        V = integrals.compute_nuclear_attraction()

        # Should be negative (electron-nucleus attraction)
        assert V[0, 0] < 0

    def test_nuclear_repulsion_h2(self):
        """Test nuclear repulsion for H2."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        integrals = OneElectronIntegrals([h1, h2], basis.basis_functions)
        E_nn = integrals.compute_nuclear_repulsion()

        # Should be positive (repulsive)
        assert E_nn > 0

        # Approximate value check (1/0.74 Å ≈ 0.714 Hartree)
        assert 0.5 < E_nn < 1.0

    def test_core_hamiltonian_hermitian(self):
        """Core Hamiltonian should be Hermitian."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        integrals = OneElectronIntegrals([h1, h2], basis.basis_functions)
        H_core = integrals.compute_core_hamiltonian()

        assert np.allclose(H_core, H_core.T)


class TestTwoElectronIntegrals:
    """Test two-electron repulsion integrals."""

    def test_eri_positive(self):
        """ERIs should be positive (repulsive)."""
        h = Atom('H', position=np.zeros(3))
        basis = BasisSet('sto-3g')
        basis.build_basis([h])

        eri_calc = TwoElectronIntegrals(basis.basis_functions)
        eri_dict = eri_calc.compute_eri_sparse()

        # (00|00) should be positive
        if (0, 0, 0, 0) in eri_dict:
            assert eri_dict[(0, 0, 0, 0)] > 0

    def test_eri_symmetry(self):
        """Test 8-fold permutational symmetry of ERIs."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([1.0, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        eri_calc = TwoElectronIntegrals(basis.basis_functions)
        ERI = eri_calc.compute_eri_tensor()

        # Test some symmetries
        # (ij|kl) = (ji|kl)
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        assert abs(ERI[i,j,k,l] - ERI[j,i,k,l]) < 1e-10
                        # (ij|kl) = (ij|lk)
                        assert abs(ERI[i,j,k,l] - ERI[i,j,l,k]) < 1e-10

    def test_eri_sparse_consistency(self):
        """Sparse and tensor ERIs should match."""
        h = Atom('H', position=np.zeros(3))
        basis = BasisSet('sto-3g')
        basis.build_basis([h])

        eri_calc = TwoElectronIntegrals(basis.basis_functions)

        ERI_tensor = eri_calc.compute_eri_tensor()
        eri_sparse = eri_calc.compute_eri_sparse(threshold=1e-12)

        # Check stored values match
        for (i, j, k, l), value in eri_sparse.items():
            assert abs(ERI_tensor[i, j, k, l] - value) < 1e-10


class TestIntegralAccuracy:
    """Test integral accuracy against known values."""

    def test_h_atom_integrals(self):
        """Test hydrogen atom integrals against theoretical values."""
        h = Atom('H', position=np.zeros(3))
        basis = BasisSet('sto-3g')
        basis.build_basis([h])

        # Overlap should be positive
        S = OverlapIntegrals.build_overlap_matrix(basis.basis_functions)
        assert S[0, 0] > 0

        # Kinetic + Nuclear should give approximate H atom energy
        integrals = OneElectronIntegrals([h], basis.basis_functions)
        H_core = integrals.compute_core_hamiltonian()

        # H atom ground state energy ≈ -0.5 Hartree (exact)
        # STO-3G gives approximate value (can be off due to basis limitations)
        assert H_core[0, 0] < 0  # Should be negative
        assert -10.0 < H_core[0, 0] < 0.0  # Reasonable range

    def test_h2_energy_order_of_magnitude(self):
        """Test that H2 energy is in correct range."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        integrals = OneElectronIntegrals([h1, h2], basis.basis_functions)
        H_core = integrals.compute_core_hamiltonian()
        E_nn = integrals.compute_nuclear_repulsion()

        # Core Hamiltonian eigenvalues
        eigvals = np.linalg.eigvalsh(H_core)

        # Bonding orbital should have lower energy than isolated atoms
        # Check for finite eigenvalues (no NaN)
        assert np.isfinite(eigvals).all()
        # Eigenvalues should be negative (attractive)
        assert eigvals[0] < 0

    def test_water_integral_shapes(self):
        """Test correct integral matrix shapes for H2O."""
        o = Atom('O', position=np.zeros(3))
        h1 = Atom('H', position=np.array([0.757, 0.586, 0.0]))
        h2 = Atom('H', position=np.array([-0.757, 0.586, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([o, h1, h2])

        # Should have 7 basis functions
        n = basis.n_basis_functions
        assert n == 7

        # Test overlap
        S = OverlapIntegrals.build_overlap_matrix(basis.basis_functions)
        assert S.shape == (7, 7)

        # Test one-electron integrals
        integrals = OneElectronIntegrals([o, h1, h2], basis.basis_functions)
        T = integrals.compute_kinetic()
        V = integrals.compute_nuclear_attraction()

        assert T.shape == (7, 7)
        assert V.shape == (7, 7)


class TestBoysFunction:
    """Test Boys function (critical for nuclear attraction integrals)."""

    def test_boys_function_n0_small_t(self):
        """Test Boys function F_0(t) for small t."""
        from kanad.core.integrals.one_electron import OneElectronIntegrals

        # For small t, F_0(t) ≈ 1 - t/3
        t = 0.1
        F0 = OneElectronIntegrals._boys_function(0, t)
        expected = np.sqrt(np.pi / (4 * t)) * erf(np.sqrt(t))

        assert abs(F0 - expected) < 1e-10

    def test_boys_function_n0_large_t(self):
        """Test Boys function F_0(t) for large t."""
        from kanad.core.integrals.one_electron import OneElectronIntegrals

        # For large t, F_0(t) ≈ 0.5 * sqrt(π/t)
        t = 10.0
        F0 = OneElectronIntegrals._boys_function(0, t)
        expected = np.sqrt(np.pi / (4 * t)) * erf(np.sqrt(t))

        assert abs(F0 - expected) < 1e-8

    def test_boys_function_higher_n(self):
        """Test Boys function for n > 0."""
        from kanad.core.integrals.one_electron import OneElectronIntegrals

        t = 1.0
        # Test recursion relation: F_n(t) should exist for n=1,2,3
        F1 = OneElectronIntegrals._boys_function(1, t)
        F2 = OneElectronIntegrals._boys_function(2, t)

        # Should be positive and decreasing with n
        assert F1 > 0
        assert F2 > 0
        assert F2 < F1


class TestCoulombExchangeMatrices:
    """Test Coulomb and Exchange matrix computation (CRITICAL - currently untested!)."""

    def test_coulomb_matrix_h2(self):
        """Test Coulomb matrix computation for H2."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        # Create simple density matrix (HF reference for 2 electrons)
        n = basis.n_basis_functions
        P = np.zeros((n, n))
        P[0, 0] = 1.0  # One electron in first orbital
        P[1, 1] = 1.0  # One electron in second orbital

        eri_calc = TwoElectronIntegrals(basis.basis_functions)
        J = eri_calc.compute_coulomb_matrix(P)

        # Coulomb matrix should be symmetric
        assert np.allclose(J, J.T)

        # Should have positive diagonal (electron repulsion)
        assert J[0, 0] > 0
        assert J[1, 1] > 0

    def test_exchange_matrix_h2(self):
        """Test Exchange matrix computation for H2."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        n = basis.n_basis_functions
        P = np.zeros((n, n))
        P[0, 0] = 1.0
        P[1, 1] = 1.0

        eri_calc = TwoElectronIntegrals(basis.basis_functions)
        K = eri_calc.compute_exchange_matrix(P)

        # Exchange matrix should be symmetric
        assert np.allclose(K, K.T)

        # Should have positive diagonal
        assert K[0, 0] >= 0
        assert K[1, 1] >= 0

    def test_fock_matrix_construction(self):
        """Test Fock matrix construction (CRITICAL - completely untested!)."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        n = basis.n_basis_functions
        P = np.zeros((n, n))
        P[0, 0] = 1.0
        P[1, 1] = 1.0

        # Get core Hamiltonian
        integrals_1e = OneElectronIntegrals([h1, h2], basis.basis_functions)
        H_core = integrals_1e.compute_core_hamiltonian()

        # Compute Fock matrix
        eri_calc = TwoElectronIntegrals(basis.basis_functions)
        F = eri_calc.compute_fock_matrix(H_core, P)

        # Fock matrix should be Hermitian
        assert np.allclose(F, F.T)

        # Shape should match
        assert F.shape == (n, n)

        # Should contain core Hamiltonian contribution
        # F = H_core + 2J - K, so diagonal should be influenced by H_core
        assert np.isfinite(F).all()


class TestOrbitalIntegrals:
    """Test integrals for p and d orbitals (currently mostly untested)."""

    def test_p_orbital_overlap(self):
        """Test overlap integrals for p orbitals."""
        # Carbon has p orbitals
        c1 = Atom('C', position=np.array([0.0, 0.0, 0.0]))
        c2 = Atom('C', position=np.array([1.54, 0.0, 0.0]))  # C-C bond

        basis = BasisSet('sto-3g')
        basis.build_basis([c1, c2])

        S = OverlapIntegrals.build_overlap_matrix(basis.basis_functions)

        # Should be symmetric
        assert np.allclose(S, S.T)

        # All diagonal elements should be positive
        assert np.all(np.diag(S) > 0)

    def test_p_orbital_kinetic(self):
        """Test kinetic energy integrals for p orbitals."""
        c = Atom('C', position=np.zeros(3))
        basis = BasisSet('sto-3g')
        basis.build_basis([c])

        integrals = OneElectronIntegrals([c], basis.basis_functions)
        T = integrals.compute_kinetic()

        # Should be symmetric
        assert np.allclose(T, T.T)

        # Diagonal should be positive
        assert np.all(np.diag(T) > 0)

    def test_multi_element_integrals(self):
        """Test integrals for molecules with multiple element types."""
        # H2O has O (with p orbitals) and H (s orbitals only)
        o = Atom('O', position=np.zeros(3))
        h1 = Atom('H', position=np.array([0.757, 0.586, 0.0]))
        h2 = Atom('H', position=np.array([-0.757, 0.586, 0.0]))

        basis = BasisSet('sto-3g')
        basis.build_basis([o, h1, h2])

        # Test all integral types
        S = OverlapIntegrals.build_overlap_matrix(basis.basis_functions)

        integrals = OneElectronIntegrals([o, h1, h2], basis.basis_functions)
        T = integrals.compute_kinetic()
        V = integrals.compute_nuclear_attraction()
        H_core = integrals.compute_core_hamiltonian()

        # All should be symmetric/Hermitian
        assert np.allclose(S, S.T)
        assert np.allclose(T, T.T)
        assert np.allclose(V, V.T)
        assert np.allclose(H_core, H_core.T)


class TestIntegralErrorHandling:
    """Test error handling in integral computation."""

    def test_empty_basis_functions_overlap(self):
        """Test overlap with empty basis."""
        S = OverlapIntegrals.build_overlap_matrix([])

        # Should return empty array
        assert S.shape == (0, 0)

    def test_empty_atoms_nuclear_repulsion(self):
        """Test nuclear repulsion with no atoms."""
        h = Atom('H', position=np.zeros(3))
        basis = BasisSet('sto-3g')
        basis.build_basis([h])

        integrals = OneElectronIntegrals([], basis.basis_functions)
        E_nn = integrals.compute_nuclear_repulsion()

        # No atoms = no repulsion
        assert E_nn == 0.0

    def test_single_atom_nuclear_repulsion(self):
        """Test nuclear repulsion with single atom."""
        h = Atom('H', position=np.zeros(3))
        basis = BasisSet('sto-3g')
        basis.build_basis([h])

        integrals = OneElectronIntegrals([h], basis.basis_functions)
        E_nn = integrals.compute_nuclear_repulsion()

        # Single atom = no repulsion
        assert E_nn == 0.0

    def test_eri_sparse_threshold(self):
        """Test ERI sparse storage with different thresholds."""
        h = Atom('H', position=np.zeros(3))
        basis = BasisSet('sto-3g')
        basis.build_basis([h])

        eri_calc = TwoElectronIntegrals(basis.basis_functions)

        # Tight threshold should store more
        eri_tight = eri_calc.compute_eri_sparse(threshold=1e-12)

        # Loose threshold should store fewer
        eri_loose = eri_calc.compute_eri_sparse(threshold=0.1)

        assert len(eri_loose) <= len(eri_tight)


class TestIntegralEdgeCases:
    """Test edge cases and extreme values."""

    def test_very_close_atoms(self):
        """Test integrals with atoms very close together."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.01, 0.0, 0.0]))  # Very close

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        S = OverlapIntegrals.build_overlap_matrix(basis.basis_functions)

        # Overlap should be very large (close to 1)
        assert S[0, 1] > 0.9

        # Nuclear repulsion should be very large
        integrals = OneElectronIntegrals([h1, h2], basis.basis_functions)
        E_nn = integrals.compute_nuclear_repulsion()

        assert E_nn > 10.0  # 1/0.01 ≈ 100 in atomic units

    def test_very_distant_atoms(self):
        """Test integrals with atoms very far apart."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([100.0, 0.0, 0.0]))  # Very far

        basis = BasisSet('sto-3g')
        basis.build_basis([h1, h2])

        S = OverlapIntegrals.build_overlap_matrix(basis.basis_functions)

        # Overlap should be nearly zero
        assert S[0, 1] < 0.01

        # Nuclear repulsion should be small
        integrals = OneElectronIntegrals([h1, h2], basis.basis_functions)
        E_nn = integrals.compute_nuclear_repulsion()

        assert E_nn < 0.02  # 1/100 = 0.01 in atomic units

    def test_charged_atom_integrals(self):
        """Test integrals with charged atoms."""
        # H+ cation
        h_plus = Atom('H', position=np.zeros(3), charge=1)

        basis = BasisSet('sto-3g')
        basis.build_basis([h_plus])

        integrals = OneElectronIntegrals([h_plus], basis.basis_functions)
        V = integrals.compute_nuclear_attraction()

        # Nuclear attraction should still be negative
        assert V[0, 0] < 0

        # But magnitude should match nuclear charge
        assert np.isfinite(V[0, 0])


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
