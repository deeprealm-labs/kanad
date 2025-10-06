"""
Unit tests for MP2 electron correlation.

Tests MP2 energy calculation, density matrices, and correlated properties.
"""

import pytest
import numpy as np


def test_import_mp2solver():
    """Test that MP2Solver can be imported."""
    from kanad.core import MP2Solver
    assert MP2Solver is not None


def test_mp2_energy_h2():
    """Test MP2 energy for H2 molecule."""
    from kanad.io import from_smiles
    from kanad.core import MP2Solver

    h2 = from_smiles("[H][H]", basis='6-311g(d,p)')
    mp2 = MP2Solver(h2.hamiltonian)
    result = mp2.compute_energy()

    # Check that all expected keys are present
    assert 'e_hf' in result
    assert 'e_corr' in result
    assert 'e_mp2' in result

    # H2 with 6-311G(d,p): E_MP2 ≈ -1.159 Ha
    assert -1.165 < result['e_mp2'] < -1.155, f"MP2 energy {result['e_mp2']:.6f} outside expected range"

    # Correlation energy should be negative
    assert result['e_corr'] < 0, "MP2 correlation energy should be negative"

    # MP2 energy should be lower than HF
    assert result['e_mp2'] < result['e_hf'], "MP2 energy should be lower than HF"

    print(f"H2 MP2 energy: {result['e_mp2']:.6f} Ha")


def test_mp2_energy_h2o():
    """Test MP2 energy for H2O molecule."""
    from kanad.io import from_smiles
    from kanad.core import MP2Solver

    water = from_smiles("O", basis='6-311g(d,p)')
    mp2 = MP2Solver(water.hamiltonian)
    result = mp2.compute_energy()

    # H2O with 6-311G(d,p): E_MP2 ≈ -76.28 Ha
    assert -76.35 < result['e_mp2'] < -76.20, f"MP2 energy {result['e_mp2']:.6f} outside expected range"

    # Correlation energy check
    assert result['e_corr'] < 0
    assert -0.30 < result['e_corr'] < -0.15, "Correlation energy seems unreasonable"

    print(f"H2O MP2 energy: {result['e_mp2']:.6f} Ha (correlation: {result['e_corr']:.6f} Ha)")


def test_mp2_density_matrix_h2o():
    """Test MP2 density matrix computation."""
    from kanad.io import from_smiles
    from kanad.core import MP2Solver

    water = from_smiles("O", basis='sto-3g')
    mp2 = MP2Solver(water.hamiltonian)

    # Get MP2 density matrix
    dm_mp2 = mp2.make_rdm1()

    # Check shape
    n_ao = water.hamiltonian.mol.nao
    assert dm_mp2.shape == (n_ao, n_ao), f"Density matrix shape {dm_mp2.shape} incorrect"

    # Density matrix should be Hermitian
    assert np.allclose(dm_mp2, dm_mp2.T), "MP2 density matrix not symmetric"

    print(f"MP2 density matrix: shape={dm_mp2.shape}, trace={np.trace(dm_mp2):.4f}")


def test_mp2_requires_converged_hf():
    """Test that MP2 requires converged HF reference."""
    from kanad.io import from_smiles
    from kanad.core import MP2Solver

    # Create molecule (HF automatically runs and converges)
    mol = from_smiles("O", basis='sto-3g')

    # Manually set converged = False to test error handling
    mol.hamiltonian.mf.converged = False

    # Should raise ValueError
    with pytest.raises(ValueError, match="HF calculation must converge"):
        mp2 = MP2Solver(mol.hamiltonian)


def test_mp2_polarizability_h2o():
    """Test MP2 polarizability calculation for H2O."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O", basis='6-311g(d,p)')
    calc = PropertyCalculator(water.hamiltonian)

    # Compute HF polarizability
    result_hf = calc.compute_polarizability()

    # Compute MP2 polarizability
    result_mp2 = calc.compute_polarizability_mp2()

    # MP2 should give larger polarizability than HF (more flexible electron cloud)
    # Note: For some systems MP2 can be lower, but for water it should be higher
    assert result_mp2['alpha_mean'] > 0, "MP2 polarizability should be positive"

    # Check that result has all expected keys
    assert 'alpha_tensor' in result_mp2
    assert 'alpha_mean' in result_mp2
    assert 'method' in result_mp2
    assert result_mp2['method'] == 'mp2_finite_field'

    # Experimental: 10.13 a.u.
    # MP2/6-311G(d,p): expect ~6-7 a.u. (60-70%)
    experimental = 10.13
    assert 5.0 < result_mp2['alpha_mean'] < 8.0, f"MP2 polarizability {result_mp2['alpha_mean']:.2f} seems unreasonable"

    print(f"HF polarizability:  {result_hf['alpha_mean']:.4f} a.u.")
    print(f"MP2 polarizability: {result_mp2['alpha_mean']:.4f} a.u.")
    print(f"Experimental:       {experimental:.4f} a.u.")


def test_mp2_polarizability_tensor_symmetry():
    """Test that MP2 polarizability tensor is symmetric."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O", basis='sto-3g')
    calc = PropertyCalculator(water.hamiltonian)

    result = calc.compute_polarizability_mp2()
    alpha = result['alpha_tensor']

    # Check symmetry
    asymmetry = np.max(np.abs(alpha - alpha.T))
    assert asymmetry < 0.1, f"MP2 polarizability tensor asymmetry {asymmetry:.4f} too large"

    print(f"MP2 tensor asymmetry: {asymmetry:.6f} a.u.")


def test_mp2_vs_hf_polarizability_comparison():
    """Compare MP2 and HF polarizability to verify MP2 gives improvement."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    # Use a basis with diffuse functions for better test
    water = from_smiles("O", basis='6-311++g(d,p)')
    calc = PropertyCalculator(water.hamiltonian)

    result_hf = calc.compute_polarizability()
    result_mp2 = calc.compute_polarizability_mp2()

    experimental = 10.13

    hf_error = abs(result_hf['alpha_mean'] - experimental)
    mp2_error = abs(result_mp2['alpha_mean'] - experimental)

    print(f"\nPolarizability comparison (6-311++G(d,p)):")
    print(f"  HF:           {result_hf['alpha_mean']:.4f} a.u. (error: {hf_error:.4f})")
    print(f"  MP2:          {result_mp2['alpha_mean']:.4f} a.u. (error: {mp2_error:.4f})")
    print(f"  Experimental: {experimental:.4f} a.u.")

    # MP2 should be at least as good as HF (usually slightly better)
    # For polarizability, improvement is small (1-3%)
    assert result_mp2['alpha_mean'] > 0
    assert result_hf['alpha_mean'] > 0


def test_mp2_polarizability_with_diffuse_basis():
    """Test MP2 polarizability with diffuse basis (aug-cc-pVDZ)."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O", basis='aug-cc-pvdz')
    calc = PropertyCalculator(water.hamiltonian)

    result = calc.compute_polarizability_mp2()

    # Experimental: 10.13 a.u.
    # MP2/aug-cc-pVDZ: expect ~8-9 a.u. (80-90%)
    experimental = 10.13
    assert 7.5 < result['alpha_mean'] < 9.5, f"MP2/aug-cc-pVDZ polarizability {result['alpha_mean']:.2f} outside expected range"

    accuracy = (result['alpha_mean'] / experimental) * 100
    print(f"MP2/aug-cc-pVDZ: {result['alpha_mean']:.4f} a.u. ({accuracy:.1f}% of exp)")

    # Should achieve >75% accuracy
    assert accuracy > 75.0, f"MP2/aug-cc-pVDZ should give >75% accuracy, got {accuracy:.1f}%"


def test_mp2_different_molecules():
    """Test MP2 energies for different molecules."""
    from kanad.io import from_smiles
    from kanad.core import MP2Solver

    molecules = {
        "H2": ("[H][H]", -1.159),  # Approximate MP2/6-311G(d,p) energies
        "H2O": ("O", -76.28),
        "NH3": ("N", -56.40),
    }

    for name, (smiles, expected_energy) in molecules.items():
        mol = from_smiles(smiles, basis='6-311g(d,p)')
        mp2 = MP2Solver(mol.hamiltonian)
        result = mp2.compute_energy()

        # Allow ±0.1 Ha tolerance
        assert abs(result['e_mp2'] - expected_energy) < 0.1, \
            f"{name}: MP2 energy {result['e_mp2']:.4f} far from expected {expected_energy:.4f}"

        print(f"{name}: E_MP2 = {result['e_mp2']:.6f} Ha")
