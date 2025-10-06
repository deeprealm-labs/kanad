"""
Unit tests for PropertyCalculator

Tests dipole moment calculations for various molecules:
- H2O (known dipole ~1.85 D experimental, ~2.1-2.3 D HF/STO-3G)
- NH3 (known dipole ~1.47 D experimental)
- CO2 (should be 0.0 - linear symmetric)
- CH4 (should be 0.0 - tetrahedral symmetric)
- H2 (should be 0.0 - homonuclear)
- Validation against PySCF
"""

import pytest
import numpy as np


def test_import():
    """Test that PropertyCalculator can be imported."""
    from kanad.analysis import PropertyCalculator
    assert PropertyCalculator is not None


def test_water_dipole():
    """Test H2O dipole moment."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O", name="Water")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_dipole_moment()

    # Experimental: 1.85 D
    # HF/STO-3G typically gives ~2.1-2.3 D (overestimates)
    assert result is not None
    assert 'dipole_magnitude' in result
    assert 'dipole_vector' in result
    assert 'components' in result

    # Check magnitude is reasonable
    dipole = result['dipole_magnitude']
    assert 1.5 < dipole < 2.5, f"Water dipole {dipole:.3f} D outside expected range"

    # Check vector has 3 components
    assert len(result['dipole_vector']) == 3
    assert len(result['dipole_au']) == 3

    # Check components sum correctly
    vec = result['dipole_vector']
    mag_from_components = np.linalg.norm(vec)
    np.testing.assert_almost_equal(mag_from_components, dipole, decimal=6)

    print(f"H2O dipole: {dipole:.3f} D (experimental: 1.85 D)")


def test_ammonia_dipole():
    """Test NH3 dipole moment."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    ammonia = from_smiles("N", name="Ammonia")
    calc = PropertyCalculator(ammonia.hamiltonian)
    result = calc.compute_dipole_moment()

    # Experimental: 1.47 D
    # HF/STO-3G: expect ~1.5-1.9 D
    dipole = result['dipole_magnitude']
    assert 1.2 < dipole < 2.2, f"NH3 dipole {dipole:.3f} D outside expected range"

    print(f"NH3 dipole: {dipole:.3f} D (experimental: 1.47 D)")


def test_symmetric_molecules_zero_dipole():
    """Test that symmetric molecules have zero dipole."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    # H2 (homonuclear diatomic)
    h2 = from_smiles("[H][H]", name="Hydrogen")
    calc_h2 = PropertyCalculator(h2.hamiltonian)
    result_h2 = calc_h2.compute_dipole_moment()
    assert result_h2['dipole_magnitude'] < 0.01, "H2 should have zero dipole"

    # CH4 (tetrahedral)
    ch4 = from_smiles("C", name="Methane")
    calc_ch4 = PropertyCalculator(ch4.hamiltonian)
    result_ch4 = calc_ch4.compute_dipole_moment()
    assert result_ch4['dipole_magnitude'] < 0.1, "CH4 should have near-zero dipole"

    print(f"H2 dipole: {result_h2['dipole_magnitude']:.6f} D (expected: 0.0)")
    print(f"CH4 dipole: {result_ch4['dipole_magnitude']:.6f} D (expected: 0.0)")


def test_dipole_components():
    """Test dipole vector components."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_dipole_moment()

    # Check components dict
    assert 'x' in result['components']
    assert 'y' in result['components']
    assert 'z' in result['components']

    # Check components match vector
    vec = result['dipole_vector']
    assert result['components']['x'] == vec[0]
    assert result['components']['y'] == vec[1]
    assert result['components']['z'] == vec[2]


def test_dipole_units():
    """Test unit conversions (a.u. and Debye)."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_dipole_moment()

    # Check that both a.u. and Debye are provided
    mu_au = result['dipole_au']
    mu_debye = result['dipole_vector']

    # Check conversion factor (1 a.u. = 2.541746 Debye)
    AU_TO_DEBYE = 2.541746
    expected_debye = mu_au * AU_TO_DEBYE

    np.testing.assert_array_almost_equal(mu_debye, expected_debye, decimal=6)


def test_pyscf_verification():
    """Test that Kanad dipole matches PySCF."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)

    verification = calc.verify_dipole_with_pyscf()

    # Check that results are close (< 0.01 D difference)
    assert verification['agree'], (
        f"Kanad ({verification['kanad_dipole']:.4f} D) differs from "
        f"PySCF ({verification['pyscf_dipole']:.4f} D) by "
        f"{verification['difference']:.4f} D"
    )

    print(f"Kanad:  {verification['kanad_dipole']:.4f} D")
    print(f"PySCF:  {verification['pyscf_dipole']:.4f} D")
    print(f"Difference: {verification['difference']:.6f} D")


def test_electronic_nuclear_contributions():
    """Test that electronic and nuclear contributions are reported."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_dipole_moment()

    assert 'electronic_contribution' in result
    assert 'nuclear_contribution' in result

    # Check that they sum to total
    mu_elec = result['electronic_contribution']
    mu_nuc = result['nuclear_contribution']
    mu_total = result['dipole_vector']

    np.testing.assert_array_almost_equal(mu_elec + mu_nuc, mu_total, decimal=6)


def test_center_of_mass():
    """Test center of mass calculation."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)

    com = calc.compute_center_of_mass()

    assert com.shape == (3,)
    # COM should be closer to oxygen (heavier)
    # Not at origin for general geometry


def test_center_of_charge():
    """Test center of nuclear charge calculation."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)

    coc = calc.compute_center_of_charge()

    assert coc.shape == (3,)


def test_dipole_custom_density():
    """Test dipole with custom density matrix."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator
    import numpy as np

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)

    # Get HF density
    dm_hf = water.hamiltonian.mf.make_rdm1()

    # Compute with explicit density
    result = calc.compute_dipole_moment(density_matrix=dm_hf)

    assert result['dipole_magnitude'] > 0


def test_multiple_molecules():
    """Test dipole for multiple molecules."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    molecules = {
        "H2O": from_smiles("O"),
        "NH3": from_smiles("N"),
        "CH4": from_smiles("C"),
    }

    results = {}
    for name, mol in molecules.items():
        calc = PropertyCalculator(mol.hamiltonian)
        result = calc.compute_dipole_moment()
        results[name] = result['dipole_magnitude']

    print("\nDipole moments:")
    for name, dipole in results.items():
        print(f"  {name}: {dipole:.3f} D")

    # Water should have largest dipole
    assert results["H2O"] > results["CH4"]
    assert results["NH3"] > results["CH4"]


# =============================================================================
# POLARIZABILITY TESTS (Feature 1.3)
# =============================================================================

def test_polarizability_h2():
    """Test H2 polarizability.

    Experimental: α = 5.4 a.u. = 0.80 Å³
    HF/STO-3G: ~0.93 a.u. (minimal basis underestimates, needs polarization functions)

    Note: Accurate polarizability requires basis sets with polarization functions (d,p).
    STO-3G only has s-orbitals and cannot describe perpendicular polarization.
    """
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    h2 = from_smiles("[H][H]", name="Hydrogen")
    calc = PropertyCalculator(h2.hamiltonian)
    result = calc.compute_polarizability()

    assert result is not None
    assert 'alpha_tensor' in result
    assert 'alpha_mean' in result
    assert 'alpha_mean_angstrom3' in result

    alpha = result['alpha_mean']
    alpha_angstrom = result['alpha_mean_angstrom3']

    # Check within reasonable range for STO-3G (0.5-2.0 a.u.)
    # STO-3G severely underestimates due to lack of polarization functions
    assert 0.5 < alpha < 2.0, f"H2 polarizability {alpha:.2f} a.u. outside expected range for STO-3G"

    # Check unit conversion (1 a.u. = 0.1482 Å³)
    expected_angstrom = alpha * 0.1482
    np.testing.assert_almost_equal(alpha_angstrom, expected_angstrom, decimal=6)

    print(f"H2 polarizability: {alpha:.2f} a.u. = {alpha_angstrom:.3f} Å³ (experimental: 5.4 a.u., STO-3G underestimates)")


def test_polarizability_h2o():
    """Test H2O polarizability.

    Experimental: α = 9.8 a.u. = 1.45 Å³
    HF/STO-3G: ~4-6 a.u. (minimal basis underestimates)

    Note: STO-3G includes p-orbitals on O, so results are better than H2.
    Still underestimates due to lack of polarization functions and correlation.
    """
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O", name="Water")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_polarizability()

    alpha = result['alpha_mean']
    alpha_angstrom = result['alpha_mean_angstrom3']

    # Check within reasonable range for STO-3G (2-8 a.u.)
    assert 2.0 < alpha < 8.0, f"H2O polarizability {alpha:.2f} a.u. outside expected range for STO-3G"

    print(f"H2O polarizability: {alpha:.2f} a.u. = {alpha_angstrom:.3f} Å³ (experimental: 9.8 a.u., STO-3G underestimates)")


def test_polarizability_tensor_shape():
    """Test that polarizability tensor has correct shape."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_polarizability()

    alpha_tensor = result['alpha_tensor']

    # Should be 3×3 matrix
    assert alpha_tensor.shape == (3, 3)

    # Should contain finite values
    assert np.all(np.isfinite(alpha_tensor))


def test_polarizability_tensor_symmetry():
    """Test that polarizability tensor is symmetric (α_ij ≈ α_ji)."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_polarizability()

    alpha = result['alpha_tensor']

    # Check symmetry: α_ij should equal α_ji
    asymmetry = np.max(np.abs(alpha - alpha.T))

    # Tolerance: finite field numerical noise
    assert asymmetry < 0.1, f"Polarizability tensor asymmetry {asymmetry:.4f} too large"

    print(f"Polarizability tensor asymmetry: {asymmetry:.6f} a.u.")


def test_polarizability_mean_calculation():
    """Test mean polarizability: ᾱ = (α_xx + α_yy + α_zz) / 3."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_polarizability()

    alpha_tensor = result['alpha_tensor']
    alpha_mean = result['alpha_mean']

    # Mean should be trace/3
    expected_mean = np.trace(alpha_tensor) / 3.0

    np.testing.assert_almost_equal(alpha_mean, expected_mean, decimal=8)


def test_polarizability_diagonal_elements():
    """Test that diagonal elements (α_xx, α_yy, α_zz) are returned."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_polarizability()

    assert 'alpha_xx' in result
    assert 'alpha_yy' in result
    assert 'alpha_zz' in result

    # Check they match tensor diagonal
    alpha_tensor = result['alpha_tensor']
    assert result['alpha_xx'] == alpha_tensor[0, 0]
    assert result['alpha_yy'] == alpha_tensor[1, 1]
    assert result['alpha_zz'] == alpha_tensor[2, 2]


def test_polarizability_eigenvalues():
    """Test that principal polarizabilities (eigenvalues) are computed."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_polarizability()

    eigenvalues = result['eigenvalues']

    # Should have 3 eigenvalues
    assert eigenvalues.shape == (3,)

    # All eigenvalues should be positive (polarizability is positive-definite)
    assert np.all(eigenvalues > 0)

    # Eigenvalues should be in ascending order (eigvalsh returns sorted)
    assert eigenvalues[0] <= eigenvalues[1] <= eigenvalues[2]


def test_polarizability_anisotropy():
    """Test polarizability anisotropy calculation."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_polarizability()

    anisotropy = result['alpha_anisotropy']

    # Anisotropy should be positive
    assert anisotropy >= 0

    # Check formula: Δα = √(3/2 ||α - ᾱI||_F)
    alpha_tensor = result['alpha_tensor']
    alpha_mean = result['alpha_mean']

    alpha_iso = alpha_mean * np.eye(3)
    alpha_aniso_tensor = alpha_tensor - alpha_iso
    expected_aniso = np.sqrt(1.5 * np.sum(alpha_aniso_tensor**2))

    np.testing.assert_almost_equal(anisotropy, expected_aniso, decimal=8)


def test_polarizability_unit_conversion():
    """Test polarizability unit conversion (a.u. ↔ Å³)."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_polarizability()

    alpha_au = result['alpha_mean']
    alpha_angstrom = result['alpha_mean_angstrom3']

    # Check conversion: 1 a.u. = 0.1482 Å³
    AU_TO_ANGSTROM3 = 0.1482
    expected_angstrom = alpha_au * AU_TO_ANGSTROM3

    np.testing.assert_almost_equal(alpha_angstrom, expected_angstrom, decimal=8)


def test_polarizability_field_strength_sensitivity():
    """Test that polarizability is relatively insensitive to field strength."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    h2 = from_smiles("[H][H]", name="Hydrogen")
    calc = PropertyCalculator(h2.hamiltonian)

    # Try different field strengths
    field_strengths = [0.0005, 0.001, 0.002]
    alphas = []

    for field in field_strengths:
        result = calc.compute_polarizability(field_strength=field)
        alphas.append(result['alpha_mean'])

    # All results should be within 5% of each other
    alpha_avg = np.mean(alphas)
    for alpha in alphas:
        rel_diff = abs(alpha - alpha_avg) / alpha_avg
        assert rel_diff < 0.05, f"Field strength sensitivity too large: {rel_diff:.1%}"

    print(f"Field strength test: {alphas[0]:.3f}, {alphas[1]:.3f}, {alphas[2]:.3f} a.u.")


def test_polarizability_metadata():
    """Test that polarizability result includes metadata."""
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    water = from_smiles("O")
    calc = PropertyCalculator(water.hamiltonian)
    result = calc.compute_polarizability(method='finite_field', field_strength=0.001)

    assert 'method' in result
    assert result['method'] == 'finite_field'

    assert 'field_strength' in result
    assert result['field_strength'] == 0.001


def test_polarizability_positive_values():
    """Test that polarizability components are positive (or zero for minimal basis).

    Note: STO-3G for H2 has no p-functions, so perpendicular components may be zero.
    This is a basis set limitation, not an error in the implementation.
    """
    from kanad.io import from_smiles
    from kanad.analysis import PropertyCalculator

    molecules = [
        from_smiles("[H][H]", name="H2"),
        from_smiles("O", name="H2O"),
        from_smiles("N", name="NH3"),
    ]

    for mol in molecules:
        calc = PropertyCalculator(mol.hamiltonian)
        result = calc.compute_polarizability()

        # Mean polarizability should be positive
        assert result['alpha_mean'] > 0, f"{mol.name}: negative mean polarizability"

        # Diagonal elements should be non-negative (can be zero for minimal basis)
        assert result['alpha_xx'] >= -1e-10, f"{mol.name}: negative alpha_xx"
        assert result['alpha_yy'] >= -1e-10, f"{mol.name}: negative alpha_yy"
        assert result['alpha_zz'] >= -1e-10, f"{mol.name}: negative alpha_zz"

        # Eigenvalues should be non-negative
        assert np.all(result['eigenvalues'] >= -1e-10), f"{mol.name}: negative eigenvalues"
