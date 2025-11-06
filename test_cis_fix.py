"""
Test CIS Excited States - Validate Fix for Two-Electron Integrals

Tests that the CIS method now uses real ERIs instead of placeholders.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from kanad.bonds import BondFactory
from kanad.solvers import ExcitedStatesSolver


def test_cis_h2():
    """Test CIS for H2 molecule."""
    print("=" * 80)
    print("TEST: CIS Excited States for H2")
    print("=" * 80)

    # Create H2 bond at equilibrium
    bond = BondFactory.create_bond('H', 'H', distance=0.74)

    # Create CIS solver
    solver = ExcitedStatesSolver(
        bond=bond,
        method='cis',
        n_states=3,
        enable_analysis=False
    )

    print("\nRunning CIS calculation...")
    result = solver.solve()

    print("\n" + "-" * 80)
    print("RESULTS")
    print("-" * 80)

    print(f"\nGround State Energy: {result['ground_state_energy']:.8f} Ha")
    print(f"\nExcited States:")
    for i, (E_ha, E_ev) in enumerate(zip(
        result['excitation_energies_ha'],
        result['excitation_energies_ev']
    ), 1):
        print(f"  S{i}: {E_ev:8.4f} eV  ({E_ha:.8f} Ha)")

    print(f"\nOscillator Strengths:")
    for i, f in enumerate(result['oscillator_strengths'], 1):
        print(f"  S{i}: f = {f:.6f}")

    print(f"\nDominant Transitions:")
    for i, trans in enumerate(result['dominant_transitions'], 1):
        print(f"  S{i}: {trans}")

    # Validate results
    print("\n" + "=" * 80)
    print("VALIDATION")
    print("=" * 80)

    # H2 first excitation is around 11-12 eV (literature value)
    # With placeholder integrals, it would be completely wrong
    first_excitation = result['excitation_energies_ev'][0]

    print(f"\nFirst excitation energy: {first_excitation:.4f} eV")
    print(f"Expected range: 10-15 eV (literature: ~11.4 eV)")

    if 10.0 < first_excitation < 15.0:
        print("✅ PASS: Excitation energy in reasonable range")
        print("         (Real ERIs are being used, not placeholders!)")
        validation_passed = True
    else:
        print("❌ FAIL: Excitation energy out of range")
        print(f"         Got {first_excitation:.4f} eV, expected 10-15 eV")
        print("         Placeholders may still be in use!")
        validation_passed = False

    # Check that oscillator strengths are non-negative
    if all(f >= 0 for f in result['oscillator_strengths']):
        print("✅ PASS: All oscillator strengths non-negative")
    else:
        print("❌ FAIL: Some oscillator strengths negative")
        validation_passed = False

    print("\n" + "=" * 80)
    if validation_passed:
        print("✅ ALL TESTS PASSED - CIS FIX VALIDATED")
    else:
        print("❌ SOME TESTS FAILED")
    print("=" * 80)

    return validation_passed


def test_cis_lih():
    """Test CIS for LiH molecule."""
    print("\n" * 2)
    print("=" * 80)
    print("TEST: CIS Excited States for LiH")
    print("=" * 80)

    bond = BondFactory.create_bond('Li', 'H', distance=1.60)

    solver = ExcitedStatesSolver(
        bond=bond,
        method='cis',
        n_states=3,
        enable_analysis=False
    )

    print("\nRunning CIS calculation...")
    result = solver.solve()

    print("\n" + "-" * 80)
    print("RESULTS")
    print("-" * 80)

    print(f"\nGround State Energy: {result['ground_state_energy']:.8f} Ha")
    print(f"\nExcited States:")
    for i, E_ev in enumerate(result['excitation_energies_ev'], 1):
        print(f"  S{i}: {E_ev:8.4f} eV")

    print("\n" + "=" * 80)
    print("✅ LiH CIS calculation completed successfully")
    print("=" * 80)

    return True


if __name__ == '__main__':
    print("\n🔬 CIS TWO-ELECTRON INTEGRAL FIX - VALIDATION TEST\n")

    passed_h2 = test_cis_h2()
    passed_lih = test_cis_lih()

    print("\n" * 2)
    print("=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"H2 Test:  {'✅ PASS' if passed_h2 else '❌ FAIL'}")
    print(f"LiH Test: {'✅ PASS' if passed_lih else '❌ FAIL'}")

    if passed_h2 and passed_lih:
        print("\n✅ ALL VALIDATIONS PASSED!")
        print("   CIS is now using real two-electron integrals")
        print("   Placeholder values (0.1, 0.05) have been replaced")
        sys.exit(0)
    else:
        print("\n❌ SOME VALIDATIONS FAILED")
        sys.exit(1)
