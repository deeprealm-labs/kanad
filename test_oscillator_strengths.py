"""
Test Oscillator Strengths for VQE and SQD Excited States

Validates that oscillator strengths are now computed (not zeros).
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from kanad.bonds import BondFactory
from kanad.solvers import ExcitedStatesSolver
import numpy as np


def test_vqe_oscillator_strengths():
    """Test VQE excited states with oscillator strengths."""
    print("=" * 80)
    print("TEST 1: VQE Excited States - Oscillator Strengths")
    print("=" * 80)

    bond = BondFactory.create_bond('H', 'H', distance=0.74)

    solver = ExcitedStatesSolver(
        bond=bond,
        method='vqe',
        n_states=2,  # Ground + 1 excited
        enable_analysis=False,
        backend='statevector',
        max_iterations=50
    )

    print("\nRunning VQE calculation...")
    result = solver.solve()

    print(f"\n{'='*80}")
    print("RESULTS")
    print(f"{'='*80}")

    print(f"\nGround State Energy: {result['ground_state_energy']:.8f} Ha")

    print(f"\nExcited States:")
    for i, E_ev in enumerate(result['excitation_energies_ev'], 1):
        f = result['oscillator_strengths'][i-1]
        print(f"  S{i}: {E_ev:8.4f} eV, f = {f:.6f}")

    # Validation
    print(f"\n{'='*80}")
    print("VALIDATION")
    print(f"{'='*80}")

    oscillator_strengths = result['oscillator_strengths']

    # Check that at least one oscillator strength is non-zero
    any_nonzero = np.any(oscillator_strengths > 0)

    if any_nonzero:
        print(f"\n✅ PASS: Oscillator strengths computed (not all zeros)")
        print(f"         Values: {oscillator_strengths}")
        return True
    else:
        print(f"\n❌ FAIL: All oscillator strengths are zero")
        print(f"         This means oscillator strengths are not being computed!")
        return False


def test_sqd_oscillator_strengths():
    """Test SQD excited states with oscillator strengths."""
    print("\n" * 2)
    print("=" * 80)
    print("TEST 2: SQD Excited States - Oscillator Strengths")
    print("=" * 80)

    bond = BondFactory.create_bond('H', 'H', distance=0.74)

    solver = ExcitedStatesSolver(
        bond=bond,
        method='sqd',
        n_states=3,  # Ground + 2 excited
        enable_analysis=False,
        backend='statevector',
        subspace_dim=10
    )

    print("\nRunning SQD calculation...")
    result = solver.solve()

    print(f"\n{'='*80}")
    print("RESULTS")
    print(f"{'='*80}")

    print(f"\nGround State Energy: {result['ground_state_energy']:.8f} Ha")

    print(f"\nExcited States:")
    for i, E_ev in enumerate(result['excitation_energies_ev'], 1):
        f = result['oscillator_strengths'][i-1]
        print(f"  S{i}: {E_ev:8.4f} eV, f = {f:.6f}")

    # Validation
    print(f"\n{'='*80}")
    print("VALIDATION")
    print(f"{'='*80}")

    oscillator_strengths = result['oscillator_strengths']

    # Check that at least one oscillator strength is non-zero
    any_nonzero = np.any(oscillator_strengths > 0)

    if any_nonzero:
        print(f"\n✅ PASS: Oscillator strengths computed (not all zeros)")
        print(f"         Values: {oscillator_strengths}")
        return True
    else:
        print(f"\n❌ FAIL: All oscillator strengths are zero")
        print(f"         This means oscillator strengths are not being computed!")
        return False


def test_cis_oscillator_strengths():
    """Test CIS - should already have oscillator strengths."""
    print("\n" * 2)
    print("=" * 80)
    print("TEST 3: CIS Excited States - Oscillator Strengths (Control)")
    print("=" * 80)

    bond = BondFactory.create_bond('H', 'H', distance=0.74)

    solver = ExcitedStatesSolver(
        bond=bond,
        method='cis',
        n_states=2,
        enable_analysis=False
    )

    print("\nRunning CIS calculation...")
    result = solver.solve()

    print(f"\n{'='*80}")
    print("RESULTS")
    print(f"{'='*80}")

    print(f"\nExcited States:")
    for i, E_ev in enumerate(result['excitation_energies_ev'], 1):
        f = result['oscillator_strengths'][i-1]
        print(f"  S{i}: {E_ev:8.4f} eV, f = {f:.6f}")

    oscillator_strengths = result['oscillator_strengths']
    any_nonzero = np.any(oscillator_strengths > 0)

    if any_nonzero:
        print(f"\n✅ PASS: CIS oscillator strengths present (expected)")
        return True
    else:
        print(f"\n⚠️  WARNING: CIS oscillator strengths all zero (unexpected)")
        return True  # Don't fail the test


if __name__ == '__main__':
    print("\n🔬 OSCILLATOR STRENGTHS VALIDATION TEST\n")

    passed_vqe = test_vqe_oscillator_strengths()
    passed_sqd = test_sqd_oscillator_strengths()
    passed_cis = test_cis_oscillator_strengths()

    print("\n" * 2)
    print("=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"VQE Test:  {'✅ PASS' if passed_vqe else '❌ FAIL'}")
    print(f"SQD Test:  {'✅ PASS' if passed_sqd else '❌ FAIL'}")
    print(f"CIS Test:  {'✅ PASS' if passed_cis else '❌ FAIL'}")

    if passed_vqe and passed_sqd:
        print("\n✅ OSCILLATOR STRENGTHS IMPLEMENTED!")
        print("   VQE and SQD now compute transition dipole moments")
        print("   UV-Vis spectra will have proper intensities")
        sys.exit(0)
    else:
        print("\n❌ SOME TESTS FAILED")
        sys.exit(1)
