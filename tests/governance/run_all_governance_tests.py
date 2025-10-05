"""
Master Governance Validation Script

Runs all three bonding type governance validations:
1. Ionic (LiH)
2. Covalent (H2)
3. Metallic (Na2)
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

# Import test modules
from tests.governance.test_ionic_governance import test_ionic_lih
from tests.governance.test_covalent_governance import test_covalent_h2
from tests.governance.test_metallic_governance import test_metallic_na2


def main():
    """Run all governance validation tests."""
    print("\n" + "=" * 80)
    print("KANAD FRAMEWORK - COMPREHENSIVE GOVERNANCE VALIDATION")
    print("=" * 80)
    print("\nTesting all three bonding types:")
    print("  1. IONIC - LiH (Lithium Hydride)")
    print("  2. COVALENT - H2 (Hydrogen)")
    print("  3. METALLIC - Na2 (Sodium Dimer)")
    print("\n" + "=" * 80)

    results = {}

    # Test 1: Ionic
    print("\n\n")
    print("▶ TEST 1/3: IONIC BONDING")
    print("=" * 80)
    try:
        results['ionic'] = test_ionic_lih()
    except Exception as e:
        print(f"✗ IONIC test failed with error: {e}")
        results['ionic'] = False

    # Test 2: Covalent
    print("\n\n")
    print("▶ TEST 2/3: COVALENT BONDING")
    print("=" * 80)
    try:
        results['covalent'] = test_covalent_h2()
    except Exception as e:
        print(f"✗ COVALENT test failed with error: {e}")
        results['covalent'] = False

    # Test 3: Metallic
    print("\n\n")
    print("▶ TEST 3/3: METALLIC BONDING")
    print("=" * 80)
    try:
        results['metallic'] = test_metallic_na2()
    except Exception as e:
        print(f"✗ METALLIC test failed with error: {e}")
        results['metallic'] = False

    # Summary
    print("\n\n")
    print("=" * 80)
    print("FINAL GOVERNANCE VALIDATION SUMMARY")
    print("=" * 80)

    print("\nTest Results:")
    for bonding_type, passed in results.items():
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"  {bonding_type.upper():<12} {status}")

    # Overall result
    print("\n" + "-" * 80)
    all_passed = all(results.values())
    passed_count = sum(results.values())
    total_count = len(results)

    if all_passed:
        print(f"🎉 ALL {total_count} GOVERNANCE TESTS PASSED")
        print("\n✓ Ionic governance: OPERATIONAL")
        print("✓ Covalent governance: OPERATIONAL")
        print("✓ Metallic governance: OPERATIONAL")
        print("\n→ Kanad governance-driven quantum chemistry VALIDATED!")
    else:
        print(f"⚠ {passed_count}/{total_count} tests passed")
        print("\nFailed tests may indicate:")
        print("  - Missing basis set implementation (covalent)")
        print("  - Pending Hamiltonian integration (metallic)")

    print("=" * 80)

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
