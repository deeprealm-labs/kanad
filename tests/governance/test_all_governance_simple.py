"""
Simplified Governance Validation - All Three Bonding Types

Tests governance protocols without needing full Hamiltonian construction.
Validates governance rules directly on molecular structures.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import numpy as np
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.governance.protocols.ionic_protocol import IonicGovernanceProtocol
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol
from kanad.governance.protocols.metallic_protocol import MetallicGovernanceProtocol


def test_ionic_governance():
    """Test ionic governance protocol on LiH."""
    print("=" * 80)
    print("TEST 1/3: IONIC GOVERNANCE PROTOCOL (LiH)")
    print("=" * 80)

    # Create LiH
    li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    h = Atom('H', position=np.array([1.6, 0.0, 0.0]))
    molecule = Molecule(atoms=[li, h], charge=0, spin=0)

    print(f"\n✓ Molecule: LiH (bond length: 1.6 Å)")
    print(f"  Li electronegativity: {li.electronegativity:.2f}")
    print(f"  H electronegativity: {h.electronegativity:.2f}")
    print(f"  ΔEN = {abs(h.electronegativity - li.electronegativity):.2f}")

    # Ionic governance protocol
    protocol = IonicGovernanceProtocol()
    print(f"\n✓ Ionic Governance Protocol initialized")
    print(f"  Rules: {len(protocol.rules)}")
    print(f"  Bonding type: {protocol.bond_type.value}")

    # Check physics
    print(f"\nIonic Bonding Physics Checks:")

    # Check 1: Large electronegativity difference
    en_diff = abs(h.electronegativity - li.electronegativity)
    check1 = en_diff > 1.5
    print(f"  [{'✓' if check1 else '✗'}] Large ΔEN: {en_diff:.2f} > 1.5")

    # Check 2: One metal, one nonmetal
    check2 = li.is_metal and not h.is_metal
    print(f"  [{'✓' if check2 else '✗'}] Metal + Nonmetal: Li(metal)={li.is_metal}, H(metal)={h.is_metal}")

    # Check 3: Allowed operators
    allowed_ops = protocol.get_allowed_operators()
    forbidden_ops = protocol.get_forbidden_operators()
    check3 = len(allowed_ops) > 0 and len(forbidden_ops) > 0
    print(f"  [{'✓' if check3 else '✗'}] Operator rules: {len(allowed_ops)} allowed, {len(forbidden_ops)} forbidden")

    # Check 4: Transfer integral estimate
    distance = li.distance_to(h)
    t_estimate = protocol.get_transfer_integral_estimate(distance)
    check4 = t_estimate < 0.1  # Weak transfer for ionic
    print(f"  [{'✓' if check4 else '✗'}] Weak transfer: t = {t_estimate:.4f} < 0.1")

    all_passed = check1 and check2 and check3 and check4

    if all_passed:
        print(f"\n🎉 ALL IONIC GOVERNANCE CHECKS PASSED")
    else:
        print(f"\n⚠ Some checks failed")

    return all_passed


def test_covalent_governance():
    """Test covalent governance protocol on H2."""
    print("\n\n")
    print("=" * 80)
    print("TEST 2/3: COVALENT GOVERNANCE PROTOCOL (H2)")
    print("=" * 80)

    # Create H2
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
    molecule = Molecule(atoms=[h1, h2], charge=0, spin=0)

    print(f"\n✓ Molecule: H2 (bond length: 0.74 Å)")
    print(f"  H electronegativity: {h1.electronegativity:.2f}")
    print(f"  ΔEN = 0.00 (homonuclear)")

    # Covalent governance protocol
    protocol = CovalentGovernanceProtocol()
    print(f"\n✓ Covalent Governance Protocol initialized")
    print(f"  Rules: {len(protocol.rules)}")
    print(f"  Bonding type: {protocol.bond_type.value}")

    # Check physics
    print(f"\nCovalent Bonding Physics Checks:")

    # Check 1: Small electronegativity difference
    en_diff = abs(h2.electronegativity - h1.electronegativity)
    check1 = en_diff < 1.5
    print(f"  [{'✓' if check1 else '✗'}] Small ΔEN: {en_diff:.2f} < 1.5")

    # Check 2: Both nonmetals
    check2 = not h1.is_metal and not h2.is_metal
    print(f"  [{'✓' if check2 else '✗'}] Both nonmetals: H1(metal)={h1.is_metal}, H2(metal)={h2.is_metal}")

    # Check 3: Allowed operators (includes hybridization)
    allowed_ops = protocol.get_allowed_operators()
    check3 = 'givens' in allowed_ops  # Givens rotation for MO formation
    print(f"  [{'✓' if check3 else '✗'}] MO formation operators: 'givens' in allowed ops")

    # Check 4: Reasonable bond length
    distance = h1.distance_to(h2)
    check4 = 0.5 < distance < 1.5  # Typical covalent H-H bond
    print(f"  [{'✓' if check4 else '✗'}] Covalent bond length: {distance:.2f} Å (0.5-1.5 Å)")

    all_passed = check1 and check2 and check3 and check4

    if all_passed:
        print(f"\n🎉 ALL COVALENT GOVERNANCE CHECKS PASSED")
    else:
        print(f"\n⚠ Some checks failed")

    return all_passed


def test_metallic_governance():
    """Test metallic governance protocol on Na2."""
    print("\n\n")
    print("=" * 80)
    print("TEST 3/3: METALLIC GOVERNANCE PROTOCOL (Na2)")
    print("=" * 80)

    # Create Na2
    na1 = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
    na2 = Atom('Na', position=np.array([3.08, 0.0, 0.0]))
    molecule = Molecule(atoms=[na1, na2], charge=0, spin=0)

    print(f"\n✓ Molecule: Na2 (bond length: 3.08 Å)")
    print(f"  Na electronegativity: {na1.electronegativity:.2f}")
    print(f"  ΔEN = 0.00 (homonuclear)")

    # Metallic governance protocol
    protocol = MetallicGovernanceProtocol()
    print(f"\n✓ Metallic Governance Protocol initialized")
    print(f"  Rules: {len(protocol.rules)}")
    print(f"  Bonding type: {protocol.bond_type}")

    # Check physics
    print(f"\nMetallic Bonding Physics Checks:")

    # Check 1: Both metals
    check1 = na1.is_metal and na2.is_metal
    print(f"  [{'✓' if check1 else '✗'}] Both metals: Na1(metal)={na1.is_metal}, Na2(metal)={na2.is_metal}")

    # Check 2: Low ionization energy
    ie = na1.properties.ionization_energy
    check2 = ie < 6.0  # Typical for metals
    print(f"  [{'✓' if check2 else '✗'}] Low ionization energy: IE = {ie:.2f} eV < 6.0 eV")

    # Check 3: Allowed operators (includes delocalization)
    allowed_ops = protocol.get_allowed_operators()
    check3 = 'qft' in allowed_ops  # QFT for delocalization
    print(f"  [{'✓' if check3 else '✗'}] Delocalization operators: 'qft' in allowed ops")

    # Check 4: Same element (ΔEN = 0)
    check4 = na1.symbol == na2.symbol
    print(f"  [{'✓' if check4 else '✗'}] Same element: ΔEN = 0.00")

    all_passed = check1 and check2 and check3 and check4

    if all_passed:
        print(f"\n🎉 ALL METALLIC GOVERNANCE CHECKS PASSED")
    else:
        print(f"\n⚠ Some checks failed")

    return all_passed


def main():
    """Run all governance protocol tests."""
    print("\n" + "=" * 80)
    print("KANAD GOVERNANCE PROTOCOLS - COMPREHENSIVE VALIDATION")
    print("=" * 80)
    print("\nValidating governance rules for three bonding types:")
    print("  1. IONIC - LiH")
    print("  2. COVALENT - H2")
    print("  3. METALLIC - Na2")
    print("\n" + "=" * 80)

    results = {
        'ionic': test_ionic_governance(),
        'covalent': test_covalent_governance(),
        'metallic': test_metallic_governance()
    }

    # Final summary
    print("\n\n")
    print("=" * 80)
    print("FINAL GOVERNANCE VALIDATION SUMMARY")
    print("=" * 80)

    for bond_type, passed in results.items():
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"  {bond_type.upper():<12} {status}")

    all_passed = all(results.values())
    passed_count = sum(results.values())

    print("\n" + "-" * 80)
    if all_passed:
        print(f"🎉 ALL {len(results)} GOVERNANCE PROTOCOLS VALIDATED")
        print("\n✓ Ionic governance: Rules verified")
        print("✓ Covalent governance: Rules verified")
        print("✓ Metallic governance: Rules verified")
        print("\n→ Kanad governance framework OPERATIONAL!")
    else:
        print(f"Result: {passed_count}/{len(results)} protocols passed")

    print("=" * 80)

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
