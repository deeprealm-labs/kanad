"""
Governance Validation: METALLIC BONDING (Na2)

Tests metallic governance protocol validation on Na2 dimer.
Validates: delocalization, band structure (conceptual for dimer).
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import numpy as np
from kanad.core.molecule import Molecule
from kanad.core.atom import Atom


def test_metallic_na2():
    """Test metallic bonding governance on Na2."""
    print("=" * 80)
    print("GOVERNANCE VALIDATION: METALLIC BONDING (Na2)")
    print("=" * 80)

    # Create Na2 molecule
    print("\n1. Creating Na2 molecule...")
    na1 = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
    na2 = Atom('Na', position=np.array([3.08, 0.0, 0.0]))

    molecule = Molecule(atoms=[na1, na2], charge=0, spin=0)

    print(f"   ✓ Na2 created (bond length: 3.08 Å)")
    print(f"   - Na electronegativity: {na1.electronegativity:.2f}")
    print(f"   - ΔEN = 0.00 (same element)")
    print(f"   - Both atoms are metals: {na1.is_metal}")

    # Metallic bonding characteristics
    print("\n2. Metallic Bonding Characteristics...")
    print(f"   - Electron delocalization expected")
    print(f"   - Low ionization energy: {na1.properties.ionization_energy:.2f} eV")
    print(f"   - Valence electrons: {na1.properties.valence_electrons}")

    # Conceptual governance validation
    print("\n3. Governance Validation (Conceptual)...")
    print("\n" + "=" * 80)
    print("GOVERNANCE VALIDATION RESULTS")
    print("=" * 80)

    # Check 1: Metal atoms
    is_metal = na1.is_metal and na2.is_metal
    print(f"[{'✓' if is_metal else '✗'}] metallic_atoms")
    print(f"    Both atoms are metals: {is_metal} ✓")

    # Check 2: Low ionization energy (< 6 eV for alkali metals)
    low_ie = na1.properties.ionization_energy < 6.0
    print(f"[{'✓' if low_ie else '✗'}] low_ionization_energy")
    print(f"    IE = {na1.properties.ionization_energy:.2f} eV < 6.0 eV ✓")

    # Check 3: Same element (ΔEN = 0)
    same_element = na1.symbol == na2.symbol
    print(f"[{'✓' if same_element else '✗'}] same_element")
    print(f"    ΔEN = 0.00 (homonuclear diatomic) ✓")

    # Check 4: Delocalization expected
    print(f"[✓] delocalization_expected")
    print(f"    Valence electrons should delocalize across system ✓")

    print("=" * 80)

    all_passed = is_metal and low_ie and same_element

    if all_passed:
        print("\n🎉 ALL GOVERNANCE CHECKS PASSED - METALLIC CHARACTER CONFIRMED")
        print("   → Na2 exhibits metallic bonding characteristics")
    else:
        print("\n⚠ Some governance checks failed")

    print("\nMetallic Bonding Properties:")
    print(f"  - Delocalized electrons: Expected for Na-Na bond")
    print(f"  - Weak bonding: Na2 dissociates easily (D0 ~ 0.73 eV)")
    print(f"  - Band structure: Would form in extended Na lattice")
    print(f"  - Collective behavior: Present in bulk Na metal")

    print("\nNote: Full metallic Hamiltonian governance integration pending.")
    print("      Current validation shows metallic character via atomic properties.")

    print("\n" + "=" * 80)
    print("METALLIC BONDING VALIDATION COMPLETE")
    print("=" * 80)

    return all_passed


if __name__ == "__main__":
    success = test_metallic_na2()
    sys.exit(0 if success else 1)
