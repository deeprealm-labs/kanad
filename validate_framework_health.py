#!/usr/bin/env python3
"""
Kanad Framework Health Validator

Tests that the framework produces ACTUAL scientific results,
not just passes tests with mock/placeholder values.
"""

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.representations.second_quantization import SecondQuantizationRepresentation


def validate_covalent_bond():
    """Test H2 covalent bond produces real quantum operators."""
    print("\n" + "="*70)
    print("VALIDATING COVALENT BOND (H2)")
    print("="*70)

    # Create H2 bond
    bond = BondFactory.create_bond('H', 'H', distance=0.74, bond_type='covalent')

    # Test 1: Energy computation returns real value
    print("\n[1] Computing H2 energy...")
    result = bond.compute_energy(method='HF')
    energy = result['energy']
    print(f"    Energy: {energy:.6f} Ha")

    if energy == 0.0 or energy > 0:
        print("    ❌ FAIL: Energy is zero or positive (should be negative binding energy)")
        return False
    print("    ✅ PASS: Energy is negative and non-zero")

    # Test 2: Representation produces qubit operators
    print("\n[2] Testing LCAO to_qubit_operator()...")
    # LCAORepresentation expects hybridization string, not hamiltonian as second arg
    repr_lcao = LCAORepresentation(bond.molecule, hybridization='sp3')
    qubit_ops = repr_lcao.to_qubit_operator()

    if qubit_ops is None:
        print("    ❌ FAIL: to_qubit_operator() returned None (PLACEHOLDER!)")
        return False

    if not isinstance(qubit_ops, dict) or len(qubit_ops) == 0:
        print(f"    ❌ FAIL: Invalid qubit operators: {qubit_ops}")
        return False

    print(f"    Qubit operators: {len(qubit_ops)} Pauli terms")
    print(f"    Sample terms: {list(qubit_ops.items())[:3]}")

    # Check for actual Pauli strings
    has_pauli_strings = any(all(c in 'IXYZ' for c in key) for key in qubit_ops.keys())
    if not has_pauli_strings:
        print("    ❌ FAIL: No valid Pauli strings found")
        return False
    print("    ✅ PASS: Valid Pauli operators generated")

    # Test 3: Bonding/antibonding split
    print("\n[3] Testing bonding/antibonding split...")
    split = repr_lcao.get_bonding_antibonding_split(bond_idx=0)

    if split is None:
        print("    ❌ FAIL: get_bonding_antibonding_split() returned None")
        return False

    bonding_e = split['bonding_energy']
    antibonding_e = split['antibonding_energy']
    splitting = split['splitting']

    print(f"    Bonding energy: {bonding_e:.6f} Ha")
    print(f"    Antibonding energy: {antibonding_e:.6f} Ha")
    print(f"    Splitting: {splitting:.6f} Ha")

    # Check for placeholder values (-1.0, 1.0, 2.0)
    if bonding_e == -1.0 and antibonding_e == 1.0 and splitting == 2.0:
        print("    ❌ FAIL: Returned placeholder values (-1.0, 1.0, 2.0)")
        return False

    # Check physical validity
    if antibonding_e <= bonding_e:
        print("    ❌ FAIL: Antibonding should have higher energy than bonding")
        return False

    if abs(splitting - (antibonding_e - bonding_e)) > 1e-6:
        print("    ❌ FAIL: Splitting doesn't match energy difference")
        return False

    print("    ✅ PASS: Valid bonding/antibonding split computed")

    return True


def validate_ionic_bond():
    """Test NaCl ionic bond produces real quantum operators."""
    print("\n" + "="*70)
    print("VALIDATING IONIC BOND (NaCl)")
    print("="*70)

    # Create NaCl bond
    bond = BondFactory.create_bond('Na', 'Cl', distance=2.36, bond_type='ionic')

    # Test 1: Energy computation returns real value
    print("\n[1] Computing NaCl energy...")
    result = bond.compute_energy(method='exact')
    energy = result['energy']
    print(f"    Energy: {energy:.6f} Ha")

    if energy == 0.0 or energy > 0:
        print("    ❌ FAIL: Energy is zero or positive (should be negative binding energy)")
        return False
    print("    ✅ PASS: Energy is negative and non-zero")

    # Test 2: Representation produces qubit operators
    print("\n[2] Testing SecondQuantization to_qubit_operator()...")
    repr_sq = SecondQuantizationRepresentation(bond.molecule)
    qubit_ops = repr_sq.to_qubit_operator()

    if qubit_ops is None:
        print("    ❌ FAIL: to_qubit_operator() returned None (PLACEHOLDER!)")
        return False

    if not isinstance(qubit_ops, dict) or len(qubit_ops) == 0:
        print(f"    ❌ FAIL: Invalid qubit operators: {qubit_ops}")
        return False

    print(f"    Qubit operators: {len(qubit_ops)} Pauli terms")
    print(f"    Sample terms: {list(qubit_ops.items())[:3]}")

    # Check for actual Pauli strings
    has_pauli_strings = any(all(c in 'IXYZ' for c in key) for key in qubit_ops.keys())
    if not has_pauli_strings:
        print("    ❌ FAIL: No valid Pauli strings found")
        return False
    print("    ✅ PASS: Valid Pauli operators generated")

    # Test 3: Site occupation computation
    print("\n[3] Testing site occupation calculation...")
    ref_state = repr_sq.get_reference_state()

    # Test occupation on first site
    occupation = repr_sq._compute_site_occupation(ref_state, site=0)
    print(f"    Site 0 occupation: {occupation:.4f}")

    if occupation < 0 or occupation > 2:
        print("    ❌ FAIL: Occupation out of physical range [0, 2]")
        return False

    if occupation == 0.0:
        print("    ⚠️  WARNING: Zero occupation (might be intentional for reference state)")

    print("    ✅ PASS: Valid occupation computed")

    return True


def validate_metallic_bond():
    """Test Na metallic bond produces real band structure."""
    print("\n" + "="*70)
    print("VALIDATING METALLIC BOND (Na chain)")
    print("="*70)

    # Create Na metallic bond
    from kanad.bonds.metallic_bond import MetallicBond
    from kanad.core.atom import Atom

    atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(6)]
    bond = MetallicBond(atoms, lattice_type='1d_chain')

    # Test 1: Energy computation
    print("\n[1] Computing metallic bond energy...")
    result = bond.compute_energy()
    energy = result['energy']
    print(f"    Energy: {energy:.6f} eV")

    if energy == 0.0:
        print("    ❌ FAIL: Energy is exactly zero (likely placeholder)")
        return False

    if energy > 0:
        print("    ❌ FAIL: Energy is positive (should be negative cohesive energy)")
        return False

    print("    ✅ PASS: Energy is negative and non-zero")

    # Test 2: Band structure
    print("\n[2] Testing band structure...")
    band_energies = result.get('band_energies')

    if band_energies is None or len(band_energies) == 0:
        print("    ❌ FAIL: No band energies computed")
        return False

    print(f"    Number of bands: {len(band_energies)}")
    print(f"    Band energies (eV): {band_energies[:5]}...")

    # Check bands are not all identical
    if len(set(band_energies)) == 1:
        print("    ❌ FAIL: All band energies are identical (placeholder)")
        return False

    print("    ✅ PASS: Valid band structure computed")

    # Test 3: Analysis doesn't return placeholder string
    print("\n[3] Testing analysis...")
    analysis = bond.analyze()
    governance = analysis.get('governance_protocol', '')

    print(f"    Governance protocol: {governance}")

    if 'placeholder' in governance.lower():
        print("    ❌ FAIL: Governance protocol contains 'placeholder' text")
        return False

    print("    ✅ PASS: No placeholder text in analysis")

    return True


def main():
    """Run all validation tests."""
    print("\n" + "#"*70)
    print("# KANAD FRAMEWORK HEALTH VALIDATION")
    print("# Testing that framework produces REAL scientific results")
    print("#"*70)

    results = {}

    # Run all validations
    try:
        results['covalent'] = validate_covalent_bond()
    except Exception as e:
        print(f"\n❌ COVALENT BOND VALIDATION CRASHED: {e}")
        import traceback
        traceback.print_exc()
        results['covalent'] = False

    try:
        results['ionic'] = validate_ionic_bond()
    except Exception as e:
        print(f"\n❌ IONIC BOND VALIDATION CRASHED: {e}")
        import traceback
        traceback.print_exc()
        results['ionic'] = False

    try:
        results['metallic'] = validate_metallic_bond()
    except Exception as e:
        print(f"\n❌ METALLIC BOND VALIDATION CRASHED: {e}")
        import traceback
        traceback.print_exc()
        results['metallic'] = False

    # Print summary
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)

    for bond_type, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {bond_type.upper():15s}: {status}")

    all_passed = all(results.values())

    print("\n" + "="*70)
    if all_passed:
        print("🎉 ALL VALIDATIONS PASSED - Framework produces REAL results!")
        print("="*70)
        return 0
    else:
        print("⚠️  SOME VALIDATIONS FAILED - Framework has placeholders/issues!")
        print("="*70)
        return 1


if __name__ == '__main__':
    exit(main())
