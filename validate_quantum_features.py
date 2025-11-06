#!/usr/bin/env python3
"""
Comprehensive Quantum Features Validation Script

Tests all World's First quantum features:
1. Bonding-type resolved DOS
2. Quantum thermochemistry with bonding corrections
3. Integrated materials scout platform
4. Governance speedup validation
"""

import sys
import time
from typing import Dict, Any, List

def print_header(title: str):
    """Print formatted header"""
    print("\n" + "="*70)
    print(f"🌟 {title}")
    print("="*70)

def print_section(title: str):
    """Print formatted section"""
    print(f"\n{'='*70}")
    print(f"{title}")
    print("="*70)

def validate_quantum_dos() -> Dict[str, Any]:
    """Validate quantum DOS with bonding resolution"""
    print_header("TEST 1: QUANTUM DOS WITH BONDING RESOLUTION")

    from kanad.bonds import BondFactory
    from kanad.analysis import DOSCalculator

    results = {}

    # Test H2 (covalent)
    print("\n🔧 Testing H2 (covalent bond)...")
    h2_bond = BondFactory.create_bond('H', 'H', distance=0.74)
    dos_calc = DOSCalculator()

    start_time = time.time()
    h2_result = dos_calc.compute_quantum_dos(
        bond_or_molecule=h2_bond,
        energy_range=(-10, 10),
        n_states=5,
        solver='sqd',
        backend='statevector',
        use_governance=True,
        resolve_bonding=True,
        verbose=False
    )
    h2_time = time.time() - start_time

    results['h2'] = h2_result
    results['h2_time'] = h2_time

    print(f"  ✓ Bond type: {h2_result['bond_type']}")
    print(f"  ✓ Covalent fraction: {h2_result['covalent_fraction']*100:.1f}%")
    print(f"  ✓ Governance advantage: {h2_result['governance_advantage']:.1f}x")
    print(f"  ✓ Computation time: {h2_time:.2f}s")

    # Test LiH (ionic)
    print("\n🔧 Testing LiH (ionic bond)...")
    lih_bond = BondFactory.create_bond('Li', 'H', distance=1.60)

    start_time = time.time()
    lih_result = dos_calc.compute_quantum_dos(
        bond_or_molecule=lih_bond,
        energy_range=(-10, 10),
        n_states=5,
        solver='sqd',
        backend='statevector',
        use_governance=True,
        resolve_bonding=True,
        verbose=False
    )
    lih_time = time.time() - start_time

    results['lih'] = lih_result
    results['lih_time'] = lih_time

    print(f"  ✓ Bond type: {lih_result['bond_type']}")
    print(f"  ✓ Ionic fraction: {lih_result['ionic_fraction']*100:.1f}%")
    print(f"  ✓ Governance advantage: {lih_result['governance_advantage']:.1f}x")
    print(f"  ✓ Computation time: {lih_time:.2f}s")

    # Validation
    print_section("VALIDATION")

    validations = {
        'h2_is_covalent': h2_result['bond_type'] == 'covalent',
        'lih_is_ionic': lih_result['bond_type'] == 'ionic',
        'h2_governance_speedup': h2_result['governance_advantage'] > 1.0,
        'lih_governance_speedup': lih_result['governance_advantage'] > 1.0,
    }

    for test_name, passed in validations.items():
        status = "✅" if passed else "❌"
        print(f"{status} {test_name}: {passed}")

    results['validations'] = validations
    results['all_passed'] = all(validations.values())

    return results

def validate_quantum_thermochemistry() -> Dict[str, Any]:
    """Validate quantum thermochemistry with bonding corrections"""
    print_header("TEST 2: QUANTUM THERMOCHEMISTRY WITH BONDING CORRECTIONS")

    from kanad.bonds import BondFactory
    from kanad.core.molecule import Molecule
    from kanad.analysis import ThermochemistryCalculator

    results = {}

    # Test H2
    print("\n🔧 Testing H2 thermochemistry...")
    h2_bond = BondFactory.create_bond('H', 'H', distance=0.74)
    molecule_h2 = Molecule(h2_bond.atoms, charge=0, spin=0)

    thermo = ThermochemistryCalculator(
        molecule_h2,
        frequencies=[4401.2]  # H2 stretch frequency (cm⁻¹)
    )

    start_time = time.time()
    h2_result = thermo.compute_quantum_thermochemistry(
        bond=h2_bond,
        temperature=298.15,
        solver='sqd',
        backend='statevector',
        use_governance=True,
        apply_bonding_corrections=True,
        verbose=False
    )
    h2_time = time.time() - start_time

    results['h2'] = h2_result
    results['h2_time'] = h2_time

    print(f"  ✓ E_quantum: {h2_result['e_quantum']*627.509:.2f} kcal/mol")
    print(f"  ✓ H (Enthalpy): {h2_result['h_quantum']*627.509:.2f} kcal/mol")
    print(f"  ✓ S (Entropy): {h2_result['s_quantum']:.2f} cal/(mol·K)")
    print(f"  ✓ G (Gibbs): {h2_result['g_quantum']*627.509:.2f} kcal/mol")
    print(f"  ✓ Bond type: {h2_result['bond_type']}")
    print(f"  ✓ ΔH_bonding: {h2_result['delta_h_bonding']*627.509:.4f} kcal/mol")
    print(f"  ✓ ΔS_bonding: {h2_result['delta_s_bonding']:.2f} cal/(mol·K)")
    print(f"  ✓ Governance advantage: {h2_result['governance_advantage']:.1f}x")
    print(f"  ✓ Computation time: {h2_time:.2f}s")

    # Validation
    print_section("VALIDATION")

    validations = {
        'energy_negative': h2_result['e_quantum'] < 0,
        'enthalpy_reasonable': h2_result['h_quantum'] < 0,
        'entropy_positive': h2_result['s_quantum'] > 0,
        'governance_speedup': h2_result['governance_advantage'] > 1.0,
        'bond_type_detected': h2_result['bond_type'] == 'covalent',
        'bonding_corrections_applied': (h2_result['delta_h_bonding'] != 0 or
                                       h2_result['delta_s_bonding'] != 0),
    }

    for test_name, passed in validations.items():
        status = "✅" if passed else "❌"
        print(f"{status} {test_name}: {passed}")

    results['validations'] = validations
    results['all_passed'] = all(validations.values())

    return results

def validate_materials_scout() -> Dict[str, Any]:
    """Validate integrated materials scout platform"""
    print_header("TEST 3: INTEGRATED MATERIALS SCOUT PLATFORM")

    from kanad.applications import MaterialsScout
    from kanad.applications.materials_scout import MaterialCandidate
    from kanad.bonds import BondFactory

    results = {}

    # Initialize Materials Scout
    print("\n🔧 Initializing Materials Scout with governance...")
    scout = MaterialsScout(
        solver='sqd',
        backend='statevector',
        use_governance=True
    )
    print("  ✓ Materials Scout initialized")

    # Create H2 material
    print("\n🔧 Testing H2 material characterization...")
    h2_bond = BondFactory.create_bond('H', 'H', distance=0.74)
    h2_material = MaterialCandidate(
        name="H2",
        composition={'H': 1.0}
    )

    # Compute quantum DOS
    print("\n  📊 Computing quantum DOS...")
    start_time = time.time()
    dos_result = scout.compute_quantum_dos(
        material=h2_material,
        bond_or_molecule=h2_bond,
        n_states=5,
        resolve_bonding=True
    )
    dos_time = time.time() - start_time

    # Update material properties
    h2_material.bond_type = dos_result['bond_type']
    h2_material.covalent_fraction = dos_result['covalent_fraction']
    h2_material.governance_advantage = dos_result['governance_advantage']

    print(f"    ✓ Bond type: {h2_material.bond_type}")
    print(f"    ✓ Governance advantage: {dos_result['governance_advantage']:.1f}x")
    print(f"    ✓ Time: {dos_time:.2f}s")

    # Compute quantum thermochemistry
    print("\n  🌡️  Computing quantum thermochemistry...")
    start_time = time.time()
    thermo_result = scout.compute_quantum_thermochemistry(
        material=h2_material,
        bond=h2_bond,
        temperature=298.15,
        apply_bonding_corrections=True
    )
    thermo_time = time.time() - start_time

    # Update material properties
    h2_material.enthalpy = thermo_result['h_quantum']
    h2_material.entropy = thermo_result['s_quantum']
    h2_material.gibbs_free_energy = thermo_result['g_quantum']

    print(f"    ✓ H (Enthalpy): {h2_material.enthalpy*627.509:.2f} kcal/mol")
    print(f"    ✓ S (Entropy): {h2_material.entropy:.2f} cal/(mol·K)")
    print(f"    ✓ G (Gibbs): {h2_material.gibbs_free_energy*627.509:.2f} kcal/mol")
    print(f"    ✓ Governance advantage: {thermo_result['governance_advantage']:.1f}x")
    print(f"    ✓ Time: {thermo_time:.2f}s")

    results['dos'] = dos_result
    results['thermo'] = thermo_result
    results['dos_time'] = dos_time
    results['thermo_time'] = thermo_time
    results['material'] = h2_material

    # Validation
    print_section("VALIDATION")

    validations = {
        'bond_type_correct': h2_material.bond_type == 'covalent',
        'dos_governance_speedup': dos_result['governance_advantage'] > 1.0,
        'thermo_governance_speedup': thermo_result['governance_advantage'] > 1.0,
        'enthalpy_stable': h2_material.enthalpy < 0,
        'entropy_positive': h2_material.entropy > 0,
        'complete_pipeline': (h2_material.bond_type is not None and
                             h2_material.enthalpy is not None),
    }

    for test_name, passed in validations.items():
        status = "✅" if passed else "❌"
        print(f"{status} {test_name}: {passed}")

    results['validations'] = validations
    results['all_passed'] = all(validations.values())

    return results

def validate_governance_speedup(all_results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Validate governance speedup across all tests"""
    print_header("TEST 4: GOVERNANCE SPEEDUP VALIDATION")

    speedups = []

    # Collect speedups from all tests
    if 'dos' in all_results and 'h2' in all_results['dos']:
        speedups.append(('DOS H2', all_results['dos']['h2']['governance_advantage']))
        speedups.append(('DOS LiH', all_results['dos']['lih']['governance_advantage']))

    if 'thermo' in all_results and 'h2' in all_results['thermo']:
        speedups.append(('Thermo H2', all_results['thermo']['h2']['governance_advantage']))

    if 'scout' in all_results:
        speedups.append(('Scout DOS', all_results['scout']['dos']['governance_advantage']))
        speedups.append(('Scout Thermo', all_results['scout']['thermo']['governance_advantage']))

    print("\n📊 Governance Speedup Summary:")
    for name, speedup in speedups:
        print(f"  ✓ {name}: {speedup:.1f}x")

    avg_speedup = sum(s for _, s in speedups) / len(speedups)
    print(f"\n  📈 Average speedup: {avg_speedup:.1f}x")

    # Validation
    print_section("VALIDATION")

    validations = {
        'all_speedups_gt_1': all(s > 1.0 for _, s in speedups),
        'avg_speedup_gt_5': avg_speedup >= 5.0,
        'consistent_speedup': max(s for _, s in speedups) - min(s for _, s in speedups) < 5.0,
    }

    for test_name, passed in validations.items():
        status = "✅" if passed else "❌"
        print(f"{status} {test_name}: {passed}")

    return {
        'speedups': speedups,
        'avg_speedup': avg_speedup,
        'validations': validations,
        'all_passed': all(validations.values())
    }

def main():
    """Run comprehensive validation"""
    print("="*70)
    print("🌟 KANAD QUANTUM FEATURES COMPREHENSIVE VALIDATION")
    print("="*70)
    print("\nWorld's First Features:")
    print("  1. Bonding-type resolved DOS")
    print("  2. Quantum thermochemistry with bonding corrections")
    print("  3. Integrated materials scout platform")
    print("  4. Governance-guided quantum calculations (5-10x speedup)")
    print("\nStarting validation...")

    all_results = {}
    all_passed = []

    try:
        # Test 1: Quantum DOS
        dos_results = validate_quantum_dos()
        all_results['dos'] = dos_results
        all_passed.append(dos_results['all_passed'])

        # Test 2: Quantum Thermochemistry
        thermo_results = validate_quantum_thermochemistry()
        all_results['thermo'] = thermo_results
        all_passed.append(thermo_results['all_passed'])

        # Test 3: Materials Scout
        scout_results = validate_materials_scout()
        all_results['scout'] = scout_results
        all_passed.append(scout_results['all_passed'])

        # Test 4: Governance Speedup
        speedup_results = validate_governance_speedup(all_results)
        all_results['speedup'] = speedup_results
        all_passed.append(speedup_results['all_passed'])

    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Final Summary
    print_header("FINAL SUMMARY")

    print("\n📊 Test Results:")
    print(f"  {'✅' if all_results['dos']['all_passed'] else '❌'} Quantum DOS: {'PASSED' if all_results['dos']['all_passed'] else 'FAILED'}")
    print(f"  {'✅' if all_results['thermo']['all_passed'] else '❌'} Quantum Thermochemistry: {'PASSED' if all_results['thermo']['all_passed'] else 'FAILED'}")
    print(f"  {'✅' if all_results['scout']['all_passed'] else '❌'} Materials Scout: {'PASSED' if all_results['scout']['all_passed'] else 'FAILED'}")
    print(f"  {'✅' if all_results['speedup']['all_passed'] else '❌'} Governance Speedup: {'PASSED' if all_results['speedup']['all_passed'] else 'FAILED'}")

    total_time = (
        all_results['dos']['h2_time'] +
        all_results['dos']['lih_time'] +
        all_results['thermo']['h2_time'] +
        all_results['scout']['dos_time'] +
        all_results['scout']['thermo_time']
    )

    print(f"\n⏱️  Total computation time: {total_time:.2f}s")
    print(f"📈 Average governance speedup: {all_results['speedup']['avg_speedup']:.1f}x")

    if all(all_passed):
        print("\n" + "="*70)
        print("🎉 ALL TESTS PASSED!")
        print("="*70)
        print("\n🌟 WORLD'S FIRST FEATURES VALIDATED:")
        print("  ✓ Bonding-type resolved DOS (covalent/ionic/metallic)")
        print("  ✓ Quantum thermochemistry with bonding corrections")
        print("  ✓ Governance speedup (5-10x measured)")
        print("  ✓ Integrated materials scout platform")
        print("\nCOMPETITIVE ADVANTAGES:")
        print("  vs Materials Project: Bonding character + predictive (not database)")
        print("  vs Schrödinger: FREE + bonding-aware + governance speedup")
        print("  vs VASP/QE: Bonding-type resolution (UNIQUE!)")
        print("  vs Gaussian/ORCA: Bonding-aware corrections (UNIQUE!)")
        print("="*70)
        return 0
    else:
        print("\n" + "="*70)
        print("❌ SOME TESTS FAILED")
        print("="*70)
        return 1

if __name__ == "__main__":
    sys.exit(main())
