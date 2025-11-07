"""
KANAD FRAMEWORK - Comprehensive Validation Test Suite
=====================================================

This script tests all major components of the Kanad framework:
1. Quantum Solvers (VQE, SQD, Krylov-SQD)
2. Application Modules (Drug Discovery, Materials Scout)
3. Analysis Modules (Properties, Spectroscopy, Thermochemistry)
4. Molecular Dynamics

Author: Framework Validation
Date: 2025-11-06
"""

import sys
import time
from pathlib import Path

# Add kanad to path
sys.path.insert(0, str(Path(__file__).parent))

from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.bonds import BondFactory

print("=" * 80)
print(" KANAD FRAMEWORK - COMPREHENSIVE VALIDATION TEST SUITE")
print("=" * 80)
print()

test_results = []
total_tests = 0
passed_tests = 0


def test_section(name):
    """Print section header"""
    print(f"\n{'=' * 80}")
    print(f" {name}")
    print(f"{'=' * 80}\n")


def test_result(name, passed, details=""):
    """Record and print test result"""
    global total_tests, passed_tests
    total_tests += 1
    if passed:
        passed_tests += 1

    status = "✓ PASSED" if passed else "✗ FAILED"
    print(f"{status}: {name}")
    if details:
        print(f"   {details}")

    test_results.append({
        'name': name,
        'passed': passed,
        'details': details
    })


# =============================================================================
# 1. QUANTUM SOLVERS TESTING
# =============================================================================

test_section("SECTION 1: QUANTUM SOLVERS")

# Test 1.1: VQE Solver with H2
print("[TEST 1.1] VQE Solver - H₂ Molecule")
try:
    start_time = time.time()

    from kanad.solvers.vqe_solver import VQESolver

    h1 = Atom('H', position=(0.0, 0.0, 0.0))
    h2 = Atom('H', position=(0.74, 0.0, 0.0))
    bond = CovalentBond(h1, h2, basis='sto-3g')

    solver = VQESolver(
        bond=bond,
        ansatz_type='hardware_efficient',
        mapper_type='jordan_wigner',
        optimizer='SLSQP',
        max_iterations=50
    )

    result = solver.solve()
    energy = result['energy']
    reference = -1.137284
    error_mHa = abs(energy - reference) * 1000

    elapsed = time.time() - start_time

    passed = error_mHa < 100  # 100 mHa tolerance
    test_result(
        "VQE Solver (H₂)",
        passed,
        f"Energy: {energy:.6f} Ha, Error: {error_mHa:.3f} mHa, Time: {elapsed:.1f}s"
    )

except Exception as e:
    test_result("VQE Solver (H₂)", False, f"Error: {str(e)[:100]}")


# Test 1.2: SQD Solver for Excited States
print("\n[TEST 1.2] SQD Solver - Excited States")
try:
    start_time = time.time()

    from kanad.solvers.sqd_solver import SQDSolver

    h1 = Atom('H', position=(0.0, 0.0, 0.0))
    h2 = Atom('H', position=(0.74, 0.0, 0.0))
    bond = CovalentBond(h1, h2, basis='sto-3g')

    solver = SQDSolver(
        bond=bond,
        mapper_type='jordan_wigner',
        num_states=3
    )

    result = solver.solve()
    energies = result.get('energies', [])

    elapsed = time.time() - start_time

    passed = len(energies) >= 1
    test_result(
        "SQD Solver (Excited States)",
        passed,
        f"Ground State: {energies[0]:.6f} Ha, States: {len(energies)}, Time: {elapsed:.1f}s"
    )

except Exception as e:
    test_result("SQD Solver (Excited States)", False, f"Error: {str(e)[:100]}")


# Test 1.3: Krylov-SQD Solver
print("\n[TEST 1.3] Krylov-SQD Solver")
try:
    start_time = time.time()

    from kanad.solvers.krylov_sqd_solver import KrylovSQDSolver

    h1 = Atom('H', position=(0.0, 0.0, 0.0))
    h2 = Atom('H', position=(0.74, 0.0, 0.0))
    bond = CovalentBond(h1, h2, basis='sto-3g')

    solver = KrylovSQDSolver(
        bond=bond,
        mapper_type='jordan_wigner',
        num_states=2,
        krylov_dimension=5
    )

    result = solver.solve()
    energies = result.get('energies', [])

    elapsed = time.time() - start_time

    passed = len(energies) >= 1
    test_result(
        "Krylov-SQD Solver",
        passed,
        f"Ground State: {energies[0]:.6f} Ha, States: {len(energies)}, Time: {elapsed:.1f}s"
    )

except Exception as e:
    test_result("Krylov-SQD Solver", False, f"Error: {str(e)[:100]}")


# =============================================================================
# 2. APPLICATION MODULES TESTING
# =============================================================================

test_section("SECTION 2: APPLICATION MODULES")

# Test 2.1: Drug Discovery Platform
print("[TEST 2.1] Drug Discovery Platform - Aspirin Analysis")
try:
    start_time = time.time()

    from kanad.applications import DrugDiscoveryPlatform

    platform = DrugDiscoveryPlatform()
    aspirin_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'

    result = platform.analyze_drug_candidate(
        smiles=aspirin_smiles,
        use_quantum=False  # Classical for speed
    )

    elapsed = time.time() - start_time

    has_adme = 'adme' in result
    has_lipinski = 'lipinski_violations' in result

    passed = has_adme and has_lipinski
    druglike = result.get('druglike', 'N/A')
    mw = result.get('adme', {}).get('molecular_weight', 'N/A')

    test_result(
        "Drug Discovery Platform",
        passed,
        f"Druglike: {druglike}, MW: {mw}, Time: {elapsed:.1f}s"
    )

except Exception as e:
    test_result("Drug Discovery Platform", False, f"Error: {str(e)[:100]}")


# Test 2.2: Materials Scout
print("\n[TEST 2.2] Materials Scout - Band Gap Calculation")
try:
    start_time = time.time()

    from kanad.applications import MaterialsScout

    scout = MaterialsScout()
    benzene_smiles = 'c1ccccc1'

    result = scout.calculate_band_gap(
        smiles=benzene_smiles,
        use_quantum=False
    )

    elapsed = time.time() - start_time

    has_band_gap = 'band_gap' in result
    band_gap = result.get('band_gap', 'N/A')

    passed = has_band_gap
    test_result(
        "Materials Scout",
        passed,
        f"Band Gap: {band_gap} eV, Time: {elapsed:.1f}s"
    )

except Exception as e:
    test_result("Materials Scout", False, f"Error: {str(e)[:100]}")


# =============================================================================
# 3. ANALYSIS MODULES TESTING
# =============================================================================

test_section("SECTION 3: ANALYSIS MODULES")

# Test 3.1: Property Calculator - Dipole Moment
print("[TEST 3.1] Property Calculator - Dipole Moment")
try:
    start_time = time.time()

    from kanad.analysis import PropertyCalculator

    # HF molecule has known dipole
    bond = BondFactory.create_bond('H', 'F', distance=0.917)
    calculator = PropertyCalculator(bond)

    dipole = calculator.calculate_dipole_moment(use_quantum=False)

    elapsed = time.time() - start_time

    passed = dipole is not None and abs(dipole) > 0
    test_result(
        "Property Calculator (Dipole)",
        passed,
        f"Dipole: {dipole:.4f} D, Time: {elapsed:.1f}s" if dipole else "No dipole calculated"
    )

except Exception as e:
    test_result("Property Calculator (Dipole)", False, f"Error: {str(e)[:100]}")


# Test 3.2: UV-Vis Spectroscopy
print("\n[TEST 3.2] UV-Vis Spectroscopy Calculator")
try:
    start_time = time.time()

    from kanad.analysis import UVVisCalculator

    bond = BondFactory.create_from_smiles('c1ccccc1')  # Benzene
    calculator = UVVisCalculator(bond)

    result = calculator.calculate_spectrum(
        use_quantum=False,
        num_states=5
    )

    elapsed = time.time() - start_time

    has_spectrum = 'wavelengths' in result and len(result['wavelengths']) > 0
    num_transitions = len(result.get('wavelengths', []))

    passed = has_spectrum
    test_result(
        "UV-Vis Spectroscopy",
        passed,
        f"Transitions: {num_transitions}, Time: {elapsed:.1f}s"
    )

except Exception as e:
    test_result("UV-Vis Spectroscopy", False, f"Error: {str(e)[:100]}")


# Test 3.3: Thermochemistry Calculator
print("\n[TEST 3.3] Thermochemistry Calculator")
try:
    start_time = time.time()

    from kanad.analysis import ThermochemistryCalculator

    bond = BondFactory.create_bond('H', 'H', distance=0.74)
    calculator = ThermochemistryCalculator(bond)

    result = calculator.calculate_thermodynamics(
        temperature=298.15,
        pressure=101325
    )

    elapsed = time.time() - start_time

    has_thermo = 'enthalpy' in result and 'entropy' in result

    passed = has_thermo
    test_result(
        "Thermochemistry Calculator",
        passed,
        f"H: {result.get('enthalpy', 'N/A')}, S: {result.get('entropy', 'N/A')}, Time: {elapsed:.1f}s"
    )

except Exception as e:
    test_result("Thermochemistry Calculator", False, f"Error: {str(e)[:100]}")


# =============================================================================
# 4. MOLECULAR DYNAMICS TESTING
# =============================================================================

test_section("SECTION 4: MOLECULAR DYNAMICS")

# Test 4.1: Classical MD
print("[TEST 4.1] Classical Molecular Dynamics")
try:
    start_time = time.time()

    from kanad.dynamics import MDSimulator

    bond = BondFactory.create_bond('H', 'H', distance=0.74)

    simulator = MDSimulator(bond, mode='classical')

    result = simulator.run(
        num_steps=50,  # Short run for testing
        timestep=0.5,
        temperature=300
    )

    elapsed = time.time() - start_time

    has_trajectory = 'trajectory' in result

    passed = has_trajectory
    test_result(
        "Classical Molecular Dynamics",
        passed,
        f"Steps: {result.get('num_steps', 'N/A')}, Time: {elapsed:.1f}s"
    )

except Exception as e:
    test_result("Classical Molecular Dynamics", False, f"Error: {str(e)[:100]}")


# =============================================================================
# SUMMARY AND FINAL REPORT
# =============================================================================

test_section("VALIDATION SUMMARY")

success_rate = (passed_tests / total_tests * 100) if total_tests > 0 else 0

print(f"Total Tests: {total_tests}")
print(f"Passed: {passed_tests}")
print(f"Failed: {total_tests - passed_tests}")
print(f"Success Rate: {success_rate:.1f}%\n")

print("Detailed Results:")
print("-" * 80)
for result in test_results:
    status = "✓" if result['passed'] else "✗"
    print(f"{status} {result['name']}")
    if result['details']:
        print(f"  {result['details']}")

print("\n" + "=" * 80)
if success_rate >= 80:
    print(" ✓✓✓ FRAMEWORK VALIDATION SUCCESSFUL ✓✓✓")
    print(f" {passed_tests}/{total_tests} tests passed ({success_rate:.1f}%)")
elif success_rate >= 50:
    print(" ⚠ FRAMEWORK PARTIALLY VALIDATED ⚠")
    print(f" {passed_tests}/{total_tests} tests passed ({success_rate:.1f}%)")
else:
    print(" ✗✗✗ FRAMEWORK VALIDATION FAILED ✗✗✗")
    print(f" Only {passed_tests}/{total_tests} tests passed ({success_rate:.1f}%)")
print("=" * 80)

# Exit with appropriate code
sys.exit(0 if success_rate >= 80 else 1)
