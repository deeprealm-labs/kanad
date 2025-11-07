"""
Comprehensive Validation Test Suite for Kanad Framework
========================================================

This test suite validates all major components:
1. Solvers (VQE, Hi-VQE, SQD, Krylov-SQD, Active Space)
2. Application Modules (Drug Discovery, Materials Scout, Catalyst Optimizer, Alloy Designer)
3. Analysis Modules (Properties, Spectroscopy, DOS, Thermochemistry)
4. Molecular Dynamics (Classical and Quantum MD)

Author: Comprehensive Framework Testing
Date: 2025-11-06
"""

import numpy as np
import time
from typing import Dict, List, Any
import json

# Import Kanad framework components
from kanad.bonds import BondFactory
from kanad.solvers import VQESolver, SQDSolver, KrylovSQDSolver, ExcitedStatesSolver
from kanad.applications import (
    DrugDiscoveryPlatform,
    MaterialsScout,
    CatalystOptimizer,
    AlloyDesigner
)
from kanad.analysis import (
    PropertyCalculator,
    UVVisCalculator,
    NMRCalculator,
    RamanIRCalculator,
    DOSCalculator,
    ThermochemistryCalculator,
    BondLengthScanner
)
from kanad.dynamics import MDSimulator


class TestResults:
    """Collect and organize all test results"""

    def __init__(self):
        self.results = {
            'solvers': {},
            'applications': {},
            'analysis': {},
            'dynamics': {},
            'summary': {}
        }
        self.start_time = time.time()

    def add_result(self, category: str, test_name: str, result: Dict[str, Any]):
        """Add a test result"""
        if category not in self.results:
            self.results[category] = {}
        self.results[category][test_name] = result

    def get_summary(self) -> Dict[str, Any]:
        """Generate summary statistics"""
        total_tests = sum(len(v) for v in self.results.values() if isinstance(v, dict))
        passed_tests = sum(
            1 for category in ['solvers', 'applications', 'analysis', 'dynamics']
            for test in self.results[category].values()
            if test.get('status') == 'PASSED'
        )

        return {
            'total_tests': total_tests - 1,  # Exclude 'summary' key
            'passed_tests': passed_tests,
            'failed_tests': total_tests - 1 - passed_tests,
            'total_time': time.time() - self.start_time
        }

    def print_report(self):
        """Print comprehensive test report"""
        print("\n" + "="*80)
        print(" KANAD FRAMEWORK COMPREHENSIVE VALIDATION REPORT")
        print("="*80 + "\n")

        # Solvers section
        print("1. QUANTUM SOLVERS VALIDATION")
        print("-" * 80)
        for test_name, result in self.results['solvers'].items():
            status = "✓ PASSED" if result['status'] == 'PASSED' else "✗ FAILED"
            print(f"  {status}: {test_name}")
            if 'energy' in result:
                print(f"    Energy: {result['energy']:.6f} Ha")
            if 'accuracy' in result:
                print(f"    Accuracy: {result['accuracy']}")
            if 'time' in result:
                print(f"    Time: {result['time']:.2f}s")
            print()

        # Applications section
        print("\n2. APPLICATION MODULES VALIDATION")
        print("-" * 80)
        for test_name, result in self.results['applications'].items():
            status = "✓ PASSED" if result['status'] == 'PASSED' else "✗ FAILED"
            print(f"  {status}: {test_name}")
            if 'metrics' in result:
                for key, value in result['metrics'].items():
                    print(f"    {key}: {value}")
            print()

        # Analysis section
        print("\n3. ANALYSIS MODULES VALIDATION")
        print("-" * 80)
        for test_name, result in self.results['analysis'].items():
            status = "✓ PASSED" if result['status'] == 'PASSED' else "✗ FAILED"
            print(f"  {status}: {test_name}")
            if 'properties' in result:
                for key, value in result['properties'].items():
                    if isinstance(value, float):
                        print(f"    {key}: {value:.6f}")
                    else:
                        print(f"    {key}: {value}")
            print()

        # Dynamics section
        print("\n4. MOLECULAR DYNAMICS VALIDATION")
        print("-" * 80)
        for test_name, result in self.results['dynamics'].items():
            status = "✓ PASSED" if result['status'] == 'PASSED' else "✗ FAILED"
            print(f"  {status}: {test_name}")
            if 'steps' in result:
                print(f"    Steps: {result['steps']}")
            if 'time' in result:
                print(f"    Time: {result['time']:.2f}s")
            print()

        # Summary
        summary = self.get_summary()
        print("\n" + "="*80)
        print(" SUMMARY")
        print("="*80)
        print(f"  Total Tests: {summary['total_tests']}")
        print(f"  Passed: {summary['passed_tests']} ({summary['passed_tests']/summary['total_tests']*100:.1f}%)")
        print(f"  Failed: {summary['failed_tests']} ({summary['failed_tests']/summary['total_tests']*100:.1f}%)")
        print(f"  Total Time: {summary['total_time']:.2f}s")
        print("="*80 + "\n")


# Initialize test results collector
test_results = TestResults()


# ============================================================================
# 1. SOLVER VALIDATION TESTS
# ============================================================================

def test_vqe_standard_h2():
    """Test standard VQE on H2 molecule"""
    print("\n[TEST] Standard VQE on H₂...")

    try:
        # Create H2 molecule
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        # Standard VQE solver
        solver = VQESolver(
            bond=bond,
            ansatz_type='ucc_singles_doubles',
            mapper_type='jordan_wigner',
            optimizer='SLSQP',
            backend='statevector',
            max_iterations=100
        )

        start_time = time.time()
        result = solver.solve()
        elapsed_time = time.time() - start_time

        # Expected energy: -1.137 Ha (literature value)
        energy = result['energy']
        expected = -1.137
        error = abs(energy - expected)

        test_results.add_result('solvers', 'VQE_Standard_H2', {
            'status': 'PASSED' if error < 0.01 else 'FAILED',
            'energy': energy,
            'expected': expected,
            'error': error,
            'accuracy': f'{error*1000:.3f} mHa',
            'time': elapsed_time,
            'converged': result.get('converged', True)
        })

        print(f"  Energy: {energy:.6f} Ha (Expected: {expected:.6f} Ha)")
        print(f"  Error: {error*1000:.3f} mHa")
        print(f"  Status: {'PASSED ✓' if error < 0.01 else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('solvers', 'VQE_Standard_H2', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_hivqe_h2():
    """Test Hi-VQE mode on H2 molecule"""
    print("\n[TEST] Hi-VQE Mode on H₂...")

    try:
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        # Hi-VQE solver
        solver = VQESolver(
            bond=bond,
            mode='hivqe',
            ansatz_type='hardware_efficient',
            mapper_type='jordan_wigner',
            optimizer='SLSQP',
            backend='statevector',
            hivqe_max_iterations=10,
            hivqe_subspace_threshold=0.05
        )

        start_time = time.time()
        result = solver.solve()
        elapsed_time = time.time() - start_time

        energy = result['energy']
        expected = -1.137
        error = abs(energy - expected)

        # Check for measurement reduction stats
        measurement_reduction = result.get('hivqe_stats', {}).get('measurement_reduction', 'N/A')

        test_results.add_result('solvers', 'HiVQE_H2', {
            'status': 'PASSED' if error < 0.01 else 'FAILED',
            'energy': energy,
            'expected': expected,
            'error': error,
            'accuracy': f'{error*1000:.3f} mHa',
            'time': elapsed_time,
            'measurement_reduction': measurement_reduction
        })

        print(f"  Energy: {energy:.6f} Ha (Expected: {expected:.6f} Ha)")
        print(f"  Error: {error*1000:.3f} mHa")
        print(f"  Measurement Reduction: {measurement_reduction}")
        print(f"  Status: {'PASSED ✓' if error < 0.01 else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('solvers', 'HiVQE_H2', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_vqe_lih():
    """Test VQE on LiH molecule"""
    print("\n[TEST] VQE on LiH...")

    try:
        bond = BondFactory.create_bond('Li', 'H', distance=1.5949)

        solver = VQESolver(
            bond=bond,
            ansatz_type='ucc_singles_doubles',
            mapper_type='jordan_wigner',
            optimizer='SLSQP',
            backend='statevector',
            max_iterations=100
        )

        start_time = time.time()
        result = solver.solve()
        elapsed_time = time.time() - start_time

        energy = result['energy']
        # Expected: -7.882 Ha (approximate)
        expected = -7.882
        error = abs(energy - expected)

        test_results.add_result('solvers', 'VQE_LiH', {
            'status': 'PASSED' if error < 0.05 else 'FAILED',
            'energy': energy,
            'expected': expected,
            'error': error,
            'accuracy': f'{error*1000:.3f} mHa',
            'time': elapsed_time
        })

        print(f"  Energy: {energy:.6f} Ha (Expected: {expected:.6f} Ha)")
        print(f"  Error: {error*1000:.3f} mHa")
        print(f"  Status: {'PASSED ✓' if error < 0.05 else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('solvers', 'VQE_LiH', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_active_space_solver():
    """Test Active Space reduction"""
    print("\n[TEST] Active Space Solver...")

    try:
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        # Regular VQE
        solver_regular = VQESolver(
            bond=bond,
            ansatz_type='ucc_singles_doubles',
            backend='statevector'
        )
        result_regular = solver_regular.solve()

        # Active space VQE
        solver_active = VQESolver(
            bond=bond,
            ansatz_type='ucc_singles_doubles',
            backend='statevector',
            use_active_space=True
        )
        result_active = solver_active.solve()

        energy_diff = abs(result_regular['energy'] - result_active['energy'])

        test_results.add_result('solvers', 'ActiveSpace_H2', {
            'status': 'PASSED' if energy_diff < 0.001 else 'FAILED',
            'energy_regular': result_regular['energy'],
            'energy_active': result_active['energy'],
            'energy_difference': energy_diff,
            'accuracy': f'{energy_diff*1000:.3f} mHa'
        })

        print(f"  Regular Energy: {result_regular['energy']:.6f} Ha")
        print(f"  Active Space Energy: {result_active['energy']:.6f} Ha")
        print(f"  Difference: {energy_diff*1000:.3f} mHa")
        print(f"  Status: {'PASSED ✓' if energy_diff < 0.001 else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('solvers', 'ActiveSpace_H2', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_sqd_solver():
    """Test SQD solver for excited states"""
    print("\n[TEST] SQD Solver for Excited States...")

    try:
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        solver = SQDSolver(
            bond=bond,
            mapper_type='jordan_wigner',
            backend='statevector',
            num_states=3  # Ground + 2 excited states
        )

        start_time = time.time()
        result = solver.solve()
        elapsed_time = time.time() - start_time

        energies = result.get('energies', [])

        test_results.add_result('solvers', 'SQD_H2_ExcitedStates', {
            'status': 'PASSED' if len(energies) >= 1 else 'FAILED',
            'ground_state': energies[0] if energies else None,
            'excited_states': energies[1:] if len(energies) > 1 else [],
            'num_states': len(energies),
            'time': elapsed_time
        })

        print(f"  Ground State: {energies[0]:.6f} Ha" if energies else "  No energies")
        if len(energies) > 1:
            for i, e in enumerate(energies[1:], 1):
                print(f"  Excited State {i}: {e:.6f} Ha")
        print(f"  Status: {'PASSED ✓' if len(energies) >= 1 else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('solvers', 'SQD_H2_ExcitedStates', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_krylov_sqd_solver():
    """Test Krylov-SQD solver"""
    print("\n[TEST] Krylov-SQD Solver...")

    try:
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        solver = KrylovSQDSolver(
            bond=bond,
            mapper_type='jordan_wigner',
            backend='statevector',
            num_states=2,
            krylov_dimension=5
        )

        start_time = time.time()
        result = solver.solve()
        elapsed_time = time.time() - start_time

        energies = result.get('energies', [])

        test_results.add_result('solvers', 'KrylovSQD_H2', {
            'status': 'PASSED' if len(energies) >= 1 else 'FAILED',
            'ground_state': energies[0] if energies else None,
            'excited_states': energies[1:] if len(energies) > 1 else [],
            'num_states': len(energies),
            'time': elapsed_time
        })

        print(f"  Ground State: {energies[0]:.6f} Ha" if energies else "  No energies")
        if len(energies) > 1:
            for i, e in enumerate(energies[1:], 1):
                print(f"  Excited State {i}: {e:.6f} Ha")
        print(f"  Status: {'PASSED ✓' if len(energies) >= 1 else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('solvers', 'KrylovSQD_H2', {
            'status': 'FAILED',
            'error': str(e)
        })


# ============================================================================
# 2. APPLICATION MODULE VALIDATION TESTS
# ============================================================================

def test_drug_discovery_platform():
    """Test Drug Discovery Platform with Aspirin"""
    print("\n[TEST] Drug Discovery Platform - Aspirin...")

    try:
        platform = DrugDiscoveryPlatform()

        # Aspirin SMILES
        aspirin_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'

        start_time = time.time()
        result = platform.analyze_drug_candidate(
            smiles=aspirin_smiles,
            use_quantum=False  # Use classical for speed
        )
        elapsed_time = time.time() - start_time

        # Check key metrics
        has_adme = 'adme' in result
        has_lipinski = 'lipinski_violations' in result
        is_druglike = result.get('druglike', False)

        test_results.add_result('applications', 'DrugDiscovery_Aspirin', {
            'status': 'PASSED' if has_adme and has_lipinski else 'FAILED',
            'metrics': {
                'druglike': is_druglike,
                'lipinski_violations': result.get('lipinski_violations', 'N/A'),
                'molecular_weight': result.get('adme', {}).get('molecular_weight', 'N/A'),
                'logP': result.get('adme', {}).get('logP', 'N/A'),
            },
            'time': elapsed_time
        })

        print(f"  Druglike: {is_druglike}")
        print(f"  Lipinski Violations: {result.get('lipinski_violations', 'N/A')}")
        print(f"  Molecular Weight: {result.get('adme', {}).get('molecular_weight', 'N/A')}")
        print(f"  Status: {'PASSED ✓' if has_adme and has_lipinski else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('applications', 'DrugDiscovery_Aspirin', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_materials_scout():
    """Test Materials Scout for band gap calculation"""
    print("\n[TEST] Materials Scout - Band Gap Calculation...")

    try:
        scout = MaterialsScout()

        # Simple molecule for testing (e.g., benzene)
        benzene_smiles = 'c1ccccc1'

        start_time = time.time()
        result = scout.calculate_band_gap(
            smiles=benzene_smiles,
            use_quantum=False
        )
        elapsed_time = time.time() - start_time

        has_band_gap = 'band_gap' in result
        band_gap = result.get('band_gap', 'N/A')

        test_results.add_result('applications', 'MaterialsScout_BandGap', {
            'status': 'PASSED' if has_band_gap else 'FAILED',
            'metrics': {
                'band_gap_eV': band_gap,
                'homo': result.get('homo', 'N/A'),
                'lumo': result.get('lumo', 'N/A'),
            },
            'time': elapsed_time
        })

        print(f"  Band Gap: {band_gap} eV")
        print(f"  HOMO: {result.get('homo', 'N/A')} eV")
        print(f"  LUMO: {result.get('lumo', 'N/A')} eV")
        print(f"  Status: {'PASSED ✓' if has_band_gap else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('applications', 'MaterialsScout_BandGap', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_catalyst_optimizer():
    """Test Catalyst Optimizer"""
    print("\n[TEST] Catalyst Optimizer - Reaction Barrier...")

    try:
        optimizer = CatalystOptimizer()

        # Simple reaction (H2 dissociation)
        reactant = BondFactory.create_bond('H', 'H', distance=0.74)
        product = BondFactory.create_bond('H', 'H', distance=2.0)

        start_time = time.time()
        result = optimizer.calculate_reaction_barrier(
            reactant=reactant,
            product=product,
            use_quantum=False
        )
        elapsed_time = time.time() - start_time

        has_barrier = 'activation_energy' in result

        test_results.add_result('applications', 'CatalystOptimizer_Barrier', {
            'status': 'PASSED' if has_barrier else 'FAILED',
            'metrics': {
                'activation_energy': result.get('activation_energy', 'N/A'),
                'reaction_energy': result.get('reaction_energy', 'N/A'),
            },
            'time': elapsed_time
        })

        print(f"  Activation Energy: {result.get('activation_energy', 'N/A')} Ha")
        print(f"  Reaction Energy: {result.get('reaction_energy', 'N/A')} Ha")
        print(f"  Status: {'PASSED ✓' if has_barrier else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('applications', 'CatalystOptimizer_Barrier', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_alloy_designer():
    """Test Alloy Designer"""
    print("\n[TEST] Alloy Designer - Composition Optimization...")

    try:
        designer = AlloyDesigner()

        # Simple binary alloy (e.g., Fe-Ni)
        start_time = time.time()
        result = designer.optimize_composition(
            elements=['Fe', 'Ni'],
            target_property='strength',
            use_quantum=False
        )
        elapsed_time = time.time() - start_time

        has_composition = 'optimal_composition' in result

        test_results.add_result('applications', 'AlloyDesigner_Composition', {
            'status': 'PASSED' if has_composition else 'FAILED',
            'metrics': {
                'optimal_composition': result.get('optimal_composition', 'N/A'),
                'predicted_property': result.get('predicted_property', 'N/A'),
            },
            'time': elapsed_time
        })

        print(f"  Optimal Composition: {result.get('optimal_composition', 'N/A')}")
        print(f"  Predicted Property: {result.get('predicted_property', 'N/A')}")
        print(f"  Status: {'PASSED ✓' if has_composition else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('applications', 'AlloyDesigner_Composition', {
            'status': 'FAILED',
            'error': str(e)
        })


# ============================================================================
# 3. ANALYSIS MODULE VALIDATION TESTS
# ============================================================================

def test_property_calculator():
    """Test Property Calculator for molecular properties"""
    print("\n[TEST] Property Calculator - Dipole Moment...")

    try:
        bond = BondFactory.create_bond('H', 'F', distance=0.917)

        calculator = PropertyCalculator(bond)

        start_time = time.time()
        dipole = calculator.calculate_dipole_moment(use_quantum=False)
        elapsed_time = time.time() - start_time

        # HF has a known dipole moment (~1.82 Debye)
        has_dipole = dipole is not None and abs(dipole) > 0

        test_results.add_result('analysis', 'PropertyCalculator_Dipole', {
            'status': 'PASSED' if has_dipole else 'FAILED',
            'properties': {
                'dipole_moment_D': dipole if dipole else 'N/A',
                'expected': '~1.82 D',
            },
            'time': elapsed_time
        })

        print(f"  Dipole Moment: {dipole:.4f} D" if dipole else "  Dipole: N/A")
        print(f"  Expected: ~1.82 D")
        print(f"  Status: {'PASSED ✓' if has_dipole else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('analysis', 'PropertyCalculator_Dipole', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_uvvis_calculator():
    """Test UV-Vis Spectroscopy Calculator"""
    print("\n[TEST] UV-Vis Spectroscopy Calculator...")

    try:
        # Benzene for UV-Vis
        bond = BondFactory.create_from_smiles('c1ccccc1')

        calculator = UVVisCalculator(bond)

        start_time = time.time()
        result = calculator.calculate_spectrum(
            use_quantum=False,
            num_states=5
        )
        elapsed_time = time.time() - start_time

        has_spectrum = 'wavelengths' in result and len(result['wavelengths']) > 0

        test_results.add_result('analysis', 'UVVis_Benzene', {
            'status': 'PASSED' if has_spectrum else 'FAILED',
            'properties': {
                'num_transitions': len(result.get('wavelengths', [])),
                'max_wavelength': max(result.get('wavelengths', [0])) if result.get('wavelengths') else 'N/A',
                'min_wavelength': min(result.get('wavelengths', [0])) if result.get('wavelengths') else 'N/A',
            },
            'time': elapsed_time
        })

        print(f"  Number of Transitions: {len(result.get('wavelengths', []))}")
        if result.get('wavelengths'):
            print(f"  Wavelength Range: {min(result['wavelengths']):.1f} - {max(result['wavelengths']):.1f} nm")
        print(f"  Status: {'PASSED ✓' if has_spectrum else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('analysis', 'UVVis_Benzene', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_nmr_calculator():
    """Test NMR Calculator"""
    print("\n[TEST] NMR Spectroscopy Calculator...")

    try:
        # Water molecule for NMR
        bond = BondFactory.create_from_smiles('O')

        calculator = NMRCalculator(bond)

        start_time = time.time()
        result = calculator.calculate_chemical_shifts(use_quantum=False)
        elapsed_time = time.time() - start_time

        has_shifts = 'chemical_shifts' in result

        test_results.add_result('analysis', 'NMR_Water', {
            'status': 'PASSED' if has_shifts else 'FAILED',
            'properties': {
                'num_atoms': len(result.get('chemical_shifts', {})),
                'shifts': result.get('chemical_shifts', 'N/A'),
            },
            'time': elapsed_time
        })

        print(f"  Number of Atoms: {len(result.get('chemical_shifts', {}))}")
        print(f"  Status: {'PASSED ✓' if has_shifts else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('analysis', 'NMR_Water', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_raman_ir_calculator():
    """Test Raman/IR Spectroscopy Calculator"""
    print("\n[TEST] Raman/IR Spectroscopy Calculator...")

    try:
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        calculator = RamanIRCalculator(bond)

        start_time = time.time()
        result = calculator.calculate_vibrational_spectrum()
        elapsed_time = time.time() - start_time

        has_frequencies = 'frequencies' in result and len(result['frequencies']) > 0

        test_results.add_result('analysis', 'RamanIR_H2', {
            'status': 'PASSED' if has_frequencies else 'FAILED',
            'properties': {
                'num_modes': len(result.get('frequencies', [])),
                'frequencies': result.get('frequencies', 'N/A'),
            },
            'time': elapsed_time
        })

        print(f"  Number of Vibrational Modes: {len(result.get('frequencies', []))}")
        if result.get('frequencies'):
            print(f"  Frequencies: {result['frequencies']} cm⁻¹")
        print(f"  Status: {'PASSED ✓' if has_frequencies else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('analysis', 'RamanIR_H2', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_dos_calculator():
    """Test DOS Calculator"""
    print("\n[TEST] DOS Calculator...")

    try:
        bond = BondFactory.create_from_smiles('c1ccccc1')  # Benzene

        calculator = DOSCalculator(bond)

        start_time = time.time()
        result = calculator.calculate_dos(
            energy_range=(-10, 10),
            num_points=100
        )
        elapsed_time = time.time() - start_time

        has_dos = 'energies' in result and 'dos' in result

        test_results.add_result('analysis', 'DOS_Benzene', {
            'status': 'PASSED' if has_dos else 'FAILED',
            'properties': {
                'num_points': len(result.get('energies', [])),
                'energy_range': f"{result.get('energies', [0])[0]:.2f} to {result.get('energies', [0])[-1]:.2f} eV" if result.get('energies') and len(result.get('energies', [])) > 0 else 'N/A',
            },
            'time': elapsed_time
        })

        print(f"  Number of Points: {len(result.get('energies', []))}")
        print(f"  Status: {'PASSED ✓' if has_dos else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('analysis', 'DOS_Benzene', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_thermochemistry_calculator():
    """Test Thermochemistry Calculator"""
    print("\n[TEST] Thermochemistry Calculator...")

    try:
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        calculator = ThermochemistryCalculator(bond)

        start_time = time.time()
        result = calculator.calculate_thermodynamics(
            temperature=298.15,  # Room temperature
            pressure=101325  # 1 atm
        )
        elapsed_time = time.time() - start_time

        has_thermo = 'enthalpy' in result and 'entropy' in result

        test_results.add_result('analysis', 'Thermochemistry_H2', {
            'status': 'PASSED' if has_thermo else 'FAILED',
            'properties': {
                'enthalpy': result.get('enthalpy', 'N/A'),
                'entropy': result.get('entropy', 'N/A'),
                'gibbs_free_energy': result.get('gibbs_free_energy', 'N/A'),
            },
            'time': elapsed_time
        })

        print(f"  Enthalpy: {result.get('enthalpy', 'N/A')}")
        print(f"  Entropy: {result.get('entropy', 'N/A')}")
        print(f"  Gibbs Free Energy: {result.get('gibbs_free_energy', 'N/A')}")
        print(f"  Status: {'PASSED ✓' if has_thermo else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('analysis', 'Thermochemistry_H2', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_bond_length_scanner():
    """Test Bond Length Scanner for PES"""
    print("\n[TEST] Bond Length Scanner - PES...")

    try:
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        scanner = BondLengthScanner(bond)

        start_time = time.time()
        result = scanner.scan(
            start_distance=0.5,
            end_distance=2.0,
            num_points=5,
            use_quantum=False
        )
        elapsed_time = time.time() - start_time

        has_pes = 'distances' in result and 'energies' in result

        test_results.add_result('analysis', 'BondLengthScanner_H2', {
            'status': 'PASSED' if has_pes else 'FAILED',
            'properties': {
                'num_points': len(result.get('distances', [])),
                'min_energy': min(result.get('energies', [0])) if result.get('energies') else 'N/A',
            },
            'time': elapsed_time
        })

        print(f"  Number of Points: {len(result.get('distances', []))}")
        if result.get('energies'):
            print(f"  Minimum Energy: {min(result['energies']):.6f} Ha")
        print(f"  Status: {'PASSED ✓' if has_pes else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('analysis', 'BondLengthScanner_H2', {
            'status': 'FAILED',
            'error': str(e)
        })


# ============================================================================
# 4. MOLECULAR DYNAMICS VALIDATION TESTS
# ============================================================================

def test_classical_md():
    """Test Classical Molecular Dynamics"""
    print("\n[TEST] Classical Molecular Dynamics...")

    try:
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        simulator = MDSimulator(bond, mode='classical')

        start_time = time.time()
        result = simulator.run(
            num_steps=100,
            timestep=0.5,  # fs
            temperature=300  # K
        )
        elapsed_time = time.time() - start_time

        has_trajectory = 'trajectory' in result

        test_results.add_result('dynamics', 'ClassicalMD_H2', {
            'status': 'PASSED' if has_trajectory else 'FAILED',
            'steps': result.get('num_steps', 'N/A'),
            'time': elapsed_time,
            'final_energy': result.get('final_energy', 'N/A')
        })

        print(f"  Steps: {result.get('num_steps', 'N/A')}")
        print(f"  Final Energy: {result.get('final_energy', 'N/A')}")
        print(f"  Status: {'PASSED ✓' if has_trajectory else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('dynamics', 'ClassicalMD_H2', {
            'status': 'FAILED',
            'error': str(e)
        })


def test_quantum_md():
    """Test Quantum Molecular Dynamics"""
    print("\n[TEST] Quantum Molecular Dynamics...")

    try:
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        simulator = MDSimulator(bond, mode='quantum')

        start_time = time.time()
        result = simulator.run(
            num_steps=10,  # Fewer steps for quantum
            timestep=0.5,
            temperature=300,
            use_quantum=True
        )
        elapsed_time = time.time() - start_time

        has_trajectory = 'trajectory' in result

        test_results.add_result('dynamics', 'QuantumMD_H2', {
            'status': 'PASSED' if has_trajectory else 'FAILED',
            'steps': result.get('num_steps', 'N/A'),
            'time': elapsed_time,
            'final_energy': result.get('final_energy', 'N/A')
        })

        print(f"  Steps: {result.get('num_steps', 'N/A')}")
        print(f"  Final Energy: {result.get('final_energy', 'N/A')}")
        print(f"  Status: {'PASSED ✓' if has_trajectory else 'FAILED ✗'}")

    except Exception as e:
        print(f"  ERROR: {str(e)}")
        test_results.add_result('dynamics', 'QuantumMD_H2', {
            'status': 'FAILED',
            'error': str(e)
        })


# ============================================================================
# MAIN TEST EXECUTION
# ============================================================================

if __name__ == '__main__':
    print("\n")
    print("="*80)
    print(" KANAD FRAMEWORK COMPREHENSIVE VALIDATION TEST SUITE")
    print("="*80)
    print("\nStarting comprehensive validation tests...\n")

    # 1. Solver Tests
    print("\n" + "="*80)
    print(" SECTION 1: QUANTUM SOLVERS")
    print("="*80)
    test_vqe_standard_h2()
    test_hivqe_h2()
    test_vqe_lih()
    test_active_space_solver()
    test_sqd_solver()
    test_krylov_sqd_solver()

    # 2. Application Module Tests
    print("\n" + "="*80)
    print(" SECTION 2: APPLICATION MODULES")
    print("="*80)
    test_drug_discovery_platform()
    test_materials_scout()
    test_catalyst_optimizer()
    test_alloy_designer()

    # 3. Analysis Module Tests
    print("\n" + "="*80)
    print(" SECTION 3: ANALYSIS MODULES")
    print("="*80)
    test_property_calculator()
    test_uvvis_calculator()
    test_nmr_calculator()
    test_raman_ir_calculator()
    test_dos_calculator()
    test_thermochemistry_calculator()
    test_bond_length_scanner()

    # 4. Molecular Dynamics Tests
    print("\n" + "="*80)
    print(" SECTION 4: MOLECULAR DYNAMICS")
    print("="*80)
    test_classical_md()
    test_quantum_md()

    # Print final report
    test_results.print_report()

    # Save results to JSON
    output_file = '/home/user/kanad/validation_results.json'
    with open(output_file, 'w') as f:
        json.dump(test_results.results, f, indent=2)
    print(f"\nResults saved to: {output_file}\n")
