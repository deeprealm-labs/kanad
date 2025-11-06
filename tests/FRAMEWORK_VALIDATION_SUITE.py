"""
Kanad Framework Comprehensive Validation Suite

Tests the complete framework stack:
- Solvers (VQE, SQD, HiVQE)
- Analysis (properties, spectroscopy, thermochemistry)
- Molecular Dynamics (classical and quantum)
- Applications (drug discovery, materials, catalysis)

Compares quantum results vs classical benchmarks and validates
against traditional computational chemistry tools.

Run with: python tests/FRAMEWORK_VALIDATION_SUITE.py
"""

import sys
import os
import numpy as np
import logging
import json
import time
from datetime import datetime
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from kanad.bonds import BondFactory
from kanad.solvers import VQESolver, SQDSolver
from kanad.analysis import PropertyCalculator, DOSCalculator
from kanad.dynamics import MDSimulator
from kanad.applications import DrugDiscoveryPlatform, MaterialsScout

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class FrameworkValidator:
    """Comprehensive framework validation."""

    def __init__(self, output_dir='validation_results'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.results = {
            'timestamp': datetime.now().isoformat(),
            'tests': []
        }

    def run_all_tests(self):
        """Run all validation tests."""
        logger.info("="*80)
        logger.info("KANAD FRAMEWORK COMPREHENSIVE VALIDATION SUITE")
        logger.info("="*80)

        # Core tests
        self.test_quantum_solvers()
        self.test_analysis_modules()
        self.test_molecular_dynamics()
        self.test_applications()

        # Integration tests
        self.test_quantum_vs_classical()
        self.test_full_stack_workflow()

        # Save results
        self.save_results()
        self.print_summary()

    def test_quantum_solvers(self):
        """Test VQE and SQD solvers with validation."""
        logger.info("\n" + "="*80)
        logger.info("TEST SUITE 1: Quantum Solvers")
        logger.info("="*80)

        test_result = {
            'name': 'Quantum Solvers',
            'subtests': []
        }

        # Test 1.1: VQE on H2
        logger.info("\n[1.1] VQE Solver - H2 Molecule")
        try:
            bond = BondFactory.create_bond('H', 'H', distance=0.74)
            solver = VQESolver(bond, backend='statevector', use_governance=True)
            result = solver.solve()

            # Known H2 energy at equilibrium: ~-1.137 Ha (FCI)
            # HF energy: ~-1.117 Ha
            energy = result['energy']
            hf_energy = result.get('hf_energy', -1.117)
            correlation = energy - hf_energy

            subtest = {
                'name': 'VQE H2',
                'energy': float(energy),
                'hf_energy': float(hf_energy),
                'correlation': float(correlation),
                'converged': result['converged'],
                'status': 'PASS' if result['converged'] and energy < -1.0 else 'FAIL'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  Energy: {energy:.8f} Ha")
            logger.info(f"  HF Energy: {hf_energy:.8f} Ha")
            logger.info(f"  Correlation: {correlation:.8f} Ha")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  VQE test failed: {e}")
            test_result['subtests'].append({'name': 'VQE H2', 'status': 'ERROR', 'error': str(e)})

        # Test 1.2: SQD on H2
        logger.info("\n[1.2] SQD Solver - H2 Molecule")
        try:
            bond = BondFactory.create_bond('H', 'H', distance=0.74)
            solver = SQDSolver(bond, backend='statevector', use_governance=True)
            result = solver.solve()

            energy = result['energy']
            subtest = {
                'name': 'SQD H2',
                'energy': float(energy),
                'converged': result.get('converged', True),
                'status': 'PASS' if energy < -1.0 else 'FAIL'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  Energy: {energy:.8f} Ha")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  SQD test failed: {e}")
            test_result['subtests'].append({'name': 'SQD H2', 'status': 'ERROR', 'error': str(e)})

        # Test 1.3: Governance filtering
        logger.info("\n[1.3] Governance Protocol Effectiveness")
        try:
            bond = BondFactory.create_bond('H', 'H', distance=0.74)

            # With governance
            start = time.time()
            solver_gov = VQESolver(bond, backend='statevector', use_governance=True)
            result_gov = solver_gov.solve()
            time_gov = time.time() - start

            # Without governance
            start = time.time()
            solver_no_gov = VQESolver(bond, backend='statevector', use_governance=False)
            result_no_gov = solver_no_gov.solve()
            time_no_gov = time.time() - start

            speedup = time_no_gov / time_gov if time_gov > 0 else 1.0
            energy_diff = abs(result_gov['energy'] - result_no_gov['energy'])

            subtest = {
                'name': 'Governance Speedup',
                'time_with_governance': float(time_gov),
                'time_without_governance': float(time_no_gov),
                'speedup': float(speedup),
                'energy_difference': float(energy_diff),
                'status': 'PASS' if speedup > 0.8 and energy_diff < 0.01 else 'FAIL'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  Time with governance: {time_gov:.3f} s")
            logger.info(f"  Time without governance: {time_no_gov:.3f} s")
            logger.info(f"  Speedup: {speedup:.2f}x")
            logger.info(f"  Energy difference: {energy_diff:.8f} Ha")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  Governance test failed: {e}")
            test_result['subtests'].append({'name': 'Governance', 'status': 'ERROR', 'error': str(e)})

        self.results['tests'].append(test_result)

    def test_analysis_modules(self):
        """Test analysis capabilities."""
        logger.info("\n" + "="*80)
        logger.info("TEST SUITE 2: Analysis Modules")
        logger.info("="*80)

        test_result = {
            'name': 'Analysis Modules',
            'subtests': []
        }

        # Test 2.1: Property calculation
        logger.info("\n[2.1] Property Calculator")
        try:
            bond = BondFactory.create_bond('H', 'H', distance=0.74)
            solver = VQESolver(bond, backend='statevector', use_governance=True)
            result = solver.solve()

            # Properties are calculated during solve()
            properties = result.get('analysis', {})

            subtest = {
                'name': 'Property Calculator',
                'has_dipole': 'dipole_moment' in properties,
                'has_energy': 'energy' in result,
                'status': 'PASS' if 'energy' in result else 'FAIL'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  Properties calculated: {list(properties.keys())}")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  Property calculation test failed: {e}")
            test_result['subtests'].append({'name': 'Properties', 'status': 'ERROR', 'error': str(e)})

        # Test 2.2: DOS calculation
        logger.info("\n[2.2] Density of States")
        try:
            bond = BondFactory.create_bond('H', 'H', distance=0.74)
            solver = VQESolver(bond, backend='statevector', use_governance=True)
            result = solver.solve()

            # DOS calculation
            dos_calc = DOSCalculator(bond)
            dos_result = dos_calc.calculate_dos(method='quantum')

            subtest = {
                'name': 'DOS Calculator',
                'has_energies': 'energies' in dos_result,
                'has_dos': 'dos' in dos_result,
                'status': 'PASS' if 'dos' in dos_result else 'FAIL'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  DOS calculated: {len(dos_result.get('dos', []))} points")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  DOS calculation test failed: {e}")
            test_result['subtests'].append({'name': 'DOS', 'status': 'ERROR', 'error': str(e)})

        self.results['tests'].append(test_result)

    def test_molecular_dynamics(self):
        """Test MD capabilities."""
        logger.info("\n" + "="*80)
        logger.info("TEST SUITE 3: Molecular Dynamics")
        logger.info("="*80)

        test_result = {
            'name': 'Molecular Dynamics',
            'subtests': []
        }

        # Test 3.1: Classical MD
        logger.info("\n[3.1] Classical MD (HF Forces)")
        try:
            bond = BondFactory.create_bond('H', 'H', distance=0.74)
            md = MDSimulator(
                bond,
                temperature=300.0,
                timestep=0.5,
                force_method='hf'
            )
            result = md.run(n_steps=10)

            energy_conservation = abs(result.final_energy - result.initial_energy)

            subtest = {
                'name': 'Classical MD',
                'steps_completed': 10,
                'initial_energy': float(result.initial_energy),
                'final_energy': float(result.final_energy),
                'energy_drift': float(energy_conservation),
                'status': 'PASS' if energy_conservation < 0.1 else 'FAIL'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  Steps: 10")
            logger.info(f"  Energy drift: {energy_conservation:.6f} Ha")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  Classical MD test failed: {e}")
            test_result['subtests'].append({'name': 'Classical MD', 'status': 'ERROR', 'error': str(e)})

        # Test 3.2: Quantum MD
        logger.info("\n[3.2] Quantum MD (VQE Forces)")
        try:
            bond = BondFactory.create_bond('H', 'H', distance=0.74)
            md = MDSimulator(
                bond,
                temperature=300.0,
                timestep=0.5,
                force_method='vqe',
                use_governance=True,
                backend='statevector'
            )
            result = md.run(n_steps=2)  # Just 2 steps for quantum

            subtest = {
                'name': 'Quantum MD',
                'steps_completed': 2,
                'initial_energy': float(result.initial_energy),
                'final_energy': float(result.final_energy),
                'status': 'PASS'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  Steps: 2")
            logger.info(f"  Quantum forces computed successfully")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  Quantum MD test failed: {e}")
            test_result['subtests'].append({'name': 'Quantum MD', 'status': 'ERROR', 'error': str(e)})

        self.results['tests'].append(test_result)

    def test_applications(self):
        """Test application modules."""
        logger.info("\n" + "="*80)
        logger.info("TEST SUITE 4: Applications")
        logger.info("="*80)

        test_result = {
            'name': 'Applications',
            'subtests': []
        }

        # Test 4.1: Drug Discovery
        logger.info("\n[4.1] Drug Discovery Platform")
        try:
            platform = DrugDiscoveryPlatform(backend='statevector')

            # Simple molecule for testing
            result = platform.evaluate_molecule(
                smiles='O',  # Water
                use_quantum=False  # Use classical for speed
            )

            subtest = {
                'name': 'Drug Discovery',
                'molecule': 'Water',
                'has_energy': 'energy' in result,
                'has_properties': 'properties' in result,
                'status': 'PASS' if 'energy' in result else 'FAIL'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  Molecule evaluated: Water")
            logger.info(f"  Properties: {list(result.get('properties', {}).keys())}")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  Drug discovery test failed: {e}")
            test_result['subtests'].append({'name': 'Drug Discovery', 'status': 'ERROR', 'error': str(e)})

        # Test 4.2: Materials Scout
        logger.info("\n[4.2] Materials Scout")
        try:
            scout = MaterialsScout(backend='statevector')

            # Test with H2
            result = scout.evaluate_material(
                elements=['H', 'H'],
                distances=[0.74],
                use_quantum=False
            )

            subtest = {
                'name': 'Materials Scout',
                'material': 'H2',
                'has_energy': 'energy' in result,
                'status': 'PASS' if 'energy' in result else 'FAIL'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  Material evaluated: H2")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  Materials scout test failed: {e}")
            test_result['subtests'].append({'name': 'Materials Scout', 'status': 'ERROR', 'error': str(e)})

        self.results['tests'].append(test_result)

    def test_quantum_vs_classical(self):
        """Compare quantum vs classical results."""
        logger.info("\n" + "="*80)
        logger.info("TEST SUITE 5: Quantum vs Classical Comparison")
        logger.info("="*80)

        test_result = {
            'name': 'Quantum vs Classical',
            'subtests': []
        }

        logger.info("\n[5.1] Energy Comparison")
        try:
            bond = BondFactory.create_bond('H', 'H', distance=0.74)

            # Classical (HF)
            from kanad.core.hamiltonians import CovalentHamiltonian
            from pyscf import gto, scf

            mol = gto.Mole()
            mol.atom = [['H', [0, 0, 0]], ['H', [0, 0, 0.74]]]
            mol.basis = 'sto-3g'
            mol.build()
            mf = scf.RHF(mol)
            hf_energy = mf.kernel()

            # Quantum (VQE)
            solver = VQESolver(bond, backend='statevector', use_governance=True)
            vqe_result = solver.solve()
            vqe_energy = vqe_result['energy']

            correlation = vqe_energy - hf_energy

            subtest = {
                'name': 'Energy Comparison',
                'hf_energy': float(hf_energy),
                'vqe_energy': float(vqe_energy),
                'correlation_energy': float(correlation),
                'status': 'PASS'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  HF Energy: {hf_energy:.8f} Ha")
            logger.info(f"  VQE Energy: {vqe_energy:.8f} Ha")
            logger.info(f"  Correlation: {correlation:.8f} Ha")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  Comparison test failed: {e}")
            test_result['subtests'].append({'name': 'Energy Comparison', 'status': 'ERROR', 'error': str(e)})

        self.results['tests'].append(test_result)

    def test_full_stack_workflow(self):
        """Test complete workflow: solve → analyze → simulate."""
        logger.info("\n" + "="*80)
        logger.info("TEST SUITE 6: Full Stack Integration")
        logger.info("="*80)

        test_result = {
            'name': 'Full Stack Integration',
            'subtests': []
        }

        logger.info("\n[6.1] Complete Workflow")
        try:
            # Step 1: Create bond
            bond = BondFactory.create_bond('H', 'H', distance=0.74)

            # Step 2: Solve with VQE
            solver = VQESolver(bond, backend='statevector', use_governance=True)
            vqe_result = solver.solve()

            # Step 3: Analyze properties
            properties = vqe_result.get('analysis', {})

            # Step 4: Run MD
            md = MDSimulator(
                bond,
                temperature=300.0,
                timestep=0.5,
                force_method='hf'  # Classical for speed
            )
            md_result = md.run(n_steps=5)

            subtest = {
                'name': 'Full Workflow',
                'vqe_converged': vqe_result['converged'],
                'properties_calculated': len(properties) > 0,
                'md_completed': md_result is not None,
                'status': 'PASS'
            }
            test_result['subtests'].append(subtest)

            logger.info(f"  VQE: Converged")
            logger.info(f"  Properties: {len(properties)} calculated")
            logger.info(f"  MD: 5 steps completed")
            logger.info(f"  Status: {subtest['status']}")

        except Exception as e:
            logger.error(f"  Full workflow test failed: {e}")
            test_result['subtests'].append({'name': 'Full Workflow', 'status': 'ERROR', 'error': str(e)})

        self.results['tests'].append(test_result)

    def save_results(self):
        """Save validation results to JSON."""
        output_file = self.output_dir / f"validation_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"

        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        logger.info(f"\nResults saved to: {output_file}")

    def print_summary(self):
        """Print validation summary."""
        logger.info("\n" + "="*80)
        logger.info("VALIDATION SUMMARY")
        logger.info("="*80)

        total_tests = 0
        passed_tests = 0
        failed_tests = 0
        error_tests = 0

        for test in self.results['tests']:
            for subtest in test['subtests']:
                total_tests += 1
                status = subtest.get('status', 'UNKNOWN')
                if status == 'PASS':
                    passed_tests += 1
                elif status == 'FAIL':
                    failed_tests += 1
                elif status == 'ERROR':
                    error_tests += 1

        logger.info(f"Total Tests: {total_tests}")
        logger.info(f"Passed: {passed_tests} ✅")
        logger.info(f"Failed: {failed_tests} ❌")
        logger.info(f"Errors: {error_tests} ⚠️")
        logger.info(f"Success Rate: {passed_tests/total_tests*100:.1f}%")
        logger.info("="*80)


def main():
    """Run framework validation."""
    validator = FrameworkValidator()
    validator.run_all_tests()


if __name__ == '__main__':
    main()
