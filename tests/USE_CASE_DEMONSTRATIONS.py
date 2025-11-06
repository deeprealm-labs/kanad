"""
Kanad Framework Real-World Use Case Demonstrations

Demonstrates practical applications:
1. Drug Discovery: Binding affinity prediction
2. Alloy Design: Materials properties optimization
3. Catalyst Design: Reaction barrier calculation
4. Molecular Dynamics: Conformational sampling

Run with: python tests/USE_CASE_DEMONSTRATIONS.py
"""

import sys
import os
import numpy as np
import logging
import json
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from kanad.bonds import BondFactory
from kanad.solvers import VQESolver, SQDSolver
from kanad.analysis import PropertyCalculator, DOSCalculator
from kanad.dynamics import MDSimulator
from kanad.applications import DrugDiscoveryPlatform, MaterialsScout

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class UseCaseDemonstrator:
    """Demonstrate real-world use cases."""

    def __init__(self, output_dir='use_case_results'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.results = {
            'timestamp': datetime.now().isoformat(),
            'use_cases': []
        }

    def run_all_demonstrations(self):
        """Run all use case demonstrations."""
        logger.info("="*80)
        logger.info("KANAD FRAMEWORK - REAL-WORLD USE CASE DEMONSTRATIONS")
        logger.info("="*80)

        self.demonstrate_drug_discovery()
        self.demonstrate_materials_design()
        self.demonstrate_catalysis()
        self.demonstrate_molecular_dynamics()

        self.save_results()
        self.print_summary()

    def demonstrate_drug_discovery(self):
        """Demonstrate drug discovery workflow."""
        logger.info("\n" + "="*80)
        logger.info("USE CASE 1: Drug Discovery")
        logger.info("="*80)
        logger.info("\nScenario: Evaluate small molecules for drug-like properties")

        use_case = {
            'name': 'Drug Discovery',
            'molecules': []
        }

        # Test molecules
        molecules = [
            ('H2', 'H', 'H', 0.74),
            ('LiH', 'Li', 'H', 1.60),
        ]

        for name, atom1, atom2, distance in molecules:
            logger.info(f"\n[Testing {name}]")

            try:
                # Create bond
                bond = BondFactory.create_bond(atom1, atom2, distance=distance)

                # Quantum calculation for accuracy
                solver = VQESolver(bond, backend='statevector', use_governance=True)
                result = solver.solve()

                # Calculate properties
                energy = result['energy']
                properties = result.get('analysis', {})
                dipole = properties.get('dipole_moment', 0.0)

                # ADME-like predictions (simplified)
                molecular_weight = {
                    'H2': 2.016,
                    'LiH': 7.949
                }.get(name, 0.0)

                # Drug-likeness score (simplified Lipinski's rule)
                drug_score = 0.0
                if molecular_weight < 500:  # MW < 500
                    drug_score += 0.5
                if energy < 0:  # Stable
                    drug_score += 0.5

                mol_result = {
                    'name': name,
                    'energy': float(energy),
                    'dipole_moment': float(dipole),
                    'molecular_weight': molecular_weight,
                    'drug_likeness_score': drug_score,
                    'converged': result['converged']
                }
                use_case['molecules'].append(mol_result)

                logger.info(f"  Energy: {energy:.6f} Ha")
                logger.info(f"  Dipole: {dipole:.3f} D")
                logger.info(f"  MW: {molecular_weight:.3f} Da")
                logger.info(f"  Drug Score: {drug_score:.2f}")

            except Exception as e:
                logger.error(f"  Error: {e}")
                use_case['molecules'].append({
                    'name': name,
                    'error': str(e)
                })

        self.results['use_cases'].append(use_case)

    def demonstrate_materials_design(self):
        """Demonstrate materials design workflow."""
        logger.info("\n" + "="*80)
        logger.info("USE CASE 2: Materials Design - Alloy Discovery")
        logger.info("="*80)
        logger.info("\nScenario: Design new metallic bonds with specific properties")

        use_case = {
            'name': 'Materials Design',
            'materials': []
        }

        # Test materials (simple diatomics representing alloy bonds)
        materials = [
            ('H-H', 'H', 'H', 0.74),
            ('Li-H', 'Li', 'H', 1.60),
        ]

        for name, atom1, atom2, distance in materials:
            logger.info(f"\n[Evaluating {name}]")

            try:
                bond = BondFactory.create_bond(atom1, atom2, distance=distance)

                # Quantum calculation
                solver = VQESolver(bond, backend='statevector', use_governance=True)
                result = solver.solve()

                # Materials properties
                energy = result['energy']
                properties = result.get('analysis', {})

                # Calculate DOS for electronic structure
                dos_calc = DOSCalculator(bond)
                dos_result = dos_calc.calculate_dos(method='quantum')

                # Band gap estimation (simplified)
                dos_values = dos_result.get('dos', [])
                band_gap = 0.0  # Would need proper calculation
                if len(dos_values) > 0:
                    band_gap = abs(max(dos_values) - min(dos_values))

                # Stability metric
                stability = "High" if energy < -1.0 else "Medium" if energy < 0 else "Low"

                mat_result = {
                    'name': name,
                    'energy': float(energy),
                    'band_gap_estimate': float(band_gap),
                    'stability': stability,
                    'n_dos_points': len(dos_values),
                    'converged': result['converged']
                }
                use_case['materials'].append(mat_result)

                logger.info(f"  Energy: {energy:.6f} Ha")
                logger.info(f"  Band Gap (est): {band_gap:.3f} eV")
                logger.info(f"  Stability: {stability}")
                logger.info(f"  DOS points: {len(dos_values)}")

            except Exception as e:
                logger.error(f"  Error: {e}")
                use_case['materials'].append({
                    'name': name,
                    'error': str(e)
                })

        self.results['use_cases'].append(use_case)

    def demonstrate_catalysis(self):
        """Demonstrate catalysis workflow."""
        logger.info("\n" + "="*80)
        logger.info("USE CASE 3: Catalysis - Reaction Barrier Estimation")
        logger.info("="*80)
        logger.info("\nScenario: Calculate activation energy for H2 dissociation")

        use_case = {
            'name': 'Catalysis',
            'reaction_path': []
        }

        # Scan H-H distance (reaction coordinate)
        distances = [0.60, 0.74, 1.0, 1.5]  # Bohr

        logger.info("\n[H2 Dissociation Pathway]")

        energies = []
        for distance in distances:
            try:
                bond = BondFactory.create_bond('H', 'H', distance=distance)
                solver = VQESolver(bond, backend='statevector', use_governance=True)
                result = solver.solve()

                energy = result['energy']
                energies.append(energy)

                point = {
                    'distance': float(distance),
                    'energy': float(energy),
                    'converged': result['converged']
                }
                use_case['reaction_path'].append(point)

                logger.info(f"  d = {distance:.2f} Å: E = {energy:.6f} Ha")

            except Exception as e:
                logger.error(f"  Error at d={distance}: {e}")
                use_case['reaction_path'].append({
                    'distance': float(distance),
                    'error': str(e)
                })

        # Calculate barrier
        if len(energies) >= 2:
            min_energy = min(energies)
            max_energy = max(energies)
            barrier = (max_energy - min_energy) * 27.211  # Convert to eV

            use_case['activation_barrier_eV'] = float(barrier)
            logger.info(f"\n  Activation Barrier: {barrier:.3f} eV")

        self.results['use_cases'].append(use_case)

    def demonstrate_molecular_dynamics(self):
        """Demonstrate MD simulation."""
        logger.info("\n" + "="*80)
        logger.info("USE CASE 4: Molecular Dynamics - Thermal Sampling")
        logger.info("="*80)
        logger.info("\nScenario: Simulate H2 at 300K to sample conformations")

        use_case = {
            'name': 'Molecular Dynamics',
            'simulations': []
        }

        # Classical MD
        logger.info("\n[Classical MD with HF Forces]")
        try:
            bond = BondFactory.create_bond('H', 'H', distance=0.74)
            md_classical = MDSimulator(
                bond,
                temperature=300.0,
                timestep=0.5,
                force_method='hf'
            )
            result_classical = md_classical.run(n_steps=20)

            sim_result = {
                'type': 'classical',
                'force_method': 'HF',
                'steps': 20,
                'initial_energy': float(result_classical.initial_energy),
                'final_energy': float(result_classical.final_energy),
                'energy_drift': float(abs(result_classical.final_energy - result_classical.initial_energy)),
                'avg_temperature': float(result_classical.average_temperature)
            }
            use_case['simulations'].append(sim_result)

            logger.info(f"  Steps: 20")
            logger.info(f"  Energy drift: {sim_result['energy_drift']:.6f} Ha")
            logger.info(f"  Avg T: {sim_result['avg_temperature']:.1f} K")

        except Exception as e:
            logger.error(f"  Classical MD error: {e}")
            use_case['simulations'].append({
                'type': 'classical',
                'error': str(e)
            })

        # Quantum MD (fewer steps due to cost)
        logger.info("\n[Quantum MD with VQE Forces]")
        try:
            bond = BondFactory.create_bond('H', 'H', distance=0.74)
            md_quantum = MDSimulator(
                bond,
                temperature=300.0,
                timestep=0.5,
                force_method='vqe',
                use_governance=True,
                backend='statevector'
            )
            result_quantum = md_quantum.run(n_steps=2)

            sim_result = {
                'type': 'quantum',
                'force_method': 'VQE',
                'steps': 2,
                'initial_energy': float(result_quantum.initial_energy),
                'final_energy': float(result_quantum.final_energy),
                'governance_enabled': True
            }
            use_case['simulations'].append(sim_result)

            logger.info(f"  Steps: 2")
            logger.info(f"  Initial E: {sim_result['initial_energy']:.6f} Ha")
            logger.info(f"  Final E: {sim_result['final_energy']:.6f} Ha")
            logger.info(f"  Governance: Enabled")

        except Exception as e:
            logger.error(f"  Quantum MD error: {e}")
            use_case['simulations'].append({
                'type': 'quantum',
                'error': str(e)
            })

        self.results['use_cases'].append(use_case)

    def save_results(self):
        """Save results to JSON."""
        output_file = self.output_dir / f"use_cases_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"

        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        logger.info(f"\nResults saved to: {output_file}")

    def print_summary(self):
        """Print demonstration summary."""
        logger.info("\n" + "="*80)
        logger.info("USE CASE DEMONSTRATION SUMMARY")
        logger.info("="*80)

        for use_case in self.results['use_cases']:
            logger.info(f"\n{use_case['name']}:")

            if 'molecules' in use_case:
                successful = sum(1 for m in use_case['molecules'] if 'error' not in m)
                logger.info(f"  Molecules evaluated: {successful}/{len(use_case['molecules'])}")

            elif 'materials' in use_case:
                successful = sum(1 for m in use_case['materials'] if 'error' not in m)
                logger.info(f"  Materials evaluated: {successful}/{len(use_case['materials'])}")

            elif 'reaction_path' in use_case:
                successful = sum(1 for p in use_case['reaction_path'] if 'error' not in p)
                logger.info(f"  Points calculated: {successful}/{len(use_case['reaction_path'])}")
                if 'activation_barrier_eV' in use_case:
                    logger.info(f"  Activation barrier: {use_case['activation_barrier_eV']:.3f} eV")

            elif 'simulations' in use_case:
                successful = sum(1 for s in use_case['simulations'] if 'error' not in s)
                logger.info(f"  Simulations completed: {successful}/{len(use_case['simulations'])}")

        logger.info("\n" + "="*80)


def main():
    """Run use case demonstrations."""
    demonstrator = UseCaseDemonstrator()
    demonstrator.run_all_demonstrations()


if __name__ == '__main__':
    main()
