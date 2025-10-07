"""
Qiskit Nature Benchmarks

Full comparison using qiskit-nature VQE implementation.
Requires Python 3.11 environment (see SETUP_GUIDE.md)

Run with:
    source venv_qiskit_nature/bin/activate
    python benchmarks/qiskit_nature_benchmarks.py
"""

import sys
import time
import json
from pathlib import Path
from datetime import datetime
import numpy as np

try:
    from qiskit_nature.units import DistanceUnit
    from qiskit_nature.second_q.drivers import PySCFDriver
    from qiskit_nature.second_q.mappers import JordanWignerMapper, ParityMapper
    from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
    from qiskit_algorithms import VQE
    from qiskit_algorithms.optimizers import SLSQP
    from qiskit.primitives import Estimator
    QISKIT_NATURE_AVAILABLE = True
except ImportError as e:
    print(f"⚠️  Qiskit Nature not available: {e}")
    print("This benchmark requires Python 3.11 and qiskit-nature==0.7.2")
    print("See benchmarks/SETUP_GUIDE.md for setup instructions")
    QISKIT_NATURE_AVAILABLE = False
    sys.exit(1)


class QiskitNatureBenchmark:
    """Qiskit Nature VQE benchmarks."""

    def __init__(self, output_dir='benchmarks/results'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {
            'framework': 'Qiskit Nature',
            'version': '0.7.2',
            'timestamp': datetime.now().isoformat(),
            'molecules': {}
        }

    def benchmark_molecule(self, name, atom1, atom2, distance, charge=0, exact_energy=None):
        """Benchmark molecule with Qiskit Nature VQE."""
        print(f"\n{'='*80}")
        print(f"Qiskit Nature: {name}")
        print(f"{'='*80}")

        mol_results = {
            'name': name,
            'atoms': [atom1, atom2],
            'distance': distance,
            'charge': charge,
            'exact_energy': exact_energy,
            'configurations': []
        }

        # Build molecule
        atom_string = f"{atom1} 0.0 0.0 0.0; {atom2} 0.0 0.0 {distance}"

        try:
            # PySCF driver
            print(f"  Building molecule: {atom1}-{atom2} at {distance} Å")
            driver = PySCFDriver(
                atom=atom_string,
                basis='sto3g',
                charge=charge,
                spin=0
            )

            problem = driver.run()

            # Get HF reference
            hf_energy = problem.reference_energy
            n_orbitals = problem.num_spatial_orbitals
            n_particles = problem.num_particles

            print(f"  Orbitals: {n_orbitals}")
            print(f"  Particles: {n_particles}")
            print(f"  HF Energy: {hf_energy:.6f} Ha")

            mol_results['hf_energy'] = float(hf_energy)
            mol_results['n_orbitals'] = int(n_orbitals)
            mol_results['n_particles'] = sum(n_particles)

            # Get Hamiltonian
            hamiltonian = problem.hamiltonian.second_q_op()

            # Mapper
            mapper = JordanWignerMapper()

            # UCCSD ansatz
            print(f"\n--- VQE with UCCSD Ansatz ---")
            config_result = self._run_uccsd_vqe(
                problem, hamiltonian, mapper,
                hf_energy, exact_energy
            )
            mol_results['configurations'].append(config_result)

        except Exception as e:
            print(f"  ✗ Error: {e}")
            import traceback
            traceback.print_exc()
            mol_results['error'] = str(e)

        self.results['molecules'][name] = mol_results
        return mol_results

    def _run_uccsd_vqe(self, problem, hamiltonian, mapper, hf_energy, exact_energy):
        """Run VQE with UCCSD ansatz."""
        config_result = {
            'description': 'Qiskit Nature VQE + UCCSD',
            'ansatz': 'UCCSD',
            'mapper': 'Jordan-Wigner',
            'success': False
        }

        try:
            t_start = time.time()

            # Map to qubits
            qubit_op = mapper.map(hamiltonian)

            # HF initial state
            hf_state = HartreeFock(
                problem.num_spatial_orbitals,
                problem.num_particles,
                mapper
            )

            # UCCSD ansatz
            ansatz = UCCSD(
                problem.num_spatial_orbitals,
                problem.num_particles,
                mapper,
                initial_state=hf_state
            )

            # VQE
            estimator = Estimator()
            optimizer = SLSQP(maxiter=1000)

            vqe = VQE(estimator, ansatz, optimizer)
            result = vqe.compute_minimum_eigenvalue(qubit_op)

            t_vqe = time.time() - t_start

            vqe_energy = result.eigenvalue
            n_iterations = result.optimizer_result.nfev if hasattr(result.optimizer_result, 'nfev') else 0

            correlation = (vqe_energy - hf_energy) * 1000  # mHa

            if exact_energy:
                error_vs_exact = abs(vqe_energy - exact_energy) * 1000
            else:
                error_vs_exact = None

            config_result.update({
                'success': True,
                'vqe_energy': float(vqe_energy),
                'correlation_energy_mha': float(correlation),
                'n_iterations': int(n_iterations),
                'time_seconds': float(t_vqe),
                'error_vs_exact_mha': error_vs_exact,
                'n_parameters': len(ansatz.parameters)
            })

            print(f"  VQE Energy: {vqe_energy:.6f} Ha")
            print(f"  Correlation: {correlation:.3f} mHa")
            print(f"  Iterations: {n_iterations}")
            print(f"  Time: {t_vqe:.2f}s")
            if error_vs_exact:
                print(f"  Error vs Exact: {error_vs_exact:.3f} mHa")

        except Exception as e:
            print(f"  ✗ Error: {e}")
            config_result['error'] = str(e)
            import traceback
            traceback.print_exc()

        return config_result

    def run_all_benchmarks(self):
        """Run all benchmarks."""
        print("="*80)
        print("QISKIT NATURE VQE BENCHMARKS")
        print("="*80)

        molecules = [
            {'name': 'H2', 'atom1': 'H', 'atom2': 'H', 'distance': 0.74, 'charge': 0, 'exact_energy': -1.1373},
            {'name': 'HeH+', 'atom1': 'He', 'atom2': 'H', 'distance': 0.772, 'charge': 1, 'exact_energy': -2.8510},
            {'name': 'LiH', 'atom1': 'Li', 'atom2': 'H', 'distance': 1.60, 'charge': 0, 'exact_energy': -7.8823}
        ]

        for mol in molecules:
            try:
                self.benchmark_molecule(
                    mol['name'],
                    mol['atom1'],
                    mol['atom2'],
                    mol['distance'],
                    mol.get('charge', 0),
                    mol.get('exact_energy')
                )
            except Exception as e:
                print(f"\n✗ Failed: {e}")

        self.save_results()
        self.print_summary()

    def save_results(self):
        """Save results."""
        output_file = self.output_dir / 'qiskit_nature_results.json'
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\n✓ Results saved to: {output_file}")

    def print_summary(self):
        """Print summary."""
        print("\n" + "="*80)
        print("QISKIT NATURE SUMMARY")
        print("="*80)

        for mol_name, mol_data in self.results['molecules'].items():
            if 'configurations' in mol_data:
                print(f"\n{mol_name}:")
                for config in mol_data['configurations']:
                    if config['success']:
                        print(f"  Energy: {config['vqe_energy']:.6f} Ha")
                        print(f"  Correlation: {config['correlation_energy_mha']:.3f} mHa")
                        print(f"  Time: {config['time_seconds']:.2f}s")


def main():
    """Run Qiskit Nature benchmarks."""
    if not QISKIT_NATURE_AVAILABLE:
        return

    benchmark = QiskitNatureBenchmark()
    benchmark.run_all_benchmarks()

    print("\n" + "="*80)
    print("✓ Qiskit Nature benchmarking complete!")
    print("="*80)


if __name__ == '__main__':
    main()
