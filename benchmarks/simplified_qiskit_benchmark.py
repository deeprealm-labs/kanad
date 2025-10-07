"""
Simplified Qiskit Benchmark

Uses qiskit-nature with current environment (Python 3.13).
Simplified to work without VQE, just gets molecular data.
"""

import sys
import time
import json
from pathlib import Path
from datetime import datetime

try:
    from qiskit_nature.second_q.drivers import PySCFDriver
    from qiskit_nature.second_q.mappers import JordanWignerMapper
    from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
    QISKIT_NATURE_AVAILABLE = True
except ImportError as e:
    print(f"⚠️  Qiskit Nature components not available: {e}")
    QISKIT_NATURE_AVAILABLE = False
    sys.exit(1)


class SimplifiedQiskitBenchmark:
    """Simplified Qiskit Nature benchmark - molecular properties only."""

    def __init__(self, output_dir='benchmarks/results'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {
            'framework': 'Qiskit Nature',
            'version': '0.7.2 (Simplified)',
            'timestamp': datetime.now().isoformat(),
            'molecules': {}
        }

    def benchmark_molecule(self, name, atom1, atom2, distance, charge=0, exact_energy=None):
        """Get molecular properties using Qiskit Nature."""
        print(f"\n{'='*80}")
        print(f"Qiskit Nature (Simplified): {name}")
        print(f"{'='*80}")

        mol_results = {
            'name': name,
            'atoms': [atom1, atom2],
            'distance': distance,
            'charge': charge,
            'exact_energy': exact_energy,
            'configurations': []
        }

        try:
            # Build molecule
            atom_string = f"{atom1} 0.0 0.0 0.0; {atom2} 0.0 0.0 {distance}"

            print(f"  Building molecule: {atom1}-{atom2} at {distance} Å")
            t_start = time.time()

            driver = PySCFDriver(
                atom=atom_string,
                basis='sto3g',
                charge=charge,
                spin=0
            )

            problem = driver.run()
            t_setup = time.time() - t_start

            # Get properties
            hf_energy = problem.reference_energy
            n_orbitals = problem.num_spatial_orbitals
            n_particles = problem.num_particles

            print(f"  Orbitals: {n_orbitals}")
            print(f"  Particles: {n_particles}")
            print(f"  HF Energy: {hf_energy:.6f} Ha")
            print(f"  Setup time: {t_setup:.3f}s")

            # Get Hamiltonian for qubit count
            hamiltonian = problem.hamiltonian.second_q_op()
            mapper = JordanWignerMapper()
            qubit_op = mapper.map(hamiltonian)

            n_qubits = qubit_op.num_qubits

            print(f"  Qubits: {n_qubits}")

            mol_results['hf_energy'] = float(hf_energy)
            mol_results['n_orbitals'] = int(n_orbitals)
            mol_results['n_particles'] = sum(n_particles)
            mol_results['n_qubits'] = int(n_qubits)

            # Add a result indicating we got molecular data
            config_result = {
                'description': 'Qiskit Nature Molecular Data',
                'method': 'PySCF Driver',
                'success': True,
                'hf_energy': float(hf_energy),
                'time_seconds': float(t_setup),
                'note': 'Molecular properties only (VQE not run due to environment constraints)'
            }

            mol_results['configurations'].append(config_result)

        except Exception as e:
            print(f"  ✗ Error: {e}")
            import traceback
            traceback.print_exc()
            mol_results['error'] = str(e)

        self.results['molecules'][name] = mol_results
        return mol_results

    def run_all_benchmarks(self):
        """Run all benchmarks."""
        print("="*80)
        print("QISKIT NATURE SIMPLIFIED BENCHMARKS")
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
            if 'hf_energy' in mol_data:
                print(f"\n{mol_name}:")
                print(f"  HF Energy: {mol_data['hf_energy']:.6f} Ha")
                print(f"  Qubits: {mol_data.get('n_qubits', 'N/A')}")


def main():
    """Run Qiskit Nature benchmarks."""
    if not QISKIT_NATURE_AVAILABLE:
        return

    benchmark = SimplifiedQiskitBenchmark()
    benchmark.run_all_benchmarks()

    print("\n" + "="*80)
    print("✓ Qiskit Nature benchmarking complete!")
    print("="*80)


if __name__ == '__main__':
    main()
