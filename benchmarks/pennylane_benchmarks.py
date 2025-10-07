"""
PennyLane QChem Benchmarks

Comparison using PennyLane's quantum chemistry module.

Install:
    pip install pennylane pennylane-qchem pyscf
"""

import sys
import time
import json
from pathlib import Path
from datetime import datetime
import numpy as np

try:
    import pennylane as qml
    from pennylane import numpy as pnp
    from pennylane import qchem
    PENNYLANE_AVAILABLE = True
except ImportError as e:
    print(f"⚠️  PennyLane not available: {e}")
    print("Install with: pip install pennylane pennylane-qchem pyscf")
    PENNYLANE_AVAILABLE = False
    sys.exit(1)


class PennyLaneBenchmark:
    """PennyLane QChem VQE benchmarks."""

    def __init__(self, output_dir='benchmarks/results'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {
            'framework': 'PennyLane QChem',
            'version': 'Latest',
            'timestamp': datetime.now().isoformat(),
            'molecules': {}
        }

    def benchmark_molecule(self, name, atom1, atom2, distance, charge=0, exact_energy=None):
        """Benchmark molecule with PennyLane."""
        print(f"\n{'='*80}")
        print(f"PennyLane QChem: {name}")
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
            symbols = [atom1, atom2]
            coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, distance]])

            print(f"  Building molecule: {atom1}-{atom2} at {distance} Å")

            # Get molecular data
            t_start = time.time()

            # Build Hamiltonian
            H, n_qubits = qchem.molecular_hamiltonian(
                symbols,
                coordinates,
                charge=charge,
                basis='sto-3g'
            )

            # Get HF state
            n_electrons = sum([self._get_atomic_number(s) for s in symbols]) - charge
            hf_state = qchem.hf_state(n_electrons, n_qubits)

            t_setup = time.time() - t_start

            print(f"  Qubits: {n_qubits}")
            print(f"  Electrons: {n_electrons}")
            print(f"  Setup time: {t_setup:.3f}s")

            mol_results['n_qubits'] = int(n_qubits)
            mol_results['n_electrons'] = int(n_electrons)

            # Run VQE
            print(f"\n--- VQE with Hardware-Efficient Ansatz ---")
            config_result = self._run_vqe(
                H, n_qubits, n_electrons, hf_state, exact_energy
            )
            mol_results['configurations'].append(config_result)

        except Exception as e:
            print(f"  ✗ Error: {e}")
            import traceback
            traceback.print_exc()
            mol_results['error'] = str(e)

        self.results['molecules'][name] = mol_results
        return mol_results

    def _get_atomic_number(self, symbol):
        """Get atomic number from symbol."""
        atomic_numbers = {
            'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5,
            'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10
        }
        return atomic_numbers.get(symbol, 0)

    def _run_vqe(self, hamiltonian, n_qubits, n_electrons, hf_state, exact_energy):
        """Run VQE with PennyLane."""
        config_result = {
            'description': 'PennyLane VQE + Hardware-Efficient',
            'ansatz': 'Hardware-Efficient',
            'mapper': 'Jordan-Wigner',
            'success': False
        }

        try:
            t_start = time.time()

            # Create device
            dev = qml.device('default.qubit', wires=n_qubits)

            # Define ansatz
            def ansatz(params):
                # HF initial state
                qml.BasisState(hf_state, wires=range(n_qubits))

                # Hardware-efficient layers (2 layers)
                n_layers = 2
                for layer in range(n_layers):
                    # Rotation layer
                    for q in range(n_qubits):
                        qml.RY(params[layer * n_qubits + q], wires=q)

                    # Entanglement layer
                    for q in range(n_qubits - 1):
                        qml.CNOT(wires=[q, q + 1])

            # VQE cost function
            @qml.qnode(dev)
            def cost_fn(params):
                ansatz(params)
                return qml.expval(hamiltonian)

            # Initial parameters
            n_params = 2 * n_qubits  # 2 layers
            params = pnp.zeros(n_params, requires_grad=True)

            # Optimize (limited iterations for benchmark)
            opt = qml.GradientDescentOptimizer(stepsize=0.4)
            max_iterations = 100  # Limited for speed

            energy = cost_fn(params)
            print(f"  Initial energy: {energy:.6f} Ha")

            for i in range(max_iterations):
                params, energy = opt.step_and_cost(cost_fn, params)

                if i % 20 == 0:
                    print(f"  Iteration {i}: {energy:.6f} Ha")

            t_vqe = time.time() - t_start

            # Get HF energy (cost at zero parameters with HF state)
            hf_params = pnp.zeros(n_params)
            hf_energy = cost_fn(hf_params)

            correlation = (energy - hf_energy) * 1000  # mHa

            if exact_energy:
                error_vs_exact = abs(energy - exact_energy) * 1000
            else:
                error_vs_exact = None

            config_result.update({
                'success': True,
                'vqe_energy': float(energy),
                'hf_energy': float(hf_energy),
                'correlation_energy_mha': float(correlation),
                'n_iterations': max_iterations,
                'time_seconds': float(t_vqe),
                'error_vs_exact_mha': error_vs_exact,
                'n_parameters': int(n_params)
            })

            print(f"  Final VQE Energy: {energy:.6f} Ha")
            print(f"  HF Energy: {hf_energy:.6f} Ha")
            print(f"  Correlation: {correlation:.3f} mHa")
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
        print("PENNYLANE QCHEM BENCHMARKS")
        print("="*80)

        molecules = [
            {'name': 'H2', 'atom1': 'H', 'atom2': 'H', 'distance': 0.74, 'charge': 0, 'exact_energy': -1.1373},
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
        output_file = self.output_dir / 'pennylane_results.json'
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\n✓ Results saved to: {output_file}")

    def print_summary(self):
        """Print summary."""
        print("\n" + "="*80)
        print("PENNYLANE SUMMARY")
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
    """Run PennyLane benchmarks."""
    if not PENNYLANE_AVAILABLE:
        return

    benchmark = PennyLaneBenchmark()
    benchmark.run_all_benchmarks()

    print("\n" + "="*80)
    print("✓ PennyLane benchmarking complete!")
    print("="*80)


if __name__ == '__main__':
    main()
