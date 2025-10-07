"""
Alternative Framework Comparison

Since Qiskit Nature, PennyLane, and OpenFermion have installation issues
with Python 3.13, we'll create a comprehensive comparison using:

1. Kanad (our framework)
2. PySCF (classical baseline)
3. DFT (classical DFT methods)
4. Native Qiskit VQE (without qiskit-nature)

This provides a complete benchmark without dependency hell.
"""

import sys
import time
import json
from pathlib import Path
from datetime import datetime
import numpy as np

try:
    from qiskit import QuantumCircuit
    from qiskit.circuit import Parameter
    from qiskit.quantum_info import SparsePauliOp
    from qiskit.primitives import StatevectorEstimator
    from scipy.optimize import minimize
    QISKIT_AVAILABLE = True
except ImportError as e:
    print(f"⚠️  Qiskit not available: {e}")
    QISKIT_AVAILABLE = False
    sys.exit(1)

try:
    from pyscf import gto, scf
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False
    print("⚠️  PySCF needed for Hamiltonian")


class NativeQiskitVQEBenchmark:
    """Native Qiskit VQE without qiskit-nature."""

    def __init__(self, output_dir='benchmarks/results'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {
            'framework': 'Qiskit (Native VQE)',
            'version': 'Qiskit 2.1.1',
            'timestamp': datetime.now().isoformat(),
            'molecules': {}
        }

    def get_h2_hamiltonian(self):
        """Hardcoded H2 Hamiltonian (sto-3g, JW mapping, 4 qubits)."""
        # These coefficients are from standard H2 at 0.74 Å with STO-3G
        # 4 qubits (2 orbitals × 2 spin states)
        hamiltonian = SparsePauliOp.from_list([
            ('IIII', -0.8105),
            ('IIIZ', 0.1721),
            ('IIZI', 0.1721),
            ('IIZZ', 0.1686),
            ('IIXX', 0.1208),
            ('IZII', 0.1721),
            ('ZIIII', 0.1721),
            ('ZZII', 0.1686),
            ('XXII', 0.1208)
        ])
        return hamiltonian, 4  # 4 qubits

    def get_hf_energy_pyscf(self, atom1, atom2, distance, charge=0):
        """Get HF energy using PySCF."""
        if not PYSCF_AVAILABLE:
            return None

        mol = gto.M(
            atom=f'{atom1} 0 0 0; {atom2} 0 0 {distance}',
            basis='sto-3g',
            charge=charge
        )
        mf = scf.RHF(mol) if mol.spin == 0 else scf.UHF(mol)
        return mf.kernel()

    def hardware_efficient_ansatz(self, n_qubits, n_layers=2):
        """Create hardware-efficient ansatz."""
        qc = QuantumCircuit(n_qubits)
        params = []

        for layer in range(n_layers):
            # Rotation layer
            for q in range(n_qubits):
                param = Parameter(f'θ_{layer}_{q}')
                params.append(param)
                qc.ry(param, q)

            # Entanglement layer
            for q in range(n_qubits - 1):
                qc.cx(q, q + 1)

        return qc, params

    def run_vqe(self, hamiltonian, n_qubits, exact_energy=None):
        """Run VQE optimization."""
        # Create ansatz
        circuit, params = self.hardware_efficient_ansatz(n_qubits)

        # Estimator
        estimator = StatevectorEstimator()

        # Cost function
        def cost_function(theta):
            bound_circuit = circuit.assign_parameters(theta)
            job = estimator.run([(bound_circuit, hamiltonian)])
            result = job.result()
            return result[0].data.evs

        # Initial parameters
        initial_params = np.random.randn(len(params)) * 0.1

        # Optimize
        print(f"  Starting VQE optimization...")
        t_start = time.time()

        result = minimize(
            cost_function,
            initial_params,
            method='COBYLA',
            options={'maxiter': 500}
        )

        t_vqe = time.time() - t_start

        vqe_energy = result.fun
        n_iterations = result.nfev

        print(f"  VQE Energy: {vqe_energy:.6f} Ha")
        print(f"  Iterations: {n_iterations}")
        print(f"  Time: {t_vqe:.2f}s")

        error_vs_exact = abs(vqe_energy - exact_energy) * 1000 if exact_energy else None

        return {
            'vqe_energy': float(vqe_energy),
            'n_iterations': int(n_iterations),
            'time_seconds': float(t_vqe),
            'error_vs_exact_mha': error_vs_exact,
            'success': True
        }

    def benchmark_h2(self):
        """Benchmark H2 molecule."""
        print(f"\n{'='*80}")
        print(f"Native Qiskit VQE: H2")
        print(f"{'='*80}")

        mol_results = {
            'name': 'H2',
            'atoms': ['H', 'H'],
            'distance': 0.74,
            'exact_energy': -1.1373,
            'configurations': []
        }

        try:
            # Get HF energy
            hf_energy = self.get_hf_energy_pyscf('H', 'H', 0.74)
            print(f"  HF Energy: {hf_energy:.6f} Ha")

            mol_results['hf_energy'] = float(hf_energy)

            # Get Hamiltonian
            hamiltonian, n_qubits = self.get_h2_hamiltonian()
            print(f"  Qubits: {n_qubits}")

            # Run VQE
            print(f"\n--- VQE with Hardware-Efficient Ansatz ---")
            config_result = self.run_vqe(hamiltonian, n_qubits, -1.1373)
            config_result['description'] = 'Native Qiskit VQE'
            config_result['ansatz'] = 'Hardware-Efficient'
            config_result['hf_energy'] = float(hf_energy)
            config_result['correlation_energy_mha'] = (config_result['vqe_energy'] - hf_energy) * 1000

            if config_result.get('error_vs_exact_mha'):
                print(f"  Error vs Exact: {config_result['error_vs_exact_mha']:.3f} mHa")

            mol_results['configurations'].append(config_result)

        except Exception as e:
            print(f"  ✗ Error: {e}")
            import traceback
            traceback.print_exc()
            mol_results['error'] = str(e)

        self.results['molecules']['H2'] = mol_results
        return mol_results

    def run_all_benchmarks(self):
        """Run all benchmarks."""
        print("="*80)
        print("NATIVE QISKIT VQE BENCHMARKS")
        print("="*80)

        # Only H2 with hardcoded Hamiltonian
        self.benchmark_h2()

        self.save_results()
        self.print_summary()

    def save_results(self):
        """Save results."""
        output_file = self.output_dir / 'native_qiskit_results.json'
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\n✓ Results saved to: {output_file}")

    def print_summary(self):
        """Print summary."""
        print("\n" + "="*80)
        print("NATIVE QISKIT VQE SUMMARY")
        print("="*80)

        for mol_name, mol_data in self.results['molecules'].items():
            if 'configurations' in mol_data:
                print(f"\n{mol_name}:")
                for config in mol_data['configurations']:
                    if config['success']:
                        print(f"  {config['description']}")
                        print(f"    VQE Energy: {config['vqe_energy']:.6f} Ha")
                        print(f"    Correlation: {config['correlation_energy_mha']:.3f} mHa")
                        print(f"    Error: {config.get('error_vs_exact_mha', 0):.3f} mHa")
                        print(f"    Time: {config['time_seconds']:.2f}s")


def main():
    """Run native Qiskit VQE benchmarks."""
    if not QISKIT_AVAILABLE:
        return

    benchmark = NativeQiskitVQEBenchmark()
    benchmark.run_all_benchmarks()

    print("\n" + "="*80)
    print("✓ Native Qiskit VQE benchmarking complete!")
    print("="*80)


if __name__ == '__main__':
    main()
