"""
cuQuantum GPU-Accelerated VQE Example

Demonstrates using NVIDIA cuQuantum for GPU-accelerated
quantum circuit simulation with significant speedup.

Requirements:
- NVIDIA GPU with CUDA support
- cuQuantum SDK: pip install cuquantum-python
- CuPy: pip install cupy-cuda11x (or cuda12x)
- Qiskit Aer GPU: pip install qiskit-aer-gpu

For installation help, see:
https://docs.nvidia.com/cuda/cuquantum/python/README.html
"""

import numpy as np
import time
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
from kanad.solvers.vqe_solver import VQESolver
from kanad.backends import check_cuquantum_available, get_gpu_info


def check_system():
    """Check if cuQuantum is available and show GPU info."""
    print("=" * 70)
    print("System Check")
    print("=" * 70)

    # Check cuQuantum
    if not check_cuquantum_available():
        print("\n❌ cuQuantum not available!")
        print("Install with:")
        print("  pip install cuquantum-python")
        print("  pip install cupy-cuda11x  # or cuda12x for CUDA 12")
        print("  pip install qiskit-aer-gpu")
        return False

    # Show GPU info
    gpu_info = get_gpu_info()
    if not gpu_info:
        print("\n❌ No GPUs detected!")
        return False

    print("\n✅ GPUs Available:")
    for gpu in gpu_info:
        print(f"  GPU {gpu['id']}: {gpu['name']}")
        print(f"    Memory: {gpu['total_memory_gb']:.2f} GB total, "
              f"{gpu['free_memory_gb']:.2f} GB free")
        print(f"    Compute: {gpu['compute_capability']}")

    return True


def run_gpu_vqe_example():
    """Run VQE with cuQuantum GPU acceleration."""

    print("\n" + "=" * 70)
    print("H2 Molecule VQE - GPU Acceleration with cuQuantum")
    print("=" * 70)

    # Create H2 molecule
    print("\n1. Setting up H2 molecule...")
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))
    molecule = Molecule([h1, h2])

    # Build Hamiltonian
    print("2. Building Hamiltonian...")
    representation = LCAORepresentation(molecule)
    hamiltonian = CovalentHamiltonian(
        molecule,
        representation,
        basis_name='sto-3g'
    )

    print(f"   Number of orbitals: {hamiltonian.n_orbitals}")
    print(f"   Number of electrons: {hamiltonian.n_electrons}")

    # Create ansatz
    print("\n3. Creating ansatz...")
    ansatz = RealAmplitudesAnsatz(
        n_qubits=hamiltonian.n_orbitals,
        n_electrons=hamiltonian.n_electrons,
        n_layers=2,
        entanglement='linear'
    )
    mapper = JordanWignerMapper()

    print(f"   Circuit qubits: {ansatz.n_qubits}")
    print(f"   Circuit parameters: {ansatz.get_num_parameters()}")

    # Run VQE with different backends
    backends_to_test = [
        ('classical', {}),
        ('cuquantum_statevector', {'precision': 'single'}),
    ]

    results = {}

    for backend_name, backend_opts in backends_to_test:
        print(f"\n{'=' * 70}")
        print(f"Running VQE with backend: {backend_name}")
        print('=' * 70)

        try:
            # Create solver
            solver = VQESolver(
                hamiltonian=hamiltonian,
                ansatz=ansatz,
                mapper=mapper,
                backend=backend_name,
                optimizer='COBYLA',
                max_iterations=50,
                **backend_opts
            )

            print(f"\n{solver}")

            # Run VQE
            print("\nOptimizing...")
            start_time = time.time()
            result = solver.solve()
            elapsed_time = time.time() - start_time

            # Store results
            results[backend_name] = {
                'energy': result['energy'],
                'time': elapsed_time,
                'iterations': result['iterations']
            }

            # Print results
            print(f"\n{'─' * 70}")
            print(f"Results:")
            print(f"  Energy: {result['energy']:.8f} Ha")
            print(f"  Iterations: {result['iterations']}")
            print(f"  Time: {elapsed_time:.2f} seconds")
            print(f"{'─' * 70}")

        except Exception as e:
            print(f"\n❌ Error with backend {backend_name}: {e}")
            import traceback
            traceback.print_exc()

    # Compare performance
    if len(results) > 1:
        print(f"\n\n{'=' * 70}")
        print("PERFORMANCE COMPARISON")
        print('=' * 70)

        cpu_time = results.get('classical', {}).get('time', 0)
        gpu_time = results.get('cuquantum_statevector', {}).get('time', 0)

        if cpu_time > 0 and gpu_time > 0:
            speedup = cpu_time / gpu_time
            print(f"CPU Time (Classical):  {cpu_time:.2f} seconds")
            print(f"GPU Time (cuQuantum):  {gpu_time:.2f} seconds")
            print(f"Speedup:               {speedup:.2f}x")

            if speedup > 1:
                print(f"\n✅ GPU is {speedup:.2f}x faster!")
            else:
                print(f"\n⚠️  For small circuits, CPU may be faster due to overhead")
                print("   GPU advantage increases with circuit size")

    return results


def benchmark_scaling():
    """Benchmark GPU vs CPU for different problem sizes."""
    print("\n" + "=" * 70)
    print("GPU vs CPU Scaling Benchmark")
    print("=" * 70)

    from kanad.backends import CuQuantumBackend

    # Test different qubit counts
    qubit_counts = [4, 6, 8, 10, 12]

    print(f"\n{'Qubits':<10} {'CPU Time':<15} {'GPU Time':<15} {'Speedup':<10}")
    print('-' * 60)

    for n_qubits in qubit_counts:
        try:
            from qiskit import QuantumCircuit

            # Create test circuit
            qc = QuantumCircuit(n_qubits)
            for i in range(n_qubits):
                qc.h(i)
            for i in range(n_qubits - 1):
                qc.cx(i, i + 1)
            qc.measure_all()

            # Initialize backends
            backend_cpu = QiskitBackend(backend_name='aer_simulator_statevector')
            backend_gpu = CuQuantumBackend(
                backend_name='cuquantum_statevector',
                simulation_method='statevector'
            )

            # Benchmark
            result = backend_gpu.benchmark_gpu_speedup(qc, n_runs=5)

            print(f"{n_qubits:<10} {result['cpu_time']*1000:>12.2f} ms "
                  f"{result['gpu_time']*1000:>12.2f} ms {result['speedup']:>8.2f}x")

        except Exception as e:
            print(f"{n_qubits:<10} Error: {e}")

    print("\n✓ Benchmark complete!")


def main():
    """Main execution."""
    print("""
╔═══════════════════════════════════════════════════════════════════╗
║                  cuQuantum GPU Acceleration Demo                  ║
║                   NVIDIA GPU-Accelerated VQE                       ║
╚═══════════════════════════════════════════════════════════════════╝
    """)

    # Check system
    if not check_system():
        print("\n⚠️  cuQuantum not available. Running in CPU-only mode.")
        print("To enable GPU acceleration, install cuQuantum:")
        print("  https://docs.nvidia.com/cuda/cuquantum/python/README.html")
        return

    # Run VQE example
    results = run_gpu_vqe_example()

    # Optional: Run scaling benchmark
    try:
        response = input("\n\nRun scaling benchmark? (y/n): ").lower()
        if response == 'y':
            benchmark_scaling()
    except:
        pass

    print("\n✅ Demo complete!")


if __name__ == "__main__":
    main()
