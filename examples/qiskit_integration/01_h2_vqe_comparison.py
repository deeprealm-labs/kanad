"""
H2 Molecule VQE: Classical vs Qiskit Backends Comparison

This example demonstrates running VQE on H2 molecule using:
1. Classical NumPy simulation (exact)
2. Qiskit Aer statevector simulator (exact)
3. Qiskit Aer QASM simulator (shot-based)

Compares accuracy and performance across backends.
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


def run_h2_vqe():
    """Run VQE on H2 molecule with different backends."""

    print("=" * 70)
    print("H2 Molecule VQE - Backend Comparison")
    print("=" * 70)

    # Create H2 molecule
    print("\n1. Setting up H2 molecule...")
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))  # 1.4 bohr ≈ 0.74 Å
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
    print(f"   Nuclear repulsion: {hamiltonian.nuclear_repulsion:.6f} Ha")

    # Create ansatz and mapper
    print("\n3. Creating ansatz and mapper...")
    n_qubits = hamiltonian.n_orbitals
    ansatz = RealAmplitudesAnsatz(
        n_qubits=n_qubits,
        n_electrons=hamiltonian.n_electrons,
        n_layers=2,
        entanglement='linear'
    )
    mapper = JordanWignerMapper()

    print(f"   Ansatz: {ansatz}")
    print(f"   Number of parameters: {ansatz.get_num_parameters()}")
    print(f"   Mapper: {mapper}")

    # Define backends to test
    backends = [
        ('classical', {}),
        ('aer_simulator_statevector', {}),
        ('aer_simulator', {'shots': 1024}),
        ('aer_simulator', {'shots': 10000}),
    ]

    results = {}

    # Run VQE with each backend
    for backend_name, backend_opts in backends:
        print(f"\n{'=' * 70}")
        print(f"Running VQE with backend: {backend_name}")
        if 'shots' in backend_opts:
            print(f"Shots: {backend_opts['shots']}")
        print('=' * 70)

        try:
            # Create VQE solver
            solver = VQESolver(
                hamiltonian=hamiltonian,
                ansatz=ansatz,
                mapper=mapper,
                backend=backend_name,
                optimizer='COBYLA',
                max_iterations=100,
                **backend_opts
            )

            print(f"\n{solver}")

            # Run optimization
            print("\nOptimizing...")
            start_time = time.time()
            result = solver.solve()
            elapsed_time = time.time() - start_time

            # Store results
            key = backend_name
            if 'shots' in backend_opts:
                key += f"_{backend_opts['shots']}"

            results[key] = {
                'energy': result['energy'],
                'time': elapsed_time,
                'iterations': result['iterations'],
                'converged': result['converged'],
                'backend': backend_name,
                'shots': backend_opts.get('shots', None)
            }

            # Print results
            print(f"\n{'─' * 70}")
            print(f"Results:")
            print(f"  Ground state energy: {result['energy']:.8f} Ha")
            print(f"  Converged: {result['converged']}")
            print(f"  Iterations: {result['iterations']}")
            print(f"  Time: {elapsed_time:.2f} seconds")
            print(f"{'─' * 70}")

        except Exception as e:
            print(f"\n❌ Error with backend {backend_name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Print comparison table
    print(f"\n\n{'=' * 70}")
    print("COMPARISON TABLE")
    print('=' * 70)
    print(f"{'Backend':<30} {'Energy (Ha)':<15} {'Time (s)':<12} {'Iterations':<12}")
    print('-' * 70)

    for key, res in results.items():
        backend_str = key
        print(f"{backend_str:<30} {res['energy']:>14.8f} {res['time']:>11.2f} {res['iterations']:>11}")

    # Energy comparison
    if 'classical' in results:
        ref_energy = results['classical']['energy']
        print(f"\n{'=' * 70}")
        print("ACCURACY ANALYSIS (vs Classical)")
        print('=' * 70)
        print(f"{'Backend':<30} {'Error (mHa)':<15}")
        print('-' * 70)

        for key, res in results.items():
            if key != 'classical':
                error_mha = (res['energy'] - ref_energy) * 1000  # Convert to mHa
                print(f"{key:<30} {error_mha:>14.4f}")

    print("\n✓ Comparison complete!")
    return results


if __name__ == "__main__":
    results = run_h2_vqe()
