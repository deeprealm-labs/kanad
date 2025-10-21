"""
Experiment 2: Ansatz Benchmarking
==================================
Objective: Compare different ansatz types on the same molecule

Test Cases:
- Use H2 as test molecule (well-studied reference)
- Compare ansatz types:
  1. UCC (Unitary Coupled Cluster)
  2. Hardware Efficient
  3. Governance-Aware (physics-driven)

Metrics:
- Energy accuracy
- Circuit depth
- Gate count
- Parameter count
- Convergence speed
- Computation time
"""

import time
import json
from datetime import datetime
from kanad import BondFactory
from kanad.solvers import VQESolver

def benchmark_ansatz(bond, ansatz_type, backend='statevector', max_iterations=200):
    """Benchmark a specific ansatz type"""
    print(f"\n{'='*60}")
    print(f"Benchmarking Ansatz: {ansatz_type.upper()}")
    print(f"{'='*60}")

    results = {
        'ansatz_type': ansatz_type,
        'backend': backend,
        'max_iterations': max_iterations,
        'timestamp': datetime.now().isoformat()
    }

    try:
        # Create solver
        print(f"Creating VQE solver with {ansatz_type} ansatz...")
        start_time = time.time()

        solver = VQESolver(
            bond,
            ansatz_type=ansatz_type,
            backend=backend,
            max_iterations=max_iterations,
            enable_analysis=True
        )

        solver_creation_time = time.time() - start_time
        results['solver_creation_time'] = solver_creation_time

        print(f"✓ Solver created in {solver_creation_time:.3f} s")

        # Get circuit info before optimization
        if hasattr(solver, 'ansatz'):
            ansatz = solver.ansatz
            if hasattr(ansatz, 'num_parameters'):
                results['num_parameters'] = ansatz.num_parameters
                print(f"  Parameters: {ansatz.num_parameters}")

        # Run VQE
        print(f"\nRunning VQE optimization...")
        start_time = time.time()
        vqe_result = solver.solve()
        computation_time = time.time() - start_time

        results.update({
            'computation_time': computation_time,
            'energy_hartree': vqe_result.get('energy'),
            'converged': vqe_result.get('converged', False),
            'n_iterations': vqe_result.get('n_iterations', 'N/A'),
            'circuit_depth': vqe_result.get('circuit_depth', 'N/A'),
            'gate_count': vqe_result.get('gate_count', 'N/A')
        })

        # Get optimization history if available
        if 'optimization_history' in vqe_result:
            history = vqe_result['optimization_history']
            results['optimization_history'] = history
            results['convergence_steps'] = len(history)

        # Analysis
        if 'analysis' in vqe_result:
            analysis = vqe_result['analysis']
            results.update({
                'charge_transfer': analysis.get('charge_transfer', 'N/A'),
                'bond_order': analysis.get('bond_order', 'N/A'),
                'correlation_energy': analysis.get('correlation_energy', 'N/A')
            })

        print(f"\n✓ VQE Completed!")
        print(f"  Energy: {results['energy_hartree']:.8f} Ha")
        print(f"  Converged: {results['converged']}")
        print(f"  Iterations: {results['n_iterations']}")
        print(f"  Time: {computation_time:.2f} s")
        print(f"  Circuit depth: {results['circuit_depth']}")
        print(f"  Gate count: {results['gate_count']}")

        results['status'] = 'success'

    except Exception as e:
        print(f"\n✗ Error: {str(e)}")
        results['status'] = 'failed'
        results['error'] = str(e)
        import traceback
        results['traceback'] = traceback.format_exc()

    return results


def main():
    """Run ansatz benchmarking experiments"""
    print("\n" + "="*60)
    print("EXPERIMENT 2: ANSATZ BENCHMARKING")
    print("="*60)
    print("\nObjective: Compare different ansatz types")
    print("Test molecule: H2 at equilibrium (0.74 Å)")

    all_results = []

    # Create H2 bond (reference molecule)
    print("\n" + "="*60)
    print("Setting up test molecule: H2")
    print("="*60)

    bond = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
    print(f"✓ Bond created: {bond.__class__.__name__}")
    print(f"  Representation: {bond.representation.__class__.__name__}")
    print(f"  Qubits: {bond.representation.num_qubits if hasattr(bond.representation, 'num_qubits') else 'N/A'}")

    # Reference energy (experimental)
    # H2 experimental: -1.174476 Ha (sto-3g basis)
    reference_energy = -1.174476
    print(f"  Reference energy: {reference_energy:.6f} Ha")

    # Test different ansatze
    ansatz_types = ['ucc', 'hardware_efficient', 'governance']

    for ansatz_type in ansatz_types:
        result = benchmark_ansatz(bond, ansatz_type, backend='statevector', max_iterations=200)

        # Calculate error if successful
        if result['status'] == 'success' and result['energy_hartree'] is not None:
            error = abs(result['energy_hartree'] - reference_energy)
            error_mhartree = error * 1000
            result['energy_error_hartree'] = error
            result['energy_error_mhartree'] = error_mhartree
            print(f"  Error: {error_mhartree:.3f} mHa")

        all_results.append(result)

    # Save results
    output_file = 'experiments/comparative_studies/ansatz_benchmark/results_ansatz_benchmark.json'
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)

    print("\n\n" + "="*60)
    print("EXPERIMENT COMPLETE!")
    print("="*60)
    print(f"\nResults saved to: {output_file}")

    # Print comparative summary
    print("\n" + "="*60)
    print("COMPARATIVE SUMMARY")
    print("="*60)

    successful = [r for r in all_results if r['status'] == 'success']

    if successful:
        print(f"\nReference energy: {reference_energy:.8f} Ha")
        print("\n" + "-"*100)
        print(f"{'Ansatz':<20} {'Energy (Ha)':<18} {'Error (mHa)':<15} {'Time (s)':<12} {'Depth':<10} {'Gates':<10}")
        print("-"*100)

        for r in successful:
            energy = r.get('energy_hartree', 'N/A')
            error = r.get('energy_error_mhartree', 'N/A')
            time_val = r.get('computation_time', 'N/A')
            depth = r.get('circuit_depth', 'N/A')
            gates = r.get('gate_count', 'N/A')

            energy_str = f"{energy:.8f}" if isinstance(energy, (int, float)) else str(energy)
            error_str = f"{error:.3f}" if isinstance(error, (int, float)) else str(error)
            time_str = f"{time_val:.2f}" if isinstance(time_val, (int, float)) else str(time_val)

            print(f"{r['ansatz_type']:<20} {energy_str:<18} {error_str:<15} {time_str:<12} {str(depth):<10} {str(gates):<10}")

        # Find best ansatz by accuracy
        if all(r.get('energy_error_mhartree') is not None for r in successful):
            best = min(successful, key=lambda x: x['energy_error_mhartree'])
            fastest = min(successful, key=lambda x: x['computation_time'])

            print("\n" + "-"*60)
            print(f"Most accurate: {best['ansatz_type']} ({best['energy_error_mhartree']:.3f} mHa)")
            print(f"Fastest: {fastest['ansatz_type']} ({fastest['computation_time']:.2f} s)")

    failed = [r for r in all_results if r['status'] == 'failed']
    if failed:
        print("\n" + "-"*60)
        print("Failed experiments:")
        print("-"*60)
        for r in failed:
            print(f"  {r['ansatz_type']}: {r.get('error', 'Unknown error')}")

    return all_results


if __name__ == "__main__":
    results = main()
