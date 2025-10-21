"""
Experiment 3: Mapper Efficiency Tests
======================================
Objective: Compare different fermion-to-qubit mapping schemes

Mappers to test:
1. Jordan-Wigner (sequential, local)
2. Bravyi-Kitaev (tree-based, logarithmic)
3. Hybrid Orbital (for covalent bonds, MO pairs)

Metrics:
- Qubit count
- Operator weight (Pauli string length)
- Energy accuracy
- Circuit complexity
- Computation time
"""

import time
import json
from datetime import datetime
from kanad import BondFactory
from kanad.solvers import VQESolver

def benchmark_mapper(bond, mapper_type, ansatz_type='ucc', backend='statevector'):
    """Benchmark a specific mapper"""
    print(f"\n{'='*60}")
    print(f"Benchmarking Mapper: {mapper_type.upper()}")
    print(f"{'='*60}")

    results = {
        'mapper_type': mapper_type,
        'ansatz_type': ansatz_type,
        'backend': backend,
        'timestamp': datetime.now().isoformat()
    }

    try:
        # Create solver with specific mapper
        print(f"Creating VQE solver with {mapper_type} mapper...")
        start_time = time.time()

        solver = VQESolver(
            bond,
            ansatz_type=ansatz_type,
            mapper_type=mapper_type,
            backend=backend,
            max_iterations=200,
            enable_analysis=True
        )

        solver_creation_time = time.time() - start_time
        results['solver_creation_time'] = solver_creation_time

        print(f"✓ Solver created in {solver_creation_time:.3f} s")

        # Get mapper info
        if hasattr(solver, 'mapper'):
            mapper = solver.mapper
            mapper_class = mapper.__class__.__name__
            results['mapper_class'] = mapper_class
            print(f"  Mapper class: {mapper_class}")

        # Get Hamiltonian operator info
        if hasattr(bond, 'hamiltonian'):
            hamiltonian = bond.hamiltonian
            if hasattr(hamiltonian, 'qubit_hamiltonian'):
                qubit_ham = hamiltonian.qubit_hamiltonian
                if qubit_ham is not None:
                    # Try to get operator statistics
                    if hasattr(qubit_ham, 'num_qubits'):
                        results['hamiltonian_qubits'] = qubit_ham.num_qubits
                        print(f"  Hamiltonian qubits: {qubit_ham.num_qubits}")

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
            'gate_count': vqe_result.get('gate_count', 'N/A'),
            'num_parameters': vqe_result.get('num_parameters', 'N/A')
        })

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
    """Run mapper efficiency experiments"""
    print("\n" + "="*60)
    print("EXPERIMENT 3: MAPPER EFFICIENCY TESTS")
    print("="*60)
    print("\nObjective: Compare different fermion-to-qubit mappers")

    all_results = []

    # Test Set 1: Covalent bond (H2) - all mappers should work
    print("\n\n" + "="*60)
    print("TEST SET 1: COVALENT BOND - H2")
    print("="*60)

    bond_h2 = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
    print(f"✓ Bond created: {bond_h2.__class__.__name__}")
    print(f"  Representation: {bond_h2.representation.__class__.__name__}")
    print(f"  Qubits: {bond_h2.representation.num_qubits if hasattr(bond_h2.representation, 'num_qubits') else 'N/A'}")

    reference_energy_h2 = -1.174476
    print(f"  Reference energy: {reference_energy_h2:.6f} Ha")

    # Test all mapper types on H2
    mapper_types = ['jordan_wigner', 'bravyi_kitaev']

    # Note: hybrid_orbital mapper is specifically for covalent bonds
    # Let's try it as well
    try:
        mapper_types.append('hybrid_orbital')
    except:
        pass

    for mapper_type in mapper_types:
        result = benchmark_mapper(bond_h2, mapper_type, ansatz_type='ucc')
        result['molecule'] = 'H2'
        result['bond_class'] = bond_h2.__class__.__name__

        # Calculate error
        if result['status'] == 'success' and result['energy_hartree'] is not None:
            error = abs(result['energy_hartree'] - reference_energy_h2)
            error_mhartree = error * 1000
            result['energy_error_hartree'] = error
            result['energy_error_mhartree'] = error_mhartree
            print(f"  Error: {error_mhartree:.3f} mHa")

        all_results.append(result)

    # Test Set 2: Ionic bond (LiH) - test mapper behavior
    print("\n\n" + "="*60)
    print("TEST SET 2: IONIC BOND - LiH")
    print("="*60)

    bond_lih = BondFactory.create_bond('Li', 'H', distance=1.6, basis='sto-3g')
    print(f"✓ Bond created: {bond_lih.__class__.__name__}")
    print(f"  Representation: {bond_lih.representation.__class__.__name__}")
    print(f"  Qubits: {bond_lih.representation.num_qubits if hasattr(bond_lih.representation, 'num_qubits') else 'N/A'}")

    # LiH reference energy (approximate, sto-3g)
    reference_energy_lih = -7.98

    for mapper_type in ['jordan_wigner', 'bravyi_kitaev']:
        result = benchmark_mapper(bond_lih, mapper_type, ansatz_type='governance')
        result['molecule'] = 'LiH'
        result['bond_class'] = bond_lih.__class__.__name__

        # Calculate error (relative to reference)
        if result['status'] == 'success' and result['energy_hartree'] is not None:
            error = abs(result['energy_hartree'] - reference_energy_lih)
            result['energy_error_hartree'] = error
            print(f"  Deviation from reference: {error:.3f} Ha")

        all_results.append(result)

    # Save results
    output_file = 'experiments/comparative_studies/mapper_efficiency/results_mapper_efficiency.json'
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
        # Group by molecule
        molecules = list(set(r['molecule'] for r in successful))

        for molecule in molecules:
            mol_results = [r for r in successful if r['molecule'] == molecule]

            print(f"\n{molecule} Results:")
            print("-"*100)
            print(f"{'Mapper':<20} {'Energy (Ha)':<18} {'Error':<15} {'Time (s)':<12} {'Depth':<10} {'Gates':<10}")
            print("-"*100)

            for r in mol_results:
                energy = r.get('energy_hartree', 'N/A')
                error = r.get('energy_error_mhartree', r.get('energy_error_hartree', 'N/A'))
                time_val = r.get('computation_time', 'N/A')
                depth = r.get('circuit_depth', 'N/A')
                gates = r.get('gate_count', 'N/A')

                energy_str = f"{energy:.8f}" if isinstance(energy, (int, float)) else str(energy)

                # Format error based on type
                if isinstance(error, (int, float)):
                    if error < 1:  # Likely in Hartree
                        error_str = f"{error:.3f} Ha"
                    else:  # Likely in mHa
                        error_str = f"{error:.3f} mHa"
                else:
                    error_str = str(error)

                time_str = f"{time_val:.2f}" if isinstance(time_val, (int, float)) else str(time_val)

                print(f"{r['mapper_type']:<20} {energy_str:<18} {error_str:<15} {time_str:<12} {str(depth):<10} {str(gates):<10}")

    failed = [r for r in all_results if r['status'] == 'failed']
    if failed:
        print("\n" + "-"*60)
        print("Failed experiments:")
        print("-"*60)
        for r in failed:
            print(f"  {r['molecule']} with {r['mapper_type']}: {r.get('error', 'Unknown error')}")

    # Key insights
    print("\n" + "="*60)
    print("KEY INSIGHTS")
    print("="*60)
    print("\nMapper Characteristics:")
    print("  Jordan-Wigner: Sequential encoding, local operators, O(n) weight")
    print("  Bravyi-Kitaev: Tree encoding, non-local operators, O(log n) weight")
    print("  Hybrid Orbital: MO pair encoding, specialized for covalent bonds")

    return all_results


if __name__ == "__main__":
    results = main()
