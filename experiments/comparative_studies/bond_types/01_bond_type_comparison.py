"""
Experiment 1: Bond Type Comparison
===================================
Objective: Compare how different bond types are represented and computed in Kanad

Test Cases:
1. H2 - Pure covalent bond
2. LiH - Polar covalent/ionic bond
3. NaCl - Ionic bond
4. Compare auto-detection vs manual specification

Metrics:
- Energy accuracy
- Qubit count
- Circuit depth
- Computation time
- Representation type used
"""

import time
import json
from datetime import datetime
from kanad import BondFactory
from kanad.solvers import VQESolver

def run_bond_experiment(atom1, atom2, distance=None, bond_type='auto', basis='sto-3g'):
    """Run a single bond experiment and collect metrics"""
    print(f"\n{'='*60}")
    print(f"Testing: {atom1}-{atom2} | Type: {bond_type} | Basis: {basis}")
    print(f"{'='*60}")

    results = {
        'molecule': f"{atom1}-{atom2}",
        'specified_bond_type': bond_type,
        'basis': basis,
        'distance': distance,
        'timestamp': datetime.now().isoformat()
    }

    try:
        # Create bond
        start_time = time.time()
        if bond_type == 'auto':
            bond = BondFactory.create_bond(atom1, atom2, distance=distance, basis=basis)
        else:
            bond = BondFactory.create_bond(atom1, atom2, bond_type=bond_type, distance=distance, basis=basis)

        bond_creation_time = time.time() - start_time

        # Get bond info
        bond_info = BondFactory.quick_bond_info(atom1, atom2)

        results.update({
            'bond_creation_time': bond_creation_time,
            'detected_bond_type': str(bond.__class__.__name__),
            'predicted_type': bond_info['predicted_type'],
            'en_difference': bond_info['electronegativity_difference'],
            'estimated_distance': bond_info.get('estimated_bond_length', 'N/A'),
            'actual_distance': bond.molecule.distance if hasattr(bond.molecule, 'distance') else None
        })

        print(f"✓ Bond created: {results['detected_bond_type']}")
        print(f"  Predicted type: {results['predicted_type']}")
        print(f"  ΔEN: {results['en_difference']:.3f}")
        print(f"  Distance: {results['actual_distance']} Å")

        # Get representation info
        repr_type = str(bond.representation.__class__.__name__)
        n_qubits = bond.representation.num_qubits if hasattr(bond.representation, 'num_qubits') else 'N/A'

        results.update({
            'representation_type': repr_type,
            'num_qubits': n_qubits
        })

        print(f"  Representation: {repr_type}")
        print(f"  Qubits needed: {n_qubits}")

        # Create solver with statevector backend (fast and exact)
        print(f"\nSetting up VQE solver...")
        solver = VQESolver(
            bond,
            ansatz_type='governance',  # Use governance-aware ansatz
            backend='statevector',
            max_iterations=100,
            enable_analysis=True
        )

        # Solve
        print(f"Running VQE computation...")
        start_time = time.time()
        vqe_result = solver.solve()
        computation_time = time.time() - start_time

        results.update({
            'energy_hartree': vqe_result.get('energy'),
            'computation_time': computation_time,
            'converged': vqe_result.get('converged', False),
            'n_iterations': vqe_result.get('n_iterations', 'N/A'),
            'circuit_depth': vqe_result.get('circuit_depth', 'N/A'),
            'gate_count': vqe_result.get('gate_count', 'N/A')
        })

        # Analysis results
        if 'analysis' in vqe_result:
            analysis = vqe_result['analysis']
            results.update({
                'charge_transfer': analysis.get('charge_transfer', 'N/A'),
                'bond_order': analysis.get('bond_order', 'N/A'),
                'correlation_energy': analysis.get('correlation_energy', 'N/A')
            })

        print(f"\n✓ VQE Completed!")
        print(f"  Energy: {results['energy_hartree']:.6f} Ha")
        print(f"  Time: {computation_time:.2f} s")
        print(f"  Iterations: {results['n_iterations']}")
        print(f"  Circuit depth: {results['circuit_depth']}")

        if 'analysis' in vqe_result:
            print(f"\n  Analysis:")
            print(f"    Charge transfer: {results.get('charge_transfer', 'N/A')}")
            print(f"    Bond order: {results.get('bond_order', 'N/A')}")

        results['status'] = 'success'

    except Exception as e:
        print(f"\n✗ Error: {str(e)}")
        results['status'] = 'failed'
        results['error'] = str(e)
        import traceback
        results['traceback'] = traceback.format_exc()

    return results


def main():
    """Run comparative bond type experiments"""
    print("\n" + "="*60)
    print("EXPERIMENT 1: BOND TYPE COMPARISON")
    print("="*60)
    print("\nObjective: Compare different bond types in Kanad framework")
    print("Using: BondFactory, VQESolver, Statevector backend")

    all_results = []

    # Experiment Set 1: Pure Covalent - H2
    print("\n\n" + "="*60)
    print("SET 1: PURE COVALENT BOND - H2")
    print("="*60)

    # Auto-detection
    result = run_bond_experiment('H', 'H', distance=0.74, bond_type='auto')
    all_results.append(result)

    # Manual specification
    result = run_bond_experiment('H', 'H', distance=0.74, bond_type='covalent')
    all_results.append(result)


    # Experiment Set 2: Polar Covalent - LiH
    print("\n\n" + "="*60)
    print("SET 2: POLAR COVALENT/IONIC - LiH")
    print("="*60)

    # Auto-detection
    result = run_bond_experiment('Li', 'H', distance=1.6, bond_type='auto')
    all_results.append(result)

    # Force as ionic
    result = run_bond_experiment('Li', 'H', distance=1.6, bond_type='ionic')
    all_results.append(result)

    # Force as covalent
    result = run_bond_experiment('Li', 'H', distance=1.6, bond_type='covalent')
    all_results.append(result)


    # Experiment Set 3: Ionic - NaCl (if computational resources allow)
    print("\n\n" + "="*60)
    print("SET 3: IONIC BOND - NaCl")
    print("="*60)

    # Auto-detection
    result = run_bond_experiment('Na', 'Cl', bond_type='auto')
    all_results.append(result)


    # Save results
    output_file = 'experiments/comparative_studies/bond_types/results_bond_comparison.json'
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)

    print("\n\n" + "="*60)
    print("EXPERIMENT COMPLETE!")
    print("="*60)
    print(f"\nResults saved to: {output_file}")

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    successful = [r for r in all_results if r['status'] == 'success']
    failed = [r for r in all_results if r['status'] == 'failed']

    print(f"\nTotal experiments: {len(all_results)}")
    print(f"Successful: {len(successful)}")
    print(f"Failed: {len(failed)}")

    if successful:
        print("\n" + "-"*60)
        print("Successful Experiments:")
        print("-"*60)
        print(f"{'Molecule':<15} {'Type':<15} {'Repr':<20} {'Qubits':<8} {'Energy (Ha)':<15}")
        print("-"*60)
        for r in successful:
            print(f"{r['molecule']:<15} {r['detected_bond_type']:<15} {r['representation_type']:<20} {str(r['num_qubits']):<8} {r['energy_hartree']:<15.6f}")

    if failed:
        print("\n" + "-"*60)
        print("Failed Experiments:")
        print("-"*60)
        for r in failed:
            print(f"  {r['molecule']} ({r['specified_bond_type']}): {r.get('error', 'Unknown error')}")

    return all_results


if __name__ == "__main__":
    results = main()
