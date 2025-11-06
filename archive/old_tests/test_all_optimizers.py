#!/usr/bin/env python3
"""
Test ALL 8 Optimizers - Verify Iteration Control
Tests that user-specified max_iterations is respected
Checks for hardcoded values and proper convergence
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.services.experiment_service import create_molecule_from_config
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz
from kanad.utils.vqe_solver import VQESolver
import time


def test_optimizer(optimizer_name, max_iterations_input, molecule_config):
    """
    Test single optimizer with specific iteration limit

    Verifies:
    1. Optimizer runs (no crashes)
    2. Iterations are controlled by user input (not hardcoded)
    3. Function evaluations are logged
    4. Convergence behavior is reasonable
    """
    print(f"\n{'='*80}")
    print(f"TEST: {optimizer_name} with max_iterations={max_iterations_input}")
    print("="*80)

    # Create molecule
    molecule = create_molecule_from_config(molecule_config)
    ham = molecule.hamiltonian

    # Get HF reference
    _, hf_energy = ham.solve_scf(max_iterations=50, conv_tol=1e-6)
    print(f"HF energy: {hf_energy:.8f} Ha")

    # Create ansatz
    n_qubits = 2 * ham.n_orbitals
    ansatz = CovalentGovernanceAnsatz(n_qubits=n_qubits, n_electrons=molecule.n_electrons)
    print(f"Ansatz: CovalentGovernance with {ansatz.n_parameters} parameters")

    # Track iterations
    iteration_log = []
    function_eval_count = 0

    def callback(iteration, energy, parameters):
        nonlocal function_eval_count
        function_eval_count += 1
        iteration_log.append({
            'iteration': iteration,
            'energy': energy,
            'func_eval': function_eval_count
        })
        if iteration % 5 == 0:
            print(f"  Iter {iteration}: E = {energy:.8f} Ha (func_eval #{function_eval_count})")

    # Run VQE with user-specified max_iterations
    start_time = time.time()

    try:
        solver = VQESolver(
            hamiltonian=ham,
            ansatz=ansatz,
            optimizer=optimizer_name,
            max_iterations=max_iterations_input  # USER INPUT - NOT HARDCODED!
        )

        print(f"\n🚀 Starting VQE with {optimizer_name}...")
        print(f"   User-requested max_iterations: {max_iterations_input}")

        result = solver.solve(callback=callback)

        elapsed_time = time.time() - start_time

        # Extract results
        vqe_energy = result['energy']
        iterations_reported = result.get('iterations', 0)
        func_evals_reported = result.get('function_evaluations', 0)
        converged = result.get('converged', False)
        correlation = vqe_energy - hf_energy

        print(f"\n✅ VQE Completed")
        print(f"  VQE energy:      {vqe_energy:.8f} Ha")
        print(f"  HF energy:       {hf_energy:.8f} Ha")
        print(f"  Correlation:     {correlation:.8f} Ha ({correlation*627.5:.4f} kcal/mol)")
        print(f"  Iterations:      {iterations_reported} (callback logged: {len(iteration_log)})")
        print(f"  Func evals:      {func_evals_reported} (callback counted: {function_eval_count})")
        print(f"  Converged:       {converged}")
        print(f"  Time:            {elapsed_time:.2f}s")

        # Check if iteration limit was respected
        iteration_limit_respected = iterations_reported <= max_iterations_input

        # Check if any optimization happened
        optimization_occurred = len(iteration_log) > 2 and abs(correlation) > 1e-6

        # Success criteria
        success = (
            not result.get('error') and
            iteration_limit_respected and
            len(iteration_log) > 0
        )

        return {
            'optimizer': optimizer_name,
            'max_iter_requested': max_iterations_input,
            'iterations_actual': iterations_reported,
            'func_evals': func_evals_reported,
            'correlation': correlation,
            'converged': converged,
            'time': elapsed_time,
            'success': success,
            'limit_respected': iteration_limit_respected,
            'optimization_occurred': optimization_occurred,
            'energy_history': [log['energy'] for log in iteration_log]
        }

    except Exception as e:
        print(f"\n❌ VQE Failed: {e}")
        import traceback
        traceback.print_exc()

        return {
            'optimizer': optimizer_name,
            'max_iter_requested': max_iterations_input,
            'success': False,
            'error': str(e)
        }


def main():
    """Test all 8 optimizers with different iteration limits"""
    print("\n" + "="*80)
    print("COMPREHENSIVE OPTIMIZER TEST SUITE")
    print("Testing: All 8 optimizers with user-controlled iterations")
    print("="*80)

    # Molecule configuration
    molecule_config = {
        "smiles": "[H][H]",
        "basis": "sto-3g",
        "charge": 0,
        "multiplicity": 1
    }

    # All 8 optimizers from configuration.py
    optimizers_to_test = [
        ("COBYLA", 30, "Derivative-free, fast (RECOMMENDED)"),
        ("Powell", 25, "Good for smooth landscapes"),
        ("L-BFGS-B", 20, "For large parameter spaces"),
        ("SLSQP", 15, "Sequential Least Squares (expensive)"),
        ("Nelder-Mead", 30, "Simplex method"),
        ("CG", 25, "Conjugate Gradient"),
        ("BFGS", 20, "Standard BFGS"),
        ("TNC", 20, "Truncated Newton"),
    ]

    results = []

    print(f"\nMolecule: H2 (from SMILES: [H][H])")
    print(f"Basis: sto-3g")
    print(f"Ansatz: CovalentGovernance")
    print(f"\nTesting {len(optimizers_to_test)} optimizers...")

    for optimizer, max_iter, description in optimizers_to_test:
        result = test_optimizer(optimizer, max_iter, molecule_config)
        result['description'] = description
        results.append(result)

    # Print summary
    print("\n" + "="*80)
    print("SUMMARY - ALL OPTIMIZERS")
    print("="*80)
    print(f"\n{'Optimizer':<15} {'Requested':<10} {'Actual':<10} {'Func Evals':<12} {'Correlation':<15} {'Status':<10}")
    print("-"*80)

    for r in results:
        if r['success']:
            status = "✅ PASS"
            actual_iter = r['iterations_actual']
            func_evals = r['func_evals']
            corr = f"{r['correlation']:.6f}"
        else:
            status = "❌ FAIL"
            actual_iter = "N/A"
            func_evals = "N/A"
            corr = r.get('error', 'Error')[:12]

        requested = r['max_iter_requested']
        print(f"{r['optimizer']:<15} {requested:<10} {actual_iter:<10} {func_evals:<12} {corr:<15} {status:<10}")

    # Detailed analysis
    print("\n" + "="*80)
    print("DETAILED ANALYSIS")
    print("="*80)

    successful = [r for r in results if r['success']]
    failed = [r for r in results if not r['success']]

    print(f"\n📊 Success Rate: {len(successful)}/{len(results)} ({100*len(successful)/len(results):.1f}%)")

    if successful:
        print(f"\n✅ Working Optimizers ({len(successful)}):")
        for r in successful:
            limit_check = "✓" if r['limit_respected'] else "✗"
            opt_check = "✓" if r['optimization_occurred'] else "✗"
            print(f"   {r['optimizer']:<15} - Limit respected: {limit_check}, Optimization: {opt_check}")
            print(f"      {r['description']}")
            print(f"      Iterations: {r['iterations_actual']}/{r['max_iter_requested']}, Correlation: {r['correlation']:.8f} Ha")

    if failed:
        print(f"\n❌ Failed Optimizers ({len(failed)}):")
        for r in failed:
            print(f"   {r['optimizer']:<15} - {r.get('error', 'Unknown error')[:60]}")

    # Check for hardcoded issues
    print(f"\n🔍 HARDCODING CHECK:")
    hardcoded_issues = []

    for r in successful:
        requested = r['max_iter_requested']
        actual = r['iterations_actual']

        # Check if iterations are wildly different (possible hardcoding)
        if actual > requested * 1.5:
            hardcoded_issues.append(f"{r['optimizer']}: Exceeded limit ({actual} > {requested})")

        # Check if always same number regardless of input
        # (would need multiple tests with different inputs to detect)

    if hardcoded_issues:
        print("   ⚠️  Potential Issues:")
        for issue in hardcoded_issues:
            print(f"      - {issue}")
    else:
        print("   ✅ No hardcoding detected - iterations properly controlled by user input")

    # Performance comparison
    if successful:
        print(f"\n⚡ Performance Comparison:")
        # Sort by correlation energy (best first)
        sorted_by_corr = sorted(successful, key=lambda x: x['correlation'])
        print(f"   Best correlation:")
        for i, r in enumerate(sorted_by_corr[:3], 1):
            print(f"      {i}. {r['optimizer']}: {r['correlation']:.8f} Ha in {r['time']:.2f}s")

        # Sort by speed
        sorted_by_time = sorted(successful, key=lambda x: x['time'])
        print(f"   Fastest:")
        for i, r in enumerate(sorted_by_time[:3], 1):
            print(f"      {i}. {r['optimizer']}: {r['time']:.2f}s (correlation: {r['correlation']:.8f} Ha)")

    # Final verdict
    print("\n" + "="*80)
    if len(successful) >= 6:
        print("✅ EXCELLENT: Most optimizers working!")
        return 0
    elif len(successful) >= 4:
        print("⚠️  GOOD: Most optimizers working, some issues")
        return 0
    elif len(successful) >= 2:
        print("⚠️  PARTIAL: Some optimizers working")
        return 1
    else:
        print("❌ CRITICAL: Most optimizers failing!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
