"""
Comprehensive Ansatz Efficiency Diagnostic

Goal: Find out WHY ansatze need 600-700+ function evaluations
and HOW to reduce this while maintaining accuracy.
"""

import sys
import numpy as np
from typing import Dict, List
import time
import logging

# Configure logging
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

def test_parameter_initialization():
    """Test 1: Does parameter initialization matter?"""
    print("="*80)
    print("TEST 1: PARAMETER INITIALIZATION IMPACT")
    print("="*80)

    from kanad.core.molecule import Molecule
    from kanad.core.hamiltonians.openfermion_jw import build_h2_hamiltonian
    from kanad.ansatze import CovalentGovernanceAnsatz
    from kanad.utils.vqe_solver import VQESolver

    # Create H2 molecule
    bond = Molecule.from_smiles("[H][H]", name="H2")
    bond.optimize_geometry()

    H = build_h2_hamiltonian(bond)
    n_qubits = 4
    n_electrons = 2

    print(f"\nMolecule: H2")
    print(f"HF energy: {bond.hf_energy:.8f} Ha")
    print(f"FCI energy: {bond.fci_energy:.8f} Ha")
    print(f"Target correlation: {bond.fci_energy - bond.hf_energy:.8f} Ha")

    ansatz = CovalentGovernanceAnsatz(n_qubits=n_qubits, n_electrons=n_electrons, n_layers=2)
    print(f"\nAnsatz: {ansatz.n_parameters} parameters")

    # Test different initializations
    inits = {
        "Zero": np.zeros(ansatz.n_parameters),
        "Small random": np.random.randn(ansatz.n_parameters) * 0.01,
        "Medium random": np.random.randn(ansatz.n_parameters) * 0.1,
        "Large random": np.random.randn(ansatz.n_parameters) * 1.0,
        "HF-like": np.zeros(ansatz.n_parameters),  # Will test if starting at HF helps
    }

    results = {}

    for name, initial_params in inits.items():
        print(f"\n{'─'*80}")
        print(f"Testing: {name} initialization")
        print(f"{'─'*80}")

        try:
            solver = VQESolver(
                hamiltonian=H,
                ansatz=ansatz,
                optimizer='SLSQP',
                max_iterations=200,
                backend='statevector',
                initial_parameters=initial_params
            )

            start_time = time.time()
            result = solver.solve()
            elapsed = time.time() - start_time

            energy = result['energy']
            n_evals = result.get('n_function_evals', 0)

            recovery = 100 * (energy - bond.hf_energy) / (bond.fci_energy - bond.hf_energy)

            print(f"✅ Energy: {energy:.8f} Ha")
            print(f"   Correlation recovery: {recovery:.1f}%")
            print(f"   Function evals: {n_evals}")
            print(f"   Time: {elapsed:.2f}s")

            results[name] = {
                'energy': energy,
                'n_evals': n_evals,
                'recovery': recovery,
                'time': elapsed
            }

        except Exception as e:
            print(f"❌ Failed: {e}")
            results[name] = {'error': str(e)}

    print(f"\n{'='*80}")
    print("INITIALIZATION SUMMARY")
    print(f"{'='*80}")
    for name, res in results.items():
        if 'error' not in res:
            print(f"{name:20s}: {res['n_evals']:4d} evals, {res['recovery']:5.1f}% recovery, {res['time']:.2f}s")
        else:
            print(f"{name:20s}: ERROR")

    return results


def test_ansatz_depth_vs_accuracy():
    """Test 2: Does circuit depth impact optimization speed?"""
    print("\n" + "="*80)
    print("TEST 2: ANSATZ DEPTH VS CONVERGENCE SPEED")
    print("="*80)

    from kanad.core.molecule import Molecule
    from kanad.core.hamiltonians.openfermion_jw import build_h2_hamiltonian
    from kanad.ansatze import CovalentGovernanceAnsatz
    from kanad.utils.vqe_solver import VQESolver

    # Create H2
    bond = Molecule.from_smiles("[H][H]", name="H2")
    bond.optimize_geometry()
    H = build_h2_hamiltonian(bond)

    n_qubits = 4
    n_electrons = 2

    print(f"\nTarget energy: {bond.fci_energy:.8f} Ha")

    # Test different layer counts
    results = {}
    for n_layers in [1, 2, 3]:
        print(f"\n{'─'*80}")
        print(f"Testing: {n_layers} layers")
        print(f"{'─'*80}")

        try:
            ansatz = CovalentGovernanceAnsatz(
                n_qubits=n_qubits,
                n_electrons=n_electrons,
                n_layers=n_layers
            )

            print(f"Parameters: {ansatz.n_parameters}")

            solver = VQESolver(
                hamiltonian=H,
                ansatz=ansatz,
                optimizer='SLSQP',
                max_iterations=200,
                backend='statevector'
            )

            start_time = time.time()
            result = solver.solve()
            elapsed = time.time() - start_time

            energy = result['energy']
            n_evals = result.get('n_function_evals', 0)
            recovery = 100 * (energy - bond.hf_energy) / (bond.fci_energy - bond.hf_energy)

            print(f"✅ Energy: {energy:.8f} Ha")
            print(f"   Correlation recovery: {recovery:.1f}%")
            print(f"   Function evals: {n_evals}")
            print(f"   Evals per parameter: {n_evals / ansatz.n_parameters:.1f}")
            print(f"   Time: {elapsed:.2f}s")

            results[n_layers] = {
                'energy': energy,
                'n_evals': n_evals,
                'n_params': ansatz.n_parameters,
                'recovery': recovery,
                'evals_per_param': n_evals / ansatz.n_parameters
            }

        except Exception as e:
            print(f"❌ Failed: {e}")
            results[n_layers] = {'error': str(e)}

    print(f"\n{'='*80}")
    print("DEPTH ANALYSIS SUMMARY")
    print(f"{'='*80}")
    for layers, res in results.items():
        if 'error' not in res:
            print(f"{layers} layers: {res['n_params']:2d} params, {res['n_evals']:4d} evals ({res['evals_per_param']:.1f}/param), {res['recovery']:5.1f}% recovery")

    return results


def test_optimizer_efficiency():
    """Test 3: Which optimizer is most efficient (evals vs accuracy)?"""
    print("\n" + "="*80)
    print("TEST 3: OPTIMIZER EFFICIENCY COMPARISON")
    print("="*80)

    from kanad.core.molecule import Molecule
    from kanad.core.hamiltonians.openfermion_jw import build_h2_hamiltonian
    from kanad.ansatze import CovalentGovernanceAnsatz
    from kanad.utils.vqe_solver import VQESolver

    # Create H2
    bond = Molecule.from_smiles("[H][H]", name="H2")
    bond.optimize_geometry()
    H = build_h2_hamiltonian(bond)

    n_qubits = 4
    n_electrons = 2
    ansatz = CovalentGovernanceAnsatz(n_qubits=n_qubits, n_electrons=n_electrons, n_layers=2)

    print(f"\nAnsatz: {ansatz.n_parameters} parameters")
    print(f"Target FCI: {bond.fci_energy:.8f} Ha")

    optimizers = ['SLSQP', 'L-BFGS-B', 'TNC', 'trust-constr']

    results = {}
    for opt in optimizers:
        print(f"\n{'─'*80}")
        print(f"Testing: {opt}")
        print(f"{'─'*80}")

        try:
            solver = VQESolver(
                hamiltonian=H,
                ansatz=ansatz,
                optimizer=opt,
                max_iterations=300,
                backend='statevector'
            )

            start_time = time.time()
            result = solver.solve()
            elapsed = time.time() - start_time

            energy = result['energy']
            n_evals = result.get('n_function_evals', 0)
            recovery = 100 * (energy - bond.hf_energy) / (bond.fci_energy - bond.hf_energy)

            print(f"✅ Energy: {energy:.8f} Ha")
            print(f"   Correlation recovery: {recovery:.1f}%")
            print(f"   Function evals: {n_evals}")
            print(f"   Time: {elapsed:.2f}s")
            print(f"   Efficiency: {recovery / n_evals:.4f} %/eval")

            results[opt] = {
                'energy': energy,
                'n_evals': n_evals,
                'recovery': recovery,
                'efficiency': recovery / n_evals if n_evals > 0 else 0
            }

        except Exception as e:
            print(f"❌ Failed: {e}")
            results[opt] = {'error': str(e)}

    print(f"\n{'='*80}")
    print("OPTIMIZER EFFICIENCY RANKING")
    print(f"{'='*80}")

    # Sort by efficiency
    sorted_opts = sorted(
        [(name, res) for name, res in results.items() if 'error' not in res],
        key=lambda x: x[1]['efficiency'],
        reverse=True
    )

    print(f"{'Optimizer':<15} {'Evals':>6} {'Recovery':>8} {'Efficiency':>12}")
    print(f"{'─'*15} {'─'*6} {'─'*8} {'─'*12}")
    for name, res in sorted_opts:
        print(f"{name:<15} {res['n_evals']:>6d} {res['recovery']:>7.1f}% {res['efficiency']:>11.4f}")

    return results


def test_gradient_based_vs_derivative_free():
    """Test 4: How much do gradients help?"""
    print("\n" + "="*80)
    print("TEST 4: GRADIENT-BASED VS DERIVATIVE-FREE")
    print("="*80)

    from kanad.core.molecule import Molecule
    from kanad.core.hamiltonians.openfermion_jw import build_h2_hamiltonian
    from kanad.ansatze import CovalentGovernanceAnsatz
    from kanad.utils.vqe_solver import VQESolver

    # Create H2
    bond = Molecule.from_smiles("[H][H]", name="H2")
    bond.optimize_geometry()
    H = build_h2_hamiltonian(bond)

    n_qubits = 4
    n_electrons = 2
    ansatz = CovalentGovernanceAnsatz(n_qubits=n_qubits, n_electrons=n_electrons, n_layers=2)

    gradient_based = ['SLSQP', 'L-BFGS-B', 'TNC', 'trust-constr']
    derivative_free = ['COBYLA', 'Powell', 'Nelder-Mead']

    results = {'gradient': {}, 'derivative_free': {}}

    print(f"\n🔹 GRADIENT-BASED OPTIMIZERS")
    for opt in gradient_based:
        print(f"\nTesting {opt}...")
        try:
            solver = VQESolver(
                hamiltonian=H,
                ansatz=ansatz,
                optimizer=opt,
                max_iterations=200,
                backend='statevector'
            )

            result = solver.solve()
            energy = result['energy']
            n_evals = result.get('n_function_evals', 0)
            recovery = 100 * (energy - bond.hf_energy) / (bond.fci_energy - bond.hf_energy)

            results['gradient'][opt] = {
                'energy': energy,
                'n_evals': n_evals,
                'recovery': recovery
            }

            print(f"  {opt}: {n_evals} evals, {recovery:.1f}% recovery")

        except Exception as e:
            print(f"  {opt}: ERROR - {e}")

    print(f"\n🔹 DERIVATIVE-FREE OPTIMIZERS")
    for opt in derivative_free:
        print(f"\nTesting {opt}...")
        try:
            solver = VQESolver(
                hamiltonian=H,
                ansatz=ansatz,
                optimizer=opt,
                max_iterations=200,
                backend='statevector'
            )

            result = solver.solve()
            energy = result['energy']
            n_evals = result.get('n_function_evals', 0)
            recovery = 100 * (energy - bond.hf_energy) / (bond.fci_energy - bond.hf_energy)

            results['derivative_free'][opt] = {
                'energy': energy,
                'n_evals': n_evals,
                'recovery': recovery
            }

            print(f"  {opt}: {n_evals} evals, {recovery:.1f}% recovery")

        except Exception as e:
            print(f"  {opt}: ERROR - {e}")

    # Calculate averages
    print(f"\n{'='*80}")
    print("GRADIENT VS DERIVATIVE-FREE COMPARISON")
    print(f"{'='*80}")

    if results['gradient']:
        avg_grad_evals = np.mean([r['n_evals'] for r in results['gradient'].values()])
        avg_grad_recovery = np.mean([r['recovery'] for r in results['gradient'].values()])
        print(f"\nGradient-based average:")
        print(f"  Evals: {avg_grad_evals:.0f}")
        print(f"  Recovery: {avg_grad_recovery:.1f}%")

    if results['derivative_free']:
        avg_df_evals = np.mean([r['n_evals'] for r in results['derivative_free'].values()])
        avg_df_recovery = np.mean([r['recovery'] for r in results['derivative_free'].values()])
        print(f"\nDerivative-free average:")
        print(f"  Evals: {avg_df_evals:.0f}")
        print(f"  Recovery: {avg_df_recovery:.1f}%")

        if results['gradient']:
            print(f"\nGradient-based is {avg_df_evals / avg_grad_evals:.2f}x more efficient (fewer evals)")

    return results


def main():
    """Run all diagnostic tests"""
    print("╔" + "="*78 + "╗")
    print("║" + " ANSATZ EFFICIENCY DIAGNOSTIC SUITE ".center(78) + "║")
    print("╚" + "="*78 + "╝")
    print()
    print("Goal: Identify why ansatze need 600-700+ function evaluations")
    print("      and find ways to reduce this while maintaining accuracy.")
    print()

    all_results = {}

    # Test 1: Parameter initialization
    try:
        all_results['initialization'] = test_parameter_initialization()
    except Exception as e:
        print(f"\n❌ Test 1 failed: {e}")
        import traceback
        traceback.print_exc()

    # Test 2: Ansatz depth
    try:
        all_results['depth'] = test_ansatz_depth_vs_accuracy()
    except Exception as e:
        print(f"\n❌ Test 2 failed: {e}")
        import traceback
        traceback.print_exc()

    # Test 3: Optimizer efficiency
    try:
        all_results['optimizers'] = test_optimizer_efficiency()
    except Exception as e:
        print(f"\n❌ Test 3 failed: {e}")
        import traceback
        traceback.print_exc()

    # Test 4: Gradients
    try:
        all_results['gradients'] = test_gradient_based_vs_derivative_free()
    except Exception as e:
        print(f"\n❌ Test 4 failed: {e}")
        import traceback
        traceback.print_exc()

    # Final summary
    print("\n" + "="*80)
    print("FINAL RECOMMENDATIONS")
    print("="*80)

    print("\nBased on the tests above, here's how to make ansatze more efficient:\n")

    print("1. OPTIMIZER CHOICE:")
    print("   - Use gradient-based optimizers (SLSQP, L-BFGS-B)")
    print("   - Avoid derivative-free (COBYLA, Powell) for > 20 parameters")
    print()

    print("2. PARAMETER INITIALIZATION:")
    print("   - Test showed which initialization strategy is best")
    print("   - Consider chemistry-informed initialization (MP2, CCSD)")
    print()

    print("3. ANSATZ DESIGN:")
    print("   - More layers = more parameters = more evals")
    print("   - Find minimum layers needed for target accuracy")
    print("   - Consider adaptive layer growth")
    print()

    print("4. CONVERGENCE CRITERIA:")
    print("   - Tighten gradient tolerance to stop early")
    print("   - Monitor energy plateau (stop if not improving)")
    print()

    return 0


if __name__ == '__main__':
    sys.exit(main())
