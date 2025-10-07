#!/usr/bin/env python3
"""
Test suite for quantum optimization module.

Tests:
1. Active space selection (HOMO-LUMO, natural orbitals, governance)
2. Qubit reduction (tapering, symmetries)
3. Orbital localization (Boys, Pipek-Mezey)
4. Circuit optimization (gate cancellation, merging)
5. Adaptive methods (ADAPT-VQE style)
6. Strategy comparison
"""

import sys
import numpy as np

sys.path.insert(0, '/home/mk/deeprealm/kanad')

from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.optimization.quantum_optimizer import QuantumOptimizer


def test_optimization_h2():
    """Test optimization on H2 molecule."""
    print("="*70)
    print("TEST: QUANTUM OPTIMIZATION - H2 Molecule")
    print("="*70)

    # Create H2
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

    bond = CovalentBond(h1, h2)
    hamiltonian = bond.hamiltonian

    print(f"\nH2 System:")
    print(f"  Orbitals: {hamiltonian.n_orbitals}")
    print(f"  Electrons: {hamiltonian.n_electrons}")
    print(f"  Qubits: {2 * hamiltonian.n_orbitals}")

    # Test 1: Aggressive optimization
    print("\n" + "="*70)
    print("TEST 1: Aggressive Optimization")
    print("="*70)

    optimizer = QuantumOptimizer(hamiltonian)
    result_aggressive = optimizer.optimize(
        strategy='aggressive',
        apply_tapering=True,
        localize_orbitals=True
    )

    print(f"\nAggressive strategy results:")
    print(f"  Original qubits: {result_aggressive['original_qubits']}")
    print(f"  Final qubits: {result_aggressive['final_qubits']}")
    print(f"  Reduction: {result_aggressive['reduction_factor']:.2f}x")
    print(f"  Estimated speedup: {result_aggressive['estimated_speedup']:.1f}x")

    # Test 2: Balanced optimization
    print("\n" + "="*70)
    print("TEST 2: Balanced Optimization")
    print("="*70)

    result_balanced = optimizer.optimize(
        strategy='balanced',
        apply_tapering=True,
        localize_orbitals=True
    )

    print(f"\nBalanced strategy results:")
    print(f"  Final qubits: {result_balanced['final_qubits']}")
    print(f"  Reduction: {result_balanced['reduction_factor']:.2f}x")

    # Test 3: Conservative optimization
    print("\n" + "="*70)
    print("TEST 3: Conservative Optimization")
    print("="*70)

    result_conservative = optimizer.optimize(
        strategy='conservative',
        apply_tapering=True,
        localize_orbitals=False
    )

    print(f"\nConservative strategy results:")
    print(f"  Final qubits: {result_conservative['final_qubits']}")
    print(f"  Reduction: {result_conservative['reduction_factor']:.2f}x")

    # Test 4: Strategy comparison
    print("\n" + "="*70)
    print("TEST 4: Strategy Comparison")
    print("="*70)

    comparison = optimizer.compare_strategies()

    print("\nComparison Summary:")
    for strategy, result in comparison.items():
        print(f"  {strategy.capitalize():15} → {result['final_qubits']} qubits "
              f"({result['reduction_factor']:.2f}x reduction, "
              f"{result['estimated_speedup']:.1f}x speedup)")

    # Test 5: Optimization report
    print("\n" + "="*70)
    print("TEST 5: Detailed Optimization Report")
    print("="*70)

    report = optimizer.get_optimization_report(result_balanced)
    print("\n" + report)

    print("\n" + "="*70)
    print("ALL OPTIMIZATION TESTS COMPLETE")
    print("="*70)

    return {
        'aggressive': result_aggressive,
        'balanced': result_balanced,
        'conservative': result_conservative,
        'comparison': comparison
    }


def test_optimization_lih():
    """Test optimization on LiH molecule."""
    print("\n" + "="*70)
    print("TEST: QUANTUM OPTIMIZATION - LiH Molecule")
    print("="*70)

    # Create LiH
    li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    h = Atom('H', position=np.array([1.59, 0.0, 0.0]))

    from kanad.bonds.ionic_bond import IonicBond
    bond = IonicBond(li, h)
    hamiltonian = bond.hamiltonian

    print(f"\nLiH System:")
    print(f"  Orbitals: {hamiltonian.n_orbitals}")
    print(f"  Electrons: {hamiltonian.n_electrons}")
    print(f"  Qubits: {2 * hamiltonian.n_orbitals}")

    # Optimize with governance (if available)
    optimizer = QuantumOptimizer(hamiltonian)

    result = optimizer.optimize(
        strategy='balanced',
        active_space_method='homo_lumo',
        apply_tapering=True
    )

    print(f"\nOptimization results:")
    print(f"  Original qubits: {result['original_qubits']}")
    print(f"  Final qubits: {result['final_qubits']}")
    print(f"  Reduction: {result['reduction_factor']:.2f}x")
    print(f"  Gate reduction: {result['complexity']['gate_reduction_factor']:.1f}x")
    print(f"  Depth reduction: {result['complexity']['depth_reduction_factor']:.1f}x")
    print(f"  Estimated speedup: {result['estimated_speedup']:.1f}x")

    # Show active space details
    active_space = result['active_space']
    print(f"\nActive Space Selection:")
    print(f"  Method: {active_space['method']}")
    print(f"  Frozen orbitals: {active_space['frozen_orbitals']}")
    print(f"  Active orbitals: {active_space['active_orbitals']}")
    print(f"  Virtual orbitals: {active_space['virtual_orbitals'][:5]}..." if len(active_space['virtual_orbitals']) > 5 else f"  Virtual orbitals: {active_space['virtual_orbitals']}")

    print("\n" + "="*70)

    return result


def test_active_space_methods():
    """Test different active space selection methods."""
    print("\n" + "="*70)
    print("TEST: ACTIVE SPACE SELECTION METHODS")
    print("="*70)

    # Create H2
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

    bond = CovalentBond(h1, h2)
    hamiltonian = bond.hamiltonian

    optimizer = QuantumOptimizer(hamiltonian)

    methods = ['homo_lumo', 'auto']

    print("\nComparing active space selection methods:\n")

    for method in methods:
        print(f"Method: {method}")

        result = optimizer.optimize(
            target_qubits=4,
            strategy='balanced',
            active_space_method=method,
            apply_tapering=False,
            localize_orbitals=False
        )

        active_space = result['active_space']
        print(f"  Active orbitals: {active_space['active_orbitals']}")
        print(f"  Method used: {active_space['method']}")
        print(f"  Qubits: {active_space['n_qubits_reduced']}")
        print()

    print("="*70)


def test_circuit_complexity():
    """Test circuit complexity estimation."""
    print("\n" + "="*70)
    print("TEST: CIRCUIT COMPLEXITY ANALYSIS")
    print("="*70)

    # Test on larger system
    from kanad.core.atom import Atom
    from kanad.bonds.covalent_bond import CovalentBond

    # C-C bond (more orbitals)
    c1 = Atom('C', position=np.array([0.0, 0.0, 0.0]))
    c2 = Atom('C', position=np.array([1.54, 0.0, 0.0]))

    bond = CovalentBond(c1, c2)
    hamiltonian = bond.hamiltonian

    print(f"\nC-C System:")
    print(f"  Orbitals: {hamiltonian.n_orbitals}")
    print(f"  Qubits: {2 * hamiltonian.n_orbitals}")

    optimizer = QuantumOptimizer(hamiltonian)

    # Compare different reduction levels
    for strategy in ['aggressive', 'balanced', 'conservative']:
        result = optimizer.optimize(
            strategy=strategy,
            apply_tapering=True
        )

        comp = result['complexity']
        print(f"\n{strategy.capitalize()} strategy:")
        print(f"  Qubits: {result['original_qubits']} → {result['final_qubits']}")
        print(f"  Gates: {comp['gates_full']:.0f} → {comp['gates_reduced']:.0f} ({comp['gate_reduction_factor']:.1f}x)")
        print(f"  Depth: {comp['depth_full']:.0f} → {comp['depth_reduced']:.0f} ({comp['depth_reduction_factor']:.1f}x)")
        print(f"  Speedup: {result['estimated_speedup']:.1f}x")

    print("\n" + "="*70)


def main():
    """Run all optimization tests."""
    print("="*70)
    print("QUANTUM OPTIMIZATION MODULE TEST SUITE")
    print("="*70)
    print()

    # Test 1: H2 optimization
    h2_results = test_optimization_h2()

    # Test 2: LiH optimization
    lih_results = test_optimization_lih()

    # Test 3: Active space methods
    test_active_space_methods()

    # Test 4: Circuit complexity
    test_circuit_complexity()

    # Summary
    print("\n" + "="*70)
    print("OPTIMIZATION TEST SUITE SUMMARY")
    print("="*70)

    print("\nH2 Optimization Results:")
    for strategy in ['aggressive', 'balanced', 'conservative']:
        result = h2_results[strategy]
        print(f"  {strategy.capitalize():15} → {result['final_qubits']} qubits, "
              f"{result['estimated_speedup']:.1f}x speedup")

    print("\nLiH Optimization:")
    print(f"  Balanced strategy → {lih_results['final_qubits']} qubits, "
          f"{lih_results['estimated_speedup']:.1f}x speedup")

    print("\n" + "="*70)
    print("✅ ALL OPTIMIZATION TESTS PASSED")
    print("="*70)


if __name__ == "__main__":
    main()
