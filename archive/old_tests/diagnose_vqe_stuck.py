#!/usr/bin/env python3
"""
Deep Diagnostic: Why does VQE get stuck at HF?
Analyze parameter evolution, energy landscape, and optimizer behavior
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.utils.vqe_solver import VQESolver

print("="*80)
print("VQE STUCK DIAGNOSTIC")
print("="*80)

# Create H2 molecule
atom1 = Atom('H', [0.0, 0.0, 0.0])
atom2 = Atom('H', [0.0, 0.0, 0.74])
bond = BondFactory.create_bond(atom1, atom2, distance=0.74, basis='sto-3g')

hf_energy = bond.hamiltonian.solve_scf(max_iterations=50)[1]
print(f"\nHF Energy: {hf_energy:.8f} Ha")

# Run VQE with detailed tracking
print("\nRunning VQE with parameter tracking...")

solver = VQESolver(
    bond=bond,
    ansatz_type='governance',
    mapper_type='jordan_wigner',
    optimizer='COBYLA',
    max_iterations=200,
    backend='statevector',
    shots=None
)

# Track parameters
param_history = []
energy_history = []

def detailed_callback(iteration, energy, parameters):
    param_history.append(parameters.copy())
    energy_history.append(energy)

result = solver.solve(callback=detailed_callback)

print(f"\n✅ VQE Complete")
print(f"  Final energy: {result['energy']:.8f} Ha")
print(f"  HF energy:    {hf_energy:.8f} Ha")
print(f"  Correlation:  {result['energy'] - hf_energy:.8f} Ha")
print(f"  Function evals: {len(energy_history)}")

# Analysis 1: Did parameters move significantly?
print("\n" + "="*80)
print("PARAMETER MOVEMENT ANALYSIS")
print("="*80)

if len(param_history) > 0:
    params_initial = np.array(param_history[0])
    params_final = np.array(param_history[-1])
    param_change = np.abs(params_final - params_initial)

    print(f"\nInitial parameters (first 8): {params_initial[:8]}")
    print(f"Final parameters (first 8):   {params_final[:8]}")
    print(f"\nParameter changes:")
    print(f"  Mean change:   {np.mean(param_change):.6f}")
    print(f"  Max change:    {np.max(param_change):.6f}")
    print(f"  Min change:    {np.min(param_change):.6f}")
    print(f"  Std of change: {np.std(param_change):.6f}")

    if np.max(param_change) < 0.01:
        print(f"\n❌ Parameters barely moved! (max change < 0.01)")
        print(f"   This indicates optimizer is stuck at initialization")
    elif np.max(param_change) < 0.1:
        print(f"\n⚠️  Parameters moved slightly (max change < 0.1)")
    else:
        print(f"\n✅ Parameters moved significantly")

    # Check if stuck near zeros
    if np.max(np.abs(params_final)) < 0.05:
        print(f"\n❌ Final parameters near zero! (max |param| < 0.05)")
        print(f"   This means ansatz is essentially at HF state")

# Analysis 2: Energy landscape exploration
print("\n" + "="*80)
print("ENERGY LANDSCAPE ANALYSIS")
print("="*80)

if len(energy_history) > 10:
    energies = np.array(energy_history)

    print(f"\nEnergy statistics:")
    print(f"  Initial:    {energies[0]:.8f} Ha")
    print(f"  Final:      {energies[-1]:.8f} Ha")
    print(f"  Min:        {np.min(energies):.8f} Ha")
    print(f"  Max:        {np.max(energies):.8f} Ha")
    print(f"  Range:      {np.max(energies) - np.min(energies):.8f} Ha")
    print(f"  Improvement: {energies[0] - energies[-1]:.8f} Ha")

    # Check if energy varied much
    energy_range = np.max(energies) - np.min(energies)
    if energy_range < 0.001:
        print(f"\n❌ Energy barely changed! (range < 0.001 Ha)")
        print(f"   Optimizer is not exploring the landscape")
    elif energy_range < 0.01:
        print(f"\n⚠️  Small energy variation (range < 0.01 Ha)")
    else:
        print(f"\n✅ Significant energy exploration")

    # Check convergence
    last_10_energies = energies[-10:]
    energy_std_final = np.std(last_10_energies)
    if energy_std_final < 1e-6:
        print(f"\n✅ Converged (last 10 iterations std < 1e-6)")
    else:
        print(f"\n⚠️  Not fully converged (std = {energy_std_final:.2e})")

# Analysis 3: Stuck at HF?
print("\n" + "="*80)
print("HF TRAP DETECTION")
print("="*80)

energy_diff_from_hf = abs(result['energy'] - hf_energy)
print(f"\nEnergy difference from HF: {energy_diff_from_hf:.8f} Ha ({energy_diff_from_hf * 627.5:.4f} kcal/mol)")

if energy_diff_from_hf < 0.001:  # < 0.6 kcal/mol
    print(f"❌ STUCK AT HF!")
    print(f"\nPossible causes:")
    print(f"  1. Optimizer not moving parameters")
    print(f"  2. Ansatz gradient flat around HF")
    print(f"  3. COBYLA step size too small")
    print(f"  4. Bad initialization landed in HF basin")

    # Check which
    if len(param_history) > 0 and np.max(param_change) < 0.01:
        print(f"\n🎯 ROOT CAUSE: Optimizer not moving parameters")
        print(f"   Initialization: {params_initial[:4]}")
        print(f"   These parameters don't create gradient away from HF")
    else:
        print(f"\n🎯 ROOT CAUSE: Ansatz landscape is flat around HF")
        print(f"   Parameters moved but energy didn't improve")
else:
    print(f"✅ Not stuck at HF")

# Create visualization
if len(param_history) > 1:
    print("\nCreating visualization...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Energy vs iteration
    ax = axes[0, 0]
    ax.plot(energy_history, 'b-', linewidth=1)
    ax.axhline(hf_energy, color='r', linestyle='--', label='HF Energy')
    ax.set_xlabel('Function Evaluation')
    ax.set_ylabel('Energy (Ha)')
    ax.set_title('VQE Energy Convergence')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Parameter evolution (first 8 params)
    ax = axes[0, 1]
    params_array = np.array(param_history)
    for i in range(min(8, params_array.shape[1])):
        ax.plot(params_array[:, i], label=f'θ{i}', alpha=0.7)
    ax.set_xlabel('Function Evaluation')
    ax.set_ylabel('Parameter Value')
    ax.set_title('Parameter Evolution (first 8)')
    ax.legend(ncol=2, fontsize=8)
    ax.grid(True, alpha=0.3)

    # Plot 3: Parameter change magnitude
    ax = axes[1, 0]
    param_changes = np.diff(params_array, axis=0)
    param_change_norms = np.linalg.norm(param_changes, axis=1)
    ax.plot(param_change_norms, 'g-', linewidth=1)
    ax.set_xlabel('Function Evaluation')
    ax.set_ylabel('||Δθ||')
    ax.set_title('Parameter Step Size')
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')

    # Plot 4: Energy improvement
    ax = axes[1, 1]
    energy_improvements = -np.diff(energy_history)  # Negative because we want decreases
    ax.plot(energy_improvements, 'r-', linewidth=1)
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel('Function Evaluation')
    ax.set_ylabel('Energy Improvement (Ha)')
    ax.set_title('Energy Change per Step')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/mk/deeprealm/kanad/vqe_diagnostic.png', dpi=150)
    print(f"✅ Saved visualization to vqe_diagnostic.png")

print("\n" + "="*80)
