#!/usr/bin/env python3
"""
Benchmark Enhanced Governance Ansatze

Compare:
1. Original Governance (baseline)
2. Governance + MP2 Init
3. Adaptive Governance
4. Hybrid Governance-UCCSD

Metrics:
- Function evaluations
- Final accuracy
- Convergence time
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import time
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.utils.vqe_solver import VQESolver
from kanad.ansatze import CovalentGovernanceAnsatz
from kanad.ansatze.governance_optimized import (
    SmartInitializer,
    AdaptiveGovernanceAnsatz,
    HybridGovernanceUCCSD
)

print("="*80)
print("GOVERNANCE ANSATZ UPGRADES BENCHMARK")
print("="*80)

# Create H2 molecule
print("\n📦 Creating H2 molecule...")
atom1 = Atom('H', [0.0, 0.0, 0.0])
atom2 = Atom('H', [0.0, 0.0, 0.74])
bond = BondFactory.create_bond(atom1, atom2, distance=0.74, basis='sto-3g')

hf_energy = bond.hamiltonian.solve_scf(max_iterations=50)[1]
fci_energy = -1.137284

print(f"✅ Molecule ready")
print(f"   HF:  {hf_energy:.8f} Ha")
print(f"   FCI: {fci_energy:.8f} Ha")
print(f"   Target correlation: {fci_energy - hf_energy:.8f} Ha")

# Benchmark configurations
configs = [
    {
        'name': 'Governance (Original)',
        'ansatz_type': 'original',
        'use_mp2': False,
        'description': 'Baseline: 2-layer governance, random init'
    },
    {
        'name': 'Governance + MP2 Init',
        'ansatz_type': 'mp2_init',
        'use_mp2': True,
        'description': 'Same ansatz, MP2 initialization'
    },
    {
        'name': 'Adaptive Governance',
        'ansatz_type': 'adaptive',
        'use_mp2': True,
        'description': 'Starts with 1 layer, grows if needed'
    },
    {
        'name': 'Hybrid Governance-UCCSD',
        'ansatz_type': 'hybrid',
        'use_mp2': True,
        'description': '1 governance layer + UCCSD corrections'
    },
]

results = []

for config in configs:
    print("\n" + "="*80)
    print(f"TEST: {config['name']}")
    print("="*80)
    print(f"Description: {config['description']}")

    try:
        # Create ansatz
        if config['ansatz_type'] == 'original':
            ansatz = CovalentGovernanceAnsatz(
                n_qubits=4,
                n_electrons=2,
                n_layers=2
            )
            initial_params = None

        elif config['ansatz_type'] == 'mp2_init':
            ansatz = CovalentGovernanceAnsatz(
                n_qubits=4,
                n_electrons=2,
                n_layers=2
            )
            # Get MP2 initialization
            initializer = SmartInitializer(hamiltonian=bond.hamiltonian)
            initial_params = initializer.get_mp2_guess(ansatz.n_parameters)

        elif config['ansatz_type'] == 'adaptive':
            ansatz = AdaptiveGovernanceAnsatz(
                n_qubits=4,
                n_electrons=2,
                max_layers=3,
                use_mp2_init=True
            )
            initial_params = ansatz.get_smart_initial_params(hamiltonian=bond.hamiltonian)

        elif config['ansatz_type'] == 'hybrid':
            ansatz = HybridGovernanceUCCSD(
                n_qubits=4,
                n_electrons=2,
                base_layers=1,
                include_singles=True,
                include_doubles=True
            )
            # Get MP2 initialization for hybrid
            initializer = SmartInitializer(hamiltonian=bond.hamiltonian)
            initial_params = initializer.get_mp2_guess(ansatz.n_parameters)

        print(f"\n✅ Ansatz created")
        print(f"   Parameters: {ansatz.n_parameters}")
        if initial_params is not None:
            print(f"   Initialization: MP2-based")
            print(f"   Initial param range: [{np.min(initial_params):.4f}, {np.max(initial_params):.4f}]")
        else:
            print(f"   Initialization: Random uniform(-0.1, 0.1)")

        # Create solver
        from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner
        from pyscf import gto

        # Build Hamiltonian
        mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g', charge=0, spin=0)
        mol.build()
        mf = mol.RHF().run(verbose=0)

        h1e = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
        h_mo = mf.mo_coeff.T @ h1e @ mf.mo_coeff

        eri_ao = mol.intor('int2e')
        eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl',
                          mf.mo_coeff, mf.mo_coeff, eri_ao, mf.mo_coeff, mf.mo_coeff,
                          optimize=True)

        hamiltonian = openfermion_jordan_wigner(
            h_mo=h_mo,
            eri_mo=eri_mo,
            nuclear_repulsion=mol.energy_nuc(),
            n_electrons=mol.nelectron
        )

        # Track function evaluations
        func_evals = [0]

        def callback(iteration, energy, parameters):
            func_evals[0] = iteration

        # Run VQE
        print(f"\n🚀 Running VQE...")
        start_time = time.time()

        solver = VQESolver(
            hamiltonian=hamiltonian,
            ansatz=ansatz,
            optimizer='SLSQP',
            max_iterations=200,
            backend='statevector',
            initial_parameters=initial_params
        )

        result = solver.solve(callback=callback)
        elapsed_time = time.time() - start_time

        # Analyze results
        energy = result['energy']
        correlation = energy - hf_energy
        expected_corr = fci_energy - hf_energy
        recovery = abs(correlation / expected_corr) * 100 if abs(expected_corr) > 1e-6 else 0
        error_kcal = abs(energy - fci_energy) * 627.5

        print(f"\n✅ VQE Complete")
        print(f"   Final energy:    {energy:.8f} Ha")
        print(f"   Correlation:     {correlation:.8f} Ha")
        print(f"   Recovery:        {recovery:.2f}%")
        print(f"   Error from FCI:  {error_kcal:.4f} kcal/mol")
        print(f"   Function evals:  {func_evals[0]}")
        print(f"   Time:            {elapsed_time:.2f} seconds")

        # Store results
        results.append({
            'name': config['name'],
            'ansatz_type': config['ansatz_type'],
            'n_params': ansatz.n_parameters,
            'func_evals': func_evals[0],
            'energy': energy,
            'correlation': correlation,
            'recovery': recovery,
            'error_kcal': error_kcal,
            'time': elapsed_time,
            'success': recovery > 90
        })

    except Exception as e:
        print(f"\n❌ Failed: {e}")
        import traceback
        traceback.print_exc()

        results.append({
            'name': config['name'],
            'ansatz_type': config['ansatz_type'],
            'n_params': 0,
            'func_evals': 0,
            'energy': 0,
            'correlation': 0,
            'recovery': 0,
            'error_kcal': 999,
            'time': 0,
            'success': False
        })

# Summary
print("\n" + "="*80)
print("BENCHMARK SUMMARY")
print("="*80)

print(f"\n{'Ansatz':<30} {'Params':<8} {'Func Evals':<12} {'Recovery':<10} {'Time (s)':<10} {'Status'}")
print("-"*80)

for r in results:
    status = "✅" if r['success'] else "❌"
    print(f"{r['name']:<30} {r['n_params']:<8} {r['func_evals']:<12} {r['recovery']:>6.1f}% {r['time']:>8.2f}s {status}")

# Compare to baseline
if len(results) > 1 and results[0]['func_evals'] > 0:
    baseline = results[0]
    print(f"\n📊 Speedup vs Baseline:")
    print("-"*80)

    for r in results[1:]:
        if r['func_evals'] > 0:
            speedup = baseline['func_evals'] / r['func_evals']
            time_speedup = baseline['time'] / r['time'] if r['time'] > 0 else 0
            print(f"   {r['name']:<30} {speedup:>5.2f}x faster ({time_speedup:.2f}x time)")

print("\n" + "="*80)
print("CONCLUSIONS")
print("="*80)

# Find best
best = max(results, key=lambda x: x['recovery'] if x['func_evals'] > 0 else 0)
fastest = min([r for r in results if r['func_evals'] > 0], key=lambda x: x['func_evals'])

print(f"\n🏆 Best accuracy: {best['name']} ({best['recovery']:.1f}% recovery)")
print(f"⚡ Fastest: {fastest['name']} ({fastest['func_evals']} function evals)")

if fastest['name'] != results[0]['name']:
    speedup = results[0]['func_evals'] / fastest['func_evals']
    print(f"\n✨ Recommended: {fastest['name']}")
    print(f"   {speedup:.1f}x faster than baseline")
    print(f"   {fastest['recovery']:.1f}% correlation recovery")
    print(f"   {fastest['error_kcal']:.4f} kcal/mol error")

print("\n" + "="*80)
