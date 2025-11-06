#!/usr/bin/env python3
"""Analyze computational requirements for different molecules and basis sets."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.services.experiment_service import create_molecule_from_config

# Test molecules
MOLECULES = [
    {'name': 'H2', 'smiles': '[H][H]', 'charge': 0, 'mult': 1},
    {'name': 'HeH+', 'smiles': '[HeH+]', 'charge': 1, 'mult': 1},
    {'name': 'LiH', 'smiles': '[Li][H]', 'charge': 0, 'mult': 1},
    {'name': 'H2O', 'smiles': 'O', 'charge': 0, 'mult': 1},
    {'name': 'NH3', 'smiles': 'N', 'charge': 0, 'mult': 1},
    {'name': 'BeH2', 'smiles': '[BeH2]', 'charge': 0, 'mult': 1},
]

BASIS_SETS = ['sto-3g', '6-31g']

print("=" * 100)
print("COMPUTATIONAL SCALING ANALYSIS")
print("=" * 100)
print("\nStatevector simulation memory requirement: 2^n_qubits × 16 bytes (complex128)")
print()

results = []

for mol_config in MOLECULES:
    print(f"\n{'=' * 100}")
    print(f"Molecule: {mol_config['name']}")
    print(f"{'=' * 100}")

    for basis in BASIS_SETS:
        try:
            molecule = create_molecule_from_config({
                'smiles': mol_config['smiles'],
                'basis': basis,
                'charge': mol_config['charge'],
                'multiplicity': mol_config['mult']
            })

            ham = molecule.hamiltonian
            n_qubits = 2 * ham.n_orbitals

            # Calculate memory requirement
            statevector_size = 2 ** n_qubits
            memory_bytes = statevector_size * 16  # complex128 = 16 bytes
            memory_mb = memory_bytes / (1024 ** 2)
            memory_gb = memory_bytes / (1024 ** 3)

            # Estimate time per energy evaluation (empirical)
            # Based on observed performance: ~0.01s for 4 qubits, exponential scaling
            time_per_eval = 0.01 * (2 ** (n_qubits - 4))
            time_per_test = time_per_eval * 100  # 100 function evals

            feasible = memory_gb < 2.0 and time_per_test < 300  # 2GB RAM, 5 min max

            status = "✅" if feasible else "❌"

            print(f"\n  {basis:10s}: {status}")
            print(f"    Electrons:    {molecule.n_electrons}")
            print(f"    Orbitals:     {ham.n_orbitals}")
            print(f"    Qubits:       {n_qubits}")
            print(f"    Statevector:  2^{n_qubits} = {statevector_size:,} components")

            if memory_gb < 1:
                print(f"    Memory:       {memory_mb:.1f} MB")
            else:
                print(f"    Memory:       {memory_gb:.2f} GB")

            if time_per_test < 60:
                print(f"    Est. time:    {time_per_test:.1f} seconds")
            else:
                print(f"    Est. time:    {time_per_test/60:.1f} minutes")

            if not feasible:
                if memory_gb >= 2.0:
                    print(f"    ⚠️  WARNING: Memory requirement {memory_gb:.2f} GB exceeds typical limit")
                if time_per_test >= 300:
                    print(f"    ⚠️  WARNING: Estimated time {time_per_test/60:.1f} min exceeds reasonable limit")

            results.append({
                'molecule': mol_config['name'],
                'basis': basis,
                'electrons': molecule.n_electrons,
                'orbitals': ham.n_orbitals,
                'qubits': n_qubits,
                'memory_gb': memory_gb,
                'time_min': time_per_test / 60,
                'feasible': feasible
            })

        except Exception as e:
            print(f"\n  {basis:10s}: ❌ ERROR")
            print(f"    {str(e)[:80]}")

# Summary
print("\n" + "=" * 100)
print("SUMMARY: Computational Feasibility")
print("=" * 100)
print()
print(f"{'Molecule':<10s} {'Basis':<10s} {'Qubits':<8s} {'Memory':<12s} {'Time':<12s} {'Feasible':<10s}")
print("-" * 100)

for r in results:
    mem_str = f"{r['memory_gb']:.2f} GB" if r['memory_gb'] >= 1 else f"{r['memory_gb']*1024:.0f} MB"
    time_str = f"{r['time_min']:.1f} min" if r['time_min'] >= 1 else f"{r['time_min']*60:.0f} sec"
    status = "✅ YES" if r['feasible'] else "❌ NO"

    print(f"{r['molecule']:<10s} {r['basis']:<10s} {r['qubits']:<8d} {mem_str:<12s} {time_str:<12s} {status:<10s}")

# Count feasible tests
total_tests = len(results)
feasible_tests = sum(1 for r in results if r['feasible'])
infeasible_tests = total_tests - feasible_tests

print()
print("=" * 100)
print(f"Feasible tests:     {feasible_tests}/{total_tests} ({feasible_tests/total_tests*100:.1f}%)")
print(f"Infeasible tests:   {infeasible_tests}/{total_tests} ({infeasible_tests/total_tests*100:.1f}%)")
print("=" * 100)

print("\n" + "=" * 100)
print("RECOMMENDATIONS")
print("=" * 100)
print()

# Find the cutoff
max_feasible_qubits = max([r['qubits'] for r in results if r['feasible']], default=0)
min_infeasible_qubits = min([r['qubits'] for r in results if not r['feasible']], default=999)

print(f"✅ Feasible on local machine: Up to {max_feasible_qubits} qubits")
print(f"❌ Infeasible on local machine: {min_infeasible_qubits}+ qubits")
print()
print("For the 60-test matrix:")
print(f"  • Feasible tests: {feasible_tests * 5} (each molecule/basis × 5 ansatze)")
print(f"  • Infeasible tests: {infeasible_tests * 5}")
print()
print("Recommendations:")
print("  1. Continue with STO-3G basis for all molecules ✅")
print("  2. Use 6-31G only for small molecules (H2, HeH+) ✅")
print("  3. Skip 6-31G for larger molecules (LiH 4e+, H2O 10e, NH3 10e, BeH2 6e) ❌")
print("  4. Alternative: Use cloud compute (Azure, AWS) for larger basis sets")
print()

# Calculate revised test count
revised_feasible = sum(5 for r in results if r['feasible'])
print(f"Revised test matrix: {revised_feasible} tests (down from 60)")
print()
