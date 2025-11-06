#!/usr/bin/env python3
"""
Phase 1B: 60-Test Comprehensive Matrix
Tests 5 working ansatze × 2 basis sets × 6 molecules = 60 tests

Ansatze:
  - HardwareEfficientAnsatz
  - RealAmplitudesAnsatz
  - EfficientSU2Ansatz
  - CovalentGovernanceAnsatz
  - IonicGovernanceAnsatz

Basis Sets:
  - sto-3g (minimal, fast)
  - 6-31g (production quality)

Molecules:
  - H2 (covalent, 2e)
  - HeH+ (ionic, 2e)
  - LiH (polar covalent, 4e)
  - H2O (bent, 10e)
  - NH3 (pyramidal, 10e)
  - BeH2 (linear, 6e)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import time
import json
from datetime import datetime
from api.services.experiment_service import create_molecule_from_config
from kanad.ansatze import (
    HardwareEfficientAnsatz,
    RealAmplitudesAnsatz,
    EfficientSU2Ansatz,
    CovalentGovernanceAnsatz,
    IonicGovernanceAnsatz
)
from kanad.utils.vqe_solver import VQESolver

# Test configuration
ANSATZE = [
    ('HardwareEfficientAnsatz', HardwareEfficientAnsatz, {'n_layers': 2}),
    ('RealAmplitudesAnsatz', RealAmplitudesAnsatz, {'n_layers': 2}),
    ('EfficientSU2Ansatz', EfficientSU2Ansatz, {'n_layers': 2}),
    ('CovalentGovernanceAnsatz', CovalentGovernanceAnsatz, {'n_layers': 2}),
    ('IonicGovernanceAnsatz', IonicGovernanceAnsatz, {'n_layers': 2}),
]

BASIS_SETS = ['sto-3g', '6-31g']

MOLECULES = [
    {'name': 'H2', 'smiles': '[H][H]', 'charge': 0, 'multiplicity': 1, 'type': 'covalent'},
    {'name': 'HeH+', 'smiles': '[HeH+]', 'charge': 1, 'multiplicity': 1, 'type': 'ionic'},
    {'name': 'LiH', 'smiles': '[Li][H]', 'charge': 0, 'multiplicity': 1, 'type': 'polar'},
    {'name': 'H2O', 'smiles': 'O', 'charge': 0, 'multiplicity': 1, 'type': 'bent'},
    {'name': 'NH3', 'smiles': 'N', 'charge': 0, 'multiplicity': 1, 'type': 'pyramidal'},
    {'name': 'BeH2', 'smiles': '[BeH2]', 'charge': 0, 'multiplicity': 1, 'type': 'linear'},
]

print("=" * 80)
print("PHASE 1B: 60-TEST COMPREHENSIVE MATRIX")
print("=" * 80)
print(f"\nStarted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"\nConfiguration:")
print(f"  Ansatze:     {len(ANSATZE)}")
print(f"  Basis Sets:  {len(BASIS_SETS)}")
print(f"  Molecules:   {len(MOLECULES)}")
print(f"  Total Tests: {len(ANSATZE) * len(BASIS_SETS) * len(MOLECULES)}")
print()

results = []
test_num = 0
total_tests = len(ANSATZE) * len(BASIS_SETS) * len(MOLECULES)

for mol_config in MOLECULES:
    for basis in BASIS_SETS:
        print("\n" + "=" * 80)
        print(f"Molecule: {mol_config['name']} ({mol_config['type']}) | Basis: {basis}")
        print("=" * 80)

        try:
            # Create molecule
            print(f"\n  Creating molecule {mol_config['name']} with {basis} basis...")
            molecule = create_molecule_from_config({
                'smiles': mol_config['smiles'],
                'basis': basis,
                'charge': mol_config['charge'],
                'multiplicity': mol_config['multiplicity']
            })

            ham = molecule.hamiltonian
            n_qubits = 2 * ham.n_orbitals

            print(f"    Electrons: {molecule.n_electrons}, Orbitals: {ham.n_orbitals}, Qubits: {n_qubits}")

            # Get HF reference
            scf_results, hf_energy = ham.solve_scf(max_iterations=50, conv_tol=1e-6)
            print(f"    HF Energy: {hf_energy:.8f} Ha")

            # Test each ansatz
            for ansatz_name, ansatz_class, ansatz_params in ANSATZE:
                test_num += 1
                print(f"\n  [{test_num}/{total_tests}] Testing {ansatz_name}...")

                try:
                    # Create ansatz
                    params = {
                        'n_qubits': n_qubits,
                        'n_electrons': molecule.n_electrons,
                        **ansatz_params
                    }
                    ansatz = ansatz_class(**params)

                    # Build circuit
                    circuit = ansatz.build_circuit()

                    # Run VQE
                    start_time = time.time()
                    vqe = VQESolver(
                        hamiltonian=ham,
                        ansatz=ansatz,
                        optimizer='COBYLA',
                        max_iterations=100
                    )
                    result = vqe.solve()
                    elapsed = time.time() - start_time

                    # Calculate correlation
                    correlation_mha = (result['energy'] - hf_energy) * 1000

                    # Assess result
                    if correlation_mha < -5:
                        status = "EXCELLENT"
                        symbol = "✅"
                    elif correlation_mha < -1:
                        status = "GOOD"
                        symbol = "✅"
                    elif correlation_mha < 0:
                        status = "FAIR"
                        symbol = "⚠️ "
                    elif abs(correlation_mha) < 0.1:
                        status = "STUCK_AT_HF"
                        symbol = "⚠️ "
                    else:
                        status = "ABOVE_HF"
                        symbol = "❌"

                    print(f"      {symbol} {status}: {correlation_mha:.3f} mHa | {result['function_evaluations']} evals | {elapsed:.2f}s")

                    results.append({
                        'test_num': test_num,
                        'molecule': mol_config['name'],
                        'molecule_type': mol_config['type'],
                        'basis': basis,
                        'ansatz': ansatz_name,
                        'n_electrons': molecule.n_electrons,
                        'n_orbitals': ham.n_orbitals,
                        'n_qubits': n_qubits,
                        'n_parameters': ansatz.n_parameters,
                        'hf_energy': hf_energy,
                        'vqe_energy': result['energy'],
                        'correlation_mha': correlation_mha,
                        'function_evals': result['function_evaluations'],
                        'time_seconds': elapsed,
                        'status': status,
                        'error': None
                    })

                except Exception as e:
                    print(f"      ❌ ERROR: {str(e)[:60]}")
                    results.append({
                        'test_num': test_num,
                        'molecule': mol_config['name'],
                        'molecule_type': mol_config['type'],
                        'basis': basis,
                        'ansatz': ansatz_name,
                        'n_electrons': molecule.n_electrons if 'molecule' in locals() else None,
                        'n_orbitals': ham.n_orbitals if 'ham' in locals() else None,
                        'n_qubits': n_qubits if 'n_qubits' in locals() else None,
                        'n_parameters': None,
                        'hf_energy': hf_energy if 'hf_energy' in locals() else None,
                        'vqe_energy': None,
                        'correlation_mha': None,
                        'function_evals': None,
                        'time_seconds': None,
                        'status': 'ERROR',
                        'error': str(e)[:100]
                    })

        except Exception as e:
            print(f"\n  ❌ Failed to create molecule: {str(e)[:100]}")
            # Skip all ansatze for this molecule/basis combination
            for ansatz_name, _, _ in ANSATZE:
                test_num += 1
                results.append({
                    'test_num': test_num,
                    'molecule': mol_config['name'],
                    'molecule_type': mol_config['type'],
                    'basis': basis,
                    'ansatz': ansatz_name,
                    'n_electrons': None,
                    'n_orbitals': None,
                    'n_qubits': None,
                    'n_parameters': None,
                    'hf_energy': None,
                    'vqe_energy': None,
                    'correlation_mha': None,
                    'function_evals': None,
                    'time_seconds': None,
                    'status': 'MOLECULE_ERROR',
                    'error': str(e)[:100]
                })

# Save results to JSON
output_file = 'tests/phase_1b_results.json'
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)

print("\n" + "=" * 80)
print("FINAL SUMMARY")
print("=" * 80)

# Count statuses
status_counts = {}
for r in results:
    status = r['status']
    status_counts[status] = status_counts.get(status, 0) + 1

print(f"\nTotal Tests: {len(results)}")
print("\nStatus Breakdown:")
for status, count in sorted(status_counts.items()):
    pct = (count / len(results)) * 100
    print(f"  {status:20s}: {count:3d} ({pct:5.1f}%)")

# Best performing ansatz per molecule
print("\n" + "=" * 80)
print("BEST ANSATZ PER MOLECULE (by correlation)")
print("=" * 80)

for mol in MOLECULES:
    mol_results = [r for r in results if r['molecule'] == mol['name'] and r['correlation_mha'] is not None and r['correlation_mha'] < 0]
    if mol_results:
        best = min(mol_results, key=lambda x: x['correlation_mha'])
        print(f"\n{mol['name']:10s} ({mol['type']:10s}): {best['ansatz']}")
        print(f"            Best correlation: {best['correlation_mha']:.3f} mHa (basis: {best['basis']})")

# Best performing basis set
print("\n" + "=" * 80)
print("BASIS SET COMPARISON")
print("=" * 80)

for basis in BASIS_SETS:
    basis_results = [r for r in results if r['basis'] == basis and r['correlation_mha'] is not None and r['correlation_mha'] < 0]
    if basis_results:
        avg_corr = sum(r['correlation_mha'] for r in basis_results) / len(basis_results)
        avg_time = sum(r['time_seconds'] for r in basis_results) / len(basis_results)
        print(f"\n{basis:10s}:")
        print(f"  Average correlation: {avg_corr:.3f} mHa")
        print(f"  Average time:        {avg_time:.2f}s")
        print(f"  Success rate:        {len(basis_results)}/{len(ANSATZE) * len(MOLECULES)} ({(len(basis_results)/(len(ANSATZE)*len(MOLECULES)))*100:.1f}%)")

print("\n" + "=" * 80)
print(f"Results saved to: {output_file}")
print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)
