#!/usr/bin/env python3
"""
Comprehensive test of all working ansatze.
Tests: RealAmplitudes, EfficientSU2, IonicGovernance, AdaptiveGovernance, TwoLocal
Molecule: H2 with STO-3G basis (simple, fast validation)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import time
from api.services.experiment_service import create_molecule_from_config
from kanad.ansatze import (
    RealAmplitudesAnsatz,
    EfficientSU2Ansatz,
    IonicGovernanceAnsatz,
    AdaptiveGovernanceAnsatz,
    TwoLocalAnsatz
)
from kanad.utils.vqe_solver import VQESolver

print("=" * 80)
print("COMPREHENSIVE ANSATZ TESTING - PHASE 1A")
print("=" * 80)
print("\nTest molecule: H2 (STO-3G basis)")
print("Goal: Validate all 5 untested ansatze work correctly\n")

# Create H2 molecule
molecule_config = {
    'smiles': '[H][H]',
    'basis': 'sto-3g',
    'charge': 0,
    'multiplicity': 1
}

print("1. Creating H2 molecule...")
h2 = create_molecule_from_config(molecule_config)
ham = h2.hamiltonian
n_qubits = 2 * ham.n_orbitals

print(f"   Molecule: {h2.n_electrons} electrons, {ham.n_orbitals} orbitals, {n_qubits} qubits")

# Get HF reference
print("\n2. Running SCF calculation...")
scf_results, hf_energy = ham.solve_scf(max_iterations=50, conv_tol=1e-6)
print(f"   HF Energy: {hf_energy:.8f} Ha")

# Define ansatze to test
ansatze_configs = [
    {
        'name': 'RealAmplitudesAnsatz',
        'class': RealAmplitudesAnsatz,
        'params': {'n_qubits': n_qubits, 'n_electrons': h2.n_electrons, 'n_layers': 2}
    },
    {
        'name': 'EfficientSU2Ansatz',
        'class': EfficientSU2Ansatz,
        'params': {'n_qubits': n_qubits, 'n_electrons': h2.n_electrons, 'n_layers': 2}
    },
    {
        'name': 'IonicGovernanceAnsatz',
        'class': IonicGovernanceAnsatz,
        'params': {'n_qubits': n_qubits, 'n_electrons': h2.n_electrons, 'n_layers': 2}
    },
    {
        'name': 'AdaptiveGovernanceAnsatz',
        'class': AdaptiveGovernanceAnsatz,
        'params': {'n_qubits': n_qubits, 'n_electrons': h2.n_electrons, 'n_layers': 2, 'bonding_type': 'auto'}
    },
    {
        'name': 'TwoLocalAnsatz',
        'class': TwoLocalAnsatz,
        'params': {
            'n_qubits': n_qubits,
            'n_electrons': h2.n_electrons,
            'n_layers': 2,
            'rotation_gates': 'ry',
            'entanglement': 'linear'
        }
    }
]

results = []

print("\n" + "=" * 80)
print("TESTING ANSATZE")
print("=" * 80)

for i, config in enumerate(ansatze_configs, 1):
    print(f"\n[{i}/5] Testing {config['name']}...")
    print("-" * 80)

    try:
        # Create ansatz
        ansatz = config['class'](**config['params'])
        print(f"   ✅ Ansatz created: {ansatz.n_parameters} parameters")

        # Build circuit
        circuit = ansatz.build_circuit()
        print(f"   ✅ Circuit built successfully")

        # Run VQE
        print(f"   🔧 Running VQE (max 100 iterations)...")
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
            status = "✅ EXCELLENT"
        elif correlation_mha < 0:
            status = "✅ GOOD"
        elif abs(correlation_mha) < 0.1:
            status = "⚠️  STUCK AT HF"
        else:
            status = "❌ ABOVE HF"

        print(f"\n   {status}")
        print(f"   VQE Energy:         {result['energy']:.8f} Ha")
        print(f"   Correlation:        {correlation_mha:.3f} mHa")
        print(f"   Function Evals:     {result['function_evaluations']}")
        print(f"   Time:               {elapsed:.2f}s")

        results.append({
            'name': config['name'],
            'status': 'PASS' if correlation_mha < -1 else 'WARN' if correlation_mha < 0 else 'FAIL',
            'energy': result['energy'],
            'correlation_mha': correlation_mha,
            'n_params': ansatz.n_parameters,
            'function_evals': result['function_evaluations'],
            'time': elapsed,
            'error': None
        })

    except Exception as e:
        print(f"   ❌ ERROR: {str(e)[:100]}")
        import traceback
        traceback.print_exc()

        results.append({
            'name': config['name'],
            'status': 'ERROR',
            'energy': None,
            'correlation_mha': None,
            'n_params': None,
            'function_evals': None,
            'time': None,
            'error': str(e)
        })

# Summary report
print("\n" + "=" * 80)
print("SUMMARY REPORT")
print("=" * 80)
print(f"\nHF Energy: {hf_energy:.8f} Ha")
print(f"Target: < -5 mHa correlation for EXCELLENT performance\n")

print(f"{'Ansatz':<30s} {'Status':<10s} {'Corr (mHa)':<12s} {'Params':<8s} {'Evals':<8s} {'Time':<8s}")
print("-" * 80)

passed = 0
warned = 0
failed = 0
errors = 0

for r in results:
    if r['status'] == 'PASS':
        passed += 1
        symbol = "✅"
    elif r['status'] == 'WARN':
        warned += 1
        symbol = "⚠️ "
    elif r['status'] == 'FAIL':
        failed += 1
        symbol = "❌"
    else:
        errors += 1
        symbol = "💥"

    corr_str = f"{r['correlation_mha']:.3f}" if r['correlation_mha'] is not None else "N/A"
    params_str = str(r['n_params']) if r['n_params'] is not None else "N/A"
    evals_str = str(r['function_evals']) if r['function_evals'] is not None else "N/A"
    time_str = f"{r['time']:.2f}s" if r['time'] is not None else "N/A"

    print(f"{symbol} {r['name']:<28s} {r['status']:<10s} {corr_str:<12s} {params_str:<8s} {evals_str:<8s} {time_str:<8s}")

    if r['error']:
        print(f"   Error: {r['error'][:60]}")

print("\n" + "=" * 80)
print("FINAL RESULTS")
print("=" * 80)
print(f"✅ PASSED:  {passed}/5 ansatze (excellent correlation recovery)")
print(f"⚠️  WARNED:  {warned}/5 ansatze (some correlation but not excellent)")
print(f"❌ FAILED:  {failed}/5 ansatze (stuck at HF or above HF)")
print(f"💥 ERRORS:  {errors}/5 ansatze (code errors)")

if passed >= 4:
    print("\n🎉 EXCELLENT! Most ansatze working perfectly!")
elif passed >= 2:
    print("\n✅ GOOD! Majority of ansatze working.")
else:
    print("\n⚠️  WARNING: Less than half of ansatze working well.")

print("\n" + "=" * 80)
print("NEXT STEPS")
print("=" * 80)

if passed >= 4:
    print("✅ Ready for Phase 1B: 60-test matrix with all molecules and basis sets")
else:
    print("⚠️  Need to investigate and fix failing ansatze before proceeding")

print("\nPhase 1A Testing Complete!")
print("=" * 80)
