#!/usr/bin/env python3
"""
VQE Optimization Fix - Comprehensive Testing
Tests governance ansatz with different optimizers
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from pyscf import gto

# Import Kanad components
from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz
from kanad.ansatze.two_local_ansatz import TwoLocalAnsatz
from kanad.utils.vqe_solver import VQESolver


def run_vqe_test(ansatz_name, optimizer_name, max_iter=50):
    """Run single VQE test"""
    print(f"\n{'='*80}")
    print(f"TEST: {ansatz_name} + {optimizer_name}")
    print("="*80)

    # Create H2 molecule
    mol = gto.M(
        atom='H 0 0 0; H 0 0 0.74',
        basis='sto-3g',
        charge=0,
        spin=0
    )
    mol.build()

    # Get HF energy
    mf = mol.RHF().run(verbose=0)
    hf_energy = mf.e_tot
    nuclear_repulsion = mol.energy_nuc()

    print(f"Molecule: H2 (0.74 Å)")
    print(f"HF energy: {hf_energy:.8f} Ha")
    print(f"Nuclear repulsion: {nuclear_repulsion:.8f} Ha")

    # Build Hamiltonian
    mo_coeff = mf.mo_coeff
    h1e = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
    h_mo = mo_coeff.T @ h1e @ mo_coeff

    eri_ao = mol.intor('int2e')
    eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', mo_coeff, mo_coeff, eri_ao, mo_coeff, mo_coeff, optimize=True)

    pauli_hamiltonian = openfermion_jordan_wigner(
        h_mo=h_mo,
        eri_mo=eri_mo,
        nuclear_repulsion=nuclear_repulsion,
        n_electrons=mol.nelectron
    )

    print(f"Hamiltonian: {len(pauli_hamiltonian)} Pauli terms")

    # Create ansatz
    n_qubits = 2 * mol.nao_nr()
    n_electrons = mol.nelectron

    if ansatz_name == "CovalentGovernance":
        ansatz = CovalentGovernanceAnsatz(n_qubits=n_qubits, n_electrons=n_electrons)
    elif ansatz_name == "TwoLocal":
        ansatz = TwoLocalAnsatz(n_qubits=n_qubits, reps=2)
    else:
        raise ValueError(f"Unknown ansatz: {ansatz_name}")

    print(f"Ansatz: {ansatz.n_parameters} parameters")

    # Track iterations
    iteration_data = []

    def callback(iteration, energy, parameters):
        iteration_data.append({'iter': iteration, 'energy': energy})
        if iteration % 10 == 0:
            print(f"  Iter {iteration}: E = {energy:.8f} Ha")

    # Create VQE solver
    try:
        solver = VQESolver(
            hamiltonian=pauli_hamiltonian,
            ansatz=ansatz,
            optimizer=optimizer_name,
            max_iterations=max_iter
        )

        print(f"\n🚀 Starting VQE optimization...")
        result = solver.solve(callback=callback)

        vqe_energy = result['energy']
        correlation = vqe_energy - hf_energy
        iterations = result.get('iterations', 0)

        print(f"\n✅ VQE completed")
        print(f"  VQE energy:  {vqe_energy:.8f} Ha")
        print(f"  HF energy:   {hf_energy:.8f} Ha")
        print(f"  Correlation: {correlation:.8f} Ha ({correlation*627.5:.4f} kcal/mol)")
        print(f"  Iterations:  {iterations}")
        print(f"  Func evals:  {len(iteration_data)}")

        # Success if recovered some correlation
        success = correlation < -0.0001 and len(iteration_data) > 3

        return {
            'success': success,
            'correlation': correlation,
            'iterations': iterations,
            'func_evals': len(iteration_data)
        }

    except Exception as e:
        print(f"\n❌ VQE failed: {e}")
        import traceback
        traceback.print_exc()
        return {'success': False, 'error': str(e)}


def main():
    """Run VQE diagnostic tests"""
    print("\n" + "="*80)
    print("VQE OPTIMIZATION FIX - DIAGNOSTIC SUITE")
    print("="*80)

    results = []

    # Test 1: CovalentGovernance + COBYLA
    r1 = run_vqe_test("CovalentGovernance", "COBYLA", max_iter=50)
    r1['config'] = "CovalentGovernance + COBYLA"
    results.append(r1)

    # Test 2: CovalentGovernance + Powell
    r2 = run_vqe_test("CovalentGovernance", "Powell", max_iter=50)
    r2['config'] = "CovalentGovernance + Powell"
    results.append(r2)

    # Test 3: TwoLocal + COBYLA
    r3 = run_vqe_test("TwoLocal", "COBYLA", max_iter=50)
    r3['config'] = "TwoLocal + COBYLA"
    results.append(r3)

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"\n{'Configuration':<35} {'Status':<10} {'Correlation (Ha)':<20}")
    print("-"*70)

    for r in results:
        status = "✅ PASS" if r['success'] else "❌ FAIL"
        corr = f"{r['correlation']:.6f}" if r['success'] else r.get('error', 'Error')[:15]
        print(f"{r['config']:<35} {status:<10} {corr:<20}")

    # Check if any worked
    working = [r for r in results if r['success']]

    if working:
        print(f"\n✅ SUCCESS! {len(working)}/{len(results)} configurations work!")
        best = max(working, key=lambda x: abs(x['correlation']))
        print(f"\nBest: {best['config']}")
        print(f"  Correlation: {best['correlation']:.8f} Ha")
        print(f"  Iterations: {best['iterations']}")
        return 0
    else:
        print(f"\n❌ ALL TESTS FAILED - VQE optimization needs debugging")
        return 1


if __name__ == "__main__":
    sys.exit(main())
