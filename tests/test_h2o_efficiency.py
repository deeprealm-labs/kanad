"""
Test H2O VQE efficiency with different optimizers
Investigate why it needs so many evaluations
"""

import sys
import numpy as np
from pyscf import gto
import time

from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz
from kanad.utils.vqe_solver import VQESolver

def test_h2o_vqe(optimizer='SLSQP', max_iterations=50):
    """Test H2O VQE with given optimizer"""

    print("="*80)
    print(f"H2O VQE Test - Optimizer: {optimizer}")
    print("="*80)

    # Create H2O molecule
    mol = gto.M(
        atom='O 0 0 0; H 0.7596 0.5845 0; H -0.7596 0.5845 0',
        basis='sto-3g',
        charge=0,
        spin=0
    )
    mol.build()

    # Get HF energy
    mf = mol.RHF().run(verbose=0)
    hf_energy = mf.e_tot
    nuclear_repulsion = mol.energy_nuc()

    print(f"\nMolecule: H2O")
    print(f"Atoms: {mol.natm}")
    print(f"Electrons: {mol.nelectron}")
    print(f"Orbitals: {mol.nao_nr()}")
    print(f"Nuclear repulsion: {nuclear_repulsion:.6f} Ha")
    print(f"HF energy: {hf_energy:.6f} Ha")

    # Build Hamiltonian
    mo_coeff = mf.mo_coeff
    h1e = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
    h_mo = mo_coeff.T @ h1e @ mo_coeff

    eri_ao = mol.intor('int2e')
    eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl',
                       mo_coeff, mo_coeff, eri_ao, mo_coeff, mo_coeff,
                       optimize=True)

    pauli_hamiltonian = openfermion_jordan_wigner(
        h_mo=h_mo,
        eri_mo=eri_mo,
        nuclear_repulsion=nuclear_repulsion,
        n_electrons=mol.nelectron
    )

    n_qubits = 2 * mol.nao_nr()
    n_electrons = mol.nelectron

    print(f"\nQubit Hamiltonian:")
    print(f"  Qubits: {n_qubits}")
    print(f"  Pauli terms: {len(pauli_hamiltonian)}")

    # Create ansatz
    ansatz = CovalentGovernanceAnsatz(
        n_qubits=n_qubits,
        n_electrons=n_electrons,
        n_layers=2
    )

    print(f"\nAnsatz: Covalent Governance")
    print(f"  Parameters: {ansatz.n_parameters}")
    print(f"  Layers: 2")

    # Track iterations
    iteration_count = 0
    eval_count = 0
    energy_history = []

    def callback(iteration, energy, parameters):
        nonlocal iteration_count, eval_count
        iteration_count = iteration
        eval_count += 1
        energy_history.append(energy)

        if iteration % 5 == 0:
            print(f"  Iter {iteration:3d}: E = {energy:.8f} Ha ({eval_count} evals)")

    # Run VQE
    print(f"\n{'─'*80}")
    print(f"Running VQE with {optimizer}...")
    print(f"{'─'*80}")

    try:
        solver = VQESolver(
            hamiltonian=pauli_hamiltonian,
            ansatz=ansatz,
            optimizer=optimizer,
            max_iterations=max_iterations,
            backend='statevector'
        )

        start_time = time.time()
        result = solver.solve(callback=callback)
        elapsed = time.time() - start_time

        vqe_energy = result['energy']
        n_func_evals = result.get('n_function_evals', eval_count)

        correlation = vqe_energy - hf_energy

        print(f"\n{'='*80}")
        print(f"RESULTS - {optimizer}")
        print(f"{'='*80}")
        print(f"VQE energy:         {vqe_energy:.8f} Ha")
        print(f"HF energy:          {hf_energy:.8f} Ha")
        print(f"Correlation:        {correlation:.8f} Ha ({correlation*627.5:.2f} kcal/mol)")
        print(f"Iterations:         {iteration_count}")
        print(f"Function evals:     {n_func_evals}")
        print(f"Evals per param:    {n_func_evals / ansatz.n_parameters:.1f}")
        print(f"Evals per iter:     {n_func_evals / max(iteration_count, 1):.1f}")
        print(f"Wall time:          {elapsed:.2f}s")
        print(f"Time per eval:      {elapsed / n_func_evals * 1000:.1f}ms")

        return {
            'optimizer': optimizer,
            'energy': vqe_energy,
            'hf_energy': hf_energy,
            'correlation': correlation,
            'iterations': iteration_count,
            'function_evals': n_func_evals,
            'time': elapsed,
            'energy_history': energy_history
        }

    except Exception as e:
        print(f"\n❌ VQE failed: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    """Test multiple optimizers on H2O"""

    print("╔" + "="*78 + "╗")
    print("║" + " H2O VQE EFFICIENCY TEST ".center(78) + "║")
    print("╚" + "="*78 + "╝")
    print()

    # Test gradient-based optimizers
    optimizers = [
        ('SLSQP', 50),
        ('L-BFGS-B', 50),
        ('BFGS', 50),
    ]

    results = []

    for opt, max_iter in optimizers:
        result = test_h2o_vqe(opt, max_iter)
        if result:
            results.append(result)
        print("\n" + "="*80 + "\n")

    # Summary comparison
    if results:
        print("="*80)
        print("COMPARISON SUMMARY")
        print("="*80)
        print(f"{'Optimizer':<15} {'Evals':>8} {'Energy':>12} {'Correlation':>12} {'Time':>8}")
        print("─"*80)
        for r in results:
            print(f"{r['optimizer']:<15} {r['function_evals']:>8d} "
                  f"{r['energy']:>12.6f} {r['correlation']:>12.6f} {r['time']:>8.2f}s")

        # Find best
        best = min(results, key=lambda x: x['function_evals'])
        print(f"\n✅ Best optimizer: {best['optimizer']}")
        print(f"   {best['function_evals']} evals, {best['time']:.2f}s")
        print(f"   Energy: {best['energy']:.6f} Ha")

    return 0

if __name__ == '__main__':
    sys.exit(main())
