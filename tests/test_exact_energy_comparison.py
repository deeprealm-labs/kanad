#!/usr/bin/env python3
"""
Compare VQE energies with exact FCI reference
Investigate why VQE converges to HF instead of exact energy
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from pyscf import gto, scf, fci

# Import Kanad components
from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz
from kanad.utils.vqe_solver import VQESolver


def get_exact_energies():
    """Get exact FCI energy for H2"""
    print("\n" + "="*80)
    print("EXACT REFERENCE ENERGIES - H2 @ 0.74 Å")
    print("="*80)

    # Create H2 molecule
    mol = gto.M(
        atom='H 0 0 0; H 0 0 0.74',
        basis='sto-3g',
        charge=0,
        spin=0
    )
    mol.build()

    # Hartree-Fock
    mf = scf.RHF(mol)
    mf.kernel()
    hf_energy = mf.e_tot

    print(f"\n1. Hartree-Fock Energy: {hf_energy:.8f} Ha")

    # Full CI (EXACT for this basis)
    cisolver = fci.FCI(mf)
    fci_energy = cisolver.kernel()[0]

    print(f"2. FCI Energy (EXACT):  {fci_energy:.8f} Ha")

    correlation = fci_energy - hf_energy
    print(f"3. Correlation Energy:  {correlation:.8f} Ha ({correlation*627.5:.4f} kcal/mol)")

    return {
        'hf': hf_energy,
        'fci': fci_energy,
        'correlation': correlation,
        'mol': mol,
        'mf': mf
    }


def test_vqe_with_different_ansatze(ref):
    """Test VQE with different ansatze"""
    print("\n" + "="*80)
    print("VQE TESTS - Comparing with Exact FCI")
    print("="*80)

    mol = ref['mol']
    mf = ref['mf']

    # Build Hamiltonian
    mo_coeff = mf.mo_coeff
    h1e = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
    h_mo = mo_coeff.T @ h1e @ mo_coeff

    eri_ao = mol.intor('int2e')
    eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', mo_coeff, mo_coeff, eri_ao, mo_coeff, mo_coeff, optimize=True)

    nuclear_repulsion = mol.energy_nuc()

    pauli_hamiltonian = openfermion_jordan_wigner(
        h_mo=h_mo,
        eri_mo=eri_mo,
        nuclear_repulsion=nuclear_repulsion,
        n_electrons=mol.nelectron
    )

    print(f"\nHamiltonian: {len(pauli_hamiltonian)} Pauli terms")
    print(f"Qubits: {pauli_hamiltonian.num_qubits}")

    # Test different configurations
    n_qubits = 2 * mol.nao_nr()  # Should be 4 for H2/STO-3G
    n_electrons = mol.nelectron  # Should be 2 for H2

    print(f"\nSystem: {n_qubits} qubits, {n_electrons} electrons")

    results = []

    # Test 1: CovalentGovernance + COBYLA
    print("\n" + "-"*80)
    print("TEST 1: CovalentGovernance + COBYLA")
    print("-"*80)

    ansatz1 = CovalentGovernanceAnsatz(n_qubits=n_qubits, n_electrons=n_electrons)
    print(f"Ansatz parameters: {ansatz1.n_parameters}")

    solver1 = VQESolver(
        hamiltonian=pauli_hamiltonian,
        ansatz=ansatz1,
        optimizer='COBYLA',
        max_iterations=100  # More iterations
    )

    result1 = solver1.solve()
    vqe_energy1 = result1['energy']

    print(f"\nResults:")
    print(f"  VQE Energy:      {vqe_energy1:.8f} Ha")
    print(f"  HF Energy:       {ref['hf']:.8f} Ha")
    print(f"  FCI Energy:      {ref['fci']:.8f} Ha")
    print(f"  VQE - HF:        {vqe_energy1 - ref['hf']:.8f} Ha")
    print(f"  VQE - FCI:       {vqe_energy1 - ref['fci']:.8f} Ha (error)")
    print(f"  FCI - HF:        {ref['correlation']:.8f} Ha (target correlation)")

    recovery_percent = abs((vqe_energy1 - ref['hf']) / ref['correlation']) * 100
    print(f"  Correlation Recovery: {recovery_percent:.2f}%")

    results.append({
        'name': 'CovalentGovernance + COBYLA',
        'energy': vqe_energy1,
        'error': vqe_energy1 - ref['fci'],
        'recovery': recovery_percent
    })

    # Test 2: CovalentGovernance + Powell
    print("\n" + "-"*80)
    print("TEST 2: CovalentGovernance + Powell")
    print("-"*80)

    ansatz2 = CovalentGovernanceAnsatz(n_qubits=n_qubits, n_electrons=n_electrons)

    solver2 = VQESolver(
        hamiltonian=pauli_hamiltonian,
        ansatz=ansatz2,
        optimizer='Powell',
        max_iterations=100
    )

    result2 = solver2.solve()
    vqe_energy2 = result2['energy']

    print(f"\nResults:")
    print(f"  VQE Energy:      {vqe_energy2:.8f} Ha")
    print(f"  HF Energy:       {ref['hf']:.8f} Ha")
    print(f"  FCI Energy:      {ref['fci']:.8f} Ha")
    print(f"  VQE - HF:        {vqe_energy2 - ref['hf']:.8f} Ha")
    print(f"  VQE - FCI:       {vqe_energy2 - ref['fci']:.8f} Ha (error)")

    recovery_percent2 = abs((vqe_energy2 - ref['hf']) / ref['correlation']) * 100
    print(f"  Correlation Recovery: {recovery_percent2:.2f}%")

    results.append({
        'name': 'CovalentGovernance + Powell',
        'energy': vqe_energy2,
        'error': vqe_energy2 - ref['fci'],
        'recovery': recovery_percent2
    })

    return results


def check_ansatz_capability():
    """Check if ansatz can represent FCI state"""
    print("\n" + "="*80)
    print("ANSATZ CAPABILITY CHECK")
    print("="*80)

    # For H2 with 4 qubits, FCI wavefunction is:
    # |Ψ⟩ = c0|0011⟩ + c1|0101⟩ + c2|0110⟩ + c3|1001⟩ + c4|1010⟩ + c5|1100⟩ + ...
    #
    # Key question: Can CovalentGovernanceAnsatz represent this?

    print("\nFor H2 (2 electrons, 4 qubits):")
    print("- Hartree-Fock state: |1100⟩ (both electrons in lowest orbital)")
    print("- FCI includes all 2-electron configurations")
    print("- Number of 2-electron configurations: C(4,2) = 6")
    print("\nConfigurations:")
    print("  |0011⟩, |0101⟩, |0110⟩, |1001⟩, |1010⟩, |1100⟩")
    print("\nFor EXACT energy, ansatz must:")
    print("  1. Start from HF state |1100⟩")
    print("  2. Mix in excited states (especially |0011⟩)")
    print("  3. Use correct amplitudes")

    # Check what CovalentGovernanceAnsatz does
    ansatz = CovalentGovernanceAnsatz(n_qubits=4, n_electrons=2)

    print(f"\nCovalentGovernanceAnsatz:")
    print(f"  Parameters: {ansatz.n_parameters}")

    # Build circuit to inspect
    params = np.zeros(ansatz.n_parameters)
    circuit = ansatz.build(params)

    print(f"  Circuit gates: {circuit.count_ops()}")
    print(f"  Circuit depth: {circuit.depth()}")

    return ansatz


def main():
    """Main analysis"""
    print("\n" + "="*80)
    print("VQE ENERGY CONVERGENCE INVESTIGATION")
    print("Why does VQE give ~-1.117 Ha instead of ~-1.137 Ha?")
    print("="*80)

    # Step 1: Get exact reference
    ref = get_exact_energies()

    # Step 2: Check ansatz capability
    ansatz = check_ansatz_capability()

    # Step 3: Test VQE
    results = test_vqe_with_different_ansatze(ref)

    # Final analysis
    print("\n" + "="*80)
    print("ANALYSIS & DIAGNOSIS")
    print("="*80)

    print("\nTarget Energies:")
    print(f"  HF:  {ref['hf']:.8f} Ha")
    print(f"  FCI: {ref['fci']:.8f} Ha")
    print(f"  Correlation: {ref['correlation']:.8f} Ha")

    print("\nVQE Results:")
    for r in results:
        print(f"  {r['name']}:")
        print(f"    Energy: {r['energy']:.8f} Ha")
        print(f"    Error from FCI: {r['error']:.8f} Ha")
        print(f"    Correlation recovery: {r['recovery']:.2f}%")

    # Diagnosis
    print("\n" + "="*80)
    print("POTENTIAL ISSUES:")
    print("="*80)

    best_recovery = max(r['recovery'] for r in results)

    if best_recovery < 50:
        print("\n❌ CRITICAL: VQE recovering < 50% of correlation energy")
        print("\nPossible causes:")
        print("  1. Ansatz not expressive enough (can't represent FCI state)")
        print("  2. Optimizer getting stuck in local minima")
        print("  3. Initial parameters poor (starting too close to HF)")
        print("  4. Hamiltonian construction issue")
        print("  5. Circuit not deep enough")
    elif best_recovery < 80:
        print("\n⚠️  WARNING: VQE recovering 50-80% of correlation")
        print("\nLikely causes:")
        print("  1. Need more optimizer iterations")
        print("  2. Better initial parameters needed")
        print("  3. Ansatz could be more expressive")
    else:
        print("\n✅ GOOD: VQE recovering > 80% of correlation")
        print("This is acceptable for variational methods")

    print("\n" + "="*80)
    print("RECOMMENDATION:")
    print("="*80)

    if best_recovery < 80:
        print("\n1. CHECK: Initial state preparation")
        print("   - Is HF state correctly prepared?")
        print("   - Are excitation operators applied?")
        print("\n2. CHECK: Ansatz depth")
        print("   - May need more layers/reps")
        print("   - UCC ansatz might work better")
        print("\n3. CHECK: Parameter initialization")
        print("   - Try small random initialization")
        print("   - Try UCCSD-like initialization")
        print("\n4. TEST: Use UCCSD ansatz as gold standard")
        print("   - UCCSD should recover ~90%+ correlation")

    return 0 if best_recovery > 50 else 1


if __name__ == "__main__":
    sys.exit(main())
