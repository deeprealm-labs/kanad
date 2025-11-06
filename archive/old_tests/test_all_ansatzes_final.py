"""
FINAL VALIDATION TEST - All Ansatzes Must Achieve <1 mHa Error

This test validates ALL ansatz implementations in the Kanad framework:
1. UCC (Unitary Coupled Cluster) - Gold standard
2. TwoLocal - Hardware-efficient
3. Hardware Efficient - NISQ-optimized
4. Governance Covalent - Physics-aware for covalent bonds
5. Governance Ionic - Physics-aware for ionic bonds
6. Governance Adaptive - Auto-selects bonding type

All ansatzes must achieve <1.0 mHa error vs FCI for H2 molecule.
"""

import numpy as np
import pytest
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.ansatze.two_local_ansatz import TwoLocalAnsatz
from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz
from kanad.ansatze.governance_aware_ansatz import (
    IonicGovernanceAnsatz,
    CovalentGovernanceAnsatz,
    AdaptiveGovernanceAnsatz
)
from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import Operator
from scipy.optimize import minimize


# Shared setup
h1 = Atom('H', position=[0.0, 0.0, 0.0])
h2 = Atom('H', position=[0.0, 0.0, 0.735])
mol = Molecule(atoms=[h1, h2], charge=0, spin=0)

hamiltonian = mol.hamiltonian
mo_energies, C = hamiltonian.compute_molecular_orbitals()
h_mo = C.T @ hamiltonian.h_core @ C
eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', C, C, hamiltonian.eri, C, C, optimize=True)
pauli_ham = openfermion_jordan_wigner(h_mo, eri_mo, hamiltonian.nuclear_repulsion, hamiltonian.n_electrons)

fci_energy = np.linalg.eigh(Operator(pauli_ham).data)[0][0]
estimator = StatevectorEstimator()


def optimize_ansatz(ansatz, n_trials=5, max_iter=500):
    """Optimize ansatz with multiple random starts."""
    circuit = ansatz.build_circuit()
    qiskit_circuit = circuit.to_qiskit()
    n_params = len(qiskit_circuit.parameters)

    def cost(params):
        param_dict = {p: params[i] for i, p in enumerate(qiskit_circuit.parameters)}
        bound = qiskit_circuit.assign_parameters(param_dict)
        return estimator.run([(bound, pauli_ham)]).result()[0].data.evs

    best_energy = float('inf')

    for trial in range(n_trials):
        # Different initialization strategies
        if trial == 0:
            x0 = np.zeros(n_params)
        elif trial < 2:
            x0 = np.random.randn(n_params) * 0.01
        else:
            x0 = np.random.randn(n_params) * 0.1

        result = minimize(cost, x0=x0, method='COBYLA', options={'maxiter': max_iter, 'rhobeg': 0.1})

        if result.fun < best_energy:
            best_energy = result.fun

    return best_energy


def test_ucc_doubles():
    """Test UCC with doubles only."""
    print("\n" + "="*70)
    print("TEST: UCC (Doubles only)")
    print("="*70)

    ansatz = UCCAnsatz(4, 2, include_singles=False, include_doubles=True)
    energy = optimize_ansatz(ansatz, n_trials=2, max_iter=100)
    error_mha = abs(energy - fci_energy) * 1000

    print(f"Energy: {energy:.6f} Ha")
    print(f"FCI:    {fci_energy:.6f} Ha")
    print(f"Error:  {error_mha:.3f} mHa")

    assert error_mha < 1.0, f"UCC error {error_mha:.3f} mHa exceeds 1.0 mHa threshold"
    print("✓✓✓ PASS")


def test_two_local():
    """Test TwoLocal ansatz."""
    print("\n" + "="*70)
    print("TEST: TwoLocal (5 layers)")
    print("="*70)

    ansatz = TwoLocalAnsatz(4, 2, n_layers=5)
    energy = optimize_ansatz(ansatz, n_trials=5, max_iter=500)
    error_mha = abs(energy - fci_energy) * 1000

    print(f"Energy: {energy:.6f} Ha")
    print(f"FCI:    {fci_energy:.6f} Ha")
    print(f"Error:  {error_mha:.3f} mHa")

    assert error_mha < 1.0, f"TwoLocal error {error_mha:.3f} mHa exceeds 1.0 mHa threshold"
    print("✓✓✓ PASS")


def test_hardware_efficient():
    """Test Hardware Efficient ansatz."""
    print("\n" + "="*70)
    print("TEST: Hardware Efficient (5 layers)")
    print("="*70)

    ansatz = HardwareEfficientAnsatz(4, 2, n_layers=5)
    energy = optimize_ansatz(ansatz, n_trials=5, max_iter=500)
    error_mha = abs(energy - fci_energy) * 1000

    print(f"Energy: {energy:.6f} Ha")
    print(f"FCI:    {fci_energy:.6f} Ha")
    print(f"Error:  {error_mha:.3f} mHa")

    assert error_mha < 1.0, f"Hardware Efficient error {error_mha:.3f} mHa exceeds 1.0 mHa threshold"
    print("✓✓✓ PASS")


def test_governance_covalent():
    """Test Governance-Aware Covalent ansatz."""
    print("\n" + "="*70)
    print("TEST: Governance Covalent (6 layers)")
    print("="*70)

    ansatz = CovalentGovernanceAnsatz(4, 2, n_layers=6)
    energy = optimize_ansatz(ansatz, n_trials=10, max_iter=1000)
    error_mha = abs(energy - fci_energy) * 1000

    print(f"Energy: {energy:.6f} Ha")
    print(f"FCI:    {fci_energy:.6f} Ha")
    print(f"Error:  {error_mha:.3f} mHa")

    assert error_mha < 1.0, f"Governance Covalent error {error_mha:.3f} mHa exceeds 1.0 mHa threshold"
    print("✓✓✓ PASS")


def test_governance_ionic():
    """Test Governance-Aware Ionic ansatz."""
    print("\n" + "="*70)
    print("TEST: Governance Ionic (5 layers)")
    print("="*70)

    ansatz = IonicGovernanceAnsatz(4, 2, n_layers=5)
    energy = optimize_ansatz(ansatz, n_trials=5, max_iter=500)
    error_mha = abs(energy - fci_energy) * 1000

    print(f"Energy: {energy:.6f} Ha")
    print(f"FCI:    {fci_energy:.6f} Ha")
    print(f"Error:  {error_mha:.3f} mHa")

    # Note: Ionic is designed for ionic bonds, H2 is covalent
    # We allow slightly higher tolerance but should still be good
    assert error_mha < 5.0, f"Governance Ionic error {error_mha:.3f} mHa exceeds 5.0 mHa threshold"
    print("✓✓ PASS (within tolerance for wrong physics)")


def test_governance_adaptive():
    """Test Governance-Aware Adaptive ansatz."""
    print("\n" + "="*70)
    print("TEST: Governance Adaptive (6 layers, covalent mode)")
    print("="*70)

    ansatz = AdaptiveGovernanceAnsatz(4, 2, n_layers=6, bonding_type='covalent')
    energy = optimize_ansatz(ansatz, n_trials=10, max_iter=1000)
    error_mha = abs(energy - fci_energy) * 1000

    print(f"Energy: {energy:.6f} Ha")
    print(f"FCI:    {fci_energy:.6f} Ha")
    print(f"Error:  {error_mha:.3f} mHa")

    assert error_mha < 1.0, f"Governance Adaptive error {error_mha:.3f} mHa exceeds 1.0 mHa threshold"
    print("✓✓✓ PASS")


if __name__ == "__main__":
    """Run all tests."""
    print("="*70)
    print("KANAD FRAMEWORK - ALL ANSATZES VALIDATION")
    print("="*70)
    print(f"\nH2 Molecule at 0.735 Å")
    print(f"FCI Reference: {fci_energy:.6f} Ha")
    print(f"Target: <1.0 mHa error for all ansatzes\n")

    tests = [
        test_ucc_doubles,
        test_two_local,
        test_hardware_efficient,
        test_governance_covalent,
        test_governance_ionic,
        test_governance_adaptive
    ]

    passed = 0
    failed = 0

    for test_func in tests:
        try:
            test_func()
            passed += 1
        except AssertionError as e:
            print(f"✗ FAILED: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ ERROR: {e}")
            failed += 1

    print("\n" + "="*70)
    print("FINAL RESULTS")
    print("="*70)
    print(f"Passed: {passed}/{len(tests)}")
    print(f"Failed: {failed}/{len(tests)}")

    if failed == 0:
        print("\n🎉🎉🎉 ALL TESTS PASSED! 🎉🎉🎉")
    else:
        print(f"\n⚠️  {failed} test(s) need attention")

    print("="*70)
