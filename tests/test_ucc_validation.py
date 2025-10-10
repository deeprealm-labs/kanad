"""
Validation test for UCC ansatz with corrected Pauli evolution gates.

This test validates that the UCC (Unitary Coupled Cluster) ansatz correctly:
1. Implements fermionic excitations via Jordan-Wigner transformation
2. Uses Pauli evolution gates to create proper quantum interference
3. Achieves FCI (Full Configuration Interaction) accuracy for H2 molecule

Key innovations:
- Uses OpenFermion's validated Jordan-Wigner transformation
- Handles complex coefficients by factoring out imaginary unit
- Negates time parameter for correct evolution direction in PauliEvolutionGate
"""

import numpy as np
import pytest
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner


def test_ucc_double_excitation_h2():
    """
    Test that UCC double excitation achieves FCI accuracy for H2.

    This is the critical test validating the UCC implementation.
    H2 at equilibrium with minimal basis requires only double excitation
    to reach FCI accuracy.
    """
    # Create H2 molecule at equilibrium
    h1 = Atom('H', position=[0.0, 0.0, 0.0])
    h2 = Atom('H', position=[0.0, 0.0, 0.735])
    mol = Molecule(atoms=[h1, h2], charge=0, spin=0)

    # Get Hamiltonian in Pauli form
    hamiltonian = mol.hamiltonian
    mo_energies, C = hamiltonian.compute_molecular_orbitals()
    h_mo = C.T @ hamiltonian.h_core @ C
    eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', C, C, hamiltonian.eri, C, C, optimize=True)
    pauli_ham = openfermion_jordan_wigner(
        h_mo, eri_mo, hamiltonian.nuclear_repulsion, hamiltonian.n_electrons
    )

    # Compute FCI reference
    from qiskit.quantum_info import Operator
    ham_matrix = Operator(pauli_ham).data
    fci_energy = np.linalg.eigh(ham_matrix)[0][0]

    # Build UCC ansatz with only double excitation
    ansatz = UCCAnsatz(
        n_qubits=4,
        n_electrons=2,
        include_singles=False,
        include_doubles=True
    )

    circuit = ansatz.build_circuit()
    qiskit_circuit = circuit.to_qiskit()

    # Run VQE optimization
    from qiskit.primitives import StatevectorEstimator
    estimator = StatevectorEstimator()

    def cost_function(params):
        bound_circuit = qiskit_circuit.assign_parameters(
            {list(qiskit_circuit.parameters)[0]: params[0]}
        )
        job = estimator.run([(bound_circuit, pauli_ham)])
        return job.result()[0].data.evs

    from scipy.optimize import minimize
    result = minimize(cost_function, x0=[0.0], method='COBYLA', options={'maxiter': 100})

    optimal_energy = result.fun

    # Validate results
    assert optimal_energy < -1.13, f"Energy {optimal_energy} not reaching correlation"
    assert abs(optimal_energy - fci_energy) < 0.001, \
        f"Error {abs(optimal_energy - fci_energy)*1000:.3f} mHa exceeds chemical accuracy"

    print(f"\n✓ UCC Double Excitation Test PASSED")
    print(f"  VQE energy: {optimal_energy:.6f} Ha")
    print(f"  FCI energy: {fci_energy:.6f} Ha")
    print(f"  Error:      {abs(optimal_energy - fci_energy)*1000:.3f} mHa")


def test_ucc_statevector_correctness():
    """
    Test that UCC gates produce correct statevector amplitudes.

    Validates the circuit creates the right superposition:
    - HF state |1100⟩ should mix with excited state |0011⟩
    - Amplitudes should match analytical Jordan-Wigner evolution
    """
    h1 = Atom('H', position=[0.0, 0.0, 0.0])
    h2 = Atom('H', position=[0.0, 0.0, 0.735])
    mol = Molecule(atoms=[h1, h2], charge=0, spin=0)

    # Build UCC circuit
    ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=False, include_doubles=True)
    circuit = ansatz.build_circuit()
    qiskit_circuit = circuit.to_qiskit()

    # Apply small rotation
    theta = 0.1
    bound_circuit = qiskit_circuit.assign_parameters({list(qiskit_circuit.parameters)[0]: theta})

    # Get statevector
    from qiskit.quantum_info import Statevector
    sv = Statevector.from_instruction(bound_circuit)

    # Expected amplitudes from direct exponentiation of Jordan-Wigner operator
    # exp(θ*T) where T is the anti-Hermitian double excitation
    expected_hf_amp = 0.995004  # cos-like behavior
    expected_excited_amp = 0.099833  # sin-like behavior

    actual_hf_amp = abs(sv.data[12])  # |1100⟩ = 12 in decimal
    actual_excited_amp = abs(sv.data[3])  # |0011⟩ = 3 in decimal

    assert abs(actual_hf_amp - expected_hf_amp) < 0.001, \
        f"HF amplitude {actual_hf_amp:.6f} != expected {expected_hf_amp:.6f}"
    assert abs(actual_excited_amp - expected_excited_amp) < 0.001, \
        f"Excited amplitude {actual_excited_amp:.6f} != expected {expected_excited_amp:.6f}"

    print(f"\n✓ UCC Statevector Test PASSED")
    print(f"  |1100⟩: {actual_hf_amp:.6f} (expected {expected_hf_amp:.6f})")
    print(f"  |0011⟩: {actual_excited_amp:.6f} (expected {expected_excited_amp:.6f})")


def test_ucc_energy_lowering():
    """
    Test that UCC excitation lowers energy below Hartree-Fock.

    This is the fundamental test: correlation energy should be negative.
    """
    h1 = Atom('H', position=[0.0, 0.0, 0.0])
    h2 = Atom('H', position=[0.0, 0.0, 0.735])
    mol = Molecule(atoms=[h1, h2], charge=0, spin=0)

    hamiltonian = mol.hamiltonian
    mo_energies, C = hamiltonian.compute_molecular_orbitals()
    h_mo = C.T @ hamiltonian.h_core @ C
    eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', C, C, hamiltonian.eri, C, C, optimize=True)
    pauli_ham = openfermion_jordan_wigner(h_mo, eri_mo, hamiltonian.nuclear_repulsion, hamiltonian.n_electrons)

    ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=False, include_doubles=True)
    circuit = ansatz.build_circuit()
    qiskit_circuit = circuit.to_qiskit()

    from qiskit.primitives import StatevectorEstimator
    estimator = StatevectorEstimator()

    # HF energy (θ=0)
    hf_circuit = qiskit_circuit.assign_parameters({list(qiskit_circuit.parameters)[0]: 0.0})
    job = estimator.run([(hf_circuit, pauli_ham)])
    energy_hf = job.result()[0].data.evs

    # Energy with excitation (θ=0.1)
    excited_circuit = qiskit_circuit.assign_parameters({list(qiskit_circuit.parameters)[0]: 0.1})
    job = estimator.run([(excited_circuit, pauli_ham)])
    energy_excited = job.result()[0].data.evs

    correlation_energy = energy_excited - energy_hf

    assert correlation_energy < -0.01, \
        f"Correlation energy {correlation_energy:.6f} not negative enough"

    print(f"\n✓ UCC Energy Lowering Test PASSED")
    print(f"  HF energy:          {energy_hf:.6f} Ha")
    print(f"  Correlated energy:  {energy_excited:.6f} Ha")
    print(f"  Correlation energy: {correlation_energy:.6f} Ha")


if __name__ == "__main__":
    """Run all tests."""
    print("="*70)
    print("UCC ANSATZ VALIDATION TESTS")
    print("="*70)

    test_ucc_statevector_correctness()
    test_ucc_energy_lowering()
    test_ucc_double_excitation_h2()

    print("\n" + "="*70)
    print("ALL TESTS PASSED! ✓✓✓")
    print("="*70)
