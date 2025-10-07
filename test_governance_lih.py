#!/usr/bin/env python3
"""
Test Governance on LiH - Larger molecule with ionic character
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.solvers.vqe_solver import VQESolver

print("="*80)
print("TEST: Governance on LiH")
print("="*80)

# Create LiH bond (ionic bonding)
print("\n📊 Creating LiH bond...")
try:
    lih = BondFactory.create_bond('Li', 'H', distance=1.60, basis='sto-3g')
    ham = lih.hamiltonian

    print(f"✓ LiH created")
    print(f"  n_orbitals: {ham.n_orbitals}")
    print(f"  n_electrons: {ham.n_electrons}")
    print(f"  n_qubits: {2 * ham.n_orbitals}")
    print(f"  use_governance: {ham.use_governance}")

    # Get HF reference
    print(f"\n📊 Computing Hartree-Fock reference...")
    dm, hf_energy = ham.solve_scf()
    print(f"  HF Energy: {hf_energy:.8f} Ha")

    # Get exact ground state via diagonalization
    print(f"\n📊 Computing exact ground state...")
    ham_matrix = ham.to_matrix(n_qubits=2*ham.n_orbitals)
    eigenvalues = np.linalg.eigvalsh(ham_matrix)
    exact_energy = eigenvalues[0]

    print(f"  Exact ground state: {exact_energy:.8f} Ha")
    print(f"  Correlation energy: {exact_energy - hf_energy:.8f} Ha ({abs(exact_energy - hf_energy)/abs(hf_energy)*100:.3f}%)")

    # Test governance ansatz
    print(f"\n📊 Testing Governance VQE...")
    vqe_gov = VQESolver(
        hamiltonian=ham,
        ansatz_type='governance',
        mapper_type='jordan_wigner',
        optimizer='SLSQP',
        max_iterations=500
    )

    print(f"  Ansatz: {type(vqe_gov.ansatz).__name__}")
    print(f"  Parameters: {vqe_gov.ansatz.n_parameters}")

    print(f"\n  Running VQE optimization...")
    result_gov = vqe_gov.solve()

    vqe_energy = result_gov['energy']
    vqe_error = abs(vqe_energy - exact_energy)

    print(f"\n  Results:")
    print(f"    VQE Energy:    {vqe_energy:.8f} Ha")
    print(f"    Exact Energy:  {exact_energy:.8f} Ha")
    print(f"    Error:         {vqe_error:.8f} Ha ({vqe_error*1000:.3f} mHa)")
    print(f"    Converged:     {result_gov.get('converged', False)}")
    print(f"    Iterations:    {result_gov.get('iterations', 'N/A')}")

    if vqe_error < 0.001:
        print(f"    ✓✓ EXCELLENT accuracy!")
    elif vqe_error < 0.01:
        print(f"    ✓ Good accuracy")
    elif vqe_error < 0.1:
        print(f"    ⚠️ Moderate accuracy")
    else:
        print(f"    ✗ Poor accuracy")

    # Compare with standard ansatz
    print(f"\n📊 Testing Hardware-Efficient VQE for comparison...")
    try:
        vqe_he = VQESolver(
            hamiltonian=ham,
            ansatz_type='hardware_efficient',
            mapper_type='jordan_wigner',
            optimizer='SLSQP',
            max_iterations=500
        )

        print(f"  Ansatz: {type(vqe_he.ansatz).__name__}")
        print(f"  Parameters: {vqe_he.ansatz.n_parameters}")

        print(f"\n  Running VQE optimization...")
        result_he = vqe_he.solve()

        he_energy = result_he['energy']
        he_error = abs(he_energy - exact_energy)

        print(f"\n  Results:")
        print(f"    VQE Energy:    {he_energy:.8f} Ha")
        print(f"    Error:         {he_error:.8f} Ha ({he_error*1000:.3f} mHa)")
        print(f"    Converged:     {result_he.get('converged', False)}")

        # Comparison
        print(f"\n📊 Comparison:")
        print(f"  Governance: {vqe_error*1000:.3f} mHa error, {vqe_gov.ansatz.n_parameters} params")
        print(f"  Hardware-Eff: {he_error*1000:.3f} mHa error, {vqe_he.ansatz.n_parameters} params")

        if vqe_error < he_error:
            print(f"  ✓ Governance is MORE accurate!")
        elif vqe_error < he_error * 1.5:
            print(f"  ≈ Similar accuracy")
        else:
            print(f"  ✗ Hardware-efficient is more accurate")

    except Exception as e:
        print(f"  ✗ Hardware-efficient failed: {e}")

    # Test with exact solver (NumPy)
    print(f"\n📊 Testing NumPy Exact Solver...")
    try:
        from kanad.solvers.numpy_solver import NumPySolver

        numpy_solver = NumPySolver(bond=lih)
        numpy_result = numpy_solver.solve()

        print(f"  NumPy Energy: {numpy_result['energy']:.8f} Ha")
        print(f"  Matches exact: {abs(numpy_result['energy'] - exact_energy) < 1e-6}")

    except Exception as e:
        print(f"  ✗ NumPy solver failed: {e}")

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

print(f"\n{'='*80}")
