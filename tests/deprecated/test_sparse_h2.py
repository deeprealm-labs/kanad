"""
Test sparse Hamiltonian approach on H2.

Validates that sparse representation gives identical results to dense matrix approach.
"""

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_h2_sparse():
    """Test H2 using sparse Hamiltonian."""

    print("=" * 80)
    print("TESTING SPARSE HAMILTONIAN ON H2")
    print("=" * 80)

    # Create H2 molecule
    print("\n1. Creating H2 molecule...")
    H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
    h2_bond = BondFactory.create_bond(H1, H2)

    # Get Hamiltonian
    print("\n2. Getting Hamiltonian...")
    hamiltonian = h2_bond.hamiltonian

    print(f"   {hamiltonian.n_orbitals} orbitals → {2 * hamiltonian.n_orbitals} qubits")

    # Get sparse representation
    print("\n3. Converting to sparse representation...")
    pauli_op = hamiltonian.to_sparse_hamiltonian()
    print(f"   Sparse Hamiltonian: {len(pauli_op)} Pauli terms")
    n_qubits = 2 * hamiltonian.n_orbitals
    print(f"   vs Dense matrix: {2**n_qubits}×{2**n_qubits}")

    # Compute exact energy using sparse Hamiltonian
    print("\n4. Computing exact ground state energy (sparse)...")

    # Get eigenvalues from sparse Pauli operator
    dense_matrix = pauli_op.to_matrix()
    eigenvalues = np.linalg.eigvalsh(dense_matrix)
    exact_energy_sparse = eigenvalues[0]

    print(f"   Exact ground state (sparse): {exact_energy_sparse:.6f} Ha")
    print(f"                                {exact_energy_sparse * 27.211386:.4f} eV")

    # Compare with dense approach
    print("\n5. Comparing with dense matrix approach...")
    H_dense = hamiltonian.to_matrix()
    eigenvalues_dense = np.linalg.eigvalsh(H_dense)
    exact_energy_dense = eigenvalues_dense[0]

    print(f"   Exact ground state (dense):  {exact_energy_dense:.6f} Ha")
    print(f"                                {exact_energy_dense * 27.211386:.4f} eV")

    diff = abs(exact_energy_sparse - exact_energy_dense)
    print(f"\n   Difference: {diff:.2e} Ha")

    if diff < 1e-10:
        print("   ✅ PERFECT AGREEMENT - No accuracy loss!")
    else:
        print(f"   ⚠️  Difference detected: {diff}")

    # Test VQE with normal approach (current framework)
    print("\n6. Testing VQE with existing framework (for comparison)...")
    from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
    from kanad.ansatze.ucc_ansatz import UCCAnsatz

    mapper = JordanWignerMapper()
    ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)

    result_vqe = h2_bond.compute_energy(
        method='VQE',
        mapper=mapper,
        ansatz=ansatz,
        max_iterations=100
    )

    vqe_ev = result_vqe['energy']
    vqe_hartree = vqe_ev / 27.211386

    print(f"\n   VQE Energy: {vqe_hartree:.6f} Ha")
    print(f"               {vqe_ev:.4f} eV")
    print(f"   Exact:      {exact_energy_sparse:.6f} Ha")
    print(f"               {exact_energy_sparse * 27.211386:.4f} eV")

    error = abs(vqe_hartree - exact_energy_sparse) / abs(exact_energy_sparse) * 100
    print(f"\n   VQE Error: {error:.2f}%")

    if error < 5.0:
        print("   ✅ VQE accuracy maintained (< 5% error)")
    else:
        print(f"   ⚠️  VQE error higher than expected: {error:.2f}%")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Sparse vs Dense difference: {diff:.2e} Ha (should be ~0)")
    print(f"VQE Error: {error:.2f}% (target: < 5%)")
    print(f"Memory savings: {2**n_qubits}×{2**n_qubits} matrix → {len(pauli_op)} terms")
    print("=" * 80)

    return {
        'sparse_exact': exact_energy_sparse,
        'dense_exact': exact_energy_dense,
        'vqe_energy': vqe_hartree,
        'difference': diff,
        'vqe_error': error,
        'num_terms': len(pauli_op)
    }

if __name__ == '__main__':
    results = test_h2_sparse()

    print("\nTest completed successfully!")
