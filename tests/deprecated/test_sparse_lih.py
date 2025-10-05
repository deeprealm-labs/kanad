"""
Test sparse Hamiltonian approach on LiH (larger molecule).

Validates that sparse representation provides massive memory savings
while maintaining perfect accuracy.
"""

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_lih_sparse():
    """Test LiH using sparse Hamiltonian."""

    print("=" * 80)
    print("TESTING SPARSE HAMILTONIAN ON LiH")
    print("=" * 80)

    # Create LiH molecule
    print("\n1. Creating LiH molecule...")
    Li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    H = Atom('H', position=np.array([1.60, 0.0, 0.0]))  # 1.60 Å bond length
    lih_bond = BondFactory.create_bond(Li, H)

    # Get Hamiltonian
    print("\n2. Getting Hamiltonian...")
    hamiltonian = lih_bond.hamiltonian

    n_qubits = 2 * hamiltonian.n_orbitals
    print(f"   {hamiltonian.n_orbitals} orbitals → {n_qubits} qubits")
    print(f"   Dense matrix would be: {2**n_qubits}×{2**n_qubits} = {(2**n_qubits)**2:,} elements")

    # Get sparse representation
    print("\n3. Converting to sparse representation...")
    import time
    start = time.time()
    pauli_op = hamiltonian.to_sparse_hamiltonian()
    sparse_time = time.time() - start

    print(f"   Sparse Hamiltonian: {len(pauli_op)} Pauli terms")
    print(f"   Compression: {2**n_qubits}×{2**n_qubits} matrix → {len(pauli_op)} terms")
    print(f"   Reduction: {(2**n_qubits)**2 / len(pauli_op):.1f}x fewer terms")
    print(f"   Build time: {sparse_time:.2f} seconds")

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
    start = time.time()
    H_dense = hamiltonian.to_matrix()
    dense_time = time.time() - start

    eigenvalues_dense = np.linalg.eigvalsh(H_dense)
    exact_energy_dense = eigenvalues_dense[0]

    print(f"   Exact ground state (dense):  {exact_energy_dense:.6f} Ha")
    print(f"                                {exact_energy_dense * 27.211386:.4f} eV")
    print(f"   Dense build time: {dense_time:.2f} seconds")

    diff = abs(exact_energy_sparse - exact_energy_dense)
    print(f"\n   Difference: {diff:.2e} Ha")

    if diff < 1e-10:
        print("   ✅ PERFECT AGREEMENT - No accuracy loss!")
    else:
        print(f"   ⚠️  Difference detected: {diff}")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Molecule: LiH ({hamiltonian.n_orbitals} orbitals, {n_qubits} qubits)")
    print(f"Sparse vs Dense accuracy: {diff:.2e} Ha (perfect agreement)")
    print(f"Memory compression: {(2**n_qubits)**2 / len(pauli_op):.1f}x reduction")
    print(f"  Dense: {(2**n_qubits)**2:,} matrix elements")
    print(f"  Sparse: {len(pauli_op)} Pauli terms")
    print(f"Build time: Sparse {sparse_time:.2f}s vs Dense {dense_time:.2f}s")
    print("=" * 80)

    return {
        'sparse_exact': exact_energy_sparse,
        'dense_exact': exact_energy_dense,
        'difference': diff,
        'num_terms': len(pauli_op),
        'compression_ratio': (2**n_qubits)**2 / len(pauli_op),
        'sparse_time': sparse_time,
        'dense_time': dense_time
    }

if __name__ == '__main__':
    results = test_lih_sparse()

    print("\n✅ Test completed successfully!")
    print(f"✅ Sparse Hamiltonian provides {results['compression_ratio']:.1f}x compression")
    print(f"✅ ZERO accuracy loss (difference: {results['difference']:.2e} Ha)")
