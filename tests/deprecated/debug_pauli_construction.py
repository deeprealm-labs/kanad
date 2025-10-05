"""Debug Pauli construction."""

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.core.hamiltonians.fast_pauli_builder import build_molecular_hamiltonian_pauli

# Create H2
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_bond = BondFactory.create_bond(H1, H2)
hamiltonian = h2_bond.hamiltonian

print("Building sparse Hamiltonian...")
sparse_op = build_molecular_hamiltonian_pauli(
    h_core=hamiltonian.h_core,
    eri=hamiltonian.eri,
    nuclear_repulsion=hamiltonian.nuclear_repulsion,
    n_orbitals=hamiltonian.n_orbitals
)

print(f"\nSparse Pauli operator: {len(sparse_op)} terms")

# Get exact energy from sparse
sparse_matrix = sparse_op.to_matrix()
sparse_eigs = np.linalg.eigvalsh(sparse_matrix)
sparse_exact = sparse_eigs[0]

# Get exact energy from dense
dense_matrix = hamiltonian.to_matrix()
dense_eigs = np.linalg.eigvalsh(dense_matrix)
dense_exact = dense_eigs[0]

print(f"\nSparse exact energy: {sparse_exact:.6f} Ha = {sparse_exact * 27.211386:.4f} eV")
print(f"Dense exact energy:  {dense_exact:.6f} Ha = {dense_exact * 27.211386:.4f} eV")
print(f"Difference: {abs(sparse_exact - dense_exact):.2e} Ha")

# Check if matrices are close
print(f"\nMatrix shape: {sparse_matrix.shape}")
print(f"Matrix close: {np.allclose(sparse_matrix, dense_matrix)}")

if not np.allclose(sparse_matrix, dense_matrix):
    diff = sparse_matrix - dense_matrix
    max_diff = np.max(np.abs(diff))
    print(f"Max difference: {max_diff}")

    # Show where the differences are
    where_diff = np.where(np.abs(diff) > 1e-6)
    print(f"Number of different elements: {len(where_diff[0])}")
    if len(where_diff[0]) < 20:
        for i, j in zip(where_diff[0], where_diff[1]):
            print(f"  [{i},{j}]: sparse={sparse_matrix[i,j]:.6f}, dense={dense_matrix[i,j]:.6f}, diff={diff[i,j]:.6f}")
