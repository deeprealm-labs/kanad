"""
Debug VQE state vs exact ground state.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.solvers.vqe_solver import VQESolver
from qiskit.quantum_info import Statevector

# Create H2
H1 = Atom('H', position=[0.0, 0.0, 0.0])
H2_atom = Atom('H', position=[0.74, 0.0, 0.0])
bond = CovalentBond(H1, H2_atom, basis='sto-3g')

# Get Hamiltonian matrix
H_matrix = bond.hamiltonian.to_matrix(n_qubits=4, use_mo_basis=True)
eigenvalues, eigenvectors = np.linalg.eigh(H_matrix)

print("Exact Solution:")
print(f"Ground state energy: {eigenvalues[0]:.8f} Ha")
print(f"Ground state vector (first 5 components):")
for i in range(min(5, len(eigenvectors[:, 0]))):
    print(f"  |{i:04b}⟩: {eigenvectors[i, 0]:.6f}")

# Run VQE with UCC
solver = VQESolver(
    bond=bond,
    ansatz_type='ucc',
    mapper_type='jordan_wigner',
    optimizer='SLSQP',
    max_iterations=100
)

result = solver.solve()
print(f"\nVQE Solution:")
print(f"Energy: {result['energy']:.8f} Ha")
print(f"Parameters: {result['parameters']}")

# Get VQE state
solver.ansatz.circuit.bind_parameters(result['parameters'])
qiskit_circuit = solver.ansatz.circuit.to_qiskit()
vqe_state = Statevector.from_instruction(qiskit_circuit)

print(f"\nVQE state vector (first 5 components):")
for i in range(min(5, len(vqe_state.data))):
    print(f"  |{i:04b}⟩: {vqe_state.data[i]:.6f}")

# Compute overlap
overlap = np.abs(np.vdot(eigenvectors[:, 0], vqe_state.data))
print(f"\nOverlap with exact ground state: {overlap:.6f}")

# Compute VQE energy manually
E_vqe_manual = np.real(np.vdot(vqe_state.data, H_matrix @ vqe_state.data))
print(f"\nManual energy calculation: {E_vqe_manual:.8f} Ha")
print(f"Difference from VQE result: {(E_vqe_manual - result['energy'])*1000:.6f} mHa")
