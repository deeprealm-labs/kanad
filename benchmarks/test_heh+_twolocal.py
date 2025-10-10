"""Test HeH+ with TwoLocal ansatz"""
import time
import numpy as np
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner
from kanad.solvers.vqe_solver import VQESolver
from kanad.ansatze.two_local_ansatz import TwoLocalAnsatz

mol = Molecule(atoms=[Atom('He', [0,0,0]), Atom('H', [0,0,0.775])], charge=1)

# FCI
mo_e, C = mol.hamiltonian.compute_molecular_orbitals()
h_mo = C.T @ mol.hamiltonian.h_core @ C
eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', C, C, mol.hamiltonian.eri, C, C)
pauli_h = openfermion_jordan_wigner(h_mo, eri_mo, mol.hamiltonian.nuclear_repulsion, 2 // 2)
E_fci = np.linalg.eigh(pauli_h.to_matrix())[0][0]

# VQE
ansatz = TwoLocalAnsatz(n_qubits=4, n_electrons=2, n_layers=2)
solver = VQESolver(hamiltonian=mol.hamiltonian, molecule=mol, ansatz=ansatz, optimizer='SLSQP', max_iterations=150)
t0 = time.time()
result = solver.solve()
dt = time.time() - t0

E_vqe = result['energy']
err = (E_vqe - E_fci) * 1000

print(f"HeH+/TwoLocal/SLSQP: E={E_vqe:.6f} Ha, error={err:+.3f} mHa, time={dt:.1f}s")
assert abs(err) < 50, f"Error too large: {err} mHa"
print("✓ PASS")
