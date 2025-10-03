#!/usr/bin/env python3
"""Compare Hamiltonian matrices from different sources."""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.solvers.vqe_solver import VQESolver

# Create H2
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
bond = CovalentBond(H1, H2)

print("="*80)
print("  HAMILTONIAN MATRIX COMPARISON")
print("="*80)

# Method 1: Direct from Hamiltonian
print("\n1. Hamiltonian.to_matrix():")
H_direct = bond.hamiltonian.to_matrix()
print(f"   Shape: {H_direct.shape}")
print(f"   Min eigenvalue: {np.linalg.eigh(H_direct)[0].min():.6f} Ha")
print(f"   Max eigenvalue: {np.linalg.eigh(H_direct)[0].max():.6f} Ha")

# Method 2: From VQE solver
print("\n2. VQESolver._build_hamiltonian_matrix():")
mapper = JordanWignerMapper()
n_orbitals = bond.hamiltonian.n_orbitals
n_qubits = 2 * n_orbitals
n_electrons = bond.molecule.n_electrons

ansatz = UCCAnsatz(n_qubits=n_qubits, n_electrons=n_electrons,
                    include_singles=True, include_doubles=True)

vqe = VQESolver(
    hamiltonian=bond.hamiltonian,
    ansatz=ansatz,
    mapper=mapper,
    backend='classical'
)

H_vqe = vqe._build_hamiltonian_matrix()
print(f"   Shape: {H_vqe.shape}")
print(f"   Min eigenvalue: {np.linalg.eigh(H_vqe)[0].min():.6f} Ha")
print(f"   Max eigenvalue: {np.linalg.eigh(H_vqe)[0].max():.6f} Ha")

# Compare
print("\n3. Comparison:")
print(f"   Are they equal? {np.allclose(H_direct, H_vqe, atol=1e-8)}")
if not np.allclose(H_direct, H_vqe, atol=1e-8):
    diff = np.abs(H_direct - H_vqe)
    print(f"   Max difference: {diff.max():.6e}")
    print(f"   Average difference: {diff.mean():.6e}")
    print(f"   Non-zero elements in diff: {np.sum(diff > 1e-10)}")

    # Check dimensions
    if H_direct.shape != H_vqe.shape:
        print(f"   ❌ Different shapes!")
    else:
        # Check if eigenvalues match
        evals_direct = np.linalg.eigh(H_direct)[0]
        evals_vqe = np.linalg.eigh(H_vqe)[0]
        print(f"   Ground state from direct: {evals_direct[0]:.6f} Ha")
        print(f"   Ground state from VQE H:  {evals_vqe[0]:.6f} Ha")
        print(f"   Difference: {abs(evals_direct[0] - evals_vqe[0]):.6e} Ha")

        if abs(evals_direct[0] - evals_vqe[0]) > 0.01:
            print(f"\n   ❌ CRITICAL: Ground state energies differ significantly!")
            print(f"   This explains the VQE error!")

# Test VQE energy at HF state (should give HF energy)
print("\n4. VQE at Hartree-Fock initial state:")
hf_params = np.zeros(vqe.n_parameters)
vqe_energy_hf = vqe.compute_energy(hf_params)
print(f"   VQE energy at θ=0: {vqe_energy_hf:.6f} Ha")

# Get actual HF energy
_, hf_energy = bond.hamiltonian.solve_scf(max_iterations=100)
print(f"   True HF energy:    {hf_energy:.6f} Ha")
print(f"   Difference:        {abs(vqe_energy_hf - hf_energy):.6f} Ha")

print("\n" + "="*80)
