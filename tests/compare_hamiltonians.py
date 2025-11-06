#!/usr/bin/env python3
"""
Compare Hamiltonians from Bond API vs Direct API
Find the exact difference causing the bug
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from pyscf import gto

print("="*80)
print("HAMILTONIAN COMPARISON")
print("="*80)

# Method 1: Bond API (broken)
print("\n📋 METHOD 1: Bond API (Dashboard Path)")
print("-"*80)

from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom

atom1 = Atom('H', [0.0, 0.0, 0.0])
atom2 = Atom('H', [0.0, 0.0, 0.74])
bond = BondFactory.create_bond(atom1, atom2, distance=0.74, basis='sto-3g')

# Get Hamiltonian from Bond
bond_hamiltonian = bond.hamiltonian
print(f"Bond Hamiltonian type: {type(bond_hamiltonian).__name__}")
print(f"Nuclear repulsion: {bond_hamiltonian.nuclear_repulsion:.8f} Ha")
print(f"N orbitals: {bond_hamiltonian.n_orbitals}")
print(f"N electrons: {bond.molecule.n_electrons}")

# Convert to sparse Pauli
bond_sparse = bond_hamiltonian.to_sparse_hamiltonian(mapper='jordan_wigner')
print(f"\nBond Sparse Hamiltonian:")
print(f"  Pauli terms: {len(bond_sparse)}")
print(f"  Num qubits: {bond_sparse.num_qubits}")

# Print first few terms
print(f"  First 5 terms:")
for i in range(min(5, len(bond_sparse))):
    pauli_str = str(bond_sparse.paulis[i])
    coeff = bond_sparse.coeffs[i]
    print(f"    {pauli_str}: {coeff.real:.8f}")

# Method 2: Direct API (working)
print("\n📋 METHOD 2: Direct API (Test Path)")
print("-"*80)

from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g', charge=0, spin=0)
mol.build()

mf = mol.RHF().run(verbose=0)
mo_coeff = mf.mo_coeff

h1e = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
h_mo = mo_coeff.T @ h1e @ mo_coeff

eri_ao = mol.intor('int2e')
eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', mo_coeff, mo_coeff, eri_ao, mo_coeff, mo_coeff, optimize=True)

nuclear_repulsion = mol.energy_nuc()

direct_sparse = openfermion_jordan_wigner(
    h_mo=h_mo,
    eri_mo=eri_mo,
    nuclear_repulsion=nuclear_repulsion,
    n_electrons=mol.nelectron
)

print(f"Direct Sparse Hamiltonian:")
print(f"  Pauli terms: {len(direct_sparse)}")
print(f"  Num qubits: {direct_sparse.num_qubits}")
print(f"  Nuclear repulsion: {nuclear_repulsion:.8f} Ha")

print(f"  First 5 terms:")
for i in range(min(5, len(direct_sparse))):
    pauli_str = str(direct_sparse.paulis[i])
    coeff = direct_sparse.coeffs[i]
    print(f"    {pauli_str}: {coeff.real:.8f}")

# Comparison
print("\n" + "="*80)
print("COMPARISON")
print("="*80)

print(f"\nNumber of terms:")
print(f"  Bond API:   {len(bond_sparse)}")
print(f"  Direct API: {len(direct_sparse)}")

if len(bond_sparse) != len(direct_sparse):
    print(f"  ❌ MISMATCH! Different number of Pauli terms")
else:
    print(f"  ✓ Same number of terms")

print(f"\nNuclear repulsion:")
print(f"  Bond API:   {bond_hamiltonian.nuclear_repulsion:.8f} Ha")
print(f"  Direct API: {nuclear_repulsion:.8f} Ha")

if abs(bond_hamiltonian.nuclear_repulsion - nuclear_repulsion) > 1e-6:
    print(f"  ❌ MISMATCH! Difference: {abs(bond_hamiltonian.nuclear_repulsion - nuclear_repulsion):.8f} Ha")
else:
    print(f"  ✓ Match")

# Compare coefficients term by term
print(f"\nTerm-by-term comparison:")
matches = 0
mismatches = 0

for i in range(min(len(bond_sparse), len(direct_sparse))):
    bond_pauli = str(bond_sparse.paulis[i])
    bond_coeff = bond_sparse.coeffs[i].real

    direct_pauli = str(direct_sparse.paulis[i])
    direct_coeff = direct_sparse.coeffs[i].real

    if bond_pauli == direct_pauli and abs(bond_coeff - direct_coeff) < 1e-8:
        matches += 1
    else:
        mismatches += 1
        if mismatches <= 10:  # Print first 10 mismatches
            print(f"  Term {i}:")
            print(f"    Bond:   {bond_pauli}: {bond_coeff:.8f}")
            print(f"    Direct: {direct_pauli}: {direct_coeff:.8f}")

print(f"\nMatches: {matches}/{min(len(bond_sparse), len(direct_sparse))}")
print(f"Mismatches: {mismatches}/{min(len(bond_sparse), len(direct_sparse))}")

if mismatches == 0:
    print(f"\n✅ Hamiltonians are IDENTICAL!")
    print(f"   Bug must be elsewhere (ansatz? initial state?)")
else:
    print(f"\n❌ Hamiltonians are DIFFERENT!")
    print(f"   This is the source of the bug")

print("\n" + "="*80)
