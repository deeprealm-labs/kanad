#!/usr/bin/env python3
"""
Debug energy calculations to find the issue.
"""

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.ucc_ansatz import UCCAnsatz

print("="*80)
print("DEBUGGING H2 ENERGY CALCULATION")
print("="*80)

# Create H2 at equilibrium bond length
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))  # Angstroms
h2_bond = BondFactory.create_bond(H1, H2)

print(f"\n1. Molecule Info:")
print(f"   Bond length: 0.74 Å")
print(f"   Electrons: {h2_bond.molecule.n_electrons}")

# Get the Hamiltonian
hamiltonian = h2_bond.hamiltonian
print(f"\n2. Hamiltonian Info:")
print(f"   Orbitals: {hamiltonian.n_orbitals}")
if hasattr(hamiltonian, 'nuclear_repulsion'):
    print(f"   Nuclear repulsion: {hamiltonian.nuclear_repulsion:.6f} Hartree")

# Compute exact energy
print(f"\n3. Computing Exact Energy...")
mapper = JordanWignerMapper()
result_exact = h2_bond.compute_energy(method='exact', mapper=mapper)

print(f"   Exact energy: {result_exact['energy']:.6f} eV")

# Convert to Hartree for comparison
HARTREE_TO_EV = 27.211386245988
exact_hartree = result_exact['energy'] / HARTREE_TO_EV

print(f"   Exact energy: {exact_hartree:.6f} Hartree")

# Known H2 ground state at 0.74 Å is approximately -1.1745 Hartree
print(f"\n   Expected H2 ground state: ~-1.1745 Hartree (~-31.96 eV)")
print(f"   Our calculation: {exact_hartree:.6f} Hartree ({result_exact['energy']:.6f} eV)")
print(f"   Difference: {abs(exact_hartree - (-1.1745)):.6f} Hartree")

# Compute VQE energy
print(f"\n4. Computing VQE Energy...")
ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)
result_vqe = h2_bond.compute_energy(
    method='VQE',
    mapper=mapper,
    ansatz=ansatz,
    max_iterations=100  # More iterations
)

print(f"   VQE energy: {result_vqe['energy']:.6f} eV")
vqe_hartree = result_vqe['energy'] / HARTREE_TO_EV
print(f"   VQE energy: {vqe_hartree:.6f} Hartree")

# Check error
error_from_exact = abs(result_vqe['energy'] - result_exact['energy'])
error_percent = (error_from_exact / abs(result_exact['energy'])) * 100

print(f"\n5. Error Analysis:")
print(f"   VQE vs Exact error: {error_from_exact:.6f} eV ({error_percent:.2f}%)")
print(f"   VQE vs Expected error: {abs(vqe_hartree - (-1.1745)):.6f} Hartree")

# Check if there's a unit conversion issue
print(f"\n6. Unit Check:")
print(f"   If VQE is in wrong units by factor:")
if result_exact['energy'] != 0:
    factor = result_vqe['energy'] / result_exact['energy']
    print(f"   Factor = {factor:.6f}")
    print(f"   VQE * {1/factor:.6f} = {result_vqe['energy'] / factor:.6f} eV")

print(f"\n" + "="*80)
