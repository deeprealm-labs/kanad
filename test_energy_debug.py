#!/usr/bin/env python3
"""Debug energy calculations to find the root cause."""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.constants.conversion_factors import ConversionFactors

# Create H2
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
bond = CovalentBond(H1, H2)

print("="*80)
print("  ENERGY DEBUG - H2 Molecule at 0.74 Å")
print("="*80)

# Direct Hamiltonian access
print("\n1. Direct Hamiltonian SCF:")
density_matrix, hf_energy_ha = bond.hamiltonian.solve_scf(max_iterations=100)
print(f"   Raw HF energy: {hf_energy_ha:.6f} Ha")
print(f"   Converted:     {hf_energy_ha * ConversionFactors.HARTREE_TO_EV:.6f} eV")

# Exact diagonalization
print("\n2. Exact Diagonalization:")
H_matrix = bond.hamiltonian.to_matrix()
eigenvalues, eigenvectors = np.linalg.eigh(H_matrix)
print(f"   Raw eigenvalue: {eigenvalues[0]:.6f} Ha")
print(f"   Converted:      {eigenvalues[0] * ConversionFactors.HARTREE_TO_EV:.6f} eV")

# Through bond interface (HF)
print("\n3. Through Bond.compute_energy (HF):")
result_hf = bond.compute_energy(method='hf', max_iterations=100)
print(f"   Returned energy: {result_hf['energy']:.6f} eV")
print(f"   Converged: {result_hf['converged']}")

# Through bond interface (Exact)
print("\n4. Through Bond.compute_energy (Exact):")
mapper = JordanWignerMapper()
result_exact = bond.compute_energy(method='exact', mapper=mapper)
print(f"   Returned energy: {result_exact['energy']:.6f} eV")

# Check nuclear repulsion
print("\n5. Nuclear Repulsion:")
print(f"   V_nn: {bond.hamiltonian.nuclear_repulsion:.6f} Ha")
print(f"   V_nn: {bond.hamiltonian.nuclear_repulsion * ConversionFactors.HARTREE_TO_EV:.6f} eV")

# Reference values
print("\n6. Reference H2 Energies (literature):")
print(f"   HF energy:    -1.133 Ha = -30.82 eV")
print(f"   Exact energy: -1.174 Ha = -31.95 eV")
print(f"   Full CI (extended basis): -2.08 Ha = -56.61 eV")

print("\n7. Analysis:")
print(f"   Our HF:    {hf_energy_ha:.6f} Ha vs -1.133 Ha (expected)")
print(f"   Our Exact: {eigenvalues[0]:.6f} Ha vs -2.08 Ha (expected for STO-3G FCI)")
print(f"   Discrepancy: HF is {abs(hf_energy_ha - (-1.133)):.6f} Ha off")
print(f"   Discrepancy: Exact is {abs(eigenvalues[0] - (-2.08)):.6f} Ha off")

print("\n8. Basis Set Info:")
print(f"   Basis: {bond.hamiltonian.basis.name}")
print(f"   # basis functions: {len(bond.hamiltonian.basis.basis_functions)}")
print(f"   # orbitals: {bond.hamiltonian.n_orbitals}")
print(f"   # electrons: {bond.molecule.n_electrons}")

# Check if Hamiltonian matrix looks reasonable
print("\n9. Hamiltonian Matrix:")
print(f"   Shape: {H_matrix.shape}")
print(f"   Is Hermitian: {np.allclose(H_matrix, H_matrix.conj().T)}")
print(f"   Min eigenvalue: {eigenvalues.min():.6f} Ha")
print(f"   Max eigenvalue: {eigenvalues.max():.6f} Ha")
print(f"   All eigenvalues (Ha): {eigenvalues}")

print("\n" + "="*80)
