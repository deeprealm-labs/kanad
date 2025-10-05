#!/usr/bin/env python3
"""
Method Comparison Validation - HF vs FCI vs VQE

Compare different quantum chemistry methods on H2 to identify bugs.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.core.constants.conversion_factors import ConversionFactors

# Create H2 molecule
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2_atom = Atom('H', position=np.array([0.74, 0.0, 0.0]))

bond = CovalentBond(H1, H2_atom)

print("="*80)
print("  H2 MOLECULE - METHOD COMPARISON")
print("="*80)
print(f"\nBond length: 0.74 Å")
print(f"Basis set: STO-3G")
print(f"n_orbitals: {bond.hamiltonian.n_orbitals}")
print(f"n_electrons: {bond.hamiltonian.n_electrons}")

# Literature reference values
print("\n" + "-"*80)
print("LITERATURE REFERENCE VALUES (H2, R=0.74 Å, STO-3G):")
print("-"*80)
print(f"  HF energy:  -30.39 eV  (-1.1167 Ha)")
print(f"  FCI energy: -30.95 eV  (-1.1373 Ha)")
print(f"  Correlation: -0.56 eV  (-0.0206 Ha)")

# 1. Hartree-Fock
print("\n" + "="*80)
print("METHOD 1: HARTREE-FOCK (HF)")
print("="*80)

hf_result = bond.compute_energy(method='hf')
print(f"\nHF Energy: {hf_result['energy']:.4f} eV")
print(f"Converged: {hf_result['converged']} ({hf_result['iterations']} iterations)")
print(f"Error vs literature: {abs(hf_result['energy'] + 30.39):.4f} eV")

# Get MO info from HF
h_core = bond.hamiltonian.h_core
eri = bond.hamiltonian.eri
V_nn = bond.hamiltonian.nuclear_repulsion

print(f"\nNuclear repulsion: {V_nn:.6f} Ha ({V_nn * ConversionFactors.HARTREE_TO_EV:.4f} eV)")
print(f"\nOne-electron Hamiltonian (h_core):")
print(h_core)
print(f"\nTwo-electron integrals (eri) shape: {eri.shape}")
print(f"  eri[0,0,0,0] = {eri[0,0,0,0]:.6f}")
print(f"  eri[0,1,0,1] = {eri[0,1,0,1]:.6f}")
print(f"  eri[0,0,1,1] = {eri[0,0,1,1]:.6f}")

# 2. Full Configuration Interaction
print("\n" + "="*80)
print("METHOD 2: FULL CONFIGURATION INTERACTION (FCI)")
print("="*80)

from kanad.solvers.fci_solver import FCISolver

fci = FCISolver(bond.hamiltonian)
fci_result = fci.solve()

E_fci_electronic = fci_result['energies'][0]
E_fci_total = E_fci_electronic + V_nn

print(f"\nFCI Determinants: {fci_result['determinants']}")
print(f"CI Dimension: {fci_result['ci_dimension']}")
print(f"\nFCI Hamiltonian matrix (electronic only):")
print(fci_result['hamiltonian_matrix'])

print(f"\nFCI Electronic Energy: {E_fci_electronic:.6f} Ha ({E_fci_electronic * ConversionFactors.HARTREE_TO_EV:.4f} eV)")
print(f"FCI Total Energy: {E_fci_total:.6f} Ha ({E_fci_total * ConversionFactors.HARTREE_TO_EV:.4f} eV)")
print(f"\nError vs literature FCI: {abs(E_fci_total * ConversionFactors.HARTREE_TO_EV + 30.95):.4f} eV")

correlation_energy = E_fci_total * ConversionFactors.HARTREE_TO_EV - hf_result['energy']
print(f"Correlation energy: {correlation_energy:.4f} eV")
print(f"Expected correlation: -0.56 eV")

# 3. VQE (if time permits)
print("\n" + "="*80)
print("METHOD 3: VARIATIONAL QUANTUM EIGENSOLVER (VQE)")
print("="*80)

try:
    vqe_result = bond.compute_energy(method='VQE', max_iterations=200)
    print(f"\nVQE Energy: {vqe_result['energy']:.4f} eV")
    print(f"Converged: {vqe_result.get('converged', 'unknown')}")
    print(f"Iterations: {vqe_result.get('iterations', 'unknown')}")
    print(f"Error vs literature FCI: {abs(vqe_result['energy'] + 30.95):.4f} eV")

    vqe_correlation = vqe_result['energy'] - hf_result['energy']
    print(f"VQE Correlation energy: {vqe_correlation:.4f} eV")
except Exception as e:
    print(f"VQE failed: {str(e)[:100]}")

# 4. "Exact" diagonalization (for comparison)
print("\n" + "="*80)
print("METHOD 4: EXACT DIAGONALIZATION")
print("="*80)

try:
    exact_result = bond.compute_energy(method='exact')
    print(f"\nExact Energy: {exact_result['energy']:.4f} eV")
    print(f"Method: {exact_result.get('method', 'unknown')}")
    print(f"Error vs literature FCI: {abs(exact_result['energy'] + 30.95):.4f} eV")
except Exception as e:
    print(f"Exact failed: {str(e)[:100]}")

# Summary
print("\n" + "="*80)
print("SUMMARY")
print("="*80)

print(f"\n{'Method':<25s} {'Energy (eV)':<15s} {'Error vs Lit.':<15s} {'Status':<10s}")
print("-"*70)
print(f"{'Literature HF':<25s} {'-30.39':<15s} {'---':<15s} {'Reference':<10s}")
print(f"{'Literature FCI':<25s} {'-30.95':<15s} {'---':<15s} {'Reference':<10s}")
print(f"{'Kanad HF':<25s} {hf_result['energy']:<15.4f} {abs(hf_result['energy'] + 30.39):<15.4f} {'✅' if abs(hf_result['energy'] + 30.39) < 0.1 else '❌':<10s}")
print(f"{'Kanad FCI':<25s} {E_fci_total * ConversionFactors.HARTREE_TO_EV:<15.4f} {abs(E_fci_total * ConversionFactors.HARTREE_TO_EV + 30.95):<15.4f} {'✅' if abs(E_fci_total * ConversionFactors.HARTREE_TO_EV + 30.95) < 0.5 else '❌':<10s}")

print("\n" + "="*80)
print("KEY FINDINGS:")
print("="*80)

if abs(hf_result['energy'] + 30.39) < 0.1:
    print("✅ HF is ACCURATE (< 0.1 eV error)")
else:
    print(f"❌ HF has {abs(hf_result['energy'] + 30.39):.2f} eV error")

if abs(E_fci_total * ConversionFactors.HARTREE_TO_EV + 30.95) < 0.5:
    print("✅ FCI is ACCURATE (< 0.5 eV error)")
else:
    print(f"❌ FCI has {abs(E_fci_total * ConversionFactors.HARTREE_TO_EV + 30.95):.2f} eV error")
    print("\nPossible FCI issues:")
    print("  1. Matrix element formulas wrong for restricted closed-shell")
    print("  2. Should transform to MO basis before CI")
    print("  3. Two-electron integral indexing issue")

print("\n" + "="*80)
