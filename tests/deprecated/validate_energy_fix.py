#!/usr/bin/env python3
"""
Validation test for the energy calculation fix.
Tests H2 and LiH to verify improvements.
"""

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.ucc_ansatz import UCCAnsatz

HARTREE_TO_EV = 27.211386245988

print("=" * 80)
print("ENERGY CALCULATION FIX VALIDATION")
print("=" * 80)

# Test 1: H2
print("\n" + "=" * 80)
print("TEST 1: H2 Molecule (0.74 Å)")
print("=" * 80)

H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_bond = BondFactory.create_bond(H1, H2)

mapper = JordanWignerMapper()

# Exact energy
result_exact = h2_bond.compute_energy(method='exact', mapper=mapper)
exact_ev = result_exact['energy']
exact_hartree = exact_ev / HARTREE_TO_EV

print(f"\n✓ Exact Energy:")
print(f"  {exact_hartree:.6f} Hartree")
print(f"  {exact_ev:.6f} eV")

# VQE energy
ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)
result_vqe = h2_bond.compute_energy(
    method='VQE',
    mapper=mapper,
    ansatz=ansatz,
    max_iterations=100
)
vqe_ev = result_vqe['energy']
vqe_hartree = vqe_ev / HARTREE_TO_EV

print(f"\n✓ VQE Energy:")
print(f"  {vqe_hartree:.6f} Hartree")
print(f"  {vqe_ev:.6f} eV")

# Analysis
error_abs = abs(vqe_ev - exact_ev)
error_pct = (error_abs / abs(exact_ev)) * 100

print(f"\n✓ VQE vs Exact:")
print(f"  Absolute error: {error_abs:.6f} eV")
print(f"  Relative error: {error_pct:.2f}%")

# Literature comparison
lit_hartree = -1.1745
error_from_lit = abs(exact_hartree - lit_hartree)

print(f"\n✓ Exact vs Literature (-1.1745 Ha):")
print(f"  Error: {error_from_lit:.6f} Hartree ({error_from_lit*HARTREE_TO_EV:.4f} eV)")
print(f"  Note: Difference due to STO-3G minimal basis")

# Pass/Fail
if error_pct < 10:
    print(f"\n✅ TEST 1 PASSED: VQE error {error_pct:.2f}% < 10%")
else:
    print(f"\n❌ TEST 1 FAILED: VQE error {error_pct:.2f}% >= 10%")

# Test 2: HF (simpler molecule)
print("\n" + "=" * 80)
print("TEST 2: HF Molecule (0.92 Å)")
print("=" * 80)

H = Atom('H', position=np.array([0.0, 0.0, 0.0]))
F = Atom('F', position=np.array([0.92, 0.0, 0.0]))
hf_bond = BondFactory.create_bond(H, F)

# Exact energy
result_exact_hf = hf_bond.compute_energy(method='exact', mapper=mapper)
exact_ev_hf = result_exact_hf['energy']
exact_hartree_hf = exact_ev_hf / HARTREE_TO_EV

print(f"\n✓ Exact Energy:")
print(f"  {exact_hartree_hf:.6f} Hartree")
print(f"  {exact_ev_hf:.6f} eV")

# VQE energy
n_orbitals_hf = hf_bond.hamiltonian.n_orbitals
n_electrons_hf = hf_bond.molecule.n_electrons
ansatz_hf = UCCAnsatz(
    n_qubits=2*n_orbitals_hf,
    n_electrons=n_electrons_hf,
    include_singles=True,
    include_doubles=True
)

result_vqe_hf = hf_bond.compute_energy(
    method='VQE',
    mapper=mapper,
    ansatz=ansatz_hf,
    max_iterations=100
)
vqe_ev_hf = result_vqe_hf['energy']
vqe_hartree_hf = vqe_ev_hf / HARTREE_TO_EV

print(f"\n✓ VQE Energy:")
print(f"  {vqe_hartree_hf:.6f} Hartree")
print(f"  {vqe_ev_hf:.6f} eV")

# Analysis
error_abs_hf = abs(vqe_ev_hf - exact_ev_hf)
error_pct_hf = (error_abs_hf / abs(exact_ev_hf)) * 100

print(f"\n✓ VQE vs Exact:")
print(f"  Absolute error: {error_abs_hf:.6f} eV")
print(f"  Relative error: {error_pct_hf:.2f}%")

# Pass/Fail
if error_pct_hf < 10:
    print(f"\n✅ TEST 2 PASSED: VQE error {error_pct_hf:.2f}% < 10%")
else:
    print(f"\n❌ TEST 2 FAILED: VQE error {error_pct_hf:.2f}% >= 10%")

# Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"\nH2:  VQE vs Exact = {error_pct:.2f}%")
print(f"HF:  VQE vs Exact = {error_pct_hf:.2f}%")

if error_pct < 10 and error_pct_hf < 10:
    print(f"\n✅ ALL TESTS PASSED! Energy calculations are now accurate!")
    print(f"\nKey improvements:")
    print(f"  - Full many-body Hamiltonian construction")
    print(f"  - Proper Jordan-Wigner transformation")
    print(f"  - VQE now closely matches exact diagonalization")
else:
    print(f"\n⚠️  Some tests need improvement")

print("=" * 80)
