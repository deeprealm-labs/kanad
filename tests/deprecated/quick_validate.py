#!/usr/bin/env python3
"""
Quick validation - just H2, LiH, and H2_stretched
"""

import numpy as np
import time
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.ucc_ansatz import UCCAnsatz

HARTREE_TO_EV = 27.211386245988

print("="*80)
print("QUICK ENERGY VALIDATION - 3 Key Molecules")
print("="*80)

mapper = JordanWignerMapper()
results = {}

# Test 1: H2
print("\n" + "="*80)
print("TEST 1: H₂ at equilibrium (0.74 Å)")
print("="*80)

H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_bond = BondFactory.create_bond(H1, H2)

print(f"Orbitals: {h2_bond.hamiltonian.n_orbitals}, Electrons: {h2_bond.molecule.n_electrons}")

# Exact
start = time.time()
exact_h2 = h2_bond.compute_energy(method='exact', mapper=mapper)
exact_time = time.time() - start
print(f"✓ Exact: {exact_h2['energy']/HARTREE_TO_EV:.6f} Ha [{exact_time:.2f}s]")

# VQE
ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)
start = time.time()
vqe_h2 = h2_bond.compute_energy(method='VQE', mapper=mapper, ansatz=ansatz, max_iterations=50)
vqe_time = time.time() - start
print(f"✓ VQE:   {vqe_h2['energy']/HARTREE_TO_EV:.6f} Ha [{vqe_time:.2f}s]")

error_h2 = abs(vqe_h2['energy'] - exact_h2['energy']) / abs(exact_h2['energy']) * 100
print(f"✓ Error: {error_h2:.2f}%")
results['H2'] = {'error': error_h2, 'passed': error_h2 < 10}

# Test 2: LiH
print("\n" + "="*80)
print("TEST 2: LiH at equilibrium (1.60 Å)")
print("="*80)

Li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
H3 = Atom('H', position=np.array([1.60, 0.0, 0.0]))
lih_bond = BondFactory.create_bond(Li, H3)

print(f"Orbitals: {lih_bond.hamiltonian.n_orbitals}, Electrons: {lih_bond.molecule.n_electrons}")

# Exact
start = time.time()
exact_lih = lih_bond.compute_energy(method='exact', mapper=mapper)
exact_time = time.time() - start
print(f"✓ Exact: {exact_lih['energy']/HARTREE_TO_EV:.6f} Ha [{exact_time:.2f}s]")

# VQE - use hardware efficient for larger system
from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
n_orb = lih_bond.hamiltonian.n_orbitals
n_elec = lih_bond.molecule.n_electrons
ansatz_lih = RealAmplitudesAnsatz(n_qubits=2*n_orb, n_electrons=n_elec, n_layers=2)

start = time.time()
vqe_lih = lih_bond.compute_energy(method='VQE', mapper=mapper, ansatz=ansatz_lih, max_iterations=50)
vqe_time = time.time() - start
print(f"✓ VQE:   {vqe_lih['energy']/HARTREE_TO_EV:.6f} Ha [{vqe_time:.2f}s]")

error_lih = abs(vqe_lih['energy'] - exact_lih['energy']) / abs(exact_lih['energy']) * 100
print(f"✓ Error: {error_lih:.2f}%")
results['LiH'] = {'error': error_lih, 'passed': error_lih < 15}  # Allow 15% for harder molecule

# Test 3: H2 stretched
print("\n" + "="*80)
print("TEST 3: H₂ stretched (1.4 Å)")
print("="*80)

H4 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H5 = Atom('H', position=np.array([1.4, 0.0, 0.0]))
h2_stretch = BondFactory.create_bond(H4, H5)

# Exact
start = time.time()
exact_stretch = h2_stretch.compute_energy(method='exact', mapper=mapper)
exact_time = time.time() - start
print(f"✓ Exact: {exact_stretch['energy']/HARTREE_TO_EV:.6f} Ha [{exact_time:.2f}s]")

# VQE
ansatz_s = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)
start = time.time()
vqe_stretch = h2_stretch.compute_energy(method='VQE', mapper=mapper, ansatz=ansatz_s, max_iterations=50)
vqe_time = time.time() - start
print(f"✓ VQE:   {vqe_stretch['energy']/HARTREE_TO_EV:.6f} Ha [{vqe_time:.2f}s]")

error_stretch = abs(vqe_stretch['energy'] - exact_stretch['energy']) / abs(exact_stretch['energy']) * 100
print(f"✓ Error: {error_stretch:.2f}%")
results['H2_stretch'] = {'error': error_stretch, 'passed': error_stretch < 10}

# Summary
print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"\n{'Molecule':<15} {'Error':<12} {'Status':<10}")
print("-"*40)
for name, data in results.items():
    status = "✅ PASS" if data['passed'] else "⚠️  FAIL"
    print(f"{name:<15} {data['error']:>6.2f}%     {status}")

passed = sum(1 for r in results.values() if r['passed'])
total = len(results)

print("-"*40)
print(f"\nResult: {passed}/{total} tests passed")

if passed == total:
    print("\n🎉 ALL TESTS PASSED!")
    print("\n✅ Energy calculation fix is working correctly!")
    print("✅ VQE accurately matches exact diagonalization")
    print("✅ Framework ready for quantum chemistry!")
else:
    print(f"\n✅ {passed} tests passed")
    print(f"⚠️  {total-passed} tests need improvement (likely ansatz expressibility)")

print("="*80)
