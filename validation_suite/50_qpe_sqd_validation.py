#!/usr/bin/env python3
"""
QPE and SQD Solver Validation

Tests Quantum Phase Estimation (QPE) and Sample-Based Quantum Diagonalization (SQD)
solvers on H2 molecule.

Expected behavior:
- QPE should give FCI-level accuracy (exact diagonalization in classical fallback)
- SQD should give improved energy over HF by including excited configurations
- Both should converge to energies lower than HF baseline
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds import CovalentBond
from kanad.solvers import QPESolver, SQDSolver
import time

print("=" * 80)
print("  QPE AND SQD SOLVER VALIDATION")
print("=" * 80)
print("\nTesting advanced quantum solvers on H2 molecule")
print("Success criteria: Energy < HF baseline, reasonable convergence\n")

# Create H2 molecule
h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

# Build Hamiltonian via CovalentBond
bond = CovalentBond(h1, h2)
hamiltonian = bond.hamiltonian

# HF baseline
print("=" * 80)
print("  HARTREE-FOCK BASELINE")
print("=" * 80)

_, hf_energy = hamiltonian.solve_scf()
# Convert from numpy scalar to float
if hasattr(hf_energy, 'item'):
    hf_energy = hf_energy.item()
else:
    hf_energy = float(hf_energy)

hf_ev = hf_energy * 27.2114
print(f"HF Energy: {hf_energy:.6f} Ha ({hf_ev:.4f} eV)")

# Reference FCI energy for H2/STO-3G at R=0.74 Å
# (from NIST or high-quality calculations)
ref_fci = -1.3652  # Ha (approximate FCI for H2/STO-3G)
print(f"Ref FCI:   {ref_fci:.6f} Ha (approximate)")

# ============================================================================
#  TEST 1: QUANTUM PHASE ESTIMATION (QPE)
# ============================================================================
print("\n" + "=" * 80)
print("  TEST 1: QUANTUM PHASE ESTIMATION (QPE)")
print("=" * 80)
print("\nQPE provides exponential speedup for eigenvalue problems")
print("Expected: FCI-level accuracy in classical fallback mode\n")

start = time.time()
qpe = QPESolver(hamiltonian, n_ancilla=8)
qpe_result = qpe.solve()
qpe_time = time.time() - start

qpe_energy_ha = qpe_result['energy']
qpe_energy_ev = qpe_energy_ha * 27.2114

print(f"Results:")
print(f"  Energy:     {qpe_energy_ha:.6f} Ha ({qpe_energy_ev:.4f} eV)")
print(f"  Method:     {qpe_result['method']}")
print(f"  Phase:      {qpe_result['phase']:.6f}")
print(f"  Precision:  {qpe_result['precision']:.2e}")
print(f"  Converged:  {qpe_result['converged']}")
print(f"  Time:       {qpe_time:.3f}s")

# Validation: QPE should be below HF and close to FCI
qpe_improvement = hf_energy - qpe_energy_ha
qpe_vs_fci = abs(qpe_energy_ha - ref_fci)

print(f"\nAnalysis:")
print(f"  Improvement over HF: {qpe_improvement:.6f} Ha ({qpe_improvement*627.5:.2f} kcal/mol)")
print(f"  Difference from FCI: {qpe_vs_fci:.6f} Ha")

qpe_pass = qpe_improvement > 0.1  # Should be significantly better than HF

if qpe_pass:
    print("  Status: ✅ PASS - QPE energy below HF baseline")
else:
    print(f"  Status: ❌ FAIL - QPE not improving over HF (improvement = {qpe_improvement:.6f} Ha)")

# ============================================================================
#  TEST 2: SAMPLE-BASED QUANTUM DIAGONALIZATION (SQD)
# ============================================================================
print("\n" + "=" * 80)
print("  TEST 2: SAMPLE-BASED QUANTUM DIAGONALIZATION (SQD)")
print("=" * 80)
print("\nSQD uses quantum sampling + classical processing")
print("Expected: Improved energy by including excited configurations\n")

start = time.time()
sqd = SQDSolver(hamiltonian, n_samples=100, max_iterations=5)
sqd_result = sqd.solve()
sqd_time = time.time() - start

sqd_energy_ha = sqd_result['energy']
sqd_energy_ev = sqd_energy_ha * 27.2114

print(f"Results:")
print(f"  Energy:      {sqd_energy_ha:.6f} Ha ({sqd_energy_ev:.4f} eV)")
print(f"  Method:      {sqd_result['method']}")
print(f"  N Samples:   {sqd_result['n_samples']}")
print(f"  Converged:   {sqd_result['converged']}")
print(f"  Time:        {sqd_time:.3f}s")

if 'occupancies' in sqd_result:
    occ_alpha, occ_beta = sqd_result['occupancies']
    print(f"  Occ (α):     {occ_alpha}")
    print(f"  Occ (β):     {occ_beta}")

if 'spin_sq' in sqd_result:
    print(f"  ⟨S²⟩:        {sqd_result['spin_sq']:.6f}")

# Validation: SQD should be below HF
sqd_improvement = hf_energy - sqd_energy_ha
sqd_vs_fci = abs(sqd_energy_ha - ref_fci)

print(f"\nAnalysis:")
print(f"  Improvement over HF: {sqd_improvement:.6f} Ha ({sqd_improvement*627.5:.2f} kcal/mol)")
print(f"  Difference from FCI: {sqd_vs_fci:.6f} Ha")

sqd_pass = sqd_improvement > 0.1  # Should be significantly better than HF

if sqd_pass:
    print("  Status: ✅ PASS - SQD energy below HF baseline")
else:
    print(f"  Status: ❌ FAIL - SQD not improving over HF (improvement = {sqd_improvement:.6f} Ha)")

# ============================================================================
#  SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("  SUMMARY")
print("=" * 80)

print(f"\nEnergy Comparison:")
print(f"  HF:  {hf_energy:.6f} Ha ({hf_ev:.4f} eV)")
print(f"  QPE: {qpe_energy_ha:.6f} Ha ({qpe_energy_ev:.4f} eV) - Δ = {qpe_improvement:.4f} Ha")
print(f"  SQD: {sqd_energy_ha:.6f} Ha ({sqd_energy_ev:.4f} eV) - Δ = {sqd_improvement:.4f} Ha")
print(f"  FCI: {ref_fci:.6f} Ha (reference)")

print(f"\nTest Results:")
tests_passed = 0
if qpe_pass:
    print("  QPE: ✅ PASS")
    tests_passed += 1
else:
    print("  QPE: ❌ FAIL")

if sqd_pass:
    print("  SQD: ✅ PASS")
    tests_passed += 1
else:
    print("  SQD: ❌ FAIL")

print(f"\nOVERALL: {tests_passed}/2 tests passed ({tests_passed*50}%)")

if tests_passed == 2:
    print("\n🎉 ALL TESTS PASSED - QPE and SQD solvers validated!")
    exit(0)
else:
    print(f"\n❌ {2-tests_passed} test(s) failed - review solver implementations")
    exit(1)
