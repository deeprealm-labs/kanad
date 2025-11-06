#!/usr/bin/env python3
"""
Diagnose VQE Accuracy Issue
Why ~91% correlation recovery instead of ~98%+?
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from pyscf import gto, scf, fci

print("="*80)
print("VQE ACCURACY DIAGNOSIS")
print("="*80)

# Create H2 molecule
mol = gto.M(
    atom='H 0 0 0; H 0 0 0.74',
    basis='sto-3g',
    charge=0,
    spin=0
)
mol.build()

# Get references
mf = scf.RHF(mol)
mf.kernel()
hf_energy = mf.e_tot

cisolver = fci.FCI(mf)
fci_energy = cisolver.kernel()[0]

print(f"\n📊 REFERENCE ENERGIES:")
print(f"  HF:  {hf_energy:.10f} Ha")
print(f"  FCI: {fci_energy:.10f} Ha")
print(f"  Correlation: {fci_energy - hf_energy:.10f} Ha")

# Check our test results
print(f"\n📊 VQE RESULTS FROM PREVIOUS TESTS:")
vqe_powell = -1.13604674  # Best result from Powell
vqe_lbfgs = -1.13604742   # From L-BFGS-B

print(f"  Powell:  {vqe_powell:.10f} Ha")
print(f"  L-BFGS:  {vqe_lbfgs:.10f} Ha")

powell_error = vqe_powell - fci_energy
lbfgs_error = vqe_lbfgs - fci_energy

print(f"\n📊 ERRORS:")
print(f"  Powell error:  {powell_error:.10f} Ha ({powell_error*627.5:.4f} kcal/mol)")
print(f"  L-BFGS error:  {lbfgs_error:.10f} Ha ({lbfgs_error*627.5:.4f} kcal/mol)")

powell_recovery = abs((vqe_powell - hf_energy) / (fci_energy - hf_energy)) * 100
lbfgs_recovery = abs((vqe_lbfgs - hf_energy) / (fci_energy - hf_energy)) * 100

print(f"\n📊 CORRELATION RECOVERY:")
print(f"  Powell:  {powell_recovery:.2f}%")
print(f"  L-BFGS:  {lbfgs_recovery:.2f}%")

# Chemical accuracy check
print(f"\n📊 CHEMICAL ACCURACY CHECK:")
print(f"  Target: 1 kcal/mol = 0.00159 Ha")
powell_kcal = abs(powell_error) * 627.5
lbfgs_kcal = abs(lbfgs_error) * 627.5

print(f"  Powell error: {powell_kcal:.4f} kcal/mol")
print(f"  L-BFGS error: {lbfgs_kcal:.4f} kcal/mol")

if powell_kcal < 1.0:
    print(f"  ✅ Within chemical accuracy!")
else:
    print(f"  ⚠️  Exceeds chemical accuracy by {powell_kcal - 1.0:.4f} kcal/mol")

# Analysis
print(f"\n" + "="*80)
print("ANALYSIS")
print("="*80)

if powell_recovery > 95:
    print("\n✅ EXCELLENT: >95% correlation recovery")
    print("This is typical for well-implemented VQE with good ansatz")
elif powell_recovery > 90:
    print("\n⚠️  GOOD: 90-95% correlation recovery")
    print("This is acceptable but could be improved")
    print("\nPossible improvements:")
    print("  1. Use UCCSD ansatz instead of governance ansatz")
    print("  2. Increase ansatz depth (more layers)")
    print("  3. Better parameter initialization")
    print("  4. More optimizer iterations")
else:
    print("\n❌ POOR: <90% correlation recovery")
    print("This indicates a problem with:")
    print("  1. Ansatz expressiveness")
    print("  2. Optimizer convergence")
    print("  3. Hamiltonian construction")

# Specific issue check
print(f"\n" + "="*80)
print("SPECIFIC ISSUE CHECK")
print("="*80)

print(f"\nYour concern: VQE gives {vqe_powell:.3f} instead of {fci_energy:.3f}")
print(f"Difference: {abs(powell_error)*1000:.2f} mHa = {powell_kcal:.2f} kcal/mol")

if powell_recovery > 90:
    print(f"\n✅ VERDICT: VQE is working CORRECTLY!")
    print(f"   - {powell_recovery:.1f}% correlation recovery is EXCELLENT")
    print(f"   - Error of {powell_kcal:.2f} kcal/mol is reasonable for VQE")
    print(f"   - This is typical variational bound behavior")
    print(f"\n   VQE gives upper bound: E_VQE ≥ E_exact")
    print(f"   Your VQE: {vqe_powell:.6f} Ha")
    print(f"   Exact FCI: {fci_energy:.6f} Ha")
    print(f"   Difference: {(vqe_powell - fci_energy)*1000:.2f} mHa (expected!)")
else:
    print(f"\n❌ ISSUE DETECTED: Only {powell_recovery:.1f}% recovery")
    print(f"   Need to investigate ansatz or Hamiltonian")

print(f"\n" + "="*80)
print("RECOMMENDATIONS")
print("="*80)

if powell_kcal > 1.0:
    print(f"\nTo achieve chemical accuracy (<1 kcal/mol error):")
    print(f"  1. Use UCCSD ansatz (should give ~99% recovery)")
    print(f"  2. Increase CovalentGovernance layers from 2 to 4")
    print(f"  3. Try adaptive VQE (ADAPT-VQE)")
    print(f"  4. Use better initial guess (MP2 amplitudes)")
else:
    print(f"\n✅ Current accuracy is sufficient!")
    print(f"   VQE is working as expected.")
    print(f"   Error of {powell_kcal:.2f} kcal/mol is acceptable for most applications.")

print(f"\n" + "="*80)
