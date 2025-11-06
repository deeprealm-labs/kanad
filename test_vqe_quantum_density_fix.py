#!/usr/bin/env python3
"""
Test VQE Quantum Density Extraction Fix

Validates that VQE extracts quantum 1-RDM from statevector (not HF).
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.solvers.vqe_solver import VQESolver

print("=" * 70)
print("VQE QUANTUM DENSITY EXTRACTION TEST")
print("=" * 70)

# Create H2 bond
h2_bond = BondFactory.create_bond('H', 'H', distance=0.74)

print(f"\nMolecule: H2 (bond length 0.74 Å)")
print(f"Basis: STO-3G")
print(f"Electrons: {h2_bond.molecule.n_electrons}")

# Create VQE solver with analysis enabled
solver = VQESolver(
    bond=h2_bond,
    backend='statevector',
    enable_analysis=True,  # CRITICAL: Must enable to compute quantum density
    optimizer='COBYLA',
    max_iterations=50
)

print(f"\nVQE Configuration:")
print(f"  Backend: statevector")
print(f"  Analysis enabled: True")
print(f"  Optimizer: COBYLA")

# Solve
print("\n" + "-" * 70)
print("Running VQE solver...")
print("-" * 70)
result = solver.solve()

# Check results
print("\n" + "=" * 70)
print("RESULTS")
print("=" * 70)

print(f"\nGround state energy: {result['energy']:.8f} Hartree")
if 'hf_energy' in result:
    print(f"HF reference energy: {result['hf_energy']:.8f} Hartree")
    print(f"Correlation energy:  {result['correlation_energy']:.8f} Hartree")

# CHECK: Is quantum_rdm1 present?
print("\n" + "-" * 70)
print("CHECK: Quantum RDM1 Computed from VQE Statevector")
print("-" * 70)

if 'quantum_rdm1' in result:
    print("✅ PASS: quantum_rdm1 found in VQE results")
    quantum_rdm = np.array(result['quantum_rdm1'])
    print(f"   Shape: {quantum_rdm.shape}")
    print(f"   Trace: {np.trace(quantum_rdm):.4f}")

    # Check properties
    is_hermitian = np.allclose(quantum_rdm, quantum_rdm.T, atol=1e-8)
    n_electrons = h2_bond.molecule.n_electrons
    trace = np.trace(quantum_rdm)
    trace_correct = np.abs(trace - n_electrons) < 0.01

    print(f"\n   Hermitian: {is_hermitian} {'✅' if is_hermitian else '❌'}")
    print(f"   Trace: {trace:.4f} (expected {n_electrons}) {'✅' if trace_correct else '❌'}")

    # Get HF density for comparison
    try:
        hf_density, _ = h2_bond.hamiltonian.solve_scf()
        difference = np.max(np.abs(quantum_rdm - hf_density))
        print(f"   Max |Quantum - HF| difference: {difference:.6e}")

        if difference > 1e-6:
            print(f"   ✅ Quantum density differs significantly from HF (includes correlation)")
        else:
            print(f"   ⚠️  Quantum density close to HF")
    except:
        print("   (Could not compare with HF)")

    print("\n   VQE Quantum Density Matrix:")
    print(f"   {quantum_rdm}")

else:
    print("❌ FAIL: quantum_rdm1 NOT found in VQE results")
    print("   VQE quantum density extraction did not work!")

print("\n" + "=" * 70)
print("VQE QUANTUM DENSITY TEST COMPLETE")
print("=" * 70)
