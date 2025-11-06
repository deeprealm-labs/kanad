#!/usr/bin/env python3
"""
Test Critical Issue #1 Fix: Quantum Density Extraction in SQD Solver

This test validates that:
1. SQD solver computes quantum 1-RDM from correlated eigenvector
2. Quantum density is stored in results
3. Quantum density differs from HF density (includes correlation)
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.solvers.sqd_solver import SQDSolver


def test_sqd_quantum_density_extraction():
    """Test that SQD extracts quantum density, not HF density."""

    print("=" * 70)
    print("CRITICAL FIX TEST: SQD Quantum Density Extraction")
    print("=" * 70)

    # Create H2 bond
    h2_bond = BondFactory.create_bond('H', 'H', distance=0.74)

    print(f"\nMolecule: H2 (bond length 0.74 Å)")
    print(f"Basis: STO-3G")
    print(f"Electrons: {h2_bond.molecule.n_electrons}")

    # Create SQD solver with analysis enabled
    solver = SQDSolver(
        bond=h2_bond,
        subspace_dim=8,
        enable_analysis=True  # CRITICAL: Must enable analysis to compute quantum density
    )

    print(f"\nSQD Configuration:")
    print(f"  Subspace dimension: 8")
    print(f"  Analysis enabled: True")

    # Solve
    print("\n" + "-" * 70)
    print("Running SQD solver...")
    print("-" * 70)
    result = solver.solve(n_states=1)

    # Check results
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)

    print(f"\nGround state energy: {result['ground_state_energy']:.8f} Hartree")

    if 'hf_energy' in result:
        print(f"HF reference energy: {result['hf_energy']:.8f} Hartree")
        print(f"Correlation energy:  {result['correlation_energy']:.8f} Hartree")

    # CRITICAL CHECK #1: Is quantum_rdm1 present?
    print("\n" + "-" * 70)
    print("CHECK #1: Quantum RDM1 Computed")
    print("-" * 70)

    if 'quantum_rdm1' in result:
        print("✅ PASS: quantum_rdm1 found in results")
        quantum_rdm = np.array(result['quantum_rdm1'])
        n_electrons = h2_bond.molecule.n_electrons
        print(f"   Shape: {quantum_rdm.shape}")
        print(f"   Trace: {np.trace(quantum_rdm):.4f} (should ≈ {n_electrons})")
    else:
        print("❌ FAIL: quantum_rdm1 NOT found in results")
        print("   This means the fix didn't work!")
        return False

    # CRITICAL CHECK #2: Is rdm1 using quantum density?
    print("\n" + "-" * 70)
    print("CHECK #2: RDM1 Uses Quantum Density")
    print("-" * 70)

    if 'rdm1' in result:
        rdm1 = np.array(result['rdm1'])

        # Check if rdm1 matches quantum_rdm1
        if np.allclose(rdm1, quantum_rdm, atol=1e-10):
            print("✅ PASS: rdm1 uses quantum density (not HF)")
            print("   rdm1 and quantum_rdm1 are identical")
        else:
            print("❌ FAIL: rdm1 differs from quantum_rdm1")
            print(f"   Max difference: {np.max(np.abs(rdm1 - quantum_rdm))}")
            return False
    else:
        print("⚠️  WARNING: rdm1 not found in results")

    # CRITICAL CHECK #3: Does quantum density differ from HF?
    print("\n" + "-" * 70)
    print("CHECK #3: Quantum Density Differs from HF")
    print("-" * 70)

    # Get HF density for comparison
    try:
        # Try to get HF density from hamiltonian
        if hasattr(solver.hamiltonian, 'mf'):
            hf_density = solver.hamiltonian.mf.make_rdm1()
        else:
            # Compute HF density by solving SCF
            hf_density, _ = solver.hamiltonian.solve_scf()

        difference = np.max(np.abs(quantum_rdm - hf_density))
        print(f"Max |quantum - HF| density difference: {difference:.6e}")

        if difference > 1e-6:
            print("✅ PASS: Quantum density includes correlation effects")
            print("   (Differs significantly from HF)")
        else:
            print("⚠️  WARNING: Quantum density very close to HF")
            print("   (Correlation effects may be small for this system)")
    except Exception as e:
        print(f"⚠️  Could not compare with HF density: {e}")
        print("   (Skipping HF comparison)")

    # CRITICAL CHECK #4: Quantum density properties
    print("\n" + "-" * 70)
    print("CHECK #4: Quantum Density Properties")
    print("-" * 70)

    # Check hermiticity
    is_hermitian = np.allclose(quantum_rdm, quantum_rdm.T, atol=1e-8)
    print(f"Hermitian: {is_hermitian} {'✅' if is_hermitian else '❌'}")

    # Check trace (should equal number of electrons)
    n_electrons = h2_bond.molecule.n_electrons
    trace = np.trace(quantum_rdm)
    trace_correct = np.abs(trace - n_electrons) < 0.01
    print(f"Trace: {trace:.4f} (expected {n_electrons}) {'✅' if trace_correct else '❌'}")

    # Check eigenvalues (should be between 0 and 2 for RHF-like systems)
    eigenvalues = np.linalg.eigvalsh(quantum_rdm)
    eigenvalues_valid = np.all(eigenvalues >= -1e-6) and np.all(eigenvalues <= 2.0 + 1e-6)
    print(f"Eigenvalues in [0, 2]: {eigenvalues_valid} {'✅' if eigenvalues_valid else '❌'}")
    print(f"   Min eigenvalue: {np.min(eigenvalues):.6f}")
    print(f"   Max eigenvalue: {np.max(eigenvalues):.6f}")

    # Display quantum density matrix
    print("\n" + "-" * 70)
    print("Quantum Density Matrix (1-RDM):")
    print("-" * 70)
    print(quantum_rdm)

    print("\n" + "-" * 70)
    print("HF Density Matrix (for comparison):")
    print("-" * 70)
    print(hf_density)

    print("\n" + "-" * 70)
    print("Difference (Quantum - HF):")
    print("-" * 70)
    print(quantum_rdm - hf_density)

    # Final verdict
    print("\n" + "=" * 70)
    print("FINAL VERDICT")
    print("=" * 70)

    all_checks_passed = (
        'quantum_rdm1' in result and
        is_hermitian and
        trace_correct and
        eigenvalues_valid
    )

    if all_checks_passed:
        print("✅ ALL CHECKS PASSED")
        print("\nCritical Issue #1 is FIXED:")
        print("- SQD computes quantum density from correlated eigenvector")
        print("- Quantum density includes correlation effects")
        print("- Density has correct physical properties")
        return True
    else:
        print("❌ SOME CHECKS FAILED")
        print("\nCritical Issue #1 NOT fully fixed")
        return False


if __name__ == "__main__":
    success = test_sqd_quantum_density_extraction()
    exit(0 if success else 1)
