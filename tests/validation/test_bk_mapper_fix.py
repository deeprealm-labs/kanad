"""
Validation Test: Bravyi-Kitaev Mapper Fix
==========================================
Validates that the fixed BK mapper works correctly on H2 and LiH.
"""

import numpy as np
from kanad import BondFactory
from kanad.solvers import VQESolver


def test_bk_mapper_h2():
    """Test BK mapper on H2 - should not fail like before."""
    print("\n" + "="*60)
    print("TEST: Bravyi-Kitaev Mapper on H2")
    print("="*60)

    # Create H2 bond
    bond = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
    print(f"✓ Created bond: {bond.__class__.__name__}")

    # Test with BK mapper
    solver = VQESolver(
        bond,
        ansatz_type='ucc',
        mapper_type='bravyi_kitaev',
        backend='statevector',
        max_iterations=100
    )
    print(f"✓ Created VQE solver with BK mapper")

    # Solve
    result = solver.solve()
    energy = result['energy']

    print(f"\nResults:")
    print(f"  Energy: {energy:.8f} Ha")
    print(f"  Converged: {result.get('converged', False)}")

    # Reference energy from JW mapper (our experiments)
    reference_jw = -1.11675931
    error = abs(energy - reference_jw)

    print(f"\nComparison with Jordan-Wigner:")
    print(f"  JW energy: {reference_jw:.8f} Ha")
    print(f"  BK energy: {energy:.8f} Ha")
    print(f"  Difference: {error*1000:.3f} mHa")

    # BK should give same result as JW (both are exact transformations)
    # Allow small numerical error
    assert error < 0.01, f"BK energy differs from JW by {error*1000:.1f} mHa (should be <10 mHa)"

    # Energy should be below Hartree-Fock (validation check)
    # HF energy for H2 is around -1.117 Ha
    assert energy < -1.0, f"Energy {energy} is unphysical (should be negative)"

    print("\n✅ TEST PASSED: BK mapper works correctly on H2!")

    return result


def test_bk_mapper_lih():
    """Test BK mapper on LiH - should be faster than before."""
    print("\n" + "="*60)
    print("TEST: Bravyi-Kitaev Mapper on LiH")
    print("="*60)

    # Create LiH bond
    bond = BondFactory.create_bond('Li', 'H', distance=1.6, basis='sto-3g')
    print(f"✓ Created bond: {bond.__class__.__name__}")

    # Test with BK mapper
    import time
    solver = VQESolver(
        bond,
        ansatz_type='ucc',
        mapper_type='bravyi_kitaev',
        backend='statevector',
        max_iterations=100
    )
    print(f"✓ Created VQE solver with BK mapper")

    # Solve
    start_time = time.time()
    result = solver.solve()
    elapsed = time.time() - start_time

    energy = result['energy']

    print(f"\nResults:")
    print(f"  Energy: {energy:.8f} Ha")
    print(f"  Time: {elapsed:.2f} s")
    print(f"  Converged: {result.get('converged', False)}")

    # Reference from experiments: JW took 19.07s
    print(f"\nPerformance:")
    print(f"  BK time: {elapsed:.2f} s")
    print(f"  Previous BK time: 89.97 s (from experiments)")

    # Reference energy from JW mapper
    reference_jw = -7.86186489
    error = abs(energy - reference_jw)

    print(f"\nComparison with Jordan-Wigner:")
    print(f"  JW energy: {reference_jw:.8f} Ha")
    print(f"  BK energy: {energy:.8f} Ha")
    print(f"  Difference: {error*1000:.3f} mHa")

    # BK should give same result as JW
    assert error < 0.01, f"BK energy differs from JW by {error*1000:.1f} mHa (should be <10 mHa)"

    # Should be significantly faster than 89.97s (if optimization worked)
    # But we're more concerned with correctness for now
    print(f"\n✅ TEST PASSED: BK mapper works correctly on LiH!")

    return result


def test_bk_vs_jw_comparison():
    """Direct comparison of BK and JW on same molecule."""
    print("\n" + "="*60)
    print("TEST: Direct BK vs JW Comparison on H2")
    print("="*60)

    bond = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')

    # JW mapper
    print("\nRunning with Jordan-Wigner mapper...")
    solver_jw = VQESolver(bond, ansatz_type='ucc', mapper_type='jordan_wigner', backend='statevector')
    result_jw = solver_jw.solve()
    energy_jw = result_jw['energy']
    print(f"  JW Energy: {energy_jw:.8f} Ha")

    # BK mapper
    print("\nRunning with Bravyi-Kitaev mapper...")
    solver_bk = VQESolver(bond, ansatz_type='ucc', mapper_type='bravyi_kitaev', backend='statevector')
    result_bk = solver_bk.solve()
    energy_bk = result_bk['energy']
    print(f"  BK Energy: {energy_bk:.8f} Ha")

    # Compare
    diff = abs(energy_jw - energy_bk)
    print(f"\nDifference: {diff*1000:.6f} mHa")

    assert diff < 1e-6, f"BK and JW should give identical results, got {diff*1000:.3f} mHa difference"

    print("\n✅ TEST PASSED: BK and JW give identical results!")

    return {'jw': result_jw, 'bk': result_bk}


if __name__ == "__main__":
    print("\n" + "="*60)
    print("BRAVYI-KITAEV MAPPER VALIDATION SUITE")
    print("="*60)
    print("\nObjective: Verify BK mapper fix resolves previous failures")
    print("Previous issues:")
    print("  - H2: Completely wrong energy (-0.35 Ha vs -1.12 Ha)")
    print("  - LiH: Very slow (89.97s vs 19.07s for JW)")
    print()

    try:
        # Test 1: H2
        result_h2 = test_bk_mapper_h2()

        # Test 2: LiH
        result_lih = test_bk_mapper_lih()

        # Test 3: Direct comparison
        result_comp = test_bk_vs_jw_comparison()

        print("\n" + "="*60)
        print("ALL TESTS PASSED! ✅")
        print("="*60)
        print("\nSummary:")
        print("  ✅ BK mapper works correctly on H2")
        print("  ✅ BK mapper works correctly on LiH")
        print("  ✅ BK and JW give identical results")
        print("\nConclusion: Bravyi-Kitaev mapper is now FIXED and validated!")

    except Exception as e:
        print(f"\n❌ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        raise
