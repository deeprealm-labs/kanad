"""
Test: Bravyi-Kitaev with SQD Solver
====================================
Verify that BK mapper works correctly with exact methods (SQD)
since BK is theoretically sound for non-variational methods.
"""

import numpy as np
from kanad import BondFactory
from kanad.solvers import VQESolver, SQDSolver

def test_bk_with_sqd():
    """Test BK mapper with SQD solver (exact method)."""
    print("\n" + "="*60)
    print("TEST: Bravyi-Kitaev with SQD Solver")
    print("="*60)

    # Create H2 bond
    bond = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
    print(f"✓ Created H2 bond")

    # Test 1: VQE with JW (baseline)
    print("\n1. VQE with Jordan-Wigner (baseline):")
    solver_vqe_jw = VQESolver(bond, ansatz_type='ucc', mapper_type='jordan_wigner', backend='statevector')
    result_vqe_jw = solver_vqe_jw.solve()
    print(f"   Energy: {result_vqe_jw['energy']:.8f} Ha")

    # Test 2: SQD with JW
    print("\n2. SQD with Jordan-Wigner:")
    try:
        solver_sqd_jw = SQDSolver(bond, mapper_type='jordan_wigner')
        result_sqd_jw = solver_sqd_jw.solve()
        print(f"   Ground state: {result_sqd_jw['ground_state_energy']:.8f} Ha")
        if 'excited_states' in result_sqd_jw:
            print(f"   Excited states: {len(result_sqd_jw['excited_states'])} found")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        result_sqd_jw = None

    # Test 3: SQD with BK (should work for exact method!)
    print("\n3. SQD with Bravyi-Kitaev:")
    try:
        solver_sqd_bk = SQDSolver(bond, mapper_type='bravyi_kitaev')
        result_sqd_bk = solver_sqd_bk.solve()
        print(f"   Ground state: {result_sqd_bk['ground_state_energy']:.8f} Ha")
        if 'excited_states' in result_sqd_bk:
            print(f"   Excited states: {len(result_sqd_bk['excited_states'])} found")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        result_sqd_bk = None

    # Compare results
    print("\n" + "="*60)
    print("COMPARISON")
    print("="*60)

    if result_sqd_jw and result_sqd_bk:
        diff = abs(result_sqd_jw['ground_state_energy'] - result_sqd_bk['ground_state_energy'])
        print(f"JW Energy:  {result_sqd_jw['ground_state_energy']:.8f} Ha")
        print(f"BK Energy:  {result_sqd_bk['ground_state_energy']:.8f} Ha")
        print(f"Difference: {diff*1000:.6f} mHa")

        if diff < 1e-6:
            print("\n✅ SUCCESS! BK works correctly with SQD solver!")
            print("   (As expected - BK is fine for exact methods)")
            return True
        else:
            print("\n❌ UNEXPECTED: BK and JW differ even with exact method!")
            return False
    else:
        print("\n⚠️  Could not complete comparison (one method failed)")
        return False


if __name__ == "__main__":
    test_bk_with_sqd()
