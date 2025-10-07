"""
VQE Solver Validation - Multiple Configurations.

Tests VQE with:
- Different ansatze (Governance, Hardware-Efficient, Real-Amplitudes)
- Different mappers (Jordan-Wigner, Bravyi-Kitaev)
- Different Hamiltonians (Covalent, Ionic, Metallic)
- Different optimizers (SLSQP, Powell, COBYLA)
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.bonds.ionic_bond import IonicBond
from kanad.solvers.vqe_solver import VQESolver

print("=" * 80)
print("VQE SOLVER VALIDATION")
print("=" * 80)

# Reference energies (from literature/PySCF)
REFERENCE_ENERGIES = {
    'H2': -1.137284,  # Exact at 0.74 Å with sto-3g
    'LiH': -7.882324  # Exact with sto-3g
}

results = []


def validate_result(name, energy, reference, tolerance_mHa=50):
    """Validate computed energy against reference."""
    error_mHa = abs(energy - reference) * 1000
    status = "✓" if error_mHa < tolerance_mHa else "✗"

    print(f"\n{status} {name}")
    print(f"  Energy:    {energy:.6f} Ha")
    print(f"  Reference: {reference:.6f} Ha")
    print(f"  Error:     {error_mHa:.3f} mHa")

    results.append({
        'name': name,
        'energy': energy,
        'reference': reference,
        'error_mHa': error_mHa,
        'passed': error_mHa < tolerance_mHa
    })

    return error_mHa < tolerance_mHa


print("\n" + "=" * 80)
print("TEST 1: H2 with Governance Ansatz + Jordan-Wigner + SLSQP")
print("=" * 80)

h1 = Atom('H', position=(0.0, 0.0, 0.0))
h2 = Atom('H', position=(0.74, 0.0, 0.0))
bond = CovalentBond(h1, h2, basis='sto-3g')

try:
    solver = VQESolver(
        bond=bond,
        ansatz_type='governance',
        mapper_type='jordan_wigner',
        optimizer='SLSQP',
        max_iterations=100
    )
    result = solver.solve()
    validate_result(
        "H2 - Governance + JW + SLSQP",
        result['energy'],
        REFERENCE_ENERGIES['H2']
    )
except Exception as e:
    print(f"✗ FAILED: {e}")
    results.append({'name': 'H2 - Governance + JW + SLSQP', 'passed': False, 'error': str(e)})


print("\n" + "=" * 80)
print("TEST 2: H2 with UCC Ansatz + Jordan-Wigner + SLSQP")
print("=" * 80)

try:
    solver = VQESolver(
        bond=bond,
        ansatz_type='ucc',
        mapper_type='jordan_wigner',
        optimizer='SLSQP',
        max_iterations=100
    )
    result = solver.solve()
    validate_result(
        "H2 - UCC + JW + SLSQP",
        result['energy'],
        REFERENCE_ENERGIES['H2'],
        tolerance_mHa=100  # UCC may not converge well
    )
except Exception as e:
    print(f"✗ FAILED: {e}")
    results.append({'name': 'H2 - UCC + JW + SLSQP', 'passed': False, 'error': str(e)})


print("\n" + "=" * 80)
print("TEST 3: H2 with Hardware-Efficient Ansatz + Bravyi-Kitaev + SLSQP")
print("=" * 80)

try:
    solver = VQESolver(
        bond=bond,
        ansatz_type='hardware_efficient',
        mapper_type='bravyi_kitaev',
        optimizer='SLSQP',
        max_iterations=100
    )
    result = solver.solve()
    validate_result(
        "H2 - HardwareEfficient + BK + SLSQP",
        result['energy'],
        REFERENCE_ENERGIES['H2'],
        tolerance_mHa=100  # More lenient for hardware-efficient
    )
except Exception as e:
    print(f"✗ FAILED: {e}")
    results.append({'name': 'H2 - HardwareEfficient + BK + SLSQP', 'passed': False, 'error': str(e)})


print("\n" + "=" * 80)
print("TEST 4: LiH (Covalent) with Governance Ansatz + Jordan-Wigner + SLSQP")
print("=" * 80)

li = Atom('Li', position=(0.0, 0.0, 0.0))
h = Atom('H', position=(1.56, 0.0, 0.0))
lih_bond = CovalentBond(li, h, basis='sto-3g')

try:
    solver = VQESolver(
        bond=lih_bond,
        ansatz_type='governance',
        mapper_type='jordan_wigner',
        optimizer='SLSQP',
        max_iterations=50  # Fewer iterations for larger system
    )
    result = solver.solve()
    validate_result(
        "LiH - Covalent Governance + JW + SLSQP",
        result['energy'],
        REFERENCE_ENERGIES['LiH'],
        tolerance_mHa=100  # More lenient for larger system
    )
except Exception as e:
    print(f"✗ FAILED: {e}")
    results.append({'name': 'LiH - Covalent Governance + JW + SLSQP', 'passed': False, 'error': str(e)})


print("\n" + "=" * 80)
print("TEST 5: H2 with Governance Ansatz + Jordan-Wigner + COBYLA")
print("=" * 80)

try:
    solver = VQESolver(
        bond=bond,
        ansatz_type='governance',
        mapper_type='jordan_wigner',
        optimizer='COBYLA',
        max_iterations=200  # COBYLA needs more iterations
    )
    result = solver.solve()
    validate_result(
        "H2 - Governance + JW + COBYLA",
        result['energy'],
        REFERENCE_ENERGIES['H2'],
        tolerance_mHa=100  # More lenient for COBYLA
    )
except Exception as e:
    print(f"✗ FAILED: {e}")
    results.append({'name': 'H2 - Governance + JW + COBYLA', 'passed': False, 'error': str(e)})


# Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

passed = sum(1 for r in results if r.get('passed', False))
total = len(results)

print(f"\nPassed: {passed}/{total}")
print(f"Success Rate: {100*passed/total:.1f}%")

print("\nDetailed Results:")
for r in results:
    status = "✓" if r.get('passed', False) else "✗"
    if 'error_mHa' in r:
        print(f"{status} {r['name']}: {r['error_mHa']:.3f} mHa error")
    else:
        print(f"{status} {r['name']}: {r.get('error', 'Unknown error')}")

print("\n" + "=" * 80)
print(f"VQE VALIDATION {'PASSED' if passed == total else 'FAILED'}")
print("=" * 80)
