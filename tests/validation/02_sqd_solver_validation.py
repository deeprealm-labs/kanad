"""
SQD (Subspace Quantum Diagonalization) Solver Validation.

Tests SQD with:
- Ground state energy
- Excited states energies
- Different subspace dimensions
- Comparison with exact diagonalization
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.solvers.sqd_solver import SQDSolver

print("=" * 80)
print("SQD SOLVER VALIDATION")
print("=" * 80)

# Reference energies for H2 at 0.74 Å with sto-3g
REFERENCE_H2 = {
    'ground_state': -1.137284,
    'first_excited': -0.475,  # Approximate
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
print("TEST 1: H2 Ground State - SQD with subspace_dim=10")
print("=" * 80)

h1 = Atom('H', position=(0.0, 0.0, 0.0))
h2 = Atom('H', position=(0.74, 0.0, 0.0))
bond = CovalentBond(h1, h2, basis='sto-3g')

try:
    solver = SQDSolver(
        bond=bond,
        subspace_dim=10,
        circuit_depth=3,
        backend='statevector',
        random_seed=42  # For reproducible results
    )
    result = solver.solve()

    print(f"\nSQD Result:")
    print(f"  Energies: {result.get('energies', [])[:5]}")

    if 'energies' in result and len(result['energies']) > 0:
        ground_energy = result['energies'][0]
        validate_result(
            "H2 Ground State (SQD, dim=10)",
            ground_energy,
            REFERENCE_H2['ground_state']
        )
    else:
        print("✗ No energies returned")
        results.append({'name': 'H2 Ground State (SQD, dim=10)', 'passed': False})

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
    results.append({'name': 'H2 Ground State (SQD, dim=10)', 'passed': False, 'error': str(e)})


print("\n" + "=" * 80)
print("TEST 2: H2 Ground State - SQD with subspace_dim=20")
print("=" * 80)

# Now uses physical basis (HF + singles + doubles) with random seed
try:
    solver2 = SQDSolver(
        bond=bond,
        subspace_dim=20,  # Larger subspace
        circuit_depth=3,
        backend='statevector',
        random_seed=42  # Deterministic
    )
    result2 = solver2.solve()

    if 'energies' in result2 and len(result2['energies']) > 0:
        ground_energy = result2['energies'][0]
        validate_result(
            "H2 Ground State (SQD, dim=20)",
            ground_energy,
            REFERENCE_H2['ground_state'],
            tolerance_mHa=1.0  # Even tighter tolerance with physical basis
        )
    else:
        print("✗ No energies returned")
        results.append({'name': 'H2 Ground State (SQD, dim=20)', 'passed': False})

except Exception as e:
    print(f"✗ FAILED: {e}")
    results.append({'name': 'H2 Ground State (SQD, dim=20)', 'passed': False, 'error': str(e)})


print("\n" + "=" * 80)
print("TEST 3: H2 Excited States - Compare multiple eigenvalues")
print("=" * 80)

try:
    solver = SQDSolver(
        bond=bond,
        subspace_dim=15,
        circuit_depth=3,
        backend='statevector',
        random_seed=42  # For reproducible results
    )
    result = solver.solve()

    if 'energies' in result and len(result['energies']) >= 2:
        energies = result['energies']
        print(f"\nFound {len(energies)} eigenvalues:")
        for i, E in enumerate(energies[:5]):
            print(f"  State {i}: {E:.6f} Ha")

        # Validate ground state
        validate_result(
            "H2 Ground State (multiple eigenvalues)",
            energies[0],
            REFERENCE_H2['ground_state']
        )

        # Check first excited state is higher
        if energies[1] > energies[0]:
            print(f"\n✓ First excited state ({energies[1]:.6f} Ha) > ground state")
            print(f"  Excitation energy: {(energies[1] - energies[0])*27.211:.3f} eV")
        else:
            print(f"\n✗ First excited state ({energies[1]:.6f} Ha) ≤ ground state")

    else:
        print("✗ Not enough eigenvalues returned")
        results.append({'name': 'H2 Excited States', 'passed': False})

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
    results.append({'name': 'H2 Excited States', 'passed': False, 'error': str(e)})


print("\n" + "=" * 80)
print("TEST 4: Compare SQD vs Exact Diagonalization")
print("=" * 80)

try:
    # Get exact energies (NumPy eigensolver)
    from kanad.solvers.vqe_solver import VQESolver
    from qiskit.quantum_info import SparsePauliOp
    from scipy.sparse.linalg import eigsh

    # Build Hamiltonian using representation's to_qubit_operator
    hamiltonian_pauli = bond.hamiltonian.representation.to_qubit_operator()

    # Convert to matrix and diagonalize
    if isinstance(hamiltonian_pauli, SparsePauliOp):
        matrix = hamiltonian_pauli.to_matrix()
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
        exact_ground = eigenvalues[0]

        print(f"\nExact diagonalization:")
        print(f"  Ground state: {exact_ground:.6f} Ha")
        print(f"  First excited: {eigenvalues[1]:.6f} Ha")

        # Compare with SQD
        solver = SQDSolver(bond=bond, subspace_dim=10, circuit_depth=3)
        result = solver.solve()

        if 'energies' in result and len(result['energies']) > 0:
            sqd_ground = result['energies'][0]
            error = abs(sqd_ground - exact_ground) * 1000

            print(f"\nSQD vs Exact:")
            print(f"  SQD:   {sqd_ground:.6f} Ha")
            print(f"  Exact: {exact_ground:.6f} Ha")
            print(f"  Error: {error:.3f} mHa")

            if error < 50:
                print(f"\n✓ SQD matches exact diagonalization (< 50 mHa)")
                results.append({'name': 'SQD vs Exact', 'passed': True, 'error_mHa': error})
            else:
                print(f"\n✗ SQD error too large (> 50 mHa)")
                results.append({'name': 'SQD vs Exact', 'passed': False, 'error_mHa': error})

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
    results.append({'name': 'SQD vs Exact', 'passed': False, 'error': str(e)})


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
print(f"SQD VALIDATION {'PASSED' if passed >= total * 0.5 else 'FAILED'}")
print("=" * 80)
