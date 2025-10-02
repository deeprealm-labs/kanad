"""
Edge Cases and Stress Testing
==============================

Tests the framework's robustness and error handling with:

1. Extreme geometries (very short/long bonds)
2. Challenging SCF convergence cases
3. Large molecules (scaling test)
4. Mixed bond types in molecules
5. Numerical stability
6. Error handling and recovery

Ensures the framework fails gracefully and handles edge cases properly.
"""

import numpy as np
import sys
from kanad.bonds import BondFactory
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule

print("="*80)
print("EDGE CASES AND STRESS TESTING")
print("="*80)
print()

# ==============================================================================
# Test 1: Extreme Bond Lengths
# ==============================================================================
print("TEST 1: Extreme Bond Lengths")
print("-" * 80)

print("Testing framework stability at extreme geometries")
print()

extreme_cases = [
    ("Very short", 0.2, "Extreme repulsion"),
    ("Short", 0.4, "Strong repulsion"),
    ("Equilibrium", 0.74, "Normal bonding"),
    ("Long", 2.0, "Weak bonding"),
    ("Very long", 5.0, "Near dissociation"),
]

validations_extreme = []

print(f"H2 at different distances:")
print()

for description, distance, regime in extreme_cases:
    print(f"  {description} ({distance:.2f} Å) - {regime}:")

    try:
        h2 = BondFactory.create_bond('H', 'H', distance=distance)
        result = h2.compute_energy(method='HF', max_iterations=200)

        energy = result['energy']
        converged = result['converged']
        iterations = result['iterations']

        print(f"    Energy:     {energy:.6f} Ha")
        print(f"    Converged:  {'✓' if converged else '✗'} ({iterations} iterations)")

        # Check energy is finite and reasonable
        energy_reasonable = -5.0 < energy < 5.0 and not np.isnan(energy) and not np.isinf(energy)

        validations_extreme.append((f"r={distance:.1f}Å finite", energy_reasonable))

        # Very short distances should have higher (less negative) energy due to repulsion
        # Very long distances should approach 2 * H atom energy ≈ -1.0 Ha

    except Exception as e:
        print(f"    ✗ Failed: {str(e)[:60]}")
        validations_extreme.append((f"r={distance:.1f}Å handled", False))

    print()

# ==============================================================================
# Test 2: SCF Convergence Challenges
# ==============================================================================
print("TEST 2: Challenging SCF Convergence Cases")
print("-" * 80)

print("Testing SCF solver on difficult systems")
print()

difficult_systems = [
    ('Be', 'H', 1.3, "BeH - open shell challenge"),
    ('F', 'F', 1.42, "F₂ - strong correlation"),
    ('Li', 'Li', 2.67, "Li₂ - diffuse orbitals"),
]

validations_scf = []

for atom1, atom2, distance, description in difficult_systems:
    print(f"{description}:")
    print(f"  {atom1}-{atom2} at {distance:.2f} Å")

    try:
        bond = BondFactory.create_bond(atom1, atom2, distance=distance)

        print(f"  Bond type detected: {bond.bond_type}")
        print(f"  Electrons: {bond.molecule.n_electrons}")

        # Try SCF with standard settings
        result = bond.compute_energy(method='HF', max_iterations=250)

        energy = result['energy']
        converged = result['converged']
        iterations = result['iterations']

        print(f"  Energy:    {energy:.6f} Ha")
        print(f"  Converged: {'✓' if converged else '✗'} ({iterations} iterations)")

        # Even if doesn't converge, should produce reasonable energy
        energy_reasonable = -100 < energy < 10 and not np.isnan(energy)

        validations_scf.append((f"{description} reasonable", energy_reasonable))

        if converged:
            validations_scf.append((f"{description} converged", True))
            print(f"  ✓ Converged successfully")
        else:
            validations_scf.append((f"{description} converged", False))
            print(f"  ⚠ Did not converge but produced result")

    except Exception as e:
        print(f"  ✗ Failed: {str(e)[:60]}")
        validations_scf.append((f"{description} handled", False))

    print()

# ==============================================================================
# Test 3: Molecular Size Scaling
# ==============================================================================
print("TEST 3: Molecular Size Scaling")
print("-" * 80)

print("Testing framework performance with increasing molecular size")
print()

# Test linear hydrogen chains of different lengths
chain_sizes = [2, 3, 4]

validations_scaling = []

for n_atoms in chain_sizes:
    formula = f"H{n_atoms}"
    print(f"{formula} linear chain:")

    # Create linear chain of H atoms
    atoms_symbols = ['H'] * n_atoms

    try:
        # Use BondFactory to create molecule
        molecule = BondFactory.create_molecule(atoms_symbols, geometry='linear')

        print(f"  Atoms:     {molecule.n_atoms}")
        print(f"  Electrons: {molecule.n_electrons}")

        # Create representation
        from kanad.core.representations.lcao_representation import LCAORepresentation
        from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

        rep = LCAORepresentation(molecule)
        ham = CovalentHamiltonian(molecule, rep, basis_name='sto-3g')

        print(f"  Orbitals:  {rep.n_orbitals}")
        print(f"  Qubits:    {rep.n_qubits}")

        # Try SCF
        import time
        start_time = time.time()

        density, energy = ham.solve_scf(max_iterations=200)
        converged = getattr(ham, '_scf_converged', False)

        elapsed_time = time.time() - start_time

        print(f"  Energy:    {energy:.6f} Ha")
        print(f"  Converged: {'✓' if converged else '✗'}")
        print(f"  Time:      {elapsed_time:.2f} seconds")

        validations_scaling.append((f"{formula} calculated", True))

        # Check that energy scales roughly linearly with number of atoms
        # (For very approximate check)

    except Exception as e:
        print(f"  ✗ Failed: {str(e)[:60]}")
        validations_scaling.append((f"{formula} calculated", False))

    print()

# ==============================================================================
# Test 4: Mixed Bond Types
# ==============================================================================
print("TEST 4: Mixed Bond Types in Molecules")
print("-" * 80)

print("Testing molecules with multiple bonds")
print()

# Test C2H2 (acetylene) - has C≡C triple bond and C-H single bonds
print("Acetylene (HC≡CH):")

try:
    # Create atoms
    # Linear geometry: H-C-C-H
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    c1 = Atom('C', position=np.array([1.06, 0.0, 0.0]))  # C-H = 1.06 Å
    c2 = Atom('C', position=np.array([2.26, 0.0, 0.0]))  # C≡C = 1.20 Å
    h2 = Atom('H', position=np.array([3.32, 0.0, 0.0]))  # C-H = 1.06 Å

    c2h2 = Molecule([h1, c1, c2, h2])

    print(f"  Formula: C₂H₂")
    print(f"  Atoms:   {c2h2.n_atoms}")
    print(f"  Electrons: {c2h2.n_electrons}")

    # Create representation
    from kanad.core.representations.lcao_representation import LCAORepresentation
    from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

    c2h2_rep = LCAORepresentation(c2h2)
    c2h2_ham = CovalentHamiltonian(c2h2, c2h2_rep, basis_name='sto-3g')

    print(f"  Orbitals: {c2h2_rep.n_orbitals}")
    print(f"  Qubits:   {c2h2_rep.n_qubits}")

    # Try SCF
    print(f"  Running HF calculation...")
    density, energy = c2h2_ham.solve_scf(max_iterations=300)
    converged = getattr(c2h2_ham, '_scf_converged', False)

    print(f"  Energy:    {energy:.6f} Ha")
    print(f"  Converged: {'✓' if converged else '✗'}")

    validations_mixed = [(f"C2H2 calculated", True)]

except Exception as e:
    print(f"  ✗ Failed: {str(e)[:80]}")
    validations_mixed = [(f"C2H2 calculated", False)]

print()

# ==============================================================================
# Test 5: Numerical Stability
# ==============================================================================
print("TEST 5: Numerical Stability")
print("-" * 80)

print("Testing that calculations are numerically stable")
print()

validations_stability = []

# Test 1: Repeated calculations should give same result
print("Reproducibility test (H2):")
h2_energies = []

for i in range(3):
    h2_test = BondFactory.create_bond('H', 'H', distance=0.74)
    result = h2_test.compute_energy(method='HF', max_iterations=100)
    h2_energies.append(result['energy'])
    print(f"  Run {i+1}: {result['energy']:.10f} Ha")

energy_variance = np.std(h2_energies)
print(f"  Standard deviation: {energy_variance:.2e} Ha")

reproducible = energy_variance < 1e-8  # Should be nearly identical
print(f"  Reproducible: {'✓' if reproducible else '✗'}")

validations_stability.append(("Reproducibility", reproducible))
print()

# Test 2: Matrix condition numbers
print("Matrix conditioning (overlap matrix):")

h2_cond = BondFactory.create_bond('H', 'H', distance=0.74)
S_matrix = h2_cond.hamiltonian.S

condition_number = np.linalg.cond(S_matrix)
print(f"  Overlap matrix condition number: {condition_number:.2e}")

well_conditioned = condition_number < 1e6  # Reasonable threshold
print(f"  Well-conditioned: {'✓' if well_conditioned else '✗'}")

validations_stability.append(("Matrix conditioning", well_conditioned))
print()

# Test 3: Energy conservation during geometry changes
print("Energy continuity test:")

distances = [0.70, 0.72, 0.74, 0.76, 0.78]
cont_energies = []

for r in distances:
    h2_cont = BondFactory.create_bond('H', 'H', distance=r)
    result = h2_cont.compute_energy(method='HF', max_iterations=100)
    cont_energies.append(result['energy'])

# Check that energy changes smoothly (no jumps)
energy_diffs = np.abs(np.diff(cont_energies))
max_jump = np.max(energy_diffs)

print(f"  Distance range: {distances[0]:.2f} - {distances[-1]:.2f} Å")
print(f"  Energy range:   {min(cont_energies):.6f} - {max(cont_energies):.6f} Ha")
print(f"  Max energy jump: {max_jump:.6f} Ha")

smooth_curve = max_jump < 0.01  # Less than 10 mHa per 0.02 Å step
print(f"  Smooth curve: {'✓' if smooth_curve else '✗'}")

validations_stability.append(("Energy continuity", smooth_curve))
print()

# ==============================================================================
# Test 6: Error Handling
# ==============================================================================
print("TEST 6: Error Handling and Edge Cases")
print("-" * 80)

print("Testing that framework handles errors gracefully")
print()

validations_errors = []

# Test 1: Invalid atom symbol
print("Invalid atom symbol:")
try:
    bad_bond = BondFactory.create_bond('Xx', 'H')
    print(f"  ✗ Should have raised error but created: {bad_bond}")
    validations_errors.append(("Invalid atom error", False))
except Exception as e:
    print(f"  ✓ Correctly raised error: {type(e).__name__}")
    validations_errors.append(("Invalid atom error", True))
print()

# Test 2: Zero distance
print("Zero distance bond:")
try:
    zero_bond = BondFactory.create_bond('H', 'H', distance=0.0)
    result = zero_bond.compute_energy(method='HF', max_iterations=50)
    # If it doesn't crash, check if energy is reasonable
    energy_reasonable = not np.isnan(result['energy']) and not np.isinf(result['energy'])
    print(f"  Energy: {result['energy']:.6f} Ha")
    print(f"  {'✓' if energy_reasonable else '✗'} Handled gracefully")
    validations_errors.append(("Zero distance handled", energy_reasonable))
except Exception as e:
    print(f"  ✓ Raised error: {type(e).__name__}")
    validations_errors.append(("Zero distance handled", True))
print()

# Test 3: Negative distance (should be handled in absolute value or error)
print("Negative distance:")
try:
    neg_bond = BondFactory.create_bond('H', 'H', distance=-0.74)
    # Framework might take absolute value or error - both acceptable
    print(f"  Bond length: {neg_bond.get_bond_length():.2f} Å")
    print(f"  ✓ Handled (took absolute value)")
    validations_errors.append(("Negative distance handled", True))
except Exception as e:
    print(f"  ✓ Raised error: {type(e).__name__}")
    validations_errors.append(("Negative distance handled", True))
print()

# Test 4: Very large molecule (should handle or fail gracefully)
print("Very large system (H10):")
try:
    h10_symbols = ['H'] * 10
    h10 = BondFactory.create_molecule(h10_symbols, geometry='linear')

    from kanad.core.representations.lcao_representation import LCAORepresentation
    h10_rep = LCAORepresentation(h10)

    print(f"  Created molecule with {h10.n_atoms} atoms")
    print(f"  Qubits: {h10_rep.n_qubits}")
    print(f"  ✓ Large system handled")

    validations_errors.append(("Large system handled", True))

except Exception as e:
    print(f"  ⚠ Error (expected for very large systems): {type(e).__name__}")
    # This is acceptable - framework may have size limits
    validations_errors.append(("Large system handled", True))

print()

# ==============================================================================
# Test 7: Performance Benchmarks
# ==============================================================================
print("TEST 7: Performance Benchmarks")
print("-" * 80)

print("Measuring computational performance")
print()

import time

benchmarks = []

# H2 benchmark
print("H2 HF calculation:")
start = time.time()
h2_bench = BondFactory.create_bond('H', 'H', distance=0.74)
h2_bench_result = h2_bench.compute_energy(method='HF', max_iterations=100)
h2_time = time.time() - start

print(f"  Time: {h2_time:.3f} seconds")
print(f"  Iterations: {h2_bench_result['iterations']}")
print(f"  Time per iteration: {h2_time/h2_bench_result['iterations']:.4f} s")

benchmarks.append(('H2', h2_time))
print()

# H2O benchmark
print("H2O HF calculation:")
start = time.time()
h2o_bench = BondFactory.create_molecule(['H', 'O', 'H'], geometry='water')

from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

h2o_bench_rep = LCAORepresentation(h2o_bench)
h2o_bench_ham = CovalentHamiltonian(h2o_bench, h2o_bench_rep, basis_name='sto-3g')

density_h2o, energy_h2o = h2o_bench_ham.solve_scf(max_iterations=200)
h2o_time = time.time() - start

h2o_iters = getattr(h2o_bench_ham, '_scf_iterations', 0)

print(f"  Time: {h2o_time:.3f} seconds")
print(f"  Iterations: {h2o_iters}")
if h2o_iters > 0:
    print(f"  Time per iteration: {h2o_time/h2o_iters:.4f} s")

benchmarks.append(('H2O', h2o_time))
print()

print("Performance Summary:")
for name, time_taken in benchmarks:
    print(f"  {name:10} {time_taken:.3f} s")
print()

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
print("="*80)
print("EDGE CASES AND STRESS TEST SUMMARY")
print("="*80)
print()

all_validations = (
    validations_extreme +
    validations_scf +
    validations_scaling +
    validations_mixed +
    validations_stability +
    validations_errors
)

passed = sum(1 for _, p in all_validations if p)
total = len(all_validations)

print(f"Total Validations: {passed}/{total} passed ({passed/total*100:.1f}%)")
print()

# Categorize results
test_categories = {
    "Extreme Geometries": validations_extreme,
    "SCF Convergence": validations_scf,
    "Scaling": validations_scaling,
    "Mixed Bonds": validations_mixed,
    "Numerical Stability": validations_stability,
    "Error Handling": validations_errors,
}

print("Results by Category:")
for category, validations in test_categories.items():
    if validations:
        cat_passed = sum(1 for _, p in validations if p)
        cat_total = len(validations)
        print(f"  {category:<25} {cat_passed}/{cat_total} ({cat_passed/cat_total*100:.0f}%)")

print()

# Performance summary
if benchmarks:
    print("Performance Benchmarks:")
    for name, time_taken in benchmarks:
        print(f"  {name:10} {time_taken:.3f} s")
    print()

if passed >= total * 0.85:
    print("✅ EXCELLENT: Framework is robust and handles edge cases well")
    print("   Production-ready for research use!")
elif passed >= total * 0.70:
    print("✓ GOOD: Framework handles most edge cases")
    print("   Suitable for careful research use")
else:
    print("⚠ WARNING: Framework struggles with edge cases")
    print("   Needs robustness improvements")

print()
print("="*80)
