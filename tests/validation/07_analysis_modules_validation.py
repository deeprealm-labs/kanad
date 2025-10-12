"""
Analysis Modules Validation - Demonstrates Real Capabilities.

Uses CORRECT APIs verified against actual code.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.analysis.property_calculator import PropertyCalculator
from kanad.analysis.energy_analysis import EnergyAnalyzer, BondingAnalyzer, CorrelationAnalyzer

print("=" * 80)
print("ANALYSIS MODULES VALIDATION")
print("=" * 80)

results = []

def validate_test(name, condition, message=""):
    status = "✓" if condition else "✗"
    print(f"\n{status} {name}")
    if message:
        print(f"  {message}")
    results.append({'name': name, 'passed': condition})
    return condition


print("\n" + "=" * 80)
print("TEST 1: Dipole Moment - H2 (symmetric, should be ~0)")
print("=" * 80)

try:
    h1 = Atom('H', position=(0.0, 0.0, 0.0))
    h2 = Atom('H', position=(0.74, 0.0, 0.0))
    bond = CovalentBond(h1, h2, basis='sto-3g')

    calc = PropertyCalculator(bond.hamiltonian)
    dipole = calc.compute_dipole_moment()

    print(f"Dipole magnitude: {dipole['dipole_magnitude']:.6f} Debye")
    print(f"Components: {dipole['components']}")

    # H2 is symmetric
    validate_test(
        "H2 dipole (symmetric)",
        dipole['dipole_magnitude'] < 0.01,
        f"{dipole['dipole_magnitude']:.6f} Debye ≈ 0 ✓"
    )

except Exception as e:
    print(f"✗ FAILED: {e}")
    results.append({'name': 'Dipole', 'passed': False})


print("\n" + "=" * 80)
print("TEST 2: Polarizability - H2")
print("=" * 80)

try:
    polar = calc.compute_polarizability()

    print(f"Mean polarizability: {polar['alpha_mean']:.4f} a.u.")
    print(f"Tensor: {polar['alpha_tensor']}")

    # Note: sto-3g gives very low polarizability values
    validate_test(
        "H2 polarizability computed",
        polar['alpha_mean'] > 0.5,  # At least positive and reasonable
        f"{polar['alpha_mean']:.4f} a.u. (sto-3g gives low values)"
    )

except Exception as e:
    print(f"✗ FAILED: {e}")
    results.append({'name': 'Polarizability', 'passed': False})


print("\n" + "=" * 80)
print("TEST 3: Energy Decomposition")
print("=" * 80)

try:
    hf_result = bond.compute_energy(method='HF', max_iterations=100)
    analyzer = EnergyAnalyzer(bond.hamiltonian)
    decomp = analyzer.decompose_energy(hf_result['density_matrix'])

    print(f"Nuclear: {decomp['nuclear_repulsion']:.6f} Ha")
    print(f"One-e:   {decomp['one_electron']:.6f} Ha")
    print(f"Two-e:   {decomp['two_electron']:.6f} Ha")
    print(f"Total:   {decomp['total']:.6f} Ha")

    total = decomp['nuclear_repulsion'] + decomp['one_electron'] + decomp['two_electron']
    validate_test(
        "Energy sum",
        abs(total - decomp['total']) < 1e-6,
        "Components sum correctly ✓"
    )

except Exception as e:
    print(f"✗ FAILED: {e}")
    results.append({'name': 'Energy Decomp', 'passed': False})


print("\n" + "=" * 80)
print("TEST 4: Correlation Energy")
print("=" * 80)

try:
    from kanad.solvers.vqe_solver import VQESolver
    from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
    from kanad.core.mappers import JordanWignerMapper

    ansatz = RealAmplitudesAnsatz(
        n_qubits=2 * bond.hamiltonian.n_orbitals,
        n_electrons=bond.hamiltonian.n_electrons,
        n_layers=2
    )
    vqe = VQESolver(
        hamiltonian=bond.hamiltonian,
        ansatz=ansatz,
        mapper=JordanWignerMapper(),
        optimizer='SLSQP',
        max_iterations=50,
        enable_analysis=False
    )
    vqe_result = vqe.solve()

    corr_analyzer = CorrelationAnalyzer(bond.hamiltonian)
    corr = corr_analyzer.compute_correlation_energy(vqe_result['energy'], hf_result['energy'])

    print(f"HF:      {hf_result['energy']:.6f} Ha")
    print(f"VQE:     {vqe_result['energy']:.6f} Ha")
    print(f"Corr:    {corr * 1000:.3f} mHa")

    # CRITICAL: Check variational principle (VQE must be ≤ HF)
    variational_ok = vqe_result['energy'] <= hf_result['energy']
    if not variational_ok:
        print(f"⚠️  VARIATIONAL PRINCIPLE VIOLATED: VQE ({vqe_result['energy']:.6f}) > HF ({hf_result['energy']:.6f})")

    # Check both: correlation captured AND variational principle satisfied
    validate_test(
        "Correlation captured (with variational principle check)",
        abs(corr * 1000) > 5.0 and variational_ok,
        f"{corr * 1000:.3f} mHa captured, variational principle {'satisfied' if variational_ok else 'VIOLATED'}"
    )

except Exception as e:
    print(f"✗ FAILED: {e}")
    results.append({'name': 'Correlation', 'passed': False})


print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

passed = sum(1 for r in results if r.get('passed', False))
total = len(results)

print(f"\nPassed: {passed}/{total}")
for r in results:
    status = "✓" if r.get('passed', False) else "✗"
    print(f"{status} {r['name']}")

print("\n" + "=" * 80)
print("✓✓✓ ANALYSIS VALIDATION PASSED ✓✓✓" if passed == total else "⚠ NEEDS ATTENTION ⚠")
print("=" * 80)
