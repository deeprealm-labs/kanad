#!/usr/bin/env python3
"""
Comprehensive Excited States Solver Validation
Tests excited states for different molecular systems.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.solvers.excited_states_solver import ExcitedStatesSolver
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.ucc_ansatz import UCCAnsatz

def print_header(title: str):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")

# Test 1: H2 Molecule Excited States
print_header("TEST 1: H2 Molecule Excited States")

H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2_atom = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_bond = CovalentBond(H1, H2_atom)

# Get Hamiltonian
hamiltonian = h2_bond.hamiltonian
n_orbitals = hamiltonian.n_orbitals
n_spin_orbitals = 2 * n_orbitals
n_electrons = h2_bond.molecule.n_electrons

print(f"  H2 System: {n_spin_orbitals} spin orbitals, {n_electrons} electrons")

# Create solver
mapper = JordanWignerMapper()
ansatz = UCCAnsatz(n_qubits=n_spin_orbitals, n_electrons=n_electrons,
                    include_singles=True, include_doubles=False)  # Singles only for speed

try:
    solver = ExcitedStatesSolver(
        hamiltonian=hamiltonian,
        ansatz=ansatz,
        mapper=mapper,
        n_states=3,
        backend='classical'
    )

    print(f"  ✅ Solver initialized: {solver.n_parameters} parameters")

    # Solve for excited states
    print(f"  Computing excited states...")
    results = solver.solve_subspace_expansion(orthogonality_weight=100.0)

    print(f"\n  Excited State Energies:")
    print(f"  {'State':<10s} {'Energy (Ha)':>15s} {'Energy (eV)':>15s} {'ΔE (eV)':>12s} {'Conv':>6s}")
    print(f"  {'-'*70}")

    for i, (energy, exc_energy, converged) in enumerate(zip(
        results['energies'],
        results['excitation_energies'],
        results['converged']
    )):
        energy_ev = energy * 27.211
        conv_str = "✅" if converged else "❌"
        print(f"  S{i:<9d} {energy:>15.6f} {energy_ev:>15.4f} {exc_energy:>12.4f} {conv_str:>6s}")

    # Get spectrum
    spectrum = solver.get_excitation_spectrum()
    print(f"\n  Excitation Spectrum:")
    print(f"  {'Transition':<15s} {'ΔE (eV)':>12s} {'λ (nm)':>12s}")
    print(f"  {'-'*45}")
    for trans, exc_e, wavelength in zip(
        spectrum['transitions'],
        spectrum['excitation_energies'],
        spectrum['wavelengths']
    ):
        print(f"  {trans:<15s} {exc_e:>12.4f} {wavelength:>12.1f}")

    print(f"\n  ✅ H2 excited states validation complete")

except Exception as e:
    print(f"  ❌ Error: {str(e)[:100]}")
    import traceback
    traceback.print_exc()

# Test 2: Compare Excitation Energies
print_header("TEST 2: Excitation Energy Analysis")

print("\nChecking physical validity of excitation energies:\n")

try:
    # For H2, typical excitations:
    # σ → σ* (HOMO → LUMO): ~11-12 eV
    # Higher excitations: >15 eV

    print(f"  Expected H2 excitations:")
    print(f"    σ → σ*:     ~11-12 eV (first excitation)")
    print(f"    Higher:     >15 eV")

    if 'results' in locals():
        first_exc = results['excitation_energies'][1] if len(results['excitation_energies']) > 1 else 0
        second_exc = results['excitation_energies'][2] if len(results['excitation_energies']) > 2 else 0

        print(f"\n  Computed:")
        print(f"    First:      {first_exc:.2f} eV")
        print(f"    Second:     {second_exc:.2f} eV")

        # Validation
        if 5 < first_exc < 20:
            print(f"  ✅ First excitation energy reasonable")
        else:
            print(f"  ⚠️  First excitation energy may be incorrect")

        if second_exc > first_exc:
            print(f"  ✅ Energy ordering correct (E_S2 > E_S1)")
        else:
            print(f"  ❌ Energy ordering incorrect")

except Exception as e:
    print(f"  ❌ Error in analysis: {str(e)[:70]}")

# Test 3: State Orthogonality
print_header("TEST 3: State Orthogonality Check")

print("\nVerifying that excited states are orthogonal to ground state:\n")

try:
    if 'solver' in locals() and len(solver.state_parameters) >= 2:
        # Check overlap between states
        ground_params = solver.state_parameters[0]
        excited_params = solver.state_parameters[1]

        overlap = solver._compute_overlap(ground_params, excited_params)

        print(f"  Overlap |⟨S0|S1⟩|: {overlap:.6f}")

        if overlap < 0.1:
            print(f"  ✅ States are approximately orthogonal")
        elif overlap < 0.3:
            print(f"  ⚠️  States have some overlap (may need higher penalty)")
        else:
            print(f"  ❌ States are not orthogonal")

        # Check all pairs
        if len(solver.state_parameters) >= 3:
            print(f"\n  All state overlaps:")
            for i in range(len(solver.state_parameters)):
                for j in range(i+1, len(solver.state_parameters)):
                    overlap_ij = solver._compute_overlap(
                        solver.state_parameters[i],
                        solver.state_parameters[j]
                    )
                    status = "✅" if overlap_ij < 0.1 else "⚠️" if overlap_ij < 0.3 else "❌"
                    print(f"    {status} |⟨S{i}|S{j}⟩|: {overlap_ij:.6f}")

except Exception as e:
    print(f"  ❌ Error: {str(e)[:70]}")

# Summary
print_header("VALIDATION SUMMARY")

print("\n✅ = Passed")
print("⚠️  = Warning/needs review")
print("❌ = Failed")

print("\nTests Performed:")
print("  1. H2 excited states (3 states)")
print("  2. Excitation energy validation")
print("  3. State orthogonality check")

print("\nKey Findings:")
if 'results' in locals():
    print(f"  - Computed {len(results['energies'])} excited states")
    print(f"  - Excitation energies: {results['excitation_energies']} eV")
    print(f"  - Ground state energy: {results['energies'][0]:.6f} Ha")
else:
    print(f"  - Could not complete excited states calculation")

print("\nPotential Issues to Check:")
print("  1. Are excitation energies physically reasonable?")
print("  2. Are states properly orthogonal?")
print("  3. Does energy ordering make sense (E_S0 < E_S1 < E_S2)?")
print("  4. Do oscillator strengths indicate allowed/forbidden transitions?")
print("  5. Are wavelengths in expected range for UV-Vis spectra?")

print("\nNote: Excited states solver uses subspace expansion with orthogonality penalty.")
print("Higher penalty weights (>100) may be needed for better orthogonality.")

print("\n✅ Excited States Solver validation complete!")
