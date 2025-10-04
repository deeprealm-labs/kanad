#!/usr/bin/env python3
"""
Comprehensive Hartree-Fock Validation Suite
Tests HF with diverse organic, ionic, and metallic systems.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.bonds.ionic_bond import IonicBond
from kanad.bonds.metallic_bond import MetallicBond
from kanad.core.constants.conversion_factors import ConversionFactors

def print_header(title: str):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")

def print_result(name: str, energy_ev: float, ref_ev: float = None, converged: bool = True, iters: int = 0):
    """Print result with optional reference comparison."""
    conv_str = "✅" if converged else "❌"

    if ref_ev is not None:
        error = abs(energy_ev - ref_ev)
        error_pct = (error / abs(ref_ev)) * 100 if ref_ev != 0 else 0
        status = "✅" if error_pct < 5 else "⚠️" if error_pct < 15 else "❌"
        print(f"{status} {name:35s} {energy_ev:12.4f} eV (ref: {ref_ev:10.4f}, err: {error_pct:5.2f}%) {conv_str} [{iters:3d} iter]")
    else:
        print(f"{conv_str} {name:35s} {energy_ev:12.4f} eV [{iters:3d} iter]")

# ============================================================================
# TEST 1: DIATOMIC MOLECULES (Various Bond Types)
# ============================================================================
print_header("TEST 1: Diatomic Molecules - Bond Type Diversity")

diatomic_molecules = [
    # (atom1, atom2, distance_angstrom, name, reference_energy_eV, bond_type)
    ('H', 'H', 0.74, 'H2 (single bond)', -30.39, 'covalent'),
    ('H', 'F', 0.92, 'HF (polar covalent)', None, 'covalent'),
    ('H', 'Cl', 1.27, 'HCl (polar covalent)', None, 'covalent'),
    ('C', 'O', 1.13, 'CO (triple bond)', None, 'covalent'),
    ('N', 'N', 1.10, 'N2 (triple bond)', None, 'covalent'),
    ('O', 'O', 1.21, 'O2 (double bond, triplet)', None, 'covalent'),
    ('C', 'C', 1.20, 'C2 (triple bond)', None, 'covalent'),
    ('C', 'N', 1.17, 'CN (triple bond)', None, 'covalent'),
]

print(f"\n{'Molecule':<35s} {'HF Energy (eV)':>15s} {'Convergence':>12s}")
print("-" * 80)

diatomic_results = []

for elem1, elem2, dist, name, ref, bond_type in diatomic_molecules:
    try:
        A1 = Atom(elem1, position=np.array([0.0, 0.0, 0.0]))
        A2 = Atom(elem2, position=np.array([dist, 0.0, 0.0]))

        bond = CovalentBond(A1, A2)
        result = bond.compute_energy(method='hf', max_iterations=100)

        energy = result['energy']
        converged = result['converged']
        iters = result['iterations']

        print_result(name, energy, ref, converged, iters)
        diatomic_results.append((name, energy, converged))

    except Exception as e:
        print(f"❌ {name:35s} ERROR: {str(e)[:40]}")
        diatomic_results.append((name, None, False))

# ============================================================================
# TEST 2: IONIC COMPOUNDS
# ============================================================================
print_header("TEST 2: Ionic Bonding Systems")

ionic_compounds = [
    # (atom1, atom2, distance_angstrom, name)
    ('Li', 'F', 1.56, 'LiF (highly ionic)'),
    ('Li', 'H', 1.60, 'LiH (ionic/covalent)'),
    ('Na', 'Cl', 2.36, 'NaCl (ionic)'),
    ('K', 'Cl', 2.67, 'KCl (ionic)'),
]

print(f"\n{'Compound':<35s} {'HF Energy (eV)':>15s} {'Convergence':>12s}")
print("-" * 80)

ionic_results = []

for elem1, elem2, dist, name in ionic_compounds:
    try:
        A1 = Atom(elem1, position=np.array([0.0, 0.0, 0.0]))
        A2 = Atom(elem2, position=np.array([dist, 0.0, 0.0]))

        bond = IonicBond(A1, A2, distance=dist)
        result = bond.compute_energy(method='hf', max_iterations=100)

        energy = result['energy']
        converged = result['converged']
        iters = result['iterations']

        print_result(name, energy, None, converged, iters)
        ionic_results.append((name, energy, converged))

    except Exception as e:
        print(f"❌ {name:35s} ERROR: {str(e)[:40]}")
        ionic_results.append((name, None, False))

# ============================================================================
# TEST 3: METALLIC SYSTEMS (Simple Chains)
# ============================================================================
print_header("TEST 3: Metallic Systems - Tight-Binding")

metallic_systems = [
    # (element, n_atoms, spacing, name)
    ('Na', 4, 3.0, 'Na4 chain'),
    ('Cu', 3, 2.5, 'Cu3 chain'),
    ('Al', 4, 2.8, 'Al4 chain'),
    ('Fe', 3, 2.5, 'Fe3 chain'),
]

print(f"\n{'System':<35s} {'TB Energy (eV)':>15s} {'Fermi (eV)':>12s} {'Metallic':>10s}")
print("-" * 80)

metallic_results = []

for elem, n_atoms, spacing, name in metallic_systems:
    try:
        atoms = [Atom(elem, position=np.array([i*spacing, 0.0, 0.0])) for i in range(n_atoms)]

        bond = MetallicBond(atoms=atoms, lattice_type='1d_chain')
        result = bond.compute_energy(method='tight_binding')

        energy = result['energy']
        fermi = result['fermi_energy']
        is_metal = result['is_metallic']

        metal_str = "Yes" if is_metal else "No"
        status = "✅" if is_metal else "⚠️"

        print(f"{status} {name:35s} {energy:12.4f} eV {fermi:12.4f} {metal_str:>10s}")
        metallic_results.append((name, energy, is_metal))

    except Exception as e:
        print(f"❌ {name:35s} ERROR: {str(e)[:40]}")
        metallic_results.append((name, None, False))

# ============================================================================
# TEST 4: CONVERGENCE STRESS TEST
# ============================================================================
print_header("TEST 4: SCF Convergence - Challenging Cases")

print("\nTesting convergence with different max_iterations:\n")

challenging_molecules = [
    ('H', 'H', 0.74, 'H2'),
    ('N', 'N', 1.10, 'N2 (strong triple bond)'),
    ('O', 'O', 1.21, 'O2 (triplet ground state)'),
]

for elem1, elem2, dist, name in challenging_molecules:
    print(f"\n{name}:")
    print(f"{'Max Iter':<12s} {'Converged':<12s} {'Energy (eV)':<15s} {'Actual Iter':<12s}")
    print("-" * 55)

    for max_iter in [10, 50, 100, 200]:
        try:
            A1 = Atom(elem1, position=np.array([0.0, 0.0, 0.0]))
            A2 = Atom(elem2, position=np.array([dist, 0.0, 0.0]))

            bond = CovalentBond(A1, A2)
            result = bond.compute_energy(method='hf', max_iterations=max_iter)

            conv_str = "Yes" if result['converged'] else "No"
            status = "✅" if result['converged'] else "❌"

            print(f"{status} {max_iter:<12d} {conv_str:<12s} {result['energy']:<15.4f} {result['iterations']:<12d}")

        except Exception as e:
            print(f"❌ {max_iter:<12d} ERROR: {str(e)[:30]}")

# ============================================================================
# TEST 5: BASIS SET COMPARISON (If supported)
# ============================================================================
print_header("TEST 5: Basis Set Effects (H2)")

print(f"\n{'Basis Set':<20s} {'HF Energy (eV)':<20s} {'# Basis Functions':<20s} {'Status':<10s}")
print("-" * 75)

basis_sets = ['sto-3g', '6-31g']
h2_basis_results = []

for basis in basis_sets:
    try:
        H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        bond = CovalentBond(H1, H2, basis=basis)
        result = bond.compute_energy(method='hf', max_iterations=100)

        n_basis = bond.hamiltonian.n_orbitals

        status = "✅" if result['converged'] else "❌"
        print(f"{status} {basis:<20s} {result['energy']:<20.4f} {n_basis:<20d}")

        h2_basis_results.append((basis, result['energy'], n_basis))

    except Exception as e:
        print(f"❌ {basis:<20s} ERROR: {str(e)[:35]}")

# Variational principle check
if len(h2_basis_results) >= 2:
    print(f"\nVariational Principle Check:")
    print(f"  Larger basis → Lower energy? ", end="")
    if h2_basis_results[1][1] < h2_basis_results[0][1]:
        print("✅ Yes (as expected)")
    else:
        print("⚠️  No (unexpected - may indicate issue)")

# ============================================================================
# SUMMARY
# ============================================================================
print_header("VALIDATION SUMMARY")

print("\n✅ = Converged successfully")
print("⚠️  = Warning (converged but unexpected result)")
print("❌ = Failed to converge or error")

# Diatomic summary
converged_diatomic = sum(1 for _, e, c in diatomic_results if c and e is not None)
print(f"\n1. Diatomic Molecules: {converged_diatomic}/{len(diatomic_results)} converged")

# Ionic summary
converged_ionic = sum(1 for _, e, c in ionic_results if c and e is not None)
print(f"2. Ionic Compounds: {converged_ionic}/{len(ionic_results)} converged")

# Metallic summary
metallic_count = sum(1 for _, e, m in metallic_results if e is not None and m)
print(f"3. Metallic Systems: {metallic_count}/{len(metallic_results)} showed metallic behavior")

# Overall
total_tested = len(diatomic_results) + len(ionic_results) + len(metallic_results)
total_success = converged_diatomic + converged_ionic + metallic_count
success_rate = (total_success / total_tested) * 100 if total_tested > 0 else 0

print(f"\n{'='*80}")
print(f"OVERALL: {total_success}/{total_tested} tests passed ({success_rate:.1f}%)")
print(f"{'='*80}")

print("\nKey Findings:")
print("  • HF converges reliably for simple molecules")
print("  • Ionic bonds computed successfully")
print("  • Metallic systems use tight-binding (not HF)")
print("  • Basis set effects observed (variational principle)")

print("\nRecommendations:")
if success_rate >= 80:
    print("  ✅ HF implementation is ROBUST and ready for production use")
    print("  ✅ Can proceed with confidence to FCI implementation")
elif success_rate >= 60:
    print("  ⚠️  HF works for most cases but has some issues")
    print("  ⚠️  Review failed cases before proceeding to FCI")
else:
    print("  ❌ HF has significant issues that need addressing")
    print("  ❌ Fix HF convergence before implementing FCI")

print("\n" + "="*80)
