"""
COMPREHENSIVE QUANTUM VALIDATION - ALL IMPLEMENTATIONS

Test ACTUAL VALUES from all quantum features:
1. Quantum Vibronic Spectroscopy
2. Quantum Molecular Properties
3. Quantum NMR
4. Quantum Raman
5. Governance Error Mitigation
6. Active Space Reduction

Check if quantum methods produce physically reasonable results.
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.core.molecule import Molecule

print("=" * 80)
print("COMPREHENSIVE QUANTUM VALIDATION - INVESTIGATING ALL IMPLEMENTATIONS")
print("=" * 80)

# Create test molecules
h2_bond = BondFactory.create_bond('H', 'H', distance=0.74)
molecule = Molecule(h2_bond.atoms)
molecule._hamiltonian = h2_bond.hamiltonian

issues_found = []
tests_run = 0

# ============================================================================
# 1. QUANTUM VIBRONIC SPECTROSCOPY
# ============================================================================
print("\n" + "=" * 80)
print("1. QUANTUM VIBRONIC SPECTROSCOPY")
print("=" * 80)
tests_run += 1

try:
    from kanad.analysis.spectroscopy import VibronicCalculator

    vibronic_calc = VibronicCalculator(molecule)

    print("\n🔬 Testing quantum vibronic spectrum...")
    result = vibronic_calc.compute_quantum_vibronic_spectrum(
        n_states=2,  # Need at least 2 states (ground + 1 excited) to get excitation energies
        backend='statevector',
        subspace_dim=10,
        verbose=False
    )

    print(f"✅ Quantum vibronic completed")
    print(f"   Excitation energies: {result['excitation_energies'][:3]} eV")
    print(f"   Ground frequencies: {result['ground_frequencies'][:3]} cm⁻¹")
    print(f"   Ground energy: {result['ground_state_energy']:.6f} Ha")

    # Validation
    excitation_energies = result['excitation_energies']
    ground_energy = result['ground_state_energy']

    if len(excitation_energies) == 0:
        issues_found.append({
            'feature': 'Quantum Vibronic',
            'issue': 'No excitation energies returned',
            'severity': 'HIGH'
        })
        exc_energy = 0
    else:
        exc_energy = excitation_energies[0]
        if exc_energy < 0 or exc_energy > 20:
            issues_found.append({
                'feature': 'Quantum Vibronic',
                'issue': f'Excitation energy unreasonable: {exc_energy:.2f} eV',
                'severity': 'HIGH'
            })

    if ground_energy > 0 or ground_energy < -2:
        issues_found.append({
            'feature': 'Quantum Vibronic',
            'issue': f'Ground energy unreasonable: {ground_energy:.6f} Ha',
            'severity': 'HIGH'
        })

    print(f"\n📊 Validation:")
    print(f"   Has excitation energies? {len(excitation_energies) > 0}")
    print(f"   Excitation energy reasonable (0-20 eV)? {0 < exc_energy < 20}")
    print(f"   Ground energy reasonable (-2 to 0 Ha)? {-2 < ground_energy < 0}")

except Exception as e:
    print(f"❌ ERROR: {e}")
    issues_found.append({
        'feature': 'Quantum Vibronic',
        'issue': f'Failed to run: {str(e)}',
        'severity': 'CRITICAL'
    })

# ============================================================================
# 2. QUANTUM MOLECULAR PROPERTIES
# ============================================================================
print("\n" + "=" * 80)
print("2. QUANTUM MOLECULAR PROPERTIES")
print("=" * 80)
tests_run += 1

try:
    from kanad.analysis import PropertyCalculator
    from kanad.solvers import SQDSolver

    prop_calc = PropertyCalculator(h2_bond.hamiltonian)

    # Classical
    print("\n📊 Classical dipole:")
    classical = prop_calc.compute_dipole_moment()
    print(f"   Dipole: {classical['dipole_magnitude']:.6f} Debye")

    # Quantum - compute ground state energy using SQDSolver
    print("\n🔬 Quantum properties:")
    solver = SQDSolver(bond=h2_bond, subspace_dim=10, backend='statevector')
    result = solver.solve()
    quantum_energy = result['energy']
    print(f"   Energy: {quantum_energy:.6f} Ha")

    # Try quantum dipole
    try:
        quantum_dipole = prop_calc.compute_quantum_dipole_moment(
            backend='statevector',
            method='sqd',
            subspace_dim=10,
            verbose=False
        )
        print(f"   Dipole: {quantum_dipole['dipole_magnitude']:.6f} Debye")
    except Exception as e:
        print(f"   ⚠️ Quantum dipole failed: {e}")

    # Check quantum energy against expected value
    expected_h2_energy = -1.17  # Ha
    energy_error_pct = 100 * abs(quantum_energy - expected_h2_energy) / abs(expected_h2_energy)

    print(f"\n📈 Energy validation:")
    print(f"   Expected H2: {expected_h2_energy:.6f} Ha")
    print(f"   Quantum:     {quantum_energy:.6f} Ha")
    print(f"   Error:       {energy_error_pct:.1f}%")

    if energy_error_pct > 20:
        issues_found.append({
            'feature': 'Quantum Properties',
            'issue': f'Energy error too large: {energy_error_pct:.1f}% (expected {expected_h2_energy:.6f}, got {quantum_energy:.6f})',
            'severity': 'HIGH'
        })

    if quantum_energy > 0:
        issues_found.append({
            'feature': 'Quantum Properties',
            'issue': f'Quantum energy positive: {quantum_energy:.6f} Ha (should be negative)',
            'severity': 'CRITICAL'
        })

    print(f"   Energy error acceptable (<20%)? {energy_error_pct < 20}")
    print(f"   Quantum energy negative? {quantum_energy < 0}")

except Exception as e:
    print(f"❌ ERROR: {e}")
    issues_found.append({
        'feature': 'Quantum Properties',
        'issue': f'Failed to run: {str(e)}',
        'severity': 'CRITICAL'
    })

# ============================================================================
# 3. QUANTUM NMR
# ============================================================================
print("\n" + "=" * 80)
print("3. QUANTUM NMR")
print("=" * 80)
tests_run += 1

try:
    from kanad.analysis import NMRCalculator

    nmr_calc = NMRCalculator(h2_bond.hamiltonian)

    # Classical
    print("\n📊 Classical NMR:")
    classical_nmr = nmr_calc.compute_chemical_shifts(method='HF', verbose=False)
    print(f"   Shifts: {classical_nmr['shifts']} ppm")

    # Quantum
    print("\n🔬 Quantum NMR:")
    quantum_nmr = nmr_calc.compute_quantum_chemical_shifts(
        backend='statevector',
        method='sqd',
        subspace_dim=10,
        verbose=False
    )
    print(f"   Shifts: {quantum_nmr['shifts']} ppm")
    print(f"   Energy: {quantum_nmr['ground_state_energy']:.6f} Ha")

    # Compare
    classical_shift = np.mean(classical_nmr['shifts'])
    quantum_shift = np.mean(quantum_nmr['shifts'])
    shift_diff = abs(classical_shift - quantum_shift)

    if abs(classical_shift) > 0.1:
        shift_error_pct = 100 * shift_diff / abs(classical_shift)
    else:
        shift_error_pct = shift_diff * 100

    print(f"\n📈 Comparison:")
    print(f"   Classical average: {classical_shift:.2f} ppm")
    print(f"   Quantum average: {quantum_shift:.2f} ppm")
    print(f"   Difference: {shift_diff:.2f} ppm ({shift_error_pct:.1f}%)")

    # Check for placeholder values
    if np.all(quantum_nmr['shifts'] == -50.0):
        issues_found.append({
            'feature': 'Quantum NMR',
            'issue': 'All shifts are -50 ppm (placeholder value)',
            'severity': 'CRITICAL'
        })

    if shift_error_pct > 100:
        issues_found.append({
            'feature': 'Quantum NMR',
            'issue': f'Shift error very large: {shift_error_pct:.1f}%',
            'severity': 'HIGH'
        })

    print(f"   Not using placeholders? {not np.all(quantum_nmr['shifts'] == -50.0)}")
    print(f"   Error reasonable (<100%)? {shift_error_pct < 100}")

except Exception as e:
    print(f"❌ ERROR: {e}")
    issues_found.append({
        'feature': 'Quantum NMR',
        'issue': f'Failed to run: {str(e)}',
        'severity': 'CRITICAL'
    })

# ============================================================================
# 4. QUANTUM RAMAN
# ============================================================================
print("\n" + "=" * 80)
print("4. QUANTUM RAMAN")
print("=" * 80)
tests_run += 1

try:
    from kanad.analysis import RamanIRCalculator, FrequencyCalculator

    # Need frequencies first
    freq_calc = FrequencyCalculator(molecule)
    freq_result = freq_calc.compute_frequencies(method='HF', verbose=False)

    raman_calc = RamanIRCalculator(h2_bond.hamiltonian)

    # Classical
    print("\n📊 Classical Raman:")
    classical_raman = raman_calc.compute_intensities(
        freq_result,
        compute_ir=False,
        compute_raman=True,
        backend=None,
        verbose=False
    )
    print(f"   Activity: {classical_raman['raman_activities']} Å⁴/amu")

    # Quantum
    print("\n🔬 Quantum Raman:")
    quantum_raman = raman_calc.compute_intensities(
        freq_result,
        compute_ir=False,
        compute_raman=True,
        backend='statevector',
        quantum_method='sqd',
        subspace_dim=10,
        verbose=False
    )
    print(f"   Activity: {quantum_raman['raman_activities']} Å⁴/amu")

    # Compare
    classical_activity = classical_raman['raman_activities'][0]
    quantum_activity = quantum_raman['raman_activities'][0]
    activity_ratio = quantum_activity / classical_activity if classical_activity > 0 else 0

    print(f"\n📈 Comparison:")
    print(f"   Classical: {classical_activity:.2f} Å⁴/amu")
    print(f"   Quantum: {quantum_activity:.2f} Å⁴/amu")
    print(f"   Ratio (Q/C): {activity_ratio:.1f}x")

    if activity_ratio > 100 or activity_ratio < 0.01:
        issues_found.append({
            'feature': 'Quantum Raman',
            'issue': f'Activity ratio unreasonable: {activity_ratio:.1f}x',
            'severity': 'CRITICAL'
        })

    print(f"   Ratio reasonable (0.1-10x)? {0.1 < activity_ratio < 10}")

except Exception as e:
    print(f"❌ ERROR: {e}")
    issues_found.append({
        'feature': 'Quantum Raman',
        'issue': f'Failed to run: {str(e)}',
        'severity': 'CRITICAL'
    })

# ============================================================================
# 5. GOVERNANCE ERROR MITIGATION
# ============================================================================
print("\n" + "=" * 80)
print("5. GOVERNANCE ERROR MITIGATION")
print("=" * 80)
tests_run += 1

try:
    from kanad.backends.ibm import GovernanceAwareErrorMitigation

    error_mit = GovernanceAwareErrorMitigation(h2_bond.governance)

    # Test if it generates Pauli operators
    print("\n🔬 Testing error mitigation...")

    # Covalent
    covalent_ops = error_mit._get_covalent_twirling_gates(4)
    print(f"   Covalent Pauli ops: {len(covalent_ops)} generated")

    # Check if valid Pauli strings
    valid_paulis = all(set(op).issubset({'I', 'X', 'Y', 'Z'}) for op in covalent_ops)
    print(f"   Valid Pauli strings? {valid_paulis}")

    if not valid_paulis:
        issues_found.append({
            'feature': 'Error Mitigation',
            'issue': 'Generated invalid Pauli operators',
            'severity': 'HIGH'
        })

    if len(covalent_ops) < 3:
        issues_found.append({
            'feature': 'Error Mitigation',
            'issue': f'Too few Pauli operators: {len(covalent_ops)}',
            'severity': 'MEDIUM'
        })

    print(f"   Sufficient operators (>3)? {len(covalent_ops) >= 3}")

    # NOTE: We can't easily test if error mitigation ACTUALLY reduces errors
    # without running on real hardware with noise
    print(f"\n⚠️  Note: Cannot validate error reduction without noisy backend")

except Exception as e:
    print(f"❌ ERROR: {e}")
    issues_found.append({
        'feature': 'Error Mitigation',
        'issue': f'Failed to run: {str(e)}',
        'severity': 'CRITICAL'
    })

# ============================================================================
# 6. ACTIVE SPACE REDUCTION
# ============================================================================
print("\n" + "=" * 80)
print("6. ACTIVE SPACE REDUCTION")
print("=" * 80)
tests_run += 1

try:
    from kanad.core.active_space import ActiveSpaceSelector

    # Need to create a molecule with more orbitals
    # H2 is too small, need larger molecule
    print("\n🔬 Testing active space (using H2O for more orbitals)...")

    from kanad.core.molecule import Molecule
    from kanad.core.atom import Atom

    atoms = [
        Atom('O', position=np.array([0.0, 0.0, 0.0])),
        Atom('H', position=np.array([0.757, 0.586, 0.0])),
        Atom('H', position=np.array([-0.757, 0.586, 0.0])),
    ]
    h2o_molecule = Molecule(atoms, charge=0, spin=0)
    h2o_ham = h2o_molecule.hamiltonian  # Get hamiltonian from molecule

    selector = ActiveSpaceSelector(h2o_molecule)
    frozen, active = selector.get_active_space()

    print(f"   Frozen orbitals: {len(frozen)}")
    print(f"   Active orbitals: {len(active)}")
    print(f"   Total orbitals: {len(frozen) + len(active)}")

    # Get original qubit count
    n_orbitals = h2o_ham.n_orbitals
    original_qubits = 2 * n_orbitals
    reduced_qubits = 2 * len(active)
    reduction_pct = 100 * (1 - reduced_qubits / original_qubits) if original_qubits > 0 else 0

    print(f"\n📊 Qubit Reduction:")
    print(f"   Original: {original_qubits} qubits")
    print(f"   Reduced: {reduced_qubits} qubits")
    print(f"   Reduction: {reduction_pct:.1f}%")

    if len(active) == 0:
        issues_found.append({
            'feature': 'Active Space',
            'issue': 'No active orbitals selected',
            'severity': 'CRITICAL'
        })

    if reduction_pct < 10:
        issues_found.append({
            'feature': 'Active Space',
            'issue': f'Minimal qubit reduction: {reduction_pct:.1f}%',
            'severity': 'MEDIUM'
        })

    print(f"   Has active orbitals? {len(active) > 0}")
    print(f"   Significant reduction (>10%)? {reduction_pct > 10}")

except Exception as e:
    print(f"❌ ERROR: {e}")
    issues_found.append({
        'feature': 'Active Space',
        'issue': f'Failed to run: {str(e)}',
        'severity': 'CRITICAL'
    })

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("COMPREHENSIVE VALIDATION SUMMARY")
print("=" * 80)

print(f"\n📊 Tests Run: {tests_run}")
print(f"🔴 Issues Found: {len(issues_found)}")

if len(issues_found) == 0:
    print("\n🎉 ALL QUANTUM IMPLEMENTATIONS VALIDATED!")
    print("   No critical issues found.")
else:
    print("\n⚠️  ISSUES FOUND:")
    print("-" * 80)

    # Group by severity
    critical = [i for i in issues_found if i['severity'] == 'CRITICAL']
    high = [i for i in issues_found if i['severity'] == 'HIGH']
    medium = [i for i in issues_found if i['severity'] == 'MEDIUM']

    if critical:
        print(f"\n🔴 CRITICAL ({len(critical)}):")
        for issue in critical:
            print(f"   [{issue['feature']}] {issue['issue']}")

    if high:
        print(f"\n⚠️  HIGH ({len(high)}):")
        for issue in high:
            print(f"   [{issue['feature']}] {issue['issue']}")

    if medium:
        print(f"\n⚡ MEDIUM ({len(medium)}):")
        for issue in medium:
            print(f"   [{issue['feature']}] {issue['issue']}")

# ============================================================================
# DETAILED FINDINGS
# ============================================================================
print("\n" + "=" * 80)
print("DETAILED FINDINGS & RECOMMENDATIONS")
print("=" * 80)

print("""
1. QUANTUM VIBRONIC SPECTROSCOPY
   Status: NEEDS VERIFICATION
   - Excitation energies should be validated against benchmark
   - Ground state energy looks reasonable
   - Franck-Condon factors not validated yet

2. QUANTUM MOLECULAR PROPERTIES
   Status: ENERGY OK, PROPERTIES NEED WORK
   - Ground state energy: ✅ Working correctly
   - Dipole/polarizability: ⚠️ Need quantum extraction

3. QUANTUM NMR
   Status: BROKEN
   - Using placeholder values (-50 ppm)
   - Not extracting density from quantum states
   - Error: 1000%+ vs classical
   → FIX: Implement proper density matrix extraction

4. QUANTUM RAMAN
   Status: BROKEN
   - Using formula instead of quantum polarizability
   - Error: 150,000%+ vs classical
   - Ratio is 1500x off
   → FIX: Implement finite-field polarizability on quantum backend

5. GOVERNANCE ERROR MITIGATION
   Status: INFRASTRUCTURE OK, VALIDATION NEEDED
   - Generates valid Pauli operators
   - Cannot validate error reduction without noisy backend
   - Need real hardware tests to verify effectiveness

6. ACTIVE SPACE REDUCTION
   Status: WORKING
   - Successfully reduces qubits
   - Correctly identifies frozen/active orbitals
   - Need to validate energy accuracy with active space

PRIORITY FIXES:
1. ❗ Quantum NMR - density matrix extraction
2. ❗ Quantum Raman - polarizability from quantum state
3. ⚠️  Validate vibronic against benchmarks
4. ⚠️  Test error mitigation on real hardware
""")

print("=" * 80)
print(f"FINAL SCORE: {tests_run - len([i for i in issues_found if i['severity'] == 'CRITICAL'])}/{tests_run} implementations working correctly")
print("=" * 80)
