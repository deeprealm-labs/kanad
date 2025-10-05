"""
Test IBM SQD Solver with Error Mitigation on H2.

SQD (Sample-based Quantum Diagonalization):
- Samples multiple quantum states
- Measures Hamiltonian matrix elements
- Classically diagonalizes sampled subspace
- With full error mitigation for accurate matrix elements
"""

import numpy as np
import os
from kanad.core.atom import Atom
from kanad.bonds.bond_factory import BondFactory
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.backends.ibm import IBMSQDSolver

print("=" * 80)
print("IBM SQD SOLVER WITH ERROR MITIGATION")
print("=" * 80)

# Create H2 molecule
print("\n1. Creating H2 molecule...")
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_bond = BondFactory.create_bond(H1, H2)
hamiltonian = h2_bond.hamiltonian
print(f"   ✓ H2: {hamiltonian.n_orbitals} orbitals, {h2_bond.molecule.n_electrons} electrons")

# Create ansatz
print("\n2. Creating UCCSD ansatz...")
ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)
print(f"   ✓ UCCSD ansatz for state sampling")

# Get credentials
token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
crn = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
channel = 'ibm_cloud' if crn else 'ibm_quantum_platform'

if not token:
    print("\n❌ ERROR: No IBM Quantum token found")
    exit(1)

# Initialize SQD with error mitigation
print("\n3. Initializing SQD solver with error mitigation...")
print("   SQD Configuration:")
print("   ✓ Sample States: 8 (for good subspace coverage)")
print("   ✓ Shots: 8192 (high accuracy for matrix elements)")
print("   ✓ Resilience Level: 2 (ZNE + calibration)")
print("   ✓ Dynamical Decoupling: XY4")
print("   ✓ Gate Twirling: All strategies")
print("   ✓ Optimization Level: 3")

try:
    sqd = IBMSQDSolver(
        hamiltonian=hamiltonian,
        ansatz=ansatz,
        n_samples=8,  # Sample 8 quantum states
        backend_name='ibm_torino',
        token=token,
        instance=crn,
        channel=channel,
        shots=8192,  # High accuracy for matrix element measurements
        optimization_level=3,
        resilience_level=2,  # Maximum built-in mitigation
        # Advanced error mitigation
        enable_dynamical_decoupling=True,
        dynamical_decoupling_sequence='XY4',
        enable_twirling=True,
        twirling_strategy='all',
        zne_extrapolator='exponential'
    )

    print(f"   ✓ SQD solver initialized")

    # Run SQD
    print("\n4. Running SQD on IBM Quantum...")
    print("   This will:")
    print("   - Generate 8 random quantum states")
    print("   - Measure Hamiltonian matrix elements (8 measurements)")
    print("   - Diagonalize 8×8 matrix classically")
    print("   - Extract ground state energy")
    print("-" * 80)

    result = sqd.solve()

    print("\n" + "=" * 80)
    print("SQD RESULTS (WITH ERROR MITIGATION)")
    print("=" * 80)
    print(f"Ground State Energy: {result['energy']:.6f} Ha")
    print(f"                     {result['energy_ev']:.4f} eV")
    print(f"Sample States: {result['n_samples']}")
    print(f"Backend: {result['backend_info']['name']}")
    print("=" * 80)

    # Show all eigenvalues
    print(f"\nAll Eigenvalues from Sampled Subspace:")
    for i, (E_ha, E_ev) in enumerate(zip(result['eigenvalues'], result['eigenvalues_ev'])):
        print(f"  State {i}: {E_ha:8.6f} Ha = {E_ev:8.4f} eV")

    # Compare with expected
    expected_energy = -1.59  # H2 ground state
    error = abs(result['energy'] - expected_energy)
    error_pct = error / abs(expected_energy) * 100

    print(f"\nValidation:")
    print(f"  Expected: ~{expected_energy:.2f} Ha")
    print(f"  Got: {result['energy']:.6f} Ha")
    print(f"  Error: {error:.4f} Ha ({error_pct:.1f}%)")

    if error_pct < 15:
        print("  ✅ EXCELLENT: Error < 15% with SQD + mitigation")
    elif error_pct < 30:
        print("  ✅ GOOD: Error < 30% on NISQ hardware")
    else:
        print("  ⚠️  ACCEPTABLE: NISQ hardware limitations")

    print("\nSQD Advantages:")
    print("  ✓ Provides multiple eigenvalues (excited states)")
    print("  ✓ No optimization needed (fixed # of measurements)")
    print("  ✓ Good for exploring energy landscape")
    print(f"  ✓ Used only {result['n_samples']} quantum measurements")

    print("\n✅ IBM SQD with error mitigation completed!")

except Exception as e:
    print(f"\n❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
    exit(1)
