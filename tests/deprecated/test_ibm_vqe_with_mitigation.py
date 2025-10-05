"""
Test IBM VQE Solver with FULL Error Mitigation on H2.

Error Mitigation Techniques:
1. Dynamical Decoupling (XY4 sequence)
2. Gate Twirling (all strategies)
3. Zero-Noise Extrapolation (resilience_level=2)
4. Measurement Error Mitigation (built-in)
"""

import numpy as np
import os
from kanad.core.atom import Atom
from kanad.bonds.bond_factory import BondFactory
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.backends.ibm import IBMVQESolver

print("=" * 80)
print("IBM VQE with FULL ERROR MITIGATION")
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
print(f"   ✓ UCCSD ansatz created")

# Get credentials
token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
crn = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
channel = 'ibm_cloud' if crn else 'ibm_quantum_platform'

if not token:
    print("\n❌ ERROR: No IBM Quantum token found")
    exit(1)

# Initialize VQE with FULL error mitigation
print("\n3. Initializing VQE with FULL error mitigation...")
print("   Error Mitigation Stack:")
print("   ✓ Resilience Level: 2 (ZNE + Measurement Calibration)")
print("   ✓ Dynamical Decoupling: XY4 sequence")
print("   ✓ Gate Twirling: All (active + passive)")
print("   ✓ ZNE Extrapolator: Exponential")
print("   ✓ Optimization Level: 3 (maximum)")

try:
    vqe = IBMVQESolver(
        hamiltonian=hamiltonian,
        ansatz=ansatz,
        backend_name='ibm_torino',
        token=token,
        instance=crn,
        channel=channel,
        # Classical optimization
        optimizer='COBYLA',
        max_iterations=30,  # More iterations for better convergence
        # Quantum execution
        shots=8192,  # More shots for better statistics
        optimization_level=3,  # Maximum transpiler optimization
        # Base error mitigation
        resilience_level=2,  # ZNE + measurement error mitigation
        # Advanced error mitigation
        enable_dynamical_decoupling=True,
        dynamical_decoupling_sequence='XY4',  # Best DD sequence
        enable_twirling=True,
        twirling_strategy='all',  # Both gate and measurement twirling
        zne_extrapolator='exponential',  # Exponential extrapolation for ZNE
    )

    print(f"   ✓ VQE solver initialized with full error mitigation")

    # Run VQE
    print("\n4. Running VQE optimization on IBM Quantum...")
    print("-" * 80)

    result = vqe.solve()

    print("\n" + "=" * 80)
    print("RESULTS (WITH FULL ERROR MITIGATION)")
    print("=" * 80)
    print(f"Final Energy: {result['energy']:.6f} Ha")
    print(f"             {result['energy_ev']:.4f} eV")
    print(f"Iterations: {result['iterations']}")
    print(f"Converged: {result['converged']}")
    print(f"Backend: {result['backend_info']['name']}")
    print("=" * 80)

    # Compare with expected
    expected_energy = -1.59  # H2 ground state
    error = abs(result['energy'] - expected_energy)
    error_pct = error / abs(expected_energy) * 100

    print(f"\nValidation:")
    print(f"  Expected: ~{expected_energy:.2f} Ha")
    print(f"  Got: {result['energy']:.6f} Ha")
    print(f"  Error: {error:.4f} Ha ({error_pct:.1f}%)")

    if error_pct < 15:
        print("  ✅ GOOD: Error < 15% with noise mitigation")
    elif error_pct < 30:
        print("  ⚠️  ACCEPTABLE: Error < 30% on NISQ hardware")
    else:
        print("  ❌ POOR: Error > 30% - may need more shots/iterations")

    print("\nError Mitigation Summary:")
    print("  This run used:")
    print("    - Dynamical Decoupling (XY4)")
    print("    - Gate Twirling (all)")
    print("    - Zero-Noise Extrapolation (exponential)")
    print("    - Measurement Error Calibration")
    print("    - Transpiler optimization (level 3)")
    print(f"    - {result.get('iterations', 0)} VQE iterations")
    print(f"    - 8192 shots per measurement")

    print("\n✅ IBM VQE with error mitigation completed!")

except Exception as e:
    print(f"\n❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
    exit(1)
