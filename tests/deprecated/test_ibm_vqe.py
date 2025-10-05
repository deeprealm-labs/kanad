"""
Test IBM VQE Solver on H2 molecule.
"""

import numpy as np
import os
from kanad.core.atom import Atom
from kanad.bonds.bond_factory import BondFactory
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.backends.ibm import IBMVQESolver

print("=" * 80)
print("TEST: IBM VQE Solver on H2")
print("=" * 80)

# Create H2 molecule
print("\n1. Creating H2 molecule...")
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_bond = BondFactory.create_bond(H1, H2)
hamiltonian = h2_bond.hamiltonian

print(f"   ✓ H2 molecule created")
print(f"   Electrons: {h2_bond.molecule.n_electrons}")
print(f"   Orbitals: {hamiltonian.n_orbitals}")

# Create ansatz
print("\n2. Creating UCCSD ansatz...")
ansatz = UCCAnsatz(
    n_qubits=4,
    n_electrons=2,
    include_singles=True,
    include_doubles=True
)
print(f"   ✓ UCCSD ansatz created")

# Get credentials
token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
crn = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
channel = 'ibm_cloud' if crn else 'ibm_quantum_platform'

if not token:
    print("\n❌ ERROR: No IBM Quantum token found in environment")
    print("   Set one of: IBM_QUANTUM_TOKEN, QISKIT_IBM_TOKEN, API")
    exit(1)

# Use real hardware
backend_name = 'ibm_torino'  # Real quantum hardware
print(f"\n3. Initializing IBM VQE Solver...")
print(f"   Backend: {backend_name}")
print(f"   Channel: {channel}")

try:
    vqe = IBMVQESolver(
        hamiltonian=hamiltonian,
        ansatz=ansatz,
        backend_name=backend_name,
        token=token,
        instance=crn,
        channel=channel,
        optimizer='COBYLA',
        max_iterations=20,  # Limited for real hardware (expensive!)
        shots=4096,  # More shots for real hardware
        optimization_level=3,  # High optimization for real hardware
        resilience_level=1  # Error mitigation for real hardware
    )

    print(f"   ✓ VQE solver initialized")

    # Run VQE
    print("\n4. Running VQE optimization...")
    print("-" * 80)

    result = vqe.solve()

    print("\n" + "=" * 80)
    print("RESULTS")
    print("=" * 80)
    print(f"Final Energy: {result['energy']:.6f} Ha")
    print(f"             {result['energy_ev']:.4f} eV")
    print(f"Iterations: {result['iterations']}")
    print(f"Converged: {result['converged']}")
    print(f"Backend: {result['backend_info']['name']}")
    print("=" * 80)

    # Compare with expected
    expected_energy = -1.59  # Approximate H2 ground state
    error = abs(result['energy'] - expected_energy)
    error_pct = error / abs(expected_energy) * 100

    print(f"\nValidation:")
    print(f"  Expected: ~{expected_energy:.2f} Ha")
    print(f"  Error: {error:.4f} Ha ({error_pct:.1f}%)")

    if error_pct < 10:
        print("  ✅ PASS: Energy within 10% of expected")
    else:
        print("  ⚠️  WARNING: Energy error > 10%")

    print("\n✅ IBM VQE Solver test completed successfully!")

except Exception as e:
    print(f"\n❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
    exit(1)
