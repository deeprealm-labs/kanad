#!/usr/bin/env python3
"""
Test using qiskit-nature's FermionicOp for proper UCC implementation
"""

import sys
sys.path.insert(0, 'kanad/external/qiskit_nature')

import numpy as np
from kanad.external.qiskit_nature.operators.fermionic_op import FermionicOp
from kanad.external.qiskit_nature.mappers.jordan_wigner_mapper import JordanWignerMapper

print("="*80)
print("TEST: Proper UCC Double Excitation using FermionicOp")
print("="*80)

# For H2: 2 spatial orbitals, 4 spin orbitals (blocked ordering)
# HF state: |0101⟩ = orbital 0 alpha occupied, orbital 0 beta occupied
# Want to excite to: |1010⟩ = orbital 1 alpha occupied, orbital 1 beta occupied

# Double excitation: move BOTH electrons from orbital 0 to orbital 1
# In spin-orbital indices (blocked: 0,1=alpha, 2,3=beta):
# - Occupied: 0 (orb 0 alpha), 2 (orb 0 beta)
# - Virtual: 1 (orb 1 alpha), 3 (orb 1 beta)

# Fermionic operator: a†_1 a†_3 a_2 a_0 (create in 1,3, annihilate from 0,2)
print("\n📊 Creating double excitation operator")
print("  Occupied: [0, 2] (orbital 0, both spins)")
print("  Virtual:  [1, 3] (orbital 1, both spins)")

# Create fermionic excitation operator
# T = a†_1 a†_3 a_2 a_0 - a†_0 a†_2 a_3 a_1 (anti-hermitian part for UCC)
excitation_op = FermionicOp({
    "+_1 +_3 -_2 -_0": 1.0,  # Excitation
    "+_0 +_2 -_3 -_1": -1.0,  # De-excitation (Hermitian conjugate)
}, num_spin_orbitals=4)

print(f"\nFermionic operator created:")
print(f"  {excitation_op}")

# Map to qubits using Jordan-Wigner
print(f"\n📊 Mapping to qubits (Jordan-Wigner)")
jw_mapper = JordanWignerMapper()

try:
    qubit_op = jw_mapper.map(excitation_op)
    print(f"  Qubit operator:")
    print(f"  {qubit_op}")
    print(f"  Number of Pauli terms: {len(qubit_op)}")

    # Check if this can be converted to a circuit
    from qiskit.synthesis import PauliEvolutionSynthesis
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit import QuantumCircuit

    # Create evolution gate
    print(f"\n📊 Creating quantum circuit from operator")
    theta = 0.1  # Small parameter

    evo_gate = PauliEvolutionGate(qubit_op, time=theta)
    qc = QuantumCircuit(4)
    qc.append(evo_gate, range(4))

    print(f"  Circuit depth: {qc.depth()}")
    print(f"  Number of gates: {qc.size()}")

    # Test on HF state
    from qiskit.quantum_info import Statevector

    # Prepare HF state |0101⟩
    qc_full = QuantumCircuit(4)
    qc_full.x(0)  # Qubit 0 = 1
    qc_full.x(2)  # Qubit 2 = 1
    # State is now |0101⟩

    # Apply double excitation
    qc_full.append(evo_gate, range(4))

    # Get statevector
    statevector = Statevector.from_instruction(qc_full)
    psi = statevector.data

    print(f"\n📊 Result after applying double excitation (theta={theta}):")
    for i, amp in enumerate(psi):
        if abs(amp) > 0.001:
            print(f"  |{i:04b}⟩: {amp:.6f}")

    # Check particle number
    n_particles = sum(abs(psi[i])**2 * bin(i).count('1') for i in range(len(psi)))
    print(f"\n  Total particles: {n_particles:.6f} (should be 2.0)")

    if abs(n_particles - 2.0) < 0.001:
        print(f"  ✓ Particle number conserved!")
    else:
        print(f"  ✗ Particle number NOT conserved!")

    # Check if we got amplitude on |1010⟩
    amp_1010 = abs(psi[0b1010])  # |1010⟩
    if amp_1010 > 0.01:
        print(f"\n  ✓ Created amplitude on |1010⟩: {amp_1010:.6f}")
    else:
        print(f"\n  ✗ No significant amplitude on |1010⟩")

except Exception as e:
    print(f"  ✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

print(f"\n{'='*80}")
