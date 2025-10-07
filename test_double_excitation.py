"""
Test if UCC double excitation circuit actually connects |5⟩ and |10⟩.
"""

import numpy as np
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from qiskit.quantum_info import Statevector

# Create UCC ansatz
ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)
ansatz.build_circuit()

# Test with different angles for the double excitation parameter
# The double excitation is parameter index 4
test_angles = [0.0, 0.1, 0.5, 1.0, np.pi/4]

print("Testing double excitation circuit:")
print("Parameter 4 controls [0,2]→[1,3] double excitation")
print()

for theta in test_angles:
    # Set all singles to zero, only vary double
    params = np.array([0.0, 0.0, 0.0, 0.0, theta])
    ansatz.circuit.bind_parameters(params)
    
    qc = ansatz.circuit.to_qiskit()
    state = Statevector.from_instruction(qc)
    
    # Check amplitudes of states 5 and 10
    amp_5 = state.data[5]
    amp_10 = state.data[10]
    
    print(f"θ = {theta:.3f}:")
    print(f"  |5⟩ amplitude: {abs(amp_5):.6f}")
    print(f"  |10⟩ amplitude: {abs(amp_10):.6f}")
    print(f"  Total probability: {abs(amp_5)**2 + abs(amp_10)**2:.6f}")

print("\nIf working correctly:")
print("  - At θ=0: should be pure |5⟩ (HF state)")
print("  - At θ>0: should have mix of |5⟩ and |10⟩")
print("  - At θ=π/4: should have significant |10⟩ component")
