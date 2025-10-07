"""
Test parameter binding in UCC ansatz.
"""

import numpy as np
from kanad.ansatze.ucc_ansatz import UCCAnsatz

# Create UCC ansatz
ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)
print(f"UCC Ansatz: {ansatz.n_parameters} parameters")

# Build circuit
circuit = ansatz.build_circuit()
print(f"\nBuilt circuit with {len(ansatz.circuit.parameters)} parameters")

# Print parameter info
for i, param in enumerate(ansatz.circuit.parameters):
    print(f"  Param {i}: {param.name}, value={param.value}")

# Bind some test parameters
test_params = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
print(f"\nBinding parameters: {test_params}")
ansatz.circuit.bind_parameters(test_params)

# Check if bound
print(f"\nAfter binding:")
for i, param in enumerate(ansatz.circuit.parameters):
    print(f"  Param {i}: {param.name}, value={param.value}")

# Convert to Qiskit
print(f"\nConverting to Qiskit...")
qiskit_circuit = ansatz.circuit.to_qiskit()
print(f"Qiskit circuit has {qiskit_circuit.num_parameters} parameters")

if qiskit_circuit.num_parameters > 0:
    print(f"\nQiskit circuit parameters:")
    for i, param in enumerate(qiskit_circuit.parameters):
        print(f"  Param {i}: {param.name}")
    print("\n✗ PROBLEM: Circuit still has symbolic parameters after binding!")
else:
    print("\n✓ SUCCESS: Circuit has no parameters (all bound)")

# Try to draw first few gates
print(f"\nFirst 10 gates in Qiskit circuit:")
for i, (instr, qargs, cargs) in enumerate(qiskit_circuit.data[:10]):
    params_str = ""
    if hasattr(instr, 'params') and instr.params:
        params_str = f" params={instr.params}"
    print(f"  {i}: {instr.name} on {[q._index for q in qargs]}{params_str}")
