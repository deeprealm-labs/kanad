"""
Test if parameter binding persists across multiple calls.
"""

import numpy as np
from kanad.ansatze.ucc_ansatz import UCCAnsatz

ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)
ansatz.build_circuit()

print("First binding:")
ansatz.circuit.bind_parameters(np.array([0.1, 0.2, 0.3, 0.4, 0.5]))
qc1 = ansatz.circuit.to_qiskit()
print(f"  Parameters in circuit: {qc1.num_parameters}")
print(f"  First RY gate param: {qc1.data[4][0].params}")

print("\nSecond binding (different values):")
ansatz.circuit.bind_parameters(np.array([1.0, 2.0, 3.0, 4.0, 5.0]))
qc2 = ansatz.circuit.to_qiskit()
print(f"  Parameters in circuit: {qc2.num_parameters}")
print(f"  First RY gate param: {qc2.data[4][0].params}")

print("\nThird call to to_qiskit() without rebinding:")
qc3 = ansatz.circuit.to_qiskit()
print(f"  Parameters in circuit: {qc3.num_parameters}")
print(f"  First RY gate param: {qc3.data[4][0].params}")

print("\nConclusion:")
if qc3.num_parameters == 0:
    print("✓ Parameter binding persists!")
else:
    print("✗ Parameter binding is reset!")
