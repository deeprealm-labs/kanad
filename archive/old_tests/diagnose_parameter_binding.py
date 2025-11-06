#!/usr/bin/env python3
"""
Diagnose Parameter Binding Issue
Check if parameters are properly bound when converting circuit
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz

print("="*80)
print("PARAMETER BINDING DIAGNOSTIC")
print("="*80)

# Create simple ansatz
ansatz = CovalentGovernanceAnsatz(n_qubits=4, n_electrons=2)
print(f"\nAnsatz: {ansatz.n_parameters} parameters")

# Build circuit
ansatz.build_circuit()
print(f"Circuit built")

# Test 1: Zero parameters
print("\n" + "-"*80)
print("TEST 1: Zero parameters (HF state)")
print("-"*80)

params_zero = np.zeros(ansatz.n_parameters)
print(f"Parameters: all zeros")

# Bind parameters
ansatz.circuit.bind_parameters(params_zero)
print("✓ Bound parameters to custom circuit")

# Convert to Qiskit
qiskit_circuit = ansatz.circuit.to_qiskit()
print(f"✓ Converted to Qiskit circuit")
print(f"  Qiskit circuit has {qiskit_circuit.num_parameters} parameters")

if qiskit_circuit.num_parameters > 0:
    print("  ⚠️  Circuit still has unbound parameters!")
    print(f"  Parameters: {[str(p) for p in list(qiskit_circuit.parameters)[:5]]}...")
else:
    print("  ✓ All parameters bound")

# Test 2: Random parameters
print("\n" + "-"*80)
print("TEST 2: Random parameters")
print("-"*80)

# Rebuild circuit
ansatz.build_circuit()

params_random = np.random.uniform(-0.1, 0.1, size=ansatz.n_parameters)
print(f"Parameters: random in [-0.1, 0.1]")
print(f"  First 5: {params_random[:5]}")

# Bind parameters
ansatz.circuit.bind_parameters(params_random)
print("✓ Bound parameters to custom circuit")

# Check if binding worked on custom circuit
print(f"Custom circuit parameters after binding:")
if hasattr(ansatz.circuit, 'parameters') and len(ansatz.circuit.parameters) > 0:
    bound_count = sum(1 for p in ansatz.circuit.parameters if not p._is_symbolic)
    print(f"  Bound: {bound_count}/{len(ansatz.circuit.parameters)}")
    if bound_count == len(ansatz.circuit.parameters):
        print(f"  ✓ All custom circuit parameters bound")
    else:
        print(f"  ⚠️  Some parameters still symbolic!")

# Convert to Qiskit
qiskit_circuit = ansatz.circuit.to_qiskit()
print(f"✓ Converted to Qiskit circuit")
print(f"  Qiskit circuit has {qiskit_circuit.num_parameters} parameters")

if qiskit_circuit.num_parameters > 0:
    print("  ⚠️  ❌ CIRCUIT STILL HAS UNBOUND PARAMETERS!")
    print("  This is the BUG! Parameters not transferred during conversion.")

    # Try binding again on Qiskit circuit
    param_dict = {qiskit_circuit.parameters[i]: params_random[i]
                  for i in range(len(params_random))}
    bound_circuit = qiskit_circuit.assign_parameters(param_dict)
    print(f"\n  After manual binding:")
    print(f"    Parameters: {bound_circuit.num_parameters}")
    if bound_circuit.num_parameters == 0:
        print(f"    ✓ Manual binding works!")
else:
    print("  ✓ All parameters bound during conversion")

# Test 3: Check if statevector changes
print("\n" + "-"*80)
print("TEST 3: Check if statevector changes with parameters")
print("-"*80)

from qiskit.quantum_info import Statevector

# Rebuild circuit
ansatz.build_circuit()

# Test with zero parameters
params_zero = np.zeros(ansatz.n_parameters)
ansatz.circuit.bind_parameters(params_zero)
qiskit_zero = ansatz.circuit.to_qiskit()

if qiskit_zero.num_parameters > 0:
    param_dict = {qiskit_zero.parameters[i]: params_zero[i] for i in range(len(params_zero))}
    qiskit_zero = qiskit_zero.assign_parameters(param_dict)

sv_zero = Statevector.from_instruction(qiskit_zero)
print(f"Statevector with zero params:")
print(f"  Amplitudes: {np.abs(sv_zero.data)[:8]}")

# Rebuild circuit
ansatz.build_circuit()

# Test with random parameters
params_random = np.random.uniform(-0.5, 0.5, size=ansatz.n_parameters)
ansatz.circuit.bind_parameters(params_random)
qiskit_random = ansatz.circuit.to_qiskit()

if qiskit_random.num_parameters > 0:
    param_dict = {qiskit_random.parameters[i]: params_random[i] for i in range(len(params_random))}
    qiskit_random = qiskit_random.assign_parameters(param_dict)

sv_random = Statevector.from_instruction(qiskit_random)
print(f"\nStatevector with random params:")
print(f"  Amplitudes: {np.abs(sv_random.data)[:8]}")

# Check if they're different
diff = np.linalg.norm(sv_zero.data - sv_random.data)
print(f"\nStatevector difference: {diff:.8f}")

if diff < 1e-6:
    print("  ❌ CRITICAL BUG: Statevectors identical!")
    print("  Parameters are NOT affecting the circuit!")
else:
    print("  ✓ Statevectors are different (parameters working)")

print("\n" + "="*80)
print("DIAGNOSIS")
print("="*80)

if qiskit_circuit.num_parameters > 0:
    print("\n❌ BUG FOUND: Parameter binding not working properly")
    print("\nThe issue:")
    print("  1. Custom circuit.bind_parameters() is called")
    print("  2. But to_qiskit() creates NEW parameter objects")
    print("  3. Bound values are NOT transferred to new parameters")
    print("  4. Result: Circuit always uses symbolic parameters")
    print("\nFix needed:")
    print("  - Modify to_qiskit() to transfer bound parameter values")
    print("  - OR skip bind_parameters() and only bind on Qiskit circuit")
else:
    print("\n✓ Parameter binding working correctly")

print("\n" + "="*80)
