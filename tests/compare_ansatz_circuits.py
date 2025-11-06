#!/usr/bin/env python3
"""
Compare Ansatz Circuits from Bond API vs Direct API
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz

print("="*80)
print("ANSATZ CIRCUIT COMPARISON")
print("="*80)

# Test parameters (will use zeros and random to check)
n_qubits = 4
n_electrons = 2

# Create ansatz
ansatz = CovalentGovernanceAnsatz(n_qubits=n_qubits, n_electrons=n_electrons)
print(f"\nAnsatz: CovalentGovernanceAnsatz")
print(f"  Qubits: {n_qubits}")
print(f"  Electrons: {n_electrons}")
print(f"  Parameters: {ansatz.n_parameters}")

# Build circuit
ansatz.build_circuit()
print(f"  Circuit built")

# Convert to Qiskit to inspect
circuit = ansatz.circuit.to_qiskit()
print(f"\nCircuit structure:")
print(f"  Qubits: {circuit.num_qubits}")
print(f"  Depth: {circuit.depth()}")
print(f"  Gates: {circuit.count_ops()}")
print(f"  Parameters: {circuit.num_parameters}")

# Test with zero parameters (should give HF state)
print("\n" + "-"*80)
print("TEST 1: Zero Parameters (HF State)")
print("-"*80)

from qiskit.quantum_info import Statevector

params_zero = np.zeros(ansatz.n_parameters)
ansatz.circuit.bind_parameters(params_zero)
qiskit_zero = ansatz.circuit.to_qiskit()

if qiskit_zero.num_parameters > 0:
    param_dict = {qiskit_zero.parameters[i]: params_zero[i] for i in range(len(params_zero))}
    qiskit_zero = qiskit_zero.assign_parameters(param_dict)

sv_zero = Statevector.from_instruction(qiskit_zero)
amplitudes_zero = np.abs(sv_zero.data)

print(f"Statevector amplitudes:")
for i, amp in enumerate(amplitudes_zero):
    if amp > 1e-10:
        bitstring = format(i, f'0{n_qubits}b')
        print(f"  |{bitstring}⟩: {amp:.8f}")

# Check if it's HF state (|1100⟩ for H2 with Jordan-Wigner)
hf_state_jw = 0b1100  # qubits 2 and 3 occupied
if amplitudes_zero[hf_state_jw] > 0.99:
    print(f"\n✓ Correctly prepared HF state |1100⟩")
else:
    print(f"\n❌ NOT HF state!")
    print(f"   Expected |1100⟩ to have amplitude ~1.0")
    print(f"   Got: {amplitudes_zero[hf_state_jw]:.8f}")

# Test with random parameters
print("\n" + "-"*80)
print("TEST 2: Random Parameters")
print("-"*80)

# Rebuild circuit
ansatz.build_circuit()

params_random = np.random.uniform(-0.1, 0.1, size=ansatz.n_parameters)
ansatz.circuit.bind_parameters(params_random)
qiskit_random = ansatz.circuit.to_qiskit()

if qiskit_random.num_parameters > 0:
    param_dict = {qiskit_random.parameters[i]: params_random[i] for i in range(len(params_random))}
    qiskit_random = qiskit_random.assign_parameters(param_dict)

sv_random = Statevector.from_instruction(qiskit_random)
amplitudes_random = np.abs(sv_random.data)

print(f"Statevector amplitudes (top 5):")
sorted_indices = np.argsort(amplitudes_random)[::-1]
for idx in sorted_indices[:5]:
    amp = amplitudes_random[idx]
    if amp > 1e-10:
        bitstring = format(idx, f'0{n_qubits}b')
        print(f"  |{bitstring}⟩: {amp:.8f}")

# Compute energy with both Hamiltonians
print("\n" + "="*80)
print("ENERGY COMPUTATION TEST")
print("="*80)

from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner
from pyscf import gto

# Create Hamiltonian
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g', charge=0, spin=0)
mol.build()

mf = mol.RHF().run(verbose=0)
mo_coeff = mf.mo_coeff

h1e = mol.intor('int1e_kin') + mol.intor('int1e_nuc')
h_mo = mo_coeff.T @ h1e @ mo_coeff

eri_ao = mol.intor('int2e')
eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', mo_coeff, mo_coeff, eri_ao, mo_coeff, mo_coeff, optimize=True)

hamiltonian = openfermion_jordan_wigner(
    h_mo=h_mo,
    eri_mo=eri_mo,
    nuclear_repulsion=mol.energy_nuc(),
    n_electrons=mol.nelectron
)

print(f"\nHamiltonian: {len(hamiltonian)} Pauli terms")

# Energy with zero parameters (should be near HF)
energy_zero = sv_zero.expectation_value(hamiltonian).real
hf_energy = mf.e_tot

print(f"\nEnergy with zero parameters:")
print(f"  E(θ=0):    {energy_zero:.8f} Ha")
print(f"  HF energy: {hf_energy:.8f} Ha")
print(f"  Difference: {abs(energy_zero - hf_energy):.8f} Ha")

if abs(energy_zero - hf_energy) < 0.001:
    print(f"  ✓ Zero parameters give HF energy (correct)")
else:
    print(f"  ❌ Zero parameters DON'T give HF energy!")
    print(f"     This could be the bug!")

# Energy with random parameters (should be lower than HF ideally)
energy_random = sv_random.expectation_value(hamiltonian).real

print(f"\nEnergy with random parameters:")
print(f"  E(θ_random): {energy_random:.8f} Ha")
print(f"  HF energy:   {hf_energy:.8f} Ha")
print(f"  Δ from HF:   {energy_random - hf_energy:.8f} Ha")

if energy_random < hf_energy:
    print(f"  ✓ Random parameters can give lower energy (good ansatz)")
else:
    print(f"  ⚠️  Random parameters gave higher energy")
    print(f"     (Not necessarily bad - just unlucky initialization)")

print("\n" + "="*80)
