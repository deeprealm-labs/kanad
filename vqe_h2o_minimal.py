#!/usr/bin/env python3
"""
H2O VQE Calculation - Minimal Implementation
Works around qiskit-nature compatibility issues
"""

import numpy as np
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import SLSQP
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import SparsePauliOp
from qiskit import QuantumCircuit

print("="*60)
print("H2O VQE Calculation - Minimal Implementation")
print("="*60)

# Define the H2O molecule geometry
print("\n1. Setting up H2O molecule...")
driver = PySCFDriver(
    atom="O 0.0 0.0 0.0; H 0.757 0.586 0.0; H -0.757 0.586 0.0",
    basis="sto3g",
    charge=0,
    spin=0,
    unit=DistanceUnit.ANGSTROM,
)

# Run the driver to get the electronic structure problem
print("   Running PySCF driver...")
es_problem = driver.run()
print(f"   ✓ H2O molecule ready ({es_problem.num_spatial_orbitals} spatial orbitals)")

# Set up the qubit mapper (Jordan-Wigner)
print("\n2. Setting up qubit mapper...")
mapper = JordanWignerMapper()
print("   ✓ Jordan-Wigner mapper ready")

# Get the second quantized Hamiltonian
print("\n3. Converting Hamiltonian to qubit form...")
second_q_hamiltonian = es_problem.hamiltonian.second_q_op()
qubit_hamiltonian = mapper.map(second_q_hamiltonian)
print(f"   ✓ Hamiltonian mapped to {qubit_hamiltonian.num_qubits} qubits")

# Create the ansatz circuit (UCCSD with Hartree-Fock initial state)
print("\n4. Creating UCCSD ansatz...")
ansatz = UCCSD(
    es_problem.num_spatial_orbitals,
    es_problem.num_particles,
    mapper,
    initial_state=HartreeFock(
        es_problem.num_spatial_orbitals,
        es_problem.num_particles,
        mapper,
    ),
)
print(f"   ✓ UCCSD ansatz ready ({ansatz.num_parameters} parameters)")

# Set up the VQE solver
print("\n5. Setting up VQE solver...")
estimator = StatevectorEstimator()
optimizer = SLSQP()
vqe_solver = VQE(estimator, ansatz, optimizer)
vqe_solver.initial_point = [0.0] * ansatz.num_parameters
print("   ✓ VQE solver ready")

# Solve for the ground state
print("\n6. Running VQE calculation...")
print("   This may take a few minutes...")
result = vqe_solver.compute_minimum_eigenvalue(qubit_hamiltonian)

# Print the results
print("\n" + "="*60)
print("RESULTS")
print("="*60)
print(f"Ground state energy: {result.eigenvalue:.6f} Hartree")
print(f"Ground state energy: {result.eigenvalue * 27.211386245988:.6f} eV")
print(f"Number of qubits: {qubit_hamiltonian.num_qubits}")
print(f"Number of electrons: {es_problem.num_particles[0] + es_problem.num_particles[1]}")
print(f"Optimization iterations: {result.cost_function_evals}")
print("="*60)
