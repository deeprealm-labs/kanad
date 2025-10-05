#!/usr/bin/env python3
"""
H2O VQE Calculation using Qiskit Nature
Fixed version that works with current Qiskit versions
"""

import numpy as np
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import SLSQP
from qiskit.primitives import StatevectorEstimator

print("="*60)
print("H2O VQE Calculation using Qiskit Nature")
print("="*60)

# Define the H2O molecule geometry
print("\n1. Setting up H2O molecule...")
driver = PySCFDriver(
    atom="""
    O 0.0 0.0 0.0;
    H 0.757 0.586 0.0;
    H -0.757 0.586 0.0
    """,
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

# Create the ansatz circuit (UCCSD with Hartree-Fock initial state)
print("\n3. Creating UCCSD ansatz...")
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
print("\n4. Setting up VQE solver...")
estimator = StatevectorEstimator()
optimizer = SLSQP()
vqe_solver = VQE(estimator, ansatz, optimizer)
vqe_solver.initial_point = [0.0] * ansatz.num_parameters
print("   ✓ VQE solver ready")

# Create the ground state eigensolver
print("\n5. Creating ground state eigensolver...")
calc = GroundStateEigensolver(mapper, vqe_solver)
print("   ✓ Ground state eigensolver ready")

# Solve for the ground state
print("\n6. Running VQE calculation...")
print("   This may take a few minutes...")
result = calc.solve(es_problem)

# Print the results
print("\n" + "="*60)
print("RESULTS")
print("="*60)
print(f"Ground state energy: {result.groundenergy:.6f} Hartree")
print(f"Ground state energy: {result.groundenergy * 27.211386245988:.6f} eV")
print(f"Number of qubits: {es_problem.num_spatial_orbitals * 2}")
print(f"Number of electrons: {es_problem.num_particles[0] + es_problem.num_particles[1]}")
print("="*60)

