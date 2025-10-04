#!/usr/bin/env python3
"""
Quick H2O VQE Performance Test

Fast benchmark focusing on:
- Setup time
- VQE convergence
- Energy accuracy
- Resource usage (qubits, gates)
"""

import numpy as np
import time
from pyscf import gto, scf as pyscf_scf
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
from qiskit.primitives import StatevectorEstimator
from qiskit_algorithms.minimum_eigensolvers import VQE
from qiskit_algorithms.optimizers import SLSQP

print("="*70)
print("H2O VQE QUICK BENCHMARK")
print("="*70)

# ============================================================================
# STEP 1: Hartree-Fock Reference
# ============================================================================
print("\n[1/4] Running Hartree-Fock Reference...")
start_hf = time.time()

mol = gto.Mole()
mol.atom = "O 0.0 0.0 0.1173; H 0.0 0.7572 -0.4692; H 0.0 -0.7572 -0.4692"
mol.basis = 'sto-3g'
mol.build()

mf = pyscf_scf.RHF(mol)
hf_energy_ha = mf.kernel()
hf_time = time.time() - start_hf

# Convert to eV
hartree_to_ev = 27.211386245988
hf_energy_ev = hf_energy_ha * hartree_to_ev

print(f"  ✓ HF Energy: {hf_energy_ha:.6f} Ha ({hf_energy_ev:.4f} eV)")
print(f"  ✓ Time: {hf_time:.3f}s")
print(f"  ✓ Converged: {mf.converged}")

# ============================================================================
# STEP 2: Setup Quantum Problem
# ============================================================================
print("\n[2/4] Setting up Quantum Problem...")
start_setup = time.time()

driver = PySCFDriver(
    atom="O 0.0 0.0 0.1173; H 0.0 0.7572 -0.4692; H 0.0 -0.7572 -0.4692",
    basis='sto-3g',
    charge=0,
    spin=0
)

problem = driver.run()
hamiltonian = problem.hamiltonian.second_q_op()
num_particles = problem.num_particles
num_spatial_orbitals = problem.num_spatial_orbitals

setup_time = time.time() - start_setup

print(f"  ✓ Spatial orbitals: {num_spatial_orbitals}")
print(f"  ✓ Electrons: {sum(num_particles)} (α={num_particles[0]}, β={num_particles[1]})")
print(f"  ✓ Setup time: {setup_time:.3f}s")

# ============================================================================
# STEP 3: Build VQE Circuit
# ============================================================================
print("\n[3/4] Building VQE Circuit...")
start_circuit = time.time()

mapper = JordanWignerMapper()
qubit_op = mapper.map(hamiltonian)
num_qubits = qubit_op.num_qubits

# Initial state (Hartree-Fock)
init_state = HartreeFock(
    num_spatial_orbitals=num_spatial_orbitals,
    num_particles=num_particles,
    qubit_mapper=mapper
)

# UCCSD ansatz
ansatz = UCCSD(
    num_spatial_orbitals=num_spatial_orbitals,
    num_particles=num_particles,
    qubit_mapper=mapper,
    initial_state=init_state
)

num_parameters = ansatz.num_parameters
circuit_time = time.time() - start_circuit

print(f"  ✓ Qubits: {num_qubits}")
print(f"  ✓ Variational parameters: {num_parameters}")
print(f"  ✓ Circuit depth: {ansatz.decompose().depth()}")
print(f"  ✓ Circuit build time: {circuit_time:.3f}s")

# ============================================================================
# STEP 4: Run VQE Optimization
# ============================================================================
print("\n[4/4] Running VQE Optimization...")
print("  (This may take 30-60 seconds...)")
start_vqe = time.time()

optimizer = SLSQP(maxiter=100)  # Limit iterations for speed
estimator = StatevectorEstimator()
vqe = VQE(estimator, ansatz, optimizer)

result = vqe.compute_minimum_eigenvalue(qubit_op)

vqe_time = time.time() - start_vqe
vqe_energy_ha = result.eigenvalue.real
vqe_energy_ev = vqe_energy_ha * hartree_to_ev

print(f"  ✓ VQE Energy: {vqe_energy_ha:.6f} Ha ({vqe_energy_ev:.4f} eV)")
print(f"  ✓ Optimizer evaluations: {result.optimizer_evals}")
print(f"  ✓ VQE time: {vqe_time:.3f}s")

# ============================================================================
# RESULTS SUMMARY
# ============================================================================
print("\n" + "="*70)
print("RESULTS SUMMARY")
print("="*70)

# Calculate errors
error_ha = vqe_energy_ha - hf_energy_ha
error_mha = error_ha * 1000  # milliHartree
error_ev = vqe_energy_ev - hf_energy_ev

print(f"\nMolecule: H2O (STO-3G basis)")
print(f"  Orbitals: {num_spatial_orbitals}")
print(f"  Electrons: {sum(num_particles)}")
print(f"  Qubits: {num_qubits}")
print(f"  Parameters: {num_parameters}")

print(f"\nEnergies:")
print(f"  HF:  {hf_energy_ha:>12.6f} Ha  ({hf_energy_ev:>10.4f} eV)")
print(f"  VQE: {vqe_energy_ha:>12.6f} Ha  ({vqe_energy_ev:>10.4f} eV)")

print(f"\nAccuracy:")
print(f"  Error vs HF: {error_mha:>8.3f} mHa  ({error_ev:>8.3f} eV)")
print(f"  Relative error: {abs(error_ha/hf_energy_ha)*100:.4f}%")

print(f"\nPerformance:")
print(f"  HF time:     {hf_time:>8.3f}s")
print(f"  Setup time:  {setup_time:>8.3f}s")
print(f"  Circuit time:{circuit_time:>8.3f}s")
print(f"  VQE time:    {vqe_time:>8.3f}s")
print(f"  Total time:  {hf_time+setup_time+circuit_time+vqe_time:>8.3f}s")

print(f"\nOptimization:")
print(f"  Iterations: {result.optimizer_evals}")
print(f"  Time per iteration: {vqe_time/result.optimizer_evals:.3f}s")

# Resource scaling estimate
print(f"\nResource Scaling (for reference):")
print(f"  Circuit gates (approx): ~{num_parameters * 5} (5 gates/param)")
print(f"  Classical memory: ~{num_qubits * 8} bytes (1 complex number per qubit state)")
print(f"  Quantum state dimension: 2^{num_qubits} = {2**num_qubits}")

print("\n" + "="*70)
print("✓ BENCHMARK COMPLETE!")
print("="*70)

# Key Insights
print("\nKey Insights:")
if abs(error_mha) < 1.0:
    print("  ✓ VQE achieved chemical accuracy (<1 mHa error)")
else:
    print(f"  ⚠ VQE error {abs(error_mha):.3f} mHa (chemical accuracy is <1 mHa)")

if vqe_time < 60:
    print(f"  ✓ Fast convergence ({vqe_time:.1f}s for {result.optimizer_evals} iterations)")
else:
    print(f"  ⚠ Slow convergence ({vqe_time:.1f}s)")

print(f"  ℹ System size: {num_qubits} qubits, {num_parameters} parameters")
print(f"  ℹ H2O requires moderate quantum resources (manageable on simulators)")

print("\n")
