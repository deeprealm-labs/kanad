#!/usr/bin/env python3
"""
Complex Molecule Test - Check if VQE is actually doing quantum simulation

This will test:
1. Larger molecule (LiH - 4 atoms)
2. Full VQE optimization
3. Timing to see if it's realistic
4. State vector size to confirm no shortcuts
"""

import time
import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.solvers.vqe_solver import VQESolver

print("="*80)
print("COMPLEX MOLECULE TEST - C-O (Carbon Monoxide)")
print("="*80)
print()

# ============================================================================
# Test 1: Create C-O bond (Carbon Monoxide)
# ============================================================================
print("1. Creating C-O triple bond...")
print("-"*80)

start_time = time.time()
bond = BondFactory.create_bond('C', 'O', distance=1.13, bond_type='covalent')
creation_time = time.time() - start_time

print(f"✓ Bond created in {creation_time:.3f}s")
print(f"  Molecule: LiH")
print(f"  Electrons: {bond.molecule.n_electrons}")
print(f"  Qubits: {bond.representation.n_qubits}")
print(f"  State vector size: 2^{bond.representation.n_qubits} = {2**bond.representation.n_qubits}")
print()

# Check if this is reasonable
n_qubits = bond.representation.n_qubits
state_size = 2**n_qubits

if state_size < 64:
    print(f"⚠️  WARNING: State vector only has {state_size} entries!")
    print(f"   This is VERY small for quantum simulation.")
    print(f"   Expected at least 256 (8 qubits) for a real molecule.")
print()

# ============================================================================
# Test 2: Build UCC Ansatz
# ============================================================================
print("2. Building UCCSD ansatz...")
print("-"*80)

start_time = time.time()
ansatz = UCCAnsatz(n_qubits, bond.molecule.n_electrons)
circuit = ansatz.build_circuit()
ansatz_time = time.time() - start_time

print(f"✓ Ansatz built in {ansatz_time:.3f}s")
print(f"  Parameters: {len(ansatz.circuit.parameters)}")
print(f"  Gates: {len(ansatz.circuit.gates)}")
print(f"  Excitations: {len(ansatz.excitations)}")
print()

# Check circuit complexity
if len(ansatz.circuit.gates) < 10:
    print(f"⚠️  WARNING: Only {len(ansatz.circuit.gates)} gates!")
    print(f"   This is VERY few for a UCC ansatz.")
    print(f"   Expected dozens or hundreds of gates.")
print()

# ============================================================================
# Test 3: Run VQE and TIME IT
# ============================================================================
print("3. Running VQE optimization (THIS IS THE KEY TEST)...")
print("-"*80)

vqe = VQESolver(
    hamiltonian=bond.hamiltonian,
    ansatz=ansatz,
    mapper=bond.mapper,
    optimizer='BFGS',
    max_iterations=50  # Limit iterations for timing
)

print(f"Starting VQE with:")
print(f"  State vector: {state_size} complex numbers")
print(f"  Gates per iteration: {len(ansatz.circuit.gates)}")
print(f"  Expected cost: O({len(ansatz.circuit.gates)} * {state_size}) per iteration")
print()

# Callback to track progress
iteration_times = []
last_time = time.time()

def callback(iter, energy, params):
    global last_time
    current_time = time.time()
    iter_time = current_time - last_time
    iteration_times.append(iter_time)
    last_time = current_time
    if iter % 5 == 0 or iter == 1:
        print(f"  Iteration {iter:3d}: E = {energy:12.8f} Ha  (took {iter_time:.4f}s)")

print("Running VQE...")
start_time = time.time()
last_time = start_time

result = vqe.solve(callback=callback)

total_time = time.time() - start_time

print()
print(f"✓ VQE completed in {total_time:.3f}s")
print(f"  Final energy: {result['energy']:.8f} Ha")
print(f"  Iterations: {result['iterations']}")
print(f"  Converged: {result['converged']}")
print(f"  Avg time per iteration: {total_time/result['iterations']:.4f}s")
print()

# ============================================================================
# Test 4: Analysis - Is this realistic?
# ============================================================================
print("4. Realism Analysis")
print("-"*80)

print(f"State vector size: {state_size}")
print(f"Total VQE time: {total_time:.3f}s")
print(f"Time per iteration: {np.mean(iteration_times):.4f}s")
print()

# Estimate expected time
# For N qubits with G gates:
# - State vector: 2^N complex numbers (16 bytes each)
# - Per gate: matrix-vector multiply O(2^N * 2^N) = O(4^N)
# - Total per iteration: O(G * 4^N)

expected_ops_per_iter = len(ansatz.circuit.gates) * (4**n_qubits)
print(f"Estimated operations per iteration: {expected_ops_per_iter:,}")
print()

# For 4 qubits (16-dim state):
# - 16 x 16 matrix multiply ≈ 4096 ops
# - With 10 gates ≈ 40,960 ops
# - Modern CPU: ~10 GFLOPS → ~0.004ms per iteration
# - Python overhead: ~100x slower → ~0.4ms per iteration

if n_qubits == 4:
    print("⚠️  ANALYSIS FOR 4-QUBIT SYSTEM:")
    print(f"  - State vector: 16 complex numbers (256 bytes)")
    print(f"  - Matrix multiply: 16x16 = 256 multiply-adds per gate")
    print(f"  - Total ops: {len(ansatz.circuit.gates)} gates * 256 ops = {len(ansatz.circuit.gates)*256:,}")
    print(f"  - Expected time (Python): ~0.001s per iteration")
    print(f"  - Actual time: {np.mean(iteration_times):.4f}s per iteration")
    print()

    if np.mean(iteration_times) < 0.0001:
        print("  🚨 SUSPICIOUSLY FAST!")
        print("     Time is < 0.1ms per iteration - might be mocked or skipped")
    elif np.mean(iteration_times) < 0.001:
        print("  ⚠️  VERY FAST!")
        print("     Might be using optimized shortcuts")
    else:
        print("  ✅ REALISTIC timing for 4-qubit system")
elif n_qubits >= 8:
    print("✅ GOOD: Testing with 8+ qubits")
    print(f"  - State vector: {state_size} complex numbers ({state_size*16} bytes)")
    print(f"  - This is a meaningful quantum simulation")

    # For 8 qubits: 256x256 matrices → 65,536 ops per gate
    # Expected: ~0.01s per iteration in Python
    expected_time = 0.01
    if np.mean(iteration_times) < expected_time / 10:
        print(f"  🚨 TOO FAST! Expected ~{expected_time:.3f}s, got {np.mean(iteration_times):.4f}s")
        print("     Might be using shortcuts or mocking")
    else:
        print(f"  ✅ Timing is realistic ({np.mean(iteration_times):.4f}s per iteration)")

print()

# ============================================================================
# Test 5: Verify actual computation is happening
# ============================================================================
print("5. Verification - Are calculations real?")
print("-"*80)

# Test that different parameters give different energies
print("Testing if VQE energy changes with parameters...")

params1 = np.zeros(len(ansatz.circuit.parameters))
params2 = np.ones(len(ansatz.circuit.parameters)) * 0.1

energy1 = vqe.compute_energy(params1)
energy2 = vqe.compute_energy(params2)

print(f"  Energy(0): {energy1:.8f} Ha")
print(f"  Energy(0.1): {energy2:.8f} Ha")
print(f"  Difference: {abs(energy2 - energy1):.8f} Ha")

if abs(energy2 - energy1) < 1e-10:
    print("  🚨 PROBLEM: Energies are identical!")
    print("     VQE might not be actually varying parameters")
else:
    print("  ✅ Energies change with parameters (good!)")

print()

# Test that energy expectation is not just returning a constant
energies_at_random = []
for _ in range(5):
    random_params = np.random.randn(len(ansatz.circuit.parameters)) * 0.5
    e = vqe.compute_energy(random_params)
    energies_at_random.append(e)

print(f"5 random parameter energies:")
for i, e in enumerate(energies_at_random):
    print(f"  {i+1}. {e:.8f} Ha")

std_dev = np.std(energies_at_random)
print(f"  Standard deviation: {std_dev:.8f} Ha")

if std_dev < 1e-8:
    print("  🚨 PROBLEM: All energies are nearly identical!")
    print("     VQE might be returning a constant")
else:
    print("  ✅ Energy varies significantly (good!)")

print()

# ============================================================================
# Test 6: Check state vector is actually evolving
# ============================================================================
print("6. State Vector Evolution Check")
print("-"*80)

# Get initial state
vqe.circuit.bind_parameters(params1)
state1 = vqe._simulate_circuit()

# Get state with different parameters
vqe.circuit.bind_parameters(params2)
state2 = vqe._simulate_circuit()

# Compute fidelity
fidelity = abs(np.vdot(state1, state2))**2
print(f"  Initial state norm: {np.linalg.norm(state1):.8f}")
print(f"  Second state norm: {np.linalg.norm(state2):.8f}")
print(f"  Fidelity: {fidelity:.8f}")
print(f"  Overlap: {abs(np.vdot(state1, state2)):.8f}")

if abs(np.linalg.norm(state1) - 1.0) > 1e-6:
    print("  🚨 PROBLEM: State is not normalized!")
elif fidelity > 0.999:
    print("  🚨 PROBLEM: States are nearly identical!")
    print("     Circuit might not be applying gates")
else:
    print("  ✅ States evolve properly")

print()

# ============================================================================
# FINAL VERDICT
# ============================================================================
print("="*80)
print("FINAL VERDICT")
print("="*80)

issues = []
goods = []

# Check 1: State vector size
if state_size < 64:
    issues.append(f"State vector too small ({state_size} entries) - molecule too simple")
else:
    goods.append(f"State vector size reasonable ({state_size} entries)")

# Check 2: Timing
if n_qubits == 4 and np.mean(iteration_times) < 0.0001:
    issues.append(f"Suspiciously fast ({np.mean(iteration_times):.6f}s/iter)")
elif n_qubits >= 8 and np.mean(iteration_times) < 0.001:
    issues.append(f"Too fast for {n_qubits} qubits ({np.mean(iteration_times):.6f}s/iter)")
else:
    goods.append(f"Timing seems realistic ({np.mean(iteration_times):.4f}s/iter)")

# Check 3: Energy variation
if std_dev < 1e-8:
    issues.append("Energy doesn't vary with parameters")
else:
    goods.append(f"Energy varies properly (σ={std_dev:.6f})")

# Check 4: State evolution
if fidelity > 0.999:
    issues.append("State vectors don't evolve")
else:
    goods.append(f"State vectors evolve (fidelity={fidelity:.4f})")

# Check 5: Convergence
if result['converged']:
    goods.append("VQE converged to solution")
else:
    issues.append("VQE did not converge")

print()
print("✅ WORKING CORRECTLY:")
for g in goods:
    print(f"  ✓ {g}")

print()
if issues:
    print("🚨 POTENTIAL ISSUES:")
    for i in issues:
        print(f"  ✗ {i}")
else:
    print("🎉 NO ISSUES FOUND - VQE appears to be working correctly!")

print()
print("="*80)
print(f"CONCLUSION: {'SUSPICIOUS' if issues else 'VQE IS REAL'}")
print("="*80)
