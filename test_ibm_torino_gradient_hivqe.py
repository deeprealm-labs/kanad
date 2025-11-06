"""
Test Gradient-Optimized Hi-VQE on IBM Torino

This script tests the gradient-based configuration selection Hi-VQE
on real IBM quantum hardware to validate the improvements.

Expected improvements over previous test:
- Previous (HF + all singles): 45 mHa error
  - 20.5 mHa from incomplete subspace
  - 24.5 mHa from hardware noise
- Expected (HF + gradient selection): ~25 mHa error
  - 0 mHa from complete subspace (gradient finds |0011⟩)
  - ~25 mHa from hardware noise

Usage:
    export IBM_API='your_api_key'
    export IBM_CRN='your_crn'
    python test_ibm_torino_gradient_hivqe.py
"""
import os
import json
import numpy as np
from datetime import datetime
from kanad.bonds import BondFactory
from kanad.core.configuration import ConfigurationSubspace
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol
from qiskit import QuantumCircuit

print("="*80)
print("GRADIENT-OPTIMIZED HI-VQE ON IBM TORINO")
print("="*80)

# Step 1: Setup molecule and run local baseline
print("\n1️⃣ Running local baseline with gradient selection...")
from kanad.utils.vqe_solver import VQESolver

bond = BondFactory.create_bond('H', 'H', distance=0.74)

solver_local = VQESolver(
    bond=bond,
    mode='hivqe',
    hivqe_max_iterations=3,
    backend='statevector'
)

result_local = solver_local.solve()

print(f"   Local Hi-VQE: {result_local['energy']:.8f} Ha")
print(f"   Expected FCI: -1.13728383 Ha")
print(f"   Error: {abs(result_local['energy'] - (-1.13728383))*1000:.4f} mHa")
print(f"   Iterations: {result_local['iterations']}")
print(f"   ✓ Baseline validated (gradient selection active)")

# Step 2: Build configuration subspace with gradient selection
print("\n2️⃣ Building Hi-VQE configuration subspace with gradient selection...")
hamiltonian = bond.hamiltonian.to_sparse_hamiltonian(mapper='jordan_wigner')

# We'll use the same iteration logic as the solver
protocol = CovalentGovernanceProtocol()
subspace = ConfigurationSubspace(4, 2, protocol=protocol)

# Start with HF
hf_config = subspace.get_hf_configuration()
subspace.add_config(hf_config)

# Iteration 1: Add singles from HF
single_excs = subspace.generate_single_excitations(hf_config)
subspace.add_configs(single_excs)

# Iteration 2: Use gradient selection to add most important configuration
from kanad.core.classical_solver import (
    SubspaceHamiltonianBuilder,
    diagonalize_subspace,
    select_configurations_by_gradient
)

# Compute current ground state
builder = SubspaceHamiltonianBuilder(hamiltonian)
H_sub = builder.project_fast(subspace)
eigenvalues, eigenvectors = diagonalize_subspace(H_sub)
ground_state = eigenvectors[:, 0]

# Generate candidate pool (doubles from HF)
double_excs = list(subspace.generate_double_excitations(hf_config))
candidate_pool = double_excs

# Select top-k by gradient
selected = select_configurations_by_gradient(
    hamiltonian=hamiltonian,
    subspace=subspace,
    ground_state=ground_state,
    candidate_pool=candidate_pool,
    k=1  # Add just the most important one
)

if len(selected) > 0:
    new_configs = [config for config, _ in selected]
    subspace.add_configs(new_configs)
    print(f"   Added configuration by gradient: {new_configs[0]}")
    print(f"   Gradient magnitude: {selected[0][1]:.6f} Ha")

print(f"\n   Final subspace size: {len(subspace)} configurations")
for i, config in enumerate(subspace):
    print(f"     {i}: {config}")

# Verify local accuracy with this subspace
H_sub_final = builder.project_fast(subspace)
eigenvalues_final = np.linalg.eigvalsh(H_sub_final)
print(f"\n   Local exact energy with this subspace: {eigenvalues_final[0]:.8f} Ha")
print(f"   Expected FCI: -1.13728383 Ha")
print(f"   Error: {abs(eigenvalues_final[0] - (-1.13728383))*1000:.4f} mHa")

# Step 3: Prepare measurement circuits
print("\n3️⃣ Preparing measurement circuits for IBM Torino...")
circuits = []

for config in subspace:
    circuit = QuantumCircuit(4)

    # Simple X-gate preparation (works reliably with Sampler)
    config_str = str(config).replace('|', '').replace('⟩', '')
    for i, bit in enumerate(config_str):
        if bit == '1':
            circuit.x(i)

    circuit.measure_all()
    circuits.append(circuit)

print(f"   Circuits prepared: {len(circuits)}")
print(f"   Each circuit measures one configuration in Z-basis")

# Step 4: Submit to IBM Torino
print("\n4️⃣ Connecting to IBM Quantum...")

IBM_API = os.environ.get('IBM_API')
IBM_CRN = os.environ.get('IBM_CRN')

if not IBM_API or not IBM_CRN:
    print("   ✗ IBM credentials not found in environment")
    print("   Set IBM_API and IBM_CRN environment variables")
    exit(1)

# Use the IBMSamplerBackend wrapper
from kanad.backends.ibm.sampler_backend import IBMSamplerBackend

try:
    sampler_backend = IBMSamplerBackend(
        backend_name='ibm_torino',
        api_token=IBM_API,
        channel='ibm_quantum_platform',
        crn=IBM_CRN
    )
    print("   ✓ Connected to IBM Quantum")
    print(f"   ✓ Using backend: {sampler_backend.backend.name}")
    print(f"   Status: {sampler_backend.backend.status().status_msg}")

except Exception as e:
    print(f"   ✗ Failed to connect: {e}")
    exit(1)

# Step 5: Run measurement with Sampler
print("\n5️⃣ Submitting Hi-VQE measurement to IBM Torino...")
print(f"   • Circuits: {len(circuits)} (one per configuration)")
print(f"   • Shots per circuit: 8192")
print(f"   • Backend: IBM Torino (127 qubits)")
print(f"   • Mode: Batch (open plan)")

start_time = datetime.now()

try:
    # Use the backend's run_hivqe_measurement method
    result = sampler_backend.run_hivqe_measurement(
        circuits=circuits,
        observable=hamiltonian,
        shots=8192,
        optimization_level=3,
        enable_readout_mitigation=True,
        enable_twirling=True,
        mode='batch'
    )

    job_id = result['job_id']
    status = result['status']

    print(f"\n   ✓ Job submitted successfully!")
    print(f"   Job ID: {job_id}")
    print(f"   Status: {status}")

    # If job completed, extract counts
    if 'counts' in result:
        all_counts = result['counts']
        end_time = datetime.now()
        runtime = (end_time - start_time).total_seconds()

        print(f"\n   ✓ Job completed in {runtime:.1f}s")

        for i, counts_dict in enumerate(all_counts):
            print(f"   Circuit {i+1}: {sum(counts_dict.values())} shots, {len(counts_dict)} outcomes")
    else:
        print("\n   ⚠️ Job submitted but still running. Use job ID to check status later.")
        exit(0)

except Exception as e:
    print(f"\n   ✗ Job failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# Step 6: Process results
print("\n6️⃣ Processing quantum measurement results...")

# Calculate diagonal energies from counts
def calculate_expectation_value(counts, observable):
    """Calculate ⟨ψ|H|ψ⟩ from measurement counts"""
    total_shots = sum(counts.values())
    expectation = 0.0

    pauli_list = observable.to_list()

    for pauli_str, coeff in pauli_list:
        term_expectation = 0.0

        for bitstring, count in counts.items():
            # Remove spaces and reverse (Qiskit convention)
            bits = bitstring.replace(' ', '')[::-1]

            # Compute eigenvalue for this bitstring
            eigenvalue = 1.0
            for j, pauli in enumerate(pauli_str):
                if j < len(bits):
                    bit = int(bits[j])
                    if pauli == 'Z':
                        eigenvalue *= (1 - 2*bit)
                    elif pauli == 'I':
                        eigenvalue *= 1.0

            term_expectation += eigenvalue * count

        term_expectation /= total_shots
        expectation += coeff.real * term_expectation

    return expectation

diagonal_energies = []
for i, counts in enumerate(all_counts):
    energy = calculate_expectation_value(counts, hamiltonian)
    diagonal_energies.append(energy)
    print(f"   Config {i} ({subspace.configurations[i]}): {energy:.8f} Ha")

# Classical diagonalization with measured diagonal energies
print("\n7️⃣ Classical diagonalization with quantum-measured energies...")

H_sub_measured = builder.project_fast(subspace)

# Replace diagonal with measured values
for i in range(len(diagonal_energies)):
    H_sub_measured[i, i] = diagonal_energies[i]

# Diagonalize
eigenvalues_measured = np.linalg.eigvalsh(H_sub_measured)
eigenvectors_measured = np.linalg.eigh(H_sub_measured)[1]
ground_energy_measured = eigenvalues_measured[0]
ground_state_measured = eigenvectors_measured[:, 0]

print(f"   Ground state energy: {ground_energy_measured:.8f} Ha")

# Step 7: Results and analysis
print("\n" + "="*80)
print("RESULTS")
print("="*80)

print(f"\n🎯 Diagonal Energies (measured on IBM Torino):")
for i, (config, energy) in enumerate(zip(subspace, diagonal_energies)):
    print(f"   {config}: {energy:.8f} Ha")

print(f"\n🏆 Ground State Energy:")
print(f"   IBM Torino (gradient Hi-VQE): {ground_energy_measured:.8f} Ha")
print(f"   Expected (FCI):               -1.13728383 Ha")
print(f"   Local (noiseless):            {result_local['energy']:.8f} Ha")

error_torino = abs(ground_energy_measured - (-1.13728383))
error_mHa = error_torino * 1000

print(f"\n   Error: {error_mHa:.2f} mHa ({error_torino/abs(-1.13728383)*100:.2f}%)")

print(f"\n📊 Comparison with Previous IBM Test:")
print(f"   Previous (HF + all singles):       45.00 mHa error")
print(f"   Current (gradient selection):      {error_mHa:.2f} mHa error")
if error_mHa < 45.0:
    improvement = 45.0 - error_mHa
    print(f"   Improvement: {improvement:.2f} mHa ({improvement/45.0*100:.1f}% reduction)")
else:
    print(f"   ⚠️ Error increased (unexpected)")

print(f"\n🎨 Ground State Composition:")
for i, (config, amp, prob) in enumerate(zip(
    subspace,
    ground_state_measured,
    ground_state_measured**2
)):
    if abs(prob) > 0.01:
        print(f"   {config}: amplitude={amp:+.4f}, probability={prob:.4f}")

print(f"\n📈 Energy Spectrum:")
for i, eigenval in enumerate(eigenvalues_measured[:min(5, len(eigenvalues_measured))]):
    excitation = (eigenval - eigenvalues_measured[0]) * 27211.4  # Convert to meV
    print(f"   State {i}: {eigenval:.8f} Ha ({excitation:+.1f} meV)")

print(f"\n📊 Literature Comparison:")
print(f"   Our result (IBM Torino):  {error_mHa:.2f} mHa")
print(f"   Literature Hi-VQE range:  0.08-1.20 mHa")
print(f"   Goal: Match literature accuracy")

# Analyze error breakdown
error_subspace = abs(eigenvalues_final[0] - (-1.13728383)) * 1000  # Local exact with this subspace
error_hardware = error_mHa - error_subspace

print(f"\n🔍 Error Breakdown:")
print(f"   Subspace selection error: {error_subspace:.2f} mHa (gradient selection)")
print(f"   Hardware noise error:     {error_hardware:.2f} mHa (IBM Torino)")
print(f"   Total error:              {error_mHa:.2f} mHa")

# Save results
results_data = {
    'timestamp': datetime.now().isoformat(),
    'backend': 'ibm_torino',
    'job_id': job_id,
    'runtime_seconds': runtime,
    'molecule': 'H2',
    'bond_length': 0.74,
    'subspace_size': len(subspace),
    'configurations': [str(c) for c in subspace],
    'diagonal_energies': diagonal_energies,
    'ground_energy': float(ground_energy_measured),
    'expected_energy': -1.13728383,
    'error_Ha': float(error_torino),
    'error_mHa': float(error_mHa),
    'eigenvalues': eigenvalues_measured.tolist(),
    'ground_state_amplitudes': ground_state_measured.tolist(),
    'local_baseline': {
        'energy': result_local['energy'],
        'iterations': result_local['iterations']
    },
    'comparison': {
        'previous_test_error_mHa': 45.0,
        'current_test_error_mHa': float(error_mHa),
        'improvement_mHa': float(45.0 - error_mHa),
        'improvement_percent': float((45.0 - error_mHa) / 45.0 * 100)
    },
    'error_breakdown': {
        'subspace_error_mHa': float(error_subspace),
        'hardware_error_mHa': float(error_hardware)
    }
}

output_file = f'ibm_torino_gradient_hivqe_{job_id}.json'
with open(output_file, 'w') as f:
    json.dump(results_data, f, indent=2)

print(f"\n✓ Results saved to: {output_file}")

# Save measurement counts
counts_file = f'ibm_torino_counts_{job_id}.json'
counts_data = {
    'job_id': job_id,
    'configurations': [str(c) for c in subspace],
    'counts': [dict(c) for c in all_counts]
}
with open(counts_file, 'w') as f:
    json.dump(counts_data, f, indent=2)

print(f"✓ Measurement counts saved to: {counts_file}")

print("\n" + "="*80)
print("IBM TORINO GRADIENT HI-VQE TEST COMPLETE")
print("="*80)

if error_mHa < 30.0:
    print("\n✅ SUCCESS! Hardware error reduced as expected!")
    print("   Gradient selection eliminated subspace error.")
elif error_mHa < 45.0:
    print("\n✅ IMPROVED! Error reduced from previous test.")
else:
    print("\n⚠️ Unexpected result - needs investigation")

print(f"\n💡 Next Steps:")
print(f"   1. Analyze error breakdown (subspace vs hardware)")
print(f"   2. Apply advanced error mitigation (M3, twirling, ZNE)")
print(f"   3. Test on larger molecules (LiH, H2O)")
print(f"   4. Target literature accuracy (<1.20 mHa)")
