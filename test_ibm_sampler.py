"""
Test Hi-VQE on IBM Quantum using Sampler

Uses Sampler (not Estimator) to get raw measurement counts,
enabling custom error mitigation and rich visualization data.

Usage:
    export IBM_API='your_token'
    export IBM_CRN='your_crn'
    python test_ibm_sampler.py
"""
import os
import json
import numpy as np
from kanad.bonds import BondFactory
from kanad.backends.ibm.sampler_backend import IBMSamplerBackend
from kanad.core.configuration import ConfigurationSubspace
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol
from qiskit import QuantumCircuit

# Credentials
IBM_API = os.getenv('IBM_API')
IBM_CRN = os.getenv('IBM_CRN')

if not IBM_API:
    print("ERROR: IBM_API not set")
    exit(1)

print("="*80)
print("HI-VQE ON IBM QUANTUM - SAMPLER APPROACH")
print("="*80)

# Step 1: Local baseline
print("\n1️⃣ Running local baseline...")
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
print(f"   Expected:     -1.13728383 Ha")
print(f"   Iterations:   {result_local['iterations']}")
print(f"   ✓ Baseline validated")

# Step 2: Build configuration subspace
print("\n2️⃣ Building Hi-VQE configuration subspace...")
hamiltonian = bond.hamiltonian.to_sparse_hamiltonian(mapper='jordan_wigner')

protocol = CovalentGovernanceProtocol()
subspace = ConfigurationSubspace(4, 2, protocol=protocol)

# Add HF configuration
hf_config = subspace.get_hf_configuration()
subspace.add_config(hf_config)

# Add single excitations
single_excs = subspace.generate_single_excitations(hf_config)
subspace.add_configs(single_excs)

print(f"   Configurations: {len(subspace)}")
for i, config in enumerate(subspace):
    print(f"     {i}: {config}")

# Step 3: Create measurement circuits (Z-basis measurements)
print("\n3️⃣ Preparing Z-basis measurement circuits...")
circuits = []

for config in subspace:
    circuit = QuantumCircuit(4)

    # Prepare state corresponding to configuration
    # For Z-basis measurement, we just need to flip bits for |1⟩ states
    config_str = str(config).replace('|', '').replace('⟩', '')  # Get binary string
    for i, bit in enumerate(config_str):
        if bit == '1':
            circuit.x(i)  # Apply X gate to flip |0⟩ → |1⟩

    # Measure all qubits
    circuit.measure_all()

    circuits.append(circuit)

print(f"   Circuits prepared: {len(circuits)}")
print(f"   All circuits measure in Z-basis (computational basis)")

# Step 4: Connect to IBM Torino
print("\n4️⃣ Connecting to IBM Torino (Sampler)...")
backend = IBMSamplerBackend(
    backend_name='ibm_torino',
    api_token=IBM_API,
    channel='ibm_cloud' if IBM_CRN else 'ibm_quantum_platform',
    crn=IBM_CRN
)

print(f"   ✓ Connected to: {backend.backend.name}")
print(f"   ✓ Qubits: {backend.backend.num_qubits}")

# Step 5: Submit job
print("\n5️⃣ Submitting Hi-VQE measurement job...")
print(f"   • Circuits: {len(circuits)} (one per configuration)")
print(f"   • Shots: 8192 per circuit")
print(f"   • Readout mitigation: Enabled")
print(f"   • Twirling: Enabled (Pauli randomization)")
print(f"   • Mode: Batch")

try:
    job_result = backend.run_hivqe_measurement(
        circuits=circuits,
        observable=hamiltonian,
        shots=8192,
        optimization_level=3,
        enable_readout_mitigation=True,
        enable_twirling=True,
        mode='batch'
    )

    print(f"\n✅ JOB SUBMITTED!")
    print(f"\n📋 Job Information:")
    print(f"   Job ID: {job_result['job_id']}")
    print(f"   Status: {job_result['status']}")
    print(f"   Backend: {job_result['backend']}")
    print(f"   Mode: {job_result['mode']}")

    # Save job info
    job_info = {
        'job_id': job_result['job_id'],
        'backend': job_result['backend'],
        'mode': job_result['mode'],
        'config': job_result['config'],
        'local_baseline': {
            'energy': result_local['energy'],
            'iterations': result_local['iterations']
        },
        'subspace': {
            'size': len(subspace),
            'configurations': [str(c) for c in subspace]
        }
    }

    with open('ibm_sampler_job.json', 'w') as f:
        json.dump(job_info, f, indent=2)

    print(f"\n✓ Job info saved to: ibm_sampler_job.json")

    print(f"\n📊 NEXT STEPS:")
    print(f"   The job is now running on IBM Torino.")
    print(f"   Expected completion: 2-5 minutes")
    print(f"\n   Check status:")
    print(f"   $ python check_ibm_jobs.py {job_result['job_id']}")
    print(f"\n   Process results when done:")
    print(f"   $ python process_sampler_results.py {job_result['job_id']}")

    print(f"\n💡 What's Different:")
    print(f"   • Using Sampler (not Estimator)")
    print(f"   • Get raw measurement counts")
    print(f"   • Custom expectation value calculation")
    print(f"   • Rich data for visualization")
    print(f"   • Better error analysis")

except Exception as e:
    print(f"\n✗ Job submission failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

print("\n" + "="*80)
print("SUBMISSION COMPLETE")
print("="*80)
