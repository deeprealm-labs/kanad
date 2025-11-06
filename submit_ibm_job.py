"""
Simple IBM Quantum Job Submission for Hi-VQE

This script submits a Hi-VQE job to IBM Brisbane and exits immediately.
Check job status later with check_ibm_jobs.py

Usage:
    export IBM_API='your_token'
    export IBM_CRN='your_crn'
    python submit_ibm_job.py
"""
import os
import numpy as np
from kanad.bonds import BondFactory
from kanad.backends.ibm.backend import IBMBackend
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
print("SUBMIT HI-VQE JOB TO IBM QUANTUM")
print("="*80)

# Step 1: Create H2 molecule and get Hamiltonian
print("\n1️⃣ Building H2 molecule...")
bond = BondFactory.create_bond('H', 'H', distance=0.74)
hamiltonian = bond.hamiltonian.to_sparse_hamiltonian(mapper='jordan_wigner')
print(f"   ✓ Hamiltonian: {len(hamiltonian)} Pauli terms, {hamiltonian.num_qubits} qubits")

# Step 2: Build configuration subspace (same as local test)
print("\n2️⃣ Building Hi-VQE configuration subspace...")
protocol = CovalentGovernanceProtocol()
subspace = ConfigurationSubspace(4, 2, protocol=protocol)

# Add HF configuration
hf_config = subspace.get_hf_configuration()
subspace.add_config(hf_config)

# Add single excitations
single_excs = subspace.generate_single_excitations(hf_config)
subspace.add_configs(single_excs)

print(f"   ✓ Configurations: {len(subspace)}")

# Step 3: Create measurement circuits
print("\n3️⃣ Preparing measurement circuits...")
circuits = []
for config in subspace:
    circuit = QuantumCircuit(4)

    # Initialize to configuration state
    state_index = config.to_int()
    state_vector = np.zeros(16)
    state_vector[state_index] = 1.0

    circuit.initialize(state_vector, range(4))
    circuit.measure_all()

    circuits.append(circuit)

print(f"   ✓ Circuits: {len(circuits)} prepared")

# Step 4: Connect to IBM Torino
print("\n4️⃣ Connecting to IBM Torino...")
backend = IBMBackend(
    backend_name='ibm_torino',
    api_token=IBM_API,
    channel='ibm_cloud' if IBM_CRN else 'ibm_quantum_platform',
    crn=IBM_CRN
)

info = backend.get_backend_info()
print(f"   ✓ Backend: {info['name']}")
print(f"   ✓ Qubits: {info['num_qubits']}")
print(f"   ✓ Operational: {info['is_operational']}")
print(f"   ✓ Queue: {info['pending_jobs']} jobs")

if not info['is_operational']:
    print("\n   ✗ Backend not operational!")
    exit(1)

# Step 5: Submit job
print("\n5️⃣ Submitting job to IBM Quantum...")
print(f"   • Shots: 8192 per circuit")
print(f"   • Error mitigation: Level 2 (ZNE + readout)")
print(f"   • Optimization: Level 3 (maximum)")
print(f"   • Mode: Batch (open plan)")

try:
    # Use batch mode (session requires premium plan)
    job_result = backend.run_batch(
        circuits=circuits,
        observables=[hamiltonian] * len(circuits),
        shots=8192,
        optimization_level=3,
        resilience_level=2
    )

    print(f"\n✅ JOB SUBMITTED SUCCESSFULLY!")
    print(f"\n📋 Job Information:")
    print(f"   Job ID: {job_result['job_id']}")
    if 'session_id' in job_result:
        print(f"   Session ID: {job_result['session_id']}")
    print(f"   Status: {job_result['status']}")
    print(f"   Backend: {job_result['backend']}")
    print(f"   Mode: {job_result['mode']}")

    # Save job info
    with open('ibm_job_info.txt', 'w') as f:
        f.write(f"Job ID: {job_result['job_id']}\n")
        if 'session_id' in job_result:
            f.write(f"Session ID: {job_result['session_id']}\n")
        f.write(f"Backend: {job_result['backend']}\n")
        f.write(f"Status: {job_result['status']}\n")
        f.write(f"Mode: {job_result['mode']}\n")

    print(f"\n✓ Job info saved to: ibm_job_info.txt")

    print(f"\n📊 NEXT STEPS:")
    print(f"   The job is now in the IBM Quantum queue.")
    print(f"   Expected completion: 5-30 minutes")
    print(f"\n   Check status:")
    print(f"   $ python check_ibm_jobs.py {job_result['job_id']}")
    print(f"\n   Or list all jobs:")
    print(f"   $ python check_ibm_jobs.py")

except Exception as e:
    print(f"\n✗ Job submission failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

print("\n" + "="*80)
print("SUBMISSION COMPLETE")
print("="*80)
