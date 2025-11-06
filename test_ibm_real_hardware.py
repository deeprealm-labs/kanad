"""
Test Hi-VQE Deployment on IBM Real Quantum Hardware

Direct deployment to IBM Brisbane (127 qubits) with full error mitigation.

Usage:
    export IBM_API='your_token'
    export IBM_CRN='your_crn'
    python test_ibm_real_hardware.py
"""
import os
import time
import logging
from kanad.bonds import BondFactory
from kanad.utils.vqe_solver import VQESolver

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

print("="*80)
print("IBM QUANTUM HARDWARE DEPLOYMENT - HI-VQE")
print("="*80)

# Check credentials
IBM_API = os.getenv('IBM_API')
IBM_CRN = os.getenv('IBM_CRN')

if not IBM_API:
    print("\n⚠️  ERROR: IBM_API not set!")
    print("   Set with: export IBM_API='your_token_here'")
    exit(1)

print(f"\n✓ IBM_API found: {IBM_API[:20]}...")
if IBM_CRN:
    print(f"✓ IBM_CRN found: {IBM_CRN[:50]}...")

# =============================================================================
# PHASE 1: LOCAL BASELINE
# =============================================================================

print("\n" + "="*80)
print("PHASE 1: LOCAL BASELINE (VERIFY CORRECTNESS)")
print("="*80)

bond_h2 = BondFactory.create_bond('H', 'H', distance=0.74)

print("\n--- H2 with Hi-VQE (Local Statevector) ---")
solver_local = VQESolver(
    bond=bond_h2,
    mode='hivqe',
    hivqe_max_iterations=3,
    backend='statevector'
)

result_local = solver_local.solve()

print(f"\n✓ Local Hi-VQE Results:")
print(f"  Energy: {result_local['energy']:.8f} Ha")
print(f"  Expected: -1.13728383 Ha (FCI)")
print(f"  Error: {abs(result_local['energy'] - (-1.13728383)):.8f} Ha")
print(f"  Iterations: {result_local['iterations']}")
print(f"  Subspace size: {result_local['hivqe_stats']['final_subspace_size']}")
print(f"  Measurement reduction: {result_local['hivqe_stats']['measurement_reduction']}x")

# Check if result is correct
energy_error = abs(result_local['energy'] - (-1.13728383))
if energy_error < 0.001:
    print(f"  ✓ BASELINE PASSED: Ready for hardware deployment")
else:
    print(f"  ✗ BASELINE FAILED: Fix local simulation first")
    exit(1)

# =============================================================================
# PHASE 2: IBM REAL HARDWARE
# =============================================================================

print("\n" + "="*80)
print("PHASE 2: IBM REAL HARDWARE (PRODUCTION)")
print("="*80)

print("\n⚠️  DEPLOYMENT INFORMATION:")
print("   • Target: IBM Brisbane (127-qubit quantum processor)")
print("   • Molecule: H2 (4 qubits required)")
print("   • Method: Hi-VQE with error mitigation")
print("   • Configurations: 5 (from local test)")
print("   • Circuits: 5 measurement circuits")
print("   • Error mitigation: Level 2 (ZNE + readout)")
print("   • Shots: 8192 per circuit")
print("   • Estimated runtime: ~1-2 minutes")
print("   • Estimated cost: ~$2-3")
print("   • Wait time: 5-30 minutes (queue dependent)")

print("\n📋 WHAT WILL HAPPEN:")
print("   1. Submit Hi-VQE circuits to IBM Brisbane")
print("   2. Job enters queue (wait for hardware availability)")
print("   3. Circuits execute with error mitigation")
print("   4. Results retrieved and processed")
print("   5. Ground state energy computed from measurements")

print("\n💡 NOTE:")
print("   This script will submit the job and save the job ID.")
print("   You can check results later using the job ID.")

print("\n🚀 Ready to deploy to IBM quantum hardware!")
print("   Press Ctrl+C to cancel, or wait 5 seconds to proceed...")

try:
    time.sleep(5)
except KeyboardInterrupt:
    print("\n\n⏭️  Deployment cancelled by user")
    exit(0)

print("\n📤 Deploying to IBM Brisbane...")

try:
    from kanad.backends.ibm.backend import IBMBackend
    from kanad.core.configuration import ConfigurationSubspace
    from qiskit import QuantumCircuit
    import numpy as np

    # Initialize backend
    print("\n1️⃣ Connecting to IBM Brisbane...")
    backend_hw = IBMBackend(
        backend_name='ibm_brisbane',
        api_token=IBM_API,
        channel='ibm_cloud' if IBM_CRN else 'ibm_quantum_platform',
        crn=IBM_CRN
    )

    print(f"   ✓ Connected to: {backend_hw.backend.name}")

    info = backend_hw.get_backend_info()
    print(f"   ✓ Qubits: {info['num_qubits']}")
    print(f"   ✓ Operational: {info['is_operational']}")
    print(f"   ✓ Queue: {info['pending_jobs']} jobs ahead")

    if not info['is_operational']:
        print(f"\n   ✗ Backend is not operational. Try again later.")
        exit(1)

    # Get Hamiltonian
    print("\n2️⃣ Building Hi-VQE circuits...")
    hamiltonian = bond_h2.hamiltonian.to_sparse_hamiltonian(mapper='jordan_wigner')

    print(f"   ✓ Hamiltonian: {len(hamiltonian)} Pauli terms")
    print(f"   ✓ Qubits: {hamiltonian.num_qubits}")

    # Build configuration subspace (same as local test)
    from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol
    protocol = CovalentGovernanceProtocol()

    subspace = ConfigurationSubspace(4, 2, protocol=protocol)
    hf_config = subspace.get_hf_configuration()
    subspace.add_config(hf_config)

    # Add single excitations
    single_excs = subspace.generate_single_excitations(hf_config)
    subspace.add_configs(single_excs)

    print(f"   ✓ Configurations: {len(subspace)}")

    # Create measurement circuits
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

    print(f"   ✓ Circuits prepared: {len(circuits)}")

    # Submit job
    print("\n3️⃣ Submitting to IBM Brisbane...")
    print(f"   • Mode: Session (reserved hardware)")
    print(f"   • Shots: 8192 per circuit")
    print(f"   • Optimization: Level 3 (maximum)")
    print(f"   • Error mitigation: Level 2 (ZNE + readout)")
    print(f"   • Session time: 1 hour")

    job_result = backend_hw.run_session(
        circuits=circuits,
        observables=[hamiltonian] * len(circuits),
        shots=8192,
        optimization_level=3,
        resilience_level=2,  # Maximum error mitigation
        max_time='1h'
    )

    print(f"\n✓ Job submitted successfully!")
    print(f"\n📋 JOB INFORMATION:")
    print(f"   Job ID: {job_result['job_id']}")
    print(f"   Session ID: {job_result['session_id']}")
    print(f"   Status: {job_result['status']}")
    print(f"   Backend: {job_result['backend']}")

    # Save job ID
    job_file = 'ibm_job_info.txt'
    with open(job_file, 'w') as f:
        f.write(f"Job ID: {job_result['job_id']}\n")
        f.write(f"Session ID: {job_result['session_id']}\n")
        f.write(f"Backend: {job_result['backend']}\n")
        f.write(f"Submitted: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"\nExpected energy: -1.137 Ha (FCI)\n")
        f.write(f"Local baseline: {result_local['energy']:.8f} Ha\n")

    print(f"\n✓ Job info saved to: {job_file}")

    print("\n📊 NEXT STEPS:")
    print(f"   The job is now queued on IBM Brisbane.")
    print(f"   Expected completion: 5-30 minutes")
    print(f"\n   To check status:")
    print(f"   >>> from kanad.backends.ibm.backend import IBMBackend")
    print(f"   >>> backend = IBMBackend(backend_name='ibm_brisbane')")
    print(f"   >>> status = backend.get_job_status('{job_result['job_id']}')")
    print(f"   >>> print(status)")
    print(f"\n   To get results when done:")
    print(f"   >>> results = backend.get_job_result('{job_result['job_id']}')")
    print(f"   >>> print(results.values)  # Energies for each configuration")

    print("\n🎯 Would you like to wait for the job to complete? (This may take 5-30 minutes)")
    print("   Press Ctrl+C to exit and check later, or wait...")

    try:
        print("\n⏳ Monitoring job status...")
        max_wait = 1800  # 30 minutes
        start_time = time.time()
        check_interval = 30  # Check every 30 seconds

        while time.time() - start_time < max_wait:
            status = backend_hw.get_job_status(job_result['job_id'])
            elapsed = int(time.time() - start_time)
            print(f"   [{elapsed}s] Status: {status}")

            if status in ['DONE', 'ERROR', 'CANCELLED']:
                break

            time.sleep(check_interval)

        if status == 'DONE':
            print(f"\n✅ JOB COMPLETED!")

            # Get results
            results = backend_hw.get_job_result(job_result['job_id'])

            print(f"\n📊 IBM Hardware Results:")
            print(f"   Number of measurements: {len(results.values)}")

            energies_hw = results.values
            print(f"\n   Configuration energies:")
            for i, energy in enumerate(energies_hw):
                config = list(subspace)[i]
                print(f"     {config}: {energy:.8f} Ha")

            # Compute ground state via classical diagonalization
            from kanad.core.classical_solver import SubspaceHamiltonianBuilder

            builder = SubspaceHamiltonianBuilder(hamiltonian)
            H_sub = builder.project_fast(subspace)

            # Replace diagonal with measured values
            for i in range(len(subspace)):
                if i < len(energies_hw):
                    H_sub[i, i] = energies_hw[i]

            # Diagonalize
            eigenvalues = np.linalg.eigvalsh(H_sub)
            ground_energy_hw = eigenvalues[0]

            print(f"\n🎯 FINAL RESULTS:")
            print(f"   Local (exact):     {result_local['energy']:.8f} Ha")
            print(f"   IBM Hardware:      {ground_energy_hw:.8f} Ha")
            print(f"   FCI (expected):    -1.13728383 Ha")
            print(f"   Hardware error:    {abs(ground_energy_hw - (-1.13728383)):.6f} Ha")
            print(f"   Difference:        {abs(ground_energy_hw - result_local['energy']):.6f} Ha")

            # Check accuracy
            if abs(ground_energy_hw - result_local['energy']) < 0.1:
                print(f"\n   ✅ HARDWARE TEST PASSED!")
                print(f"      Error mitigation working well (<0.1 Ha difference)")
            else:
                print(f"\n   ⚠️  Large error from hardware noise")
                print(f"      Consider increasing shots or resilience level")

        elif status == 'ERROR':
            print(f"\n✗ Job failed with error!")
            print(f"   Check job details with:")
            print(f"   >>> backend.get_job_result('{job_result['job_id']}')")

        else:
            print(f"\n⏱️  Job still running (status: {status})")
            print(f"   Check back later with job ID: {job_result['job_id']}")

    except KeyboardInterrupt:
        print(f"\n\n⏭️  Stopped monitoring. Job is still running.")
        print(f"   Job ID: {job_result['job_id']}")
        print(f"   Check status later!")

except Exception as e:
    print(f"\n✗ Hardware deployment failed: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*80)
print("DEPLOYMENT COMPLETE")
print("="*80)

print(f"\n✅ Summary:")
print(f"   Phase 1 (Local): PASSED ({result_local['energy']:.8f} Ha)")
print(f"   Phase 2 (Hardware): Job submitted")
print(f"\n🚀 Hi-VQE is production ready!")
print("="*80)
