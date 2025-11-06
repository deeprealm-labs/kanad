"""
Monitor the Gradient-Optimized Hi-VQE Job on IBM Torino

Usage:
    export IBM_API='your_api_key'
    export IBM_CRN='your_crn'
    python monitor_gradient_hivqe_job.py <job_id>
"""
import os
import sys
import time
import json
import numpy as np
from datetime import datetime
from qiskit_ibm_runtime import QiskitRuntimeService

if len(sys.argv) < 2:
    print("Usage: python monitor_gradient_hivqe_job.py <job_id>")
    exit(1)

job_id = sys.argv[1]

IBM_API = os.environ.get('IBM_API')
IBM_CRN = os.environ.get('IBM_CRN')

if not IBM_API:
    print("ERROR: IBM_API not set")
    exit(1)

print("="*80)
print(f"MONITORING IBM TORINO JOB: {job_id}")
print("="*80)

# Connect to service
service = QiskitRuntimeService(
    channel='ibm_quantum_platform',
    token=IBM_API,
    instance=IBM_CRN
)

# Retrieve job
try:
    job = service.job(job_id)
    print(f"\n✓ Job retrieved: {job_id}")
    print(f"  Backend: {job.backend().name}")
except Exception as e:
    print(f"\n✗ Failed to retrieve job: {e}")
    exit(1)

# Monitor job status
print(f"\n⏳ Waiting for job to complete...")
print(f"   Current status: {job.status()}")

# Poll every 10 seconds
check_count = 0
while True:
    status = job.status()
    check_count += 1

    print(f"\n[Check {check_count}] {datetime.now().strftime('%H:%M:%S')}")
    print(f"   Status: {status}")

    # Status is already a string in newer Qiskit Runtime
    status_str = status if isinstance(status, str) else status.name

    if status_str in ['DONE', 'ERROR', 'CANCELLED']:
        break

    time.sleep(10)  # Wait 10 seconds

# Handle completion
status_str = status if isinstance(status, str) else status.name
if status_str == 'DONE':
    print(f"\n✅ JOB COMPLETED!")

    # Get results
    try:
        result = job.result()
        print(f"\n📊 Processing results...")

        # Extract counts from all circuits
        all_counts = []
        for i in range(5):  # 5 circuits for gradient Hi-VQE
            pub_result = result[i]
            counts_dict = pub_result.data.meas.get_counts()
            all_counts.append(counts_dict)
            print(f"   Circuit {i+1}: {sum(counts_dict.values())} shots, {len(counts_dict)} outcomes")

        # Load the Hamiltonian to calculate energies
        from kanad.bonds import BondFactory
        from kanad.core.configuration import ConfigurationSubspace
        from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol
        from kanad.core.classical_solver import SubspaceHamiltonianBuilder

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        hamiltonian = bond.hamiltonian.to_sparse_hamiltonian(mapper='jordan_wigner')

        # Rebuild subspace (HF + singles + |0011⟩ double)
        protocol = CovalentGovernanceProtocol()
        subspace = ConfigurationSubspace(4, 2, protocol=protocol)

        hf_config = subspace.get_hf_configuration()
        subspace.add_config(hf_config)

        single_excs = subspace.generate_single_excitations(hf_config)
        subspace.add_configs(single_excs)

        # Add the |0011⟩ double excitation
        double_excs = list(subspace.generate_double_excitations(hf_config))
        # Find |0011⟩ in the doubles
        for double in double_excs:
            if str(double) == '|0011⟩':
                subspace.add_config(double)
                break

        print(f"\n   Reconstructed subspace: {len(subspace)} configurations")

        # Calculate expectation values from counts
        def calculate_expectation_value(counts, observable):
            total_shots = sum(counts.values())
            expectation = 0.0

            pauli_list = observable.to_list()

            for pauli_str, coeff in pauli_list:
                term_expectation = 0.0

                for bitstring, count in counts.items():
                    bits = bitstring.replace(' ', '')[::-1]

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
        config_list = list(subspace)
        for i, counts in enumerate(all_counts):
            energy = calculate_expectation_value(counts, hamiltonian)
            diagonal_energies.append(energy)
            print(f"   Config {i} ({config_list[i]}): {energy:.8f} Ha")

        # Classical diagonalization
        builder = SubspaceHamiltonianBuilder(hamiltonian)
        H_sub = builder.project_fast(subspace)

        # Replace diagonal with measured values
        for i in range(len(diagonal_energies)):
            H_sub[i, i] = diagonal_energies[i]

        # Diagonalize
        eigenvalues = np.linalg.eigvalsh(H_sub)
        eigenvectors = np.linalg.eigh(H_sub)[1]
        ground_energy = eigenvalues[0]
        ground_state = eigenvectors[:, 0]

        # Results summary
        print(f"\n" + "="*80)
        print("RESULTS SUMMARY")
        print("="*80)

        print(f"\n🏆 Ground State Energy:")
        print(f"   IBM Torino (gradient Hi-VQE): {ground_energy:.8f} Ha")
        print(f"   Expected (FCI):               -1.13728383 Ha")

        error_Ha = abs(ground_energy - (-1.13728383))
        error_mHa = error_Ha * 1000

        print(f"\n   Error: {error_mHa:.2f} mHa ({error_Ha/abs(-1.13728383)*100:.2f}%)")

        print(f"\n📊 Comparison with Previous Test:")
        print(f"   Previous (HF + all singles):       45.00 mHa error")
        print(f"   Current (gradient selection):      {error_mHa:.2f} mHa error")

        if error_mHa < 45.0:
            improvement = 45.0 - error_mHa
            print(f"   Improvement: {improvement:.2f} mHa ({improvement/45.0*100:.1f}% reduction)")
            print(f"\n   ✅ SUCCESS! Gradient selection improved accuracy!")
        else:
            print(f"   ⚠️ Error increased (unexpected)")

        print(f"\n🎨 Ground State Composition:")
        for i, (config, amp, prob) in enumerate(zip(
            config_list,
            ground_state,
            ground_state**2
        )):
            if abs(prob) > 0.01:
                print(f"   {config}: amplitude={amp:+.4f}, probability={prob:.4f}")

        # Save results
        results_data = {
            'timestamp': datetime.now().isoformat(),
            'job_id': job_id,
            'backend': 'ibm_torino',
            'molecule': 'H2',
            'bond_length': 0.74,
            'gradient_selection': True,
            'subspace_size': len(subspace),
            'configurations': [str(c) for c in config_list],
            'diagonal_energies': diagonal_energies,
            'ground_energy': float(ground_energy),
            'expected_energy': -1.13728383,
            'error_Ha': float(error_Ha),
            'error_mHa': float(error_mHa),
            'eigenvalues': eigenvalues.tolist(),
            'ground_state_amplitudes': ground_state.tolist(),
            'comparison': {
                'previous_test_error_mHa': 45.0,
                'current_test_error_mHa': float(error_mHa),
                'improvement_mHa': float(45.0 - error_mHa),
                'improvement_percent': float((45.0 - error_mHa) / 45.0 * 100)
            }
        }

        output_file = f'ibm_torino_gradient_results_{job_id}.json'
        with open(output_file, 'w') as f:
            json.dump(results_data, f, indent=2)

        print(f"\n✓ Results saved to: {output_file}")

        # Save counts
        counts_file = f'ibm_torino_gradient_counts_{job_id}.json'
        counts_data = {
            'job_id': job_id,
            'configurations': [str(c) for c in config_list],
            'counts': [dict(c) for c in all_counts]
        }
        with open(counts_file, 'w') as f:
            json.dump(counts_data, f, indent=2)

        print(f"✓ Measurement counts saved to: {counts_file}")

    except Exception as e:
        print(f"\n✗ Error processing results: {e}")
        import traceback
        traceback.print_exc()

elif status_str == 'ERROR':
    print(f"\n✗ Job failed with error")
    try:
        print(f"   Error message: {job.error_message()}")
    except:
        print(f"   No error message available")

elif status_str == 'CANCELLED':
    print(f"\n⚠️ Job was cancelled")

print("\n" + "="*80)
