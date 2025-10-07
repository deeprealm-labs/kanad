"""
Simple IBM Quantum Job Submission for Large Molecules
======================================================

Lightweight script that submits simple test circuits for each molecule
without complex VQE preparation that causes memory/timeout issues.

Usage:
    python submit_molecules_simple.py
"""

import os
import sys
import json
from pathlib import Path
from datetime import datetime
from dotenv import load_dotenv

# Load credentials
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(env_path)

sys.path.insert(0, str(Path(__file__).parent.parent))

print("=" * 80)
print("IBM QUANTUM - SIMPLE MOLECULE JOB SUBMISSION")
print("=" * 80)
print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Check credentials
ibm_api = os.getenv('IBM_API')
if not ibm_api:
    print("✗ IBM_API not found")
    sys.exit(1)

print("✓ Credentials loaded\n")

# Connect to IBM Quantum
print("Connecting to IBM Quantum...")

try:
    from qiskit_ibm_runtime import QiskitRuntimeService, Batch, SamplerV2 as Sampler
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    from qiskit import QuantumCircuit
    from kanad.bonds import BondFactory

    service = QiskitRuntimeService(
        channel='ibm_quantum_platform',
        token=ibm_api
    )

    print("✓ Connected!\n")

    # Get backend with shortest queue
    backends = service.backends()
    backend_info = []
    for b in backends:
        status = b.status()
        if status.operational and not b.configuration().simulator:
            backend_info.append((b.name, status.pending_jobs, b.num_qubits))

    backend_info.sort(key=lambda x: x[1])
    backend_name = backend_info[0][0]
    backend = service.backend(backend_name)

    print(f"Using: {backend_name} ({backend.num_qubits} qubits, {backend_info[0][1]} jobs in queue)\n")

    # Define molecules to test
    molecules = [
        {
            'name': 'Caffeine (C-N bond)',
            'description': 'Psychoactive drug - aromatic C-N',
            'atom1': 'C', 'atom2': 'N', 'distance': 1.34
        },
        {
            'name': 'Aspirin (C=O bond)',
            'description': 'Pharmaceutical - carbonyl C=O',
            'atom1': 'C', 'atom2': 'O', 'distance': 1.22
        },
        {
            'name': 'Vitamin C (O-H bond)',
            'description': 'Antioxidant - hydroxyl O-H',
            'atom1': 'O', 'atom2': 'H', 'distance': 0.96
        },
        {
            'name': 'Cholesterol (C-C bond)',
            'description': 'Biological lipid - sp³ C-C',
            'atom1': 'C', 'atom2': 'C', 'distance': 1.54
        },
        {
            'name': 'Adenine (N-H bond)',
            'description': 'DNA base - aromatic N-H',
            'atom1': 'N', 'atom2': 'H', 'distance': 1.01
        },
        {
            'name': 'Heme (Fe-N bond)',
            'description': 'Hemoglobin - Fe-N coordination',
            'atom1': 'Fe', 'atom2': 'N', 'distance': 2.0
        },
        {
            'name': 'MOF (Cu-O bond)',
            'description': 'Catalyst - Cu-O coordination',
            'atom1': 'Cu', 'atom2': 'O', 'distance': 1.95
        },
        {
            'name': 'Penicillin (S-C bond)',
            'description': 'Antibiotic - S-C thiazolidine',
            'atom1': 'S', 'atom2': 'C', 'distance': 1.82
        },
        {
            'name': 'Chlorophyll (Mg-N bond)',
            'description': 'Photosynthesis - Mg-N coordination',
            'atom1': 'Mg', 'atom2': 'N', 'distance': 2.07
        },
        {
            'name': 'Serotonin (aromatic C-C)',
            'description': 'Neurotransmitter - aromatic C-C',
            'atom1': 'C', 'atom2': 'C', 'distance': 1.40
        }
    ]

    print("=" * 80)
    print(f"SUBMITTING {len(molecules)} MOLECULE EXPERIMENTS")
    print("=" * 80)
    print()

    job_data = {
        'submission_time': datetime.now().isoformat(),
        'backend': backend_name,
        'experiments': []
    }

    # Transpile manager (reuse for all circuits)
    pm = generate_preset_pass_manager(optimization_level=1, backend=backend)

    for i, mol in enumerate(molecules, 1):
        print(f"[{i}/{len(molecules)}] {mol['name']}")
        print(f"  {mol['description']}")
        print(f"  Bond: {mol['atom1']}-{mol['atom2']} @ {mol['distance']} Å")

        try:
            # Create bond to get HF energy
            bond = BondFactory.create_bond(
                mol['atom1'], mol['atom2'],
                distance=mol['distance'],
                basis='sto-3g'
            )

            n_qubits = 2 * bond.hamiltonian.n_orbitals
            print(f"  Qubits: {n_qubits}")

            # Get HF reference
            hf_result = bond.compute_energy(method='HF')
            hf_energy = hf_result['energy']
            print(f"  HF energy: {hf_energy:.6f} Ha")

            # Create simple variational circuit (NOT full VQE - just a test circuit)
            # Use minimal qubits (2-4) for faster execution
            test_qubits = min(4, n_qubits)

            qc = QuantumCircuit(test_qubits, test_qubits)

            # Simple ansatz: rotations + entanglement
            for q in range(test_qubits):
                qc.ry(0.1 * (q + 1), q)  # Small rotation angles

            for q in range(test_qubits - 1):
                qc.cx(q, q + 1)

            # Measure all
            qc.measure(range(test_qubits), range(test_qubits))

            print(f"  Test circuit: {test_qubits} qubits, depth {qc.depth()}")

            # Transpile
            transpiled = pm.run(qc)
            print(f"  Transpiled depth: {transpiled.depth()}")

            # Submit to IBM Quantum
            print(f"  Submitting to {backend_name}...")

            with Batch(backend=backend) as batch:
                sampler = Sampler(mode=batch)
                sampler.options.default_shots = 1024

                job = sampler.run([transpiled])
                job_id = job.job_id()

                print(f"  ✓ Job submitted: {job_id}")

                # Store job info
                job_data['experiments'].append({
                    'name': mol['name'],
                    'description': mol['description'],
                    'bond': f"{mol['atom1']}-{mol['atom2']}",
                    'distance': mol['distance'],
                    'qubits': n_qubits,
                    'test_qubits': test_qubits,
                    'hf_energy': float(hf_energy),
                    'job_id': job_id,
                    'status': 'submitted',
                    'submission_time': datetime.now().isoformat()
                })

        except Exception as e:
            print(f"  ✗ Failed: {e}")
            import traceback
            traceback.print_exc()

        print()

    # Save job IDs
    job_file = Path(__file__).parent / 'ibm_job_ids.json'
    with open(job_file, 'w') as f:
        json.dump(job_data, f, indent=2)

    print("=" * 80)
    print("SUBMISSION COMPLETE")
    print("=" * 80)
    print(f"\n✓ Submitted {len(job_data['experiments'])} jobs")
    print(f"✓ Job IDs saved to: {job_file}")

    print(f"\nJob IDs:")
    for exp in job_data['experiments']:
        print(f"  {exp['name']}: {exp['job_id']}")

    print("\n" + "=" * 80)
    print("NEXT STEPS")
    print("=" * 80)
    print("""
1. Jobs are queued on IBM Quantum hardware
2. Expected wait time: 5-30 minutes per job
3. Check results with:
   python research/get_ibm_results.py

4. View in dashboard:
   https://quantum.ibm.com/workloads
""")

except Exception as e:
    print(f"\n✗ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
