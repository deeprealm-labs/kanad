"""
Submit Single Molecule with Full VQE to IBM Quantum
====================================================

Uses lightweight preparation to avoid memory issues.
Submits actual VQE circuit with Hamiltonian observable.

Usage:
    python submit_vqe_single_molecule.py --molecule caffeine
"""

import os
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime
from dotenv import load_dotenv

# Load credentials
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(env_path)

sys.path.insert(0, str(Path(__file__).parent.parent))

print("=" * 80)
print("IBM QUANTUM - VQE SINGLE MOLECULE SUBMISSION")
print("=" * 80)
print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--molecule', default='caffeine',
                   choices=['caffeine', 'aspirin', 'vitamin_c', 'h2'],
                   help='Molecule to test')
args = parser.parse_args()

# Check credentials
ibm_api = os.getenv('IBM_API')
if not ibm_api:
    print("✗ IBM_API not found")
    sys.exit(1)

print("✓ Credentials loaded\n")

# Define molecule
molecules = {
    'caffeine': ('C', 'N', 1.34, 'Caffeine - aromatic C-N bond'),
    'aspirin': ('C', 'O', 1.22, 'Aspirin - carbonyl C=O bond'),
    'vitamin_c': ('O', 'H', 0.96, 'Vitamin C - hydroxyl O-H bond'),
    'h2': ('H', 'H', 0.74, 'Hydrogen molecule - simplest test')
}

atom1, atom2, distance, description = molecules[args.molecule]

print(f"Molecule: {args.molecule.replace('_', ' ').title()}")
print(f"Description: {description}")
print(f"Bond: {atom1}-{atom2} @ {distance} Å\n")

try:
    from qiskit_ibm_runtime import QiskitRuntimeService, Batch, EstimatorV2 as Estimator
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    from kanad.bonds import BondFactory
    from kanad.backends.ibm.lightweight_preparation import LightweightIBMPreparation

    # Connect to IBM
    print("Connecting to IBM Quantum...")
    service = QiskitRuntimeService(
        channel='ibm_quantum_platform',
        token=ibm_api
    )
    print("✓ Connected!\n")

    # Get best backend
    backends = service.backends()
    backend_info = []
    for b in backends:
        status = b.status()
        if status.operational and not b.configuration().simulator:
            backend_info.append((b.name, status.pending_jobs, b.num_qubits))

    backend_info.sort(key=lambda x: x[1])
    backend_name = backend_info[0][0]
    backend = service.backend(backend_name)

    print(f"Backend: {backend_name}")
    print(f"  Qubits: {backend.num_qubits}")
    print(f"  Queue: {backend_info[0][1]} jobs\n")

    # Create molecule
    print("=" * 80)
    print("CREATING MOLECULE")
    print("=" * 80)

    bond = BondFactory.create_bond(atom1, atom2, distance=distance, basis='sto-3g')
    n_qubits = 2 * bond.hamiltonian.n_orbitals

    print(f"Bond type: {bond.bond_type}")
    print(f"Orbitals: {bond.hamiltonian.n_orbitals}")
    print(f"Qubits: {n_qubits}\n")

    # Prepare VQE with lightweight preparation
    print("=" * 80)
    print("PREPARING VQE CIRCUIT (Lightweight)")
    print("=" * 80)

    prep = LightweightIBMPreparation(bond, ansatz_type='hardware_efficient', n_layers=2)

    summary = prep.get_preparation_summary()
    print(f"\nPreparation Summary:")
    print(f"  HF Energy: {summary['hf_energy']:.6f} Ha")
    print(f"  Ansatz: {summary['ansatz_type']}")
    print(f"  Layers: {summary['n_layers']}")
    print(f"  Parameters: {summary['n_parameters']}")
    print(f"  Circuit depth: {summary['circuit_depth']}\n")

    # Get VQE circuits
    circuits, observables = prep.prepare_vqe_circuits()

    print(f"VQE circuit prepared:")
    print(f"  Qubits: {circuits[0].num_qubits}")
    print(f"  Depth: {circuits[0].depth()}")
    print(f"  Gates: {circuits[0].size()}")
    print(f"  Observable terms: {len(observables[0])}\n")

    # Transpile with initial layout to preserve qubit count
    print("=" * 80)
    print("TRANSPILING FOR HARDWARE")
    print("=" * 80)

    # Use first N qubits of hardware
    initial_layout = list(range(n_qubits))

    pm = generate_preset_pass_manager(
        optimization_level=1,
        backend=backend,
        initial_layout=initial_layout
    )
    transpiled_circuit = pm.run(circuits[0])

    print(f"\nTranspilation complete:")
    print(f"  Original qubits: {circuits[0].num_qubits}")
    print(f"  Transpiled qubits: {transpiled_circuit.num_qubits}")
    print(f"  Original depth: {circuits[0].depth()}")
    print(f"  Transpiled depth: {transpiled_circuit.depth()}")
    print(f"  Hardware gates: {transpiled_circuit.num_nonlocal_gates()}\n")

    # Submit to IBM Quantum
    print("=" * 80)
    print("SUBMITTING TO IBM QUANTUM")
    print("=" * 80)

    print(f"\nBackend: {backend_name}")
    print(f"Shots: 2048")
    print(f"Primitive: Estimator (for energy)")
    print(f"\nSubmitting...\n")

    with Batch(backend=backend) as batch:
        estimator = Estimator(mode=batch)
        estimator.options.default_shots = 2048

        # Submit as (circuit, observable) pub
        job = estimator.run([(transpiled_circuit, observables[0])])

        job_id = job.job_id()

        print(f"✓ JOB SUBMITTED!")
        print(f"\nJob ID: {job_id}")
        print(f"Status: {job.status()}")
        print(f"Backend: {backend_name}")

    # Save job info
    job_file = Path(__file__).parent / f'vqe_{args.molecule}_job.json'
    with open(job_file, 'w') as f:
        json.dump({
            'molecule': args.molecule,
            'description': description,
            'bond': f"{atom1}-{atom2}",
            'distance': distance,
            'qubits': n_qubits,
            'hf_energy': summary['hf_energy'],
            'job_id': job_id,
            'backend': backend_name,
            'shots': 2048,
            'submission_time': datetime.now().isoformat(),
            'preparation': summary
        }, f, indent=2)

    print(f"\nJob info saved to: {job_file}")

    print("\n" + "=" * 80)
    print("SUCCESS!")
    print("=" * 80)
    print(f"""
VQE job submitted for {args.molecule.replace('_', ' ').title()}!

HF Reference Energy: {summary['hf_energy']:.6f} Ha
Expected VQE Energy: < {summary['hf_energy']:.6f} Ha (lower due to correlation)

To check results:
  python research/get_ibm_results.py --job-id {job_id}

Dashboard:
  https://quantum.ibm.com/workloads

Expected wait time: 5-30 minutes
""")

except Exception as e:
    print(f"\n✗ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
