"""
Submit Single Test Job to IBM Quantum
======================================

Quickly submits ONE simple job to test the workflow.
"""

import os
import sys
from pathlib import Path
from datetime import datetime
from dotenv import load_dotenv

# Load credentials
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(env_path)

print("=" * 80)
print("IBM QUANTUM - SINGLE TEST JOB SUBMISSION")
print("=" * 80)
print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Add kanad to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from kanad.bonds import BondFactory

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
    from qiskit.circuit.library import EfficientSU2

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

    backend_info.sort(key=lambda x: x[1])  # Sort by queue length

    print("Available backends (sorted by queue):")
    for name, queue, qubits in backend_info[:3]:
        print(f"  {name}: {qubits} qubits, {queue} jobs in queue")

    # Use the one with shortest queue
    backend_name = backend_info[0][0]
    backend = service.backend(backend_name)

    print(f"\nUsing: {backend_name}\n")

    # Create simple H2 molecule
    print("=" * 80)
    print("CREATING H2 MOLECULE (TEST)")
    print("=" * 80)

    bond = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
    n_qubits = 2 * bond.hamiltonian.n_orbitals

    print(f"Bond: H-H @ 0.74 Å")
    print(f"Qubits: {n_qubits}")
    print(f"Bond type: {bond.bond_type}")

    # Get HF energy
    hf_result = bond.compute_energy(method='HF')
    print(f"HF energy: {hf_result['energy']:.6f} Ha\n")

    # Create simple test circuit
    print("Creating test circuit...")

    from qiskit import QuantumCircuit

    # Simple Bell state circuit (2 qubits)
    qc = QuantumCircuit(2, 2)
    qc.h(0)
    qc.cx(0, 1)
    qc.measure([0, 1], [0, 1])

    print(f"Circuit: Bell state (2 qubits)")
    print(f"Circuit depth: {qc.depth()}")
    print(f"Gates: {qc.num_nonlocal_gates()}\n")

    # Transpile for hardware
    print("Transpiling for hardware...")
    pm = generate_preset_pass_manager(optimization_level=1, backend=backend)
    transpiled_circuit = pm.run(qc)

    print(f"Transpiled depth: {transpiled_circuit.depth()}")
    print(f"Transpiled gates: {transpiled_circuit.num_nonlocal_gates()}\n")

    # Submit job
    print("=" * 80)
    print("SUBMITTING JOB TO IBM QUANTUM")
    print("=" * 80)
    print(f"Backend: {backend_name}")
    print(f"Shots: 1024")
    print(f"Submitting...\n")

    with Batch(backend=backend) as batch:
        sampler = Sampler(mode=batch)
        sampler.options.default_shots = 1024

        job = sampler.run([transpiled_circuit])

        job_id = job.job_id()

        print(f"✓ JOB SUBMITTED!")
        print(f"\nJob ID: {job_id}")
        print(f"Backend: {backend_name}")
        print(f"Status: {job.status()}")

        # Save job ID
        import json
        job_file = Path(__file__).parent / 'test_job_id.json'
        with open(job_file, 'w') as f:
            json.dump({
                'job_id': job_id,
                'backend': backend_name,
                'molecule': 'H2',
                'hf_energy': float(hf_result['energy']),
                'submission_time': datetime.now().isoformat(),
                'qubits': n_qubits
            }, f, indent=2)

        print(f"\nJob info saved to: {job_file}")

    print("\n" + "=" * 80)
    print("SUCCESS! CHECK YOUR IBM QUANTUM DASHBOARD")
    print("=" * 80)
    print(f"\nWorkload should appear at: https://quantum.ibm.com/workloads")
    print(f"\nTo check status later, run:")
    print(f"  python -c \"")
    print(f"from qiskit_ibm_runtime import QiskitRuntimeService;")
    print(f"s = QiskitRuntimeService(channel='ibm_quantum_platform');")
    print(f"j = s.job('{job_id}');")
    print(f"print(f'Status: {{j.status()}}')\"")

except Exception as e:
    print(f"\n✗ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
