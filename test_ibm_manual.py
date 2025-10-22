#!/usr/bin/env python3
"""
Manual test script for IBM Quantum backend.
Tests the full integration with real credentials.
"""

import os
import sys

# Set credentials from known values
os.environ['IBM_API'] = 'fPgicVCXNNj_Q26nqeEGUPpcaVET5dbjA6cv-cldgVpz'
os.environ['IBM_CRN'] = 'crn:v1:bluemix:public:quantum-computing:us-east:a/2b8172b414734f9c921f8173f2d420c4:d601ae59-5aaf-40b5-8240-a4e4c71521d3::'

print("=" * 80)
print("IBM QUANTUM BACKEND TEST")
print("=" * 80)

# Test 1: Backend initialization
print("\n[Test 1] Initializing IBM backend...")
try:
    from kanad.backends.ibm.backend import IBMBackend

    backend = IBMBackend(
        backend_name='ibm_brisbane',
        api_token=os.environ['IBM_API'],
        instance=os.environ['IBM_CRN']
    )

    print(f"✓ Backend initialized successfully")
    print(f"  Name: {backend.backend.name}")
    print(f"  Qubits: {backend.backend.num_qubits}")
    print(f"  Simulator: {backend.backend.simulator}")

except Exception as e:
    print(f"✗ Failed: {e}")
    sys.exit(1)

# Test 2: Get backend info
print("\n[Test 2] Getting backend information...")
try:
    info = backend.get_backend_info()
    print(f"✓ Backend info retrieved:")
    print(f"  Operational: {info['is_operational']}")
    print(f"  Pending jobs: {info['pending_jobs']}")
    print(f"  Basis gates: {info['basis_gates'][:5] if info['basis_gates'] else 'N/A'}...")

except Exception as e:
    print(f"✗ Failed: {e}")
    sys.exit(1)

# Test 3: Submit a simple circuit
print("\n[Test 3] Submitting simple test circuit...")
try:
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import SparsePauliOp

    # Create Bell state circuit
    circuit = QuantumCircuit(2)
    circuit.h(0)
    circuit.cx(0, 1)

    # Observable: measure Z on first qubit
    observable = SparsePauliOp.from_list([("ZI", 1.0)])

    print(f"  Circuit: {circuit.num_qubits} qubits, depth {circuit.depth()}")
    print(f"  Submitting to {backend.backend.name}...")

    result = backend.run_batch(
        circuits=[circuit],
        observables=[observable],
        shots=1024
    )

    print(f"✓ Job submitted successfully")
    print(f"  Job ID: {result['job_id']}")
    print(f"  Status: {result['status']}")
    print(f"  Backend: {result['backend']}")

    job_id = result['job_id']

except Exception as e:
    print(f"✗ Failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 4: Check job status
print("\n[Test 4] Checking job status...")
try:
    import time

    for i in range(5):
        status = backend.get_job_status(job_id)
        print(f"  Attempt {i+1}/5: {status}")

        if status in ['DONE', 'ERROR', 'CANCELLED']:
            break

        time.sleep(2)

    print(f"✓ Job status check working")

except Exception as e:
    print(f"✗ Failed: {e}")

# Test 5: Test with VQE solver
print("\n[Test 5] Testing VQE integration with IBM backend...")
print("  Creating H2 molecule...")

try:
    from kanad.core.molecule import Molecule
    from kanad.core.atom import Atom
    from kanad.solvers.vqe_solver import VQESolver

    # Create H2 molecule
    h1 = Atom('H', position=(0.0, 0.0, 0.0))
    h2 = Atom('H', position=(0.74, 0.0, 0.0))

    molecule = Molecule(
        atoms=[h1, h2],
        charge=0,
        spin=0,
        basis='sto-3g'
    )

    print(f"✓ Molecule created: H2")
    print(f"  Electrons: {molecule.n_electrons}")
    print(f"  Orbitals: {molecule.n_orbitals}")

    # Get Hamiltonian from molecule (framework automatically creates the right type)
    hamiltonian = molecule.hamiltonian
    print(f"✓ Hamiltonian obtained from molecule")
    print(f"  Type: {type(hamiltonian).__name__}")

    # Create VQE solver with IBM backend
    print(f"\n  Creating VQE solver with IBM backend...")
    print(f"  WARNING: This will submit jobs to IBM Quantum!")
    print(f"  Using max_iterations=2 to minimize IBM usage")

    solver = VQESolver(
        hamiltonian=hamiltonian,
        ansatz_type='ucc',
        mapper_type='jordan_wigner',
        molecule=molecule,
        optimizer='SLSQP',
        max_iterations=2,  # Minimal iterations for testing
        backend='ibm',  # Pass string 'ibm', not the object
        shots=1024,
        enable_analysis=False,  # Disable to keep it simple
        # Pass IBM credentials via kwargs
        backend_name='ibm_brisbane',
        api_token=os.environ['IBM_API'],
        instance=os.environ['IBM_CRN']
    )

    print(f"✓ VQE solver created")
    print(f"  Ansatz: {solver.ansatz_type}")
    print(f"  Backend: {solver.backend}")
    print(f"  _use_statevector: {solver._use_statevector}")

    # Check if backend is properly set
    if solver._use_statevector:
        print(f"✗ ERROR: VQE solver is using statevector instead of IBM backend!")
        sys.exit(1)

    if not hasattr(solver, '_ibm_backend'):
        print(f"✗ ERROR: VQE solver does not have _ibm_backend attribute!")
        sys.exit(1)

    print(f"✓ IBM backend properly configured in VQE solver")
    print(f"  IBM Backend object: {solver._ibm_backend}")

    # Run ONE iteration to test
    print(f"\n  Running VQE with 2 iterations...")
    print(f"  This will submit 2 jobs to IBM Quantum")

    response = input("  Continue? (yes/no): ")
    if response.lower() != 'yes':
        print("  Test cancelled by user")
        sys.exit(0)

    result = solver.solve()

    print(f"\n✓ VQE completed!")
    print(f"  Energy: {result['energy']:.8f} Ha")
    print(f"  HF Energy: {result['hf_energy']:.8f} Ha")
    print(f"  Iterations: {result['iterations']}")
    print(f"  Converged: {result['converged']}")

except Exception as e:
    print(f"✗ Failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 80)
print("ALL TESTS PASSED!")
print("=" * 80)
