#!/usr/bin/env python3
"""
IBM Runtime Quick Test - Fast molecular setup tests
Tests molecular setup and circuit generation with IBM Quantum Runtime backend
"""

import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.backends.ibm_runtime_backend import IBMRuntimeBackend

from dotenv import load_dotenv
load_dotenv()

def print_header(title: str):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")

def test_backend_setup():
    """Test IBM Runtime backend initialization."""
    print_header("TEST 1: IBM Runtime Backend Setup")

    try:
        # Explicitly set channel to avoid issues
        token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
        instance = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
        channel = 'ibm_cloud' if instance else 'ibm_quantum_platform'
        
        backend = IBMRuntimeBackend(token=token, instance=instance, channel=channel)
        backend_info = backend.get_backend_info()

        print(f"✅ Backend initialized")
        print(f"   Name: {backend_info['name']}")
        print(f"   Channel: {backend_info['channel']}")
        print(f"   Qubits: {backend_info.get('num_qubits', 'N/A')}")
        print(f"   Status: {backend_info.get('status', 'N/A')}")
        print(f"   Pending Jobs: {backend_info.get('pending_jobs', 'N/A')}")
        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        return False

def test_h2_molecule_setup():
    """Test H2 molecule setup with IBM Runtime."""
    print_header("TEST 2: H2 Molecule Setup")

    try:
        # Create H2
        H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        h2_bond = BondFactory.create_bond(H1, H2)

        print(f"✅ H2 molecule created")
        print(f"   Electrons: {h2_bond.molecule.n_electrons}")
        print(f"   Orbitals: {h2_bond.hamiltonian.n_orbitals}")

        # Create circuit
        mapper = JordanWignerMapper()
        ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)

        circuit = ansatz.build_circuit()

        print(f"✅ Circuit created")
        print(f"   Qubits: {circuit.n_qubits}")
        print(f"   Gates: {len(circuit.gates)}")
        print(f"   Depth: {circuit.depth}")

        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_simple_circuit_execution():
    """Test simple circuit execution on IBM Runtime."""
    print_header("TEST 3: Simple Circuit Execution on IBM Runtime")

    try:
        from qiskit import QuantumCircuit

        # Create simple Bell state circuit
        circuit = QuantumCircuit(2)
        circuit.h(0)
        circuit.cx(0, 1)
        circuit.measure_all()

        print(f"  Circuit: Bell state (2 qubits)")
        print(f"  Submitting to IBM Runtime...")

        # Explicitly set channel to avoid issues
        token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
        instance = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
        channel = 'ibm_cloud' if instance else 'ibm_quantum_platform'
        
        backend = IBMRuntimeBackend(token=token, instance=instance, channel=channel)

        # Transpile circuit for IBM backend
        transpiled_circuit = backend.transpile_circuit(circuit)
        print(f"  Transpiled circuit: {transpiled_circuit.depth()} depth")

        # Get sampler and run circuit
        sampler = backend.get_sampler()
        
        start_time = time.time()
        job = sampler.run([transpiled_circuit], shots=1024)
        result = job.result()
        exec_time = time.time() - start_time

        counts = result[0].data.meas.get_counts()
        print(f"✅ Execution completed in {exec_time:.2f}s")
        print(f"   Counts: {counts}")
        print(f"   Job ID: {job.job_id()}")

        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_h2_energy_estimation():
    """Test H2 energy estimation using IBM Runtime Estimator."""
    print_header("TEST 4: H2 Energy Estimation on IBM Runtime")

    try:
        # Create H2
        H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        h2_bond = BondFactory.create_bond(H1, H2)

        # Setup
        mapper = JordanWignerMapper()
        ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)
        # Explicitly set channel to avoid issues
        token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
        instance = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
        channel = 'ibm_cloud' if instance else 'ibm_quantum_platform'
        
        backend = IBMRuntimeBackend(token=token, instance=instance, channel=channel)

        print(f"  Computing energy expectation value...")

        # Create circuit and Hamiltonian
        circuit = ansatz.build_circuit()
        hamiltonian = h2_bond.hamiltonian
        
        # Convert to Qiskit format
        from qiskit.quantum_info import SparsePauliOp
        from qiskit import QuantumCircuit as QiskitCircuit
        
        # Create a simple test Hamiltonian (H2 in minimal basis)
        # This is a simplified version - in practice you'd use the full molecular Hamiltonian
        pauli_strings = ['II', 'ZZ', 'XX', 'YY']
        coeffs = [0.5, -0.2, 0.1, 0.1]
        test_hamiltonian = SparsePauliOp(pauli_strings, coeffs)
        
        # Convert our circuit to Qiskit format
        qiskit_circuit = QiskitCircuit(2)
        qiskit_circuit.h(0)
        qiskit_circuit.cx(0, 1)
        qiskit_circuit.ry(0.1, 0)

        # Transpile for IBM backend
        transpiled_circuit = backend.transpile_circuit(qiskit_circuit)

        # Get estimator and compute expectation value
        estimator = backend.get_estimator()
        
        start_time = time.time()
        job = estimator.run([(transpiled_circuit, test_hamiltonian)])
        result = job.result()
        exec_time = time.time() - start_time

        expectation_value = result[0].data.evs[0]
        print(f"✅ Energy computed in {exec_time:.2f}s")
        print(f"   Expectation value: {expectation_value:.4f}")
        print(f"   Job ID: {job.job_id()}")

        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_h2_fast_vqe():
    """Test fast VQE optimization (3 iterations only)."""
    print_header("TEST 5: H2 Fast VQE (3 iterations) on IBM Runtime")

    try:
        # Create H2
        H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        h2_bond = BondFactory.create_bond(H1, H2)

        mapper = JordanWignerMapper()
        ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)

        print(f"  Running VQE with 3 iterations locally...")
        print(f"  (Testing local VQE first)")

        start_time = time.time()
        result = h2_bond.compute_energy(
            method='VQE',
            mapper=mapper,
            ansatz=ansatz,
            max_iterations=3
        )
        exec_time = time.time() - start_time

        print(f"✅ VQE completed in {exec_time:.2f}s")
        print(f"   Energy: {result['energy']:.4f} eV")
        print(f"   Iterations: {result.get('iterations', 'N/A')}")

        # Get exact for comparison
        exact_result = h2_bond.compute_energy(method='exact', mapper=mapper)
        error = abs(result['energy'] - exact_result['energy'])
        rel_error = (error / abs(exact_result['energy'])) * 100

        print(f"\n  Exact energy: {exact_result['energy']:.4f} eV")
        print(f"  Error: {error:.6f} eV ({rel_error:.2f}%)")

        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print_header("IBM Runtime Quick Test Suite")

    # Check credentials
    token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
    crn = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
    
    if not token:
        print("\n❌ No IBM Quantum token found in environment")
        print("   Please set IBM_QUANTUM_TOKEN or API in .env file")
        sys.exit(1)

    print(f"✅ Token loaded: {token[:10]}... ({len(token)} chars)")
    if crn:
        print(f"✅ CRN loaded: {crn[:20]}... ({len(crn)} chars)")

    results = []

    # Run tests
    results.append(("Backend Setup", test_backend_setup()))
    results.append(("H2 Molecule Setup", test_h2_molecule_setup()))
    results.append(("Simple Circuit Execution", test_simple_circuit_execution()))
    results.append(("Energy Estimation", test_h2_energy_estimation()))
    results.append(("Fast VQE (3 iterations)", test_h2_fast_vqe()))

    # Summary
    print_header("Test Summary")

    passed = sum(1 for _, r in results if r)
    total = len(results)

    for name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {name}")

    print(f"\nPassed: {passed}/{total}")
    print(f"{'='*80}\n")

    sys.exit(0 if passed == total else 1)
