#!/usr/bin/env python3
"""
BlueQubit Quick Test - Fast molecular setup tests
Tests molecular setup and circuit generation without full VQE optimization
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
from kanad.backends.bluequbit_backend import get_bluequbit_backend

from dotenv import load_dotenv
load_dotenv()

def print_header(title: str):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")

def test_backend_setup():
    """Test BlueQubit backend initialization."""
    print_header("TEST 1: Backend Setup")

    try:
        backend = get_bluequbit_backend(device='cpu')
        device_info = backend.get_device_info()

        print(f"✅ Backend initialized")
        print(f"   Device: {device_info['current_device']}")
        print(f"   Qubits: {device_info['device_info'].get('qubits', 'N/A')}")
        print(f"   Type: {device_info['device_info'].get('type', 'N/A')}")
        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        return False

def test_h2_molecule_setup():
    """Test H2 molecule setup with BlueQubit."""
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

        from qiskit import QuantumCircuit
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
    """Test simple circuit execution on BlueQubit."""
    print_header("TEST 3: Simple Circuit Execution on Cloud")

    try:
        from qiskit import QuantumCircuit

        # Create simple Bell state circuit
        circuit = QuantumCircuit(2)
        circuit.h(0)
        circuit.cx(0, 1)
        circuit.measure_all()

        print(f"  Circuit: Bell state (2 qubits)")
        print(f"  Submitting to BlueQubit cloud...")

        backend = get_bluequbit_backend(device='cpu')

        start_time = time.time()
        result = backend.run_circuit(circuit, shots=1024)
        exec_time = time.time() - start_time

        print(f"✅ Execution completed in {exec_time:.2f}s")
        print(f"   Counts: {result['counts']}")
        print(f"   Job ID: {result.get('job_id', 'N/A')}")

        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_h2_single_vqe_step():
    """Test single VQE energy evaluation on cloud."""
    print_header("TEST 4: H2 Single Energy Evaluation on Cloud")

    try:
        # Create H2
        H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        h2_bond = BondFactory.create_bond(H1, H2)

        # Setup
        mapper = JordanWignerMapper()
        ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)
        backend = get_bluequbit_backend(device='cpu')

        print(f"  Evaluating energy with initial parameters...")

        # Create a simple test circuit
        from qiskit import QuantumCircuit
        test_circuit = QuantumCircuit(2)
        test_circuit.h(0)
        test_circuit.cx(0, 1)
        test_circuit.ry(0.1, 0)
        test_circuit.measure_all()

        start_time = time.time()
        result = backend.run_circuit(test_circuit, shots=1024)
        exec_time = time.time() - start_time

        print(f"✅ Circuit executed in {exec_time:.2f}s")
        print(f"   Counts: {result['counts']}")
        print(f"   Job ID: {result.get('job_id', 'N/A')}")

        # Compare with local
        print(f"\n  Running same calculation locally...")
        from kanad.solvers.vqe_solver import VQESolver

        vqe_local = VQESolver(
            hamiltonian=h2_bond.hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            backend='classical'
        )

        params = np.zeros(ansatz.get_num_parameters())
        start_time = time.time()
        energy_local = vqe_local.compute_energy(params)
        local_time = time.time() - start_time

        print(f"✅ Local energy: {energy_local:.4f} eV ({local_time:.2f}s)")

        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_h2_fast_vqe():
    """Test fast VQE optimization (3 iterations only)."""
    print_header("TEST 5: H2 Fast VQE (3 iterations) on Cloud")

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
    print_header("BlueQubit Quick Test Suite")

    token = os.getenv('BLUEQUBIT_API_TOKEN') or os.getenv('TOKEN')
    if not token:
        print("\n❌ No BlueQubit token found in environment")
        sys.exit(1)

    print(f"✅ Token loaded: {token[:10]}... ({len(token)} chars)")

    results = []

    # Run tests
    results.append(("Backend Setup", test_backend_setup()))
    results.append(("H2 Molecule Setup", test_h2_molecule_setup()))
    results.append(("Simple Circuit Execution", test_simple_circuit_execution()))
    results.append(("Single Energy Evaluation", test_h2_single_vqe_step()))
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
