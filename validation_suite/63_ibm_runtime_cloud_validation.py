#!/usr/bin/env python3
"""
IBM Runtime Cloud Backend Validation Suite
Tests molecular calculations on IBM Quantum Runtime cloud platform.

This script runs real molecular calculations on IBM Quantum's cloud infrastructure
using the ibm_torino backend (133 qubits).
"""

import os
import sys
import time
import numpy as np
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
from kanad.backends.ibm_runtime_backend import IBMRuntimeBackend

# Load environment variables
from dotenv import load_dotenv
load_dotenv()

def print_header(title: str):
    """Print a formatted header."""
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")

def print_subheader(title: str):
    """Print a formatted subheader."""
    print(f"\n{'-'*80}")
    print(f"  {title}")
    print(f"{'-'*80}")

def print_result(name: str, energy: float, time_sec: float = None, status: str = "✅"):
    """Print formatted result."""
    time_str = f" | Time: {time_sec:.2f}s" if time_sec else ""
    print(f"{status} {name:50s} Energy: {energy:10.4f} eV{time_str}")

def print_comparison(name: str, cloud_energy: float, local_energy: float,
                    cloud_time: float = None, local_time: float = None):
    """Print comparison between cloud and local results."""
    error = abs(cloud_energy - local_energy)
    rel_error = (error / abs(local_energy)) * 100 if local_energy != 0 else 0
    status = "✅" if rel_error < 5 else "⚠️" if rel_error < 30 else "❌"

    time_str = ""
    if cloud_time and local_time:
        speedup = local_time / cloud_time if cloud_time > 0 else 0
        time_str = f" | Cloud: {cloud_time:.2f}s, Local: {local_time:.2f}s (Speedup: {speedup:.2f}x)"

    print(f"{status} {name:40s} Cloud: {cloud_energy:10.4f} eV | Local: {local_energy:10.4f} eV | Error: {rel_error:6.2f}%{time_str}")

def check_ibm_setup():
    """Check if IBM Runtime is properly set up."""
    print_header("IBM Runtime Cloud Setup Check")

    # Check token
    token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
    if not token:
        print("❌ ERROR: IBM Quantum token not found!")
        print("   Please set IBM_QUANTUM_TOKEN or API in .env file")
        return False

    print(f"✅ API Token loaded: {token[:10]}... ({len(token)} chars)")

    # Check CRN
    crn = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
    if crn:
        print(f"✅ CRN loaded: {crn[:20]}... ({len(crn)} chars)")
        channel = 'ibm_cloud'
    else:
        print("ℹ️  No CRN found, using IBM Quantum Platform")
        channel = 'ibm_quantum_platform'

    # Test backend initialization
    try:
        backend = IBMRuntimeBackend(token=token, instance=crn, channel=channel)
        backend_info = backend.get_backend_info()
        print(f"✅ Backend initialized successfully")
        print(f"   Backend: {backend_info['name']}")
        print(f"   Channel: {backend_info['channel']}")
        print(f"   Qubits: {backend_info.get('num_qubits', 'N/A')}")
        print(f"   Status: {backend_info.get('status', 'N/A')}")
        print(f"   Pending Jobs: {backend_info.get('pending_jobs', 'N/A')}")
        backend.close()
        return True
    except Exception as e:
        print(f"❌ ERROR: Failed to initialize backend: {e}")
        import traceback
        traceback.print_exc()
        return False

# ============================================================================
# MAIN VALIDATION TESTS
# ============================================================================

def test_h2_molecule_cloud():
    """Test H2 molecule on IBM Runtime cloud."""
    print_header("TEST 1: H2 Molecule on IBM Runtime Cloud")

    # Create H2 molecule
    H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
    h2_bond = BondFactory.create_bond(H1, H2)

    print(f"  Molecule: H2")
    print(f"  Bond length: 0.74 Å")
    print(f"  Electrons: {h2_bond.molecule.n_electrons}")

    # Get backend with optimized settings for real hardware
    token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
    crn = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
    channel = 'ibm_cloud' if crn else 'ibm_quantum_platform'

    backend = IBMRuntimeBackend(
        token=token,
        instance=crn,
        channel=channel,
        shots=4096,  # Increased for better statistics on real hardware
        optimization_level=3,  # Maximum transpiler optimization
        resilience_level=1  # Error mitigation enabled
    )

    # Log backend info
    backend_info = backend.get_backend_info()
    print(f"  Backend: {backend_info['name']} ({backend_info.get('num_qubits', 'N/A')} qubits)")
    print(f"  Shots: {backend.shots}, Resilience: {backend.resilience_level}")

    # Test with different ansätze
    mapper = JordanWignerMapper()

    results = {}

    # 1. UCCSD Ansatz
    print_subheader("1.1 Running VQE with UCCSD Ansatz on Cloud")
    ansatz_uccsd = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)

    # Show circuit info before execution
    print(f"  VQE iterations: 50")
    print(f"  Submitting to {backend_info['name']}...")

    start_time = time.time()
    try:
        result_cloud = h2_bond.compute_energy(
            method='VQE',
            mapper=mapper,
            ansatz=ansatz_uccsd,
            max_iterations=50,
            backend=backend
        )
        cloud_time = time.time() - start_time
        results['uccsd_cloud'] = result_cloud['energy']
        results['uccsd_cloud_time'] = cloud_time

        # Show additional details if available
        if 'iterations' in result_cloud:
            print(f"  Converged in {result_cloud['iterations']} iterations")

        print_result("UCCSD on IBM Runtime Cloud", result_cloud['energy'], cloud_time, "✅")
    except Exception as e:
        print(f"❌ UCCSD Cloud execution failed: {e}")
        import traceback
        traceback.print_exc()
        results['uccsd_cloud'] = None

    # 2. Local comparison
    print_subheader("1.2 Running VQE with UCCSD Ansatz Locally (for comparison)")
    start_time = time.time()
    result_local = h2_bond.compute_energy(
        method='VQE',
        mapper=mapper,
        ansatz=ansatz_uccsd,
        max_iterations=50
    )
    local_time = time.time() - start_time
    results['uccsd_local'] = result_local['energy']
    results['uccsd_local_time'] = local_time
    print_result("UCCSD Locally", result_local['energy'], local_time, "✅")

    # 3. Exact energy
    print_subheader("1.3 Computing Exact Energy")
    result_exact = h2_bond.compute_energy(method='exact', mapper=mapper)
    results['exact'] = result_exact['energy']
    print_result("Exact Energy", result_exact['energy'], status="✅")

    # Comparison
    print_subheader("1.4 Results Comparison")
    if results.get('uccsd_cloud'):
        print_comparison(
            "Cloud vs Local",
            results['uccsd_cloud'],
            results['uccsd_local'],
            results.get('uccsd_cloud_time'),
            results.get('uccsd_local_time')
        )
        error_cloud = abs(results['uccsd_cloud'] - results['exact'])
        error_local = abs(results['uccsd_local'] - results['exact'])
        print(f"  Cloud error from exact: {error_cloud:.6f} eV ({(error_cloud/abs(results['exact']))*100:.2f}%)")
        print(f"  Local error from exact: {error_local:.6f} eV ({(error_local/abs(results['exact']))*100:.2f}%)")

    backend.close()
    return results

def test_lih_molecule_cloud():
    """Test LiH molecule on IBM Runtime cloud."""
    print_header("TEST 2: LiH Ionic Molecule on IBM Runtime Cloud")

    # Create LiH molecule
    Li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    H = Atom('H', position=np.array([1.60, 0.0, 0.0]))
    lih_bond = BondFactory.create_bond(Li, H)

    # Get system info
    hamiltonian = lih_bond.hamiltonian
    n_orbitals = hamiltonian.n_orbitals
    n_spin_orbitals = 2 * n_orbitals
    n_electrons = lih_bond.molecule.n_electrons

    print(f"  Molecule: LiH")
    print(f"  Bond length: 1.60 Å")
    print(f"  Spin orbitals: {n_spin_orbitals}")
    print(f"  Electrons: {n_electrons}")

    # Get backend with optimized settings for real hardware
    token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
    crn = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
    channel = 'ibm_cloud' if crn else 'ibm_quantum_platform'

    backend = IBMRuntimeBackend(
        token=token,
        instance=crn,
        channel=channel,
        shots=4096,  # Increased for better statistics on real hardware
        optimization_level=3,  # Maximum transpiler optimization
        resilience_level=1  # Error mitigation enabled
    )

    # Log backend info
    backend_info = backend.get_backend_info()
    print(f"  Backend: {backend_info['name']} ({backend_info.get('num_qubits', 'N/A')} qubits)")
    print(f"  Shots: {backend.shots}, Resilience: {backend.resilience_level}")

    mapper = JordanWignerMapper()
    results = {}

    # 1. Hardware-efficient ansatz (lighter than UCCSD for larger systems)
    print_subheader("2.1 Running VQE with Hardware-Efficient Ansatz on Cloud")
    ansatz_he = RealAmplitudesAnsatz(n_qubits=n_spin_orbitals, n_electrons=n_electrons, n_layers=2)

    start_time = time.time()
    try:
        result_cloud = lih_bond.compute_energy(
            method='VQE',
            mapper=mapper,
            ansatz=ansatz_he,
            max_iterations=50,
            backend=backend
        )
        cloud_time = time.time() - start_time
        results['he_cloud'] = result_cloud['energy']
        results['he_cloud_time'] = cloud_time
        print_result("Hardware-Efficient on Cloud", result_cloud['energy'], cloud_time, "✅")
    except Exception as e:
        print(f"❌ Cloud execution failed: {e}")
        results['he_cloud'] = None

    # 2. Local comparison
    print_subheader("2.2 Running VQE Locally (for comparison)")
    start_time = time.time()
    result_local = lih_bond.compute_energy(
        method='VQE',
        mapper=mapper,
        ansatz=ansatz_he,
        max_iterations=50
    )
    local_time = time.time() - start_time
    results['he_local'] = result_local['energy']
    results['he_local_time'] = local_time
    print_result("Hardware-Efficient Locally", result_local['energy'], local_time, "✅")

    # 3. Exact energy
    print_subheader("2.3 Computing Exact Energy")
    result_exact = lih_bond.compute_energy(method='exact', mapper=mapper)
    results['exact'] = result_exact['energy']
    print_result("Exact Energy", result_exact['energy'], status="✅")

    # Comparison
    print_subheader("2.4 Results Comparison")
    if results.get('he_cloud'):
        print_comparison(
            "Cloud vs Local",
            results['he_cloud'],
            results['he_local'],
            results.get('he_cloud_time'),
            results.get('he_local_time')
        )
        error_cloud = abs(results['he_cloud'] - results['exact'])
        error_local = abs(results['he_local'] - results['exact'])
        print(f"  Cloud error from exact: {error_cloud:.6f} eV ({(error_cloud/abs(results['exact']))*100:.2f}%)")
        print(f"  Local error from exact: {error_local:.6f} eV ({(error_local/abs(results['exact']))*100:.2f}%)")

    backend.close()
    return results

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

        # Get backend
        token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
        crn = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
        channel = 'ibm_cloud' if crn else 'ibm_quantum_platform'

        backend = IBMRuntimeBackend(
            token=token,
            instance=crn,
            channel=channel,
            shots=4096,  # Increased for better statistics
            optimization_level=3
        )

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

        backend.close()
        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print_header(f"IBM Runtime Cloud Validation Suite - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Check setup
    if not check_ibm_setup():
        print("\n❌ Setup check failed. Please fix the issues above and try again.")
        sys.exit(1)

    print("\n✅ All checks passed! Starting molecular validation tests...")
    print("\n⚠️  NOTE: These tests will submit jobs to IBM Quantum cloud and may take several minutes.")
    print("         Make sure you have a valid API token with available credits.\n")

    all_results = {}

    # Run tests
    try:
        # Test 1: H2
        h2_results = test_h2_molecule_cloud()
        all_results['H2'] = h2_results

        # Test 2: LiH
        lih_results = test_lih_molecule_cloud()
        all_results['LiH'] = lih_results

        # Test 3: Simple circuit execution
        circuit_result = test_simple_circuit_execution()
        all_results['circuit'] = circuit_result

    except KeyboardInterrupt:
        print("\n\n⚠️  Tests interrupted by user")
    except Exception as e:
        print(f"\n\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()

    # Final summary
    print_header("VALIDATION SUMMARY")

    success_count = 0
    total_count = 0

    for molecule, results in all_results.items():
        if molecule == 'circuit':
            if results:
                success_count += 1
                print(f"✅ Circuit execution: Cloud execution successful")
            else:
                print(f"❌ Circuit execution: Cloud execution failed")
            total_count += 1
        else:
            if results and any(v is not None for k, v in results.items() if 'cloud' in k):
                success_count += 1
                print(f"✅ {molecule}: Cloud execution successful")
            elif results:
                print(f"⚠️  {molecule}: Local execution only (cloud failed)")
            total_count += 1

    print(f"\nTests completed: {success_count}/{total_count} tests ran successfully on IBM Runtime cloud")
    print(f"\n{'='*80}")
    print("Validation suite complete!")
    print(f"{'='*80}\n")
