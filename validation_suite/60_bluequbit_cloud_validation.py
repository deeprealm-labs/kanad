#!/usr/bin/env python3
"""
BlueQubit Cloud Backend Validation Suite
Tests molecular calculations on BlueQubit cloud quantum computing platform.

This script runs real molecular calculations on BlueQubit's cloud infrastructure
using the free 34-qubit CPU simulator.
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
from kanad.backends.bluequbit_backend import BlueQubitBackend, get_bluequbit_backend

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

def check_bluequbit_setup():
    """Check if BlueQubit is properly set up."""
    print_header("BlueQubit Cloud Setup Check")

    # Check token
    token = os.getenv('BLUEQUBIT_API_TOKEN') or os.getenv('TOKEN')
    if not token:
        print("❌ ERROR: BlueQubit API token not found!")
        print("   Please set BLUEQUBIT_API_TOKEN or TOKEN in .env file")
        return False

    print(f"✅ API Token loaded: {token[:10]}... ({len(token)} chars)")

    # Check SDK
    try:
        import bluequbit
        print(f"✅ BlueQubit SDK installed: version {bluequbit.__version__}")
    except ImportError:
        print("❌ ERROR: BlueQubit SDK not installed!")
        print("   Run: pip install bluequbit")
        return False

    # Test backend initialization
    try:
        backend = get_bluequbit_backend(device='cpu')
        device_info = backend.get_device_info()
        print(f"✅ Backend initialized successfully")
        print(f"   Device: {device_info['current_device']}")
        current_info = device_info['device_info']
        print(f"   Qubits: {current_info.get('qubits', 'N/A')}")
        print(f"   Type: {current_info.get('type', 'N/A')}")
        print(f"   Cost: {current_info.get('cost', 'N/A')}")
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
    """Test H2 molecule on BlueQubit cloud."""
    print_header("TEST 1: H2 Molecule on BlueQubit Cloud (34-qubit CPU)")

    # Create H2 molecule
    H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
    h2_bond = BondFactory.create_bond(H1, H2)

    print(f"  Molecule: H2")
    print(f"  Bond length: 0.74 Å")
    print(f"  Electrons: {h2_bond.molecule.n_electrons}")

    # Get backend
    backend = get_bluequbit_backend(device='cpu')

    # Test with different ansätze
    mapper = JordanWignerMapper()

    results = {}

    # 1. UCCSD Ansatz
    print_subheader("1.1 Running VQE with UCCSD Ansatz on Cloud")
    ansatz_uccsd = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)

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
        print_result("UCCSD on BlueQubit Cloud", result_cloud['energy'], cloud_time, "✅")
    except Exception as e:
        print(f"❌ UCCSD Cloud execution failed: {e}")
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

    return results

def test_lih_molecule_cloud():
    """Test LiH molecule on BlueQubit cloud."""
    print_header("TEST 2: LiH Ionic Molecule on BlueQubit Cloud")

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

    # Get backend
    backend = get_bluequbit_backend(device='cpu')

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

    return results

def test_h2o_molecule_cloud():
    """Test H2O molecule on BlueQubit cloud (more complex system)."""
    print_header("TEST 3: H2O Molecule on BlueQubit Cloud (Complex System)")

    # Create H2O molecule
    O = Atom('O', position=np.array([0.0, 0.0, 0.0]))
    H1 = Atom('H', position=np.array([0.757, 0.586, 0.0]))
    H2 = Atom('H', position=np.array([-0.757, 0.586, 0.0]))

    from kanad.core.molecule import Molecule
    h2o = Molecule([O, H1, H2])

    print(f"  Molecule: H2O")
    print(f"  Atoms: O, H, H")
    print(f"  Electrons: {h2o.n_electrons}")
    print(f"  Charge: {h2o.charge}")

    # Build Hamiltonian
    from kanad.core.hamiltonian_builder import HamiltonianBuilder
    ham_builder = HamiltonianBuilder()
    hamiltonian = ham_builder.build(h2o)

    n_orbitals = hamiltonian.n_orbitals
    n_spin_orbitals = 2 * n_orbitals

    print(f"  Spin orbitals: {n_spin_orbitals}")

    # Check if within BlueQubit CPU limits (34 qubits)
    if n_spin_orbitals > 34:
        print(f"⚠️  WARNING: System requires {n_spin_orbitals} qubits, but BlueQubit CPU supports 34 qubits")
        print(f"   Skipping this test. Consider using active space reduction.")
        return None

    # Get backend
    backend = get_bluequbit_backend(device='cpu')

    mapper = JordanWignerMapper()
    results = {}

    # Use lightweight ansatz for complex system
    print_subheader("3.1 Running VQE with Hardware-Efficient Ansatz on Cloud")
    ansatz = RealAmplitudesAnsatz(n_qubits=n_spin_orbitals, n_electrons=h2o.n_electrons, n_layers=1)

    start_time = time.time()
    try:
        # Note: This might take longer due to system complexity
        print("  ⏳ This may take several minutes for H2O...")

        from kanad.solvers.vqe_solver import VQESolver
        vqe_solver = VQESolver(
            hamiltonian=hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            backend=backend
        )
        result_cloud = vqe_solver.solve(max_iterations=30)

        cloud_time = time.time() - start_time
        results['cloud'] = result_cloud['energy']
        results['cloud_time'] = cloud_time
        print_result("H2O on BlueQubit Cloud", result_cloud['energy'], cloud_time, "✅")
    except Exception as e:
        print(f"❌ Cloud execution failed: {e}")
        import traceback
        traceback.print_exc()
        results['cloud'] = None

    return results

def test_device_comparison():
    """Test same molecule on different BlueQubit devices."""
    print_header("TEST 4: Device Comparison (H2 on different backends)")

    # Create H2 molecule
    H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
    h2_bond = BondFactory.create_bond(H1, H2)

    mapper = JordanWignerMapper()
    ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)

    devices = ['cpu']  # Only free tier for now
    # If user has paid access, can test: ['cpu', 'gpu', 'mps.cpu']

    results = {}

    for device in devices:
        print_subheader(f"4.{devices.index(device)+1} Testing on {device.upper()} device")
        try:
            backend = get_bluequbit_backend(device=device)

            start_time = time.time()
            result = h2_bond.compute_energy(
                method='VQE',
                mapper=mapper,
                ansatz=ansatz,
                max_iterations=50,
                backend=backend
            )
            exec_time = time.time() - start_time

            results[device] = {
                'energy': result['energy'],
                'time': exec_time
            }
            print_result(f"{device.upper()} device", result['energy'], exec_time, "✅")
        except Exception as e:
            print(f"❌ {device.upper()} execution failed: {e}")
            results[device] = None

    return results

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print_header(f"BlueQubit Cloud Validation Suite - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Check setup
    if not check_bluequbit_setup():
        print("\n❌ Setup check failed. Please fix the issues above and try again.")
        sys.exit(1)

    print("\n✅ All checks passed! Starting molecular validation tests...")
    print("\n⚠️  NOTE: These tests will submit jobs to BlueQubit cloud and may take several minutes.")
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

        # Test 3: H2O (complex system - may take longer)
        print("\n⚠️  H2O test may take 5-10 minutes on cloud. Skip? (y/N): ", end='')
        skip_h2o = input().strip().lower() == 'y'
        if not skip_h2o:
            h2o_results = test_h2o_molecule_cloud()
            all_results['H2O'] = h2o_results
        else:
            print("  Skipping H2O test")

        # Test 4: Device comparison
        device_results = test_device_comparison()
        all_results['devices'] = device_results

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
        if molecule == 'devices':
            continue
        if results and any(v is not None for k, v in results.items() if 'cloud' in k):
            success_count += 1
            print(f"✅ {molecule}: Cloud execution successful")
        elif results:
            print(f"⚠️  {molecule}: Local execution only (cloud failed)")
        total_count += 1

    print(f"\nTests completed: {success_count}/{total_count} molecules ran successfully on BlueQubit cloud")
    print(f"\n{'='*80}")
    print("Validation suite complete!")
    print(f"{'='*80}\n")
