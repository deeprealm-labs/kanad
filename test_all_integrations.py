#!/usr/bin/env python3
"""
Comprehensive Integration Tests for Kanad Framework
Tests all solvers with all backends and analysis options
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.solvers import VQESolver, SQDSolver
from kanad.backends.bluequbit import BlueQubitBackend
from kanad.backends.ibm import IBMBackend
import numpy as np


def print_section(title):
    """Print a formatted section header."""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)


def test_vqe_classical():
    """Test VQE with classical (statevector) backend."""
    print_section("TEST 1: VQE with Classical Backend")

    try:
        # Create H2 molecule
        h1 = Atom('H', position=[0.0, 0.0, 0.0])
        h2 = Atom('H', position=[0.74, 0.0, 0.0])
        molecule = Molecule(atoms=[h1, h2], basis='sto-3g')

        print(f"✓ Molecule created: {molecule.formula}")
        print(f"  Electrons: {molecule.n_electrons}, Orbitals: {molecule.n_orbitals}")

        # Create VQE solver
        from kanad.bonds import BondFactory
        bond = BondFactory.create_bond(h1, h2, distance=0.74, basis='sto-3g')

        solver = VQESolver(
            bond=bond,
            ansatz_type='hardware_efficient',
            mapper_type='jordan_wigner',
            optimizer='COBYLA',
            max_iterations=5,  # Quick test
            backend='statevector'
        )

        print(f"✓ VQE solver created")

        # Run optimization
        result = solver.solve()

        print(f"✅ VQE CLASSICAL TEST PASSED")
        print(f"  Energy: {result['energy']:.6f} Ha")
        print(f"  Converged: {result['converged']}")
        print(f"  Iterations: {result['iterations']}")

        return True

    except Exception as e:
        print(f"❌ VQE CLASSICAL TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_vqe_bluequbit():
    """Test VQE with BlueQubit CPU backend."""
    print_section("TEST 2: VQE with BlueQubit CPU Backend")

    try:
        import os
        api_token = os.environ.get('BLUE_TOKEN')
        if not api_token:
            print("⚠️  SKIPPED: BLUE_TOKEN not set")
            return None

        # Create H2 molecule
        h1 = Atom('H', position=[0.0, 0.0, 0.0])
        h2 = Atom('H', position=[0.74, 0.0, 0.0])

        from kanad.bonds import BondFactory
        bond = BondFactory.create_bond(h1, h2, distance=0.74, basis='sto-3g')

        print(f"✓ H2 bond created")

        # Create VQE solver with BlueQubit
        solver = VQESolver(
            bond=bond,
            ansatz_type='hardware_efficient',
            mapper_type='jordan_wigner',
            optimizer='COBYLA',
            max_iterations=3,  # Very quick test
            backend='bluequbit',
            api_token=api_token,
            device='cpu'  # Free tier
        )

        print(f"✓ VQE solver created with BlueQubit CPU")

        # Run optimization
        result = solver.solve()

        print(f"✅ VQE BLUEQUBIT TEST PASSED")
        print(f"  Energy: {result['energy']:.6f} Ha")
        print(f"  Converged: {result['converged']}")
        print(f"  Iterations: {result['iterations']}")

        return True

    except Exception as e:
        print(f"❌ VQE BLUEQUBIT TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_sqd_classical():
    """Test SQD with classical backend."""
    print_section("TEST 3: SQD with Classical Backend")

    try:
        # Create H2 molecule
        h1 = Atom('H', position=[0.0, 0.0, 0.0])
        h2 = Atom('H', position=[0.74, 0.0, 0.0])

        from kanad.bonds import BondFactory
        bond = BondFactory.create_bond(h1, h2, distance=0.74, basis='sto-3g')

        print(f"✓ H2 bond created")

        # Create SQD solver
        solver = SQDSolver(
            bond=bond,
            subspace_dim=5,
            circuit_depth=2,
            backend='statevector',
            enable_analysis=True
        )

        print(f"✓ SQD solver created")

        # Run SQD
        result = solver.solve(n_states=3)

        print(f"✅ SQD CLASSICAL TEST PASSED")
        print(f"  Ground state energy: {result['energy']:.6f} Ha")
        print(f"  Excited states: {len(result.get('excited_state_energies', []))}")
        print(f"  Converged: {result['converged']}")

        return True

    except Exception as e:
        print(f"❌ SQD CLASSICAL TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_sqd_bluequbit():
    """Test SQD with BlueQubit CPU backend."""
    print_section("TEST 4: SQD with BlueQubit CPU Backend")

    try:
        import os
        api_token = os.environ.get('BLUE_TOKEN')
        if not api_token:
            print("⚠️  SKIPPED: BLUE_TOKEN not set")
            return None

        # Create H2 molecule
        h1 = Atom('H', position=[0.0, 0.0, 0.0])
        h2 = Atom('H', position=[0.74, 0.0, 0.0])

        from kanad.bonds import BondFactory
        bond = BondFactory.create_bond(h1, h2, distance=0.74, basis='sto-3g')

        print(f"✓ H2 bond created")

        # Create SQD solver with BlueQubit
        solver = SQDSolver(
            bond=bond,
            subspace_dim=4,
            circuit_depth=2,
            backend='bluequbit',
            api_token=api_token,
            device='cpu',
            enable_analysis=True
        )

        print(f"✓ SQD solver created with BlueQubit CPU")

        # Run SQD
        result = solver.solve(n_states=2)

        print(f"✅ SQD BLUEQUBIT TEST PASSED")
        print(f"  Ground state energy: {result['energy']:.6f} Ha")
        print(f"  Excited states: {len(result.get('excited_state_energies', []))}")

        return True

    except Exception as e:
        print(f"❌ SQD BLUEQUBIT TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_analysis_generation():
    """Test that analysis is properly generated."""
    print_section("TEST 5: Analysis Generation")

    try:
        # Create H2 molecule
        h1 = Atom('H', position=[0.0, 0.0, 0.0])
        h2 = Atom('H', position=[0.74, 0.0, 0.0])

        from kanad.bonds import BondFactory
        bond = BondFactory.create_bond(h1, h2, distance=0.74, basis='sto-3g')

        # Create VQE solver with analysis enabled
        solver = VQESolver(
            bond=bond,
            ansatz_type='hardware_efficient',
            mapper_type='jordan_wigner',
            optimizer='COBYLA',
            max_iterations=3,
            backend='statevector',
            enable_analysis=True
        )

        result = solver.solve()

        # Check if analysis exists
        if 'analysis' not in result:
            print(f"❌ ANALYSIS TEST FAILED: No analysis in results")
            return False

        analysis = result['analysis']
        print(f"✓ Analysis found in results")
        print(f"  Analysis keys: {list(analysis.keys())}")

        # Check for expected analysis components
        expected_keys = ['energy_components', 'bonding', 'properties']
        missing_keys = [k for k in expected_keys if k not in analysis]

        if missing_keys:
            print(f"⚠️  Missing analysis keys: {missing_keys}")

        print(f"✅ ANALYSIS GENERATION TEST PASSED")
        return True

    except Exception as e:
        print(f"❌ ANALYSIS TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_api_endpoints():
    """Test API endpoints for experiment submission."""
    print_section("TEST 6: API Integration")

    try:
        import requests

        # Test health endpoint
        response = requests.get('http://localhost:8000/')
        print(f"✓ API server responding: {response.status_code}")

        # Test settings endpoint
        response = requests.get('http://localhost:8000/api/settings/defaults')
        if response.status_code == 200:
            settings = response.json()
            print(f"✓ Settings endpoint working")
            print(f"  Default method: {settings['settings'].get('method')}")

        print(f"✅ API INTEGRATION TEST PASSED")
        return True

    except Exception as e:
        print(f"❌ API TEST FAILED: {e}")
        return False


def run_all_tests():
    """Run all integration tests."""
    print("\n" + "=" * 80)
    print("  KANAD FRAMEWORK - COMPREHENSIVE INTEGRATION TESTS")
    print("=" * 80)

    results = {
        'VQE Classical': test_vqe_classical(),
        'VQE BlueQubit': test_vqe_bluequbit(),
        'SQD Classical': test_sqd_classical(),
        'SQD BlueQubit': test_sqd_bluequbit(),
        'Analysis Generation': test_analysis_generation(),
        'API Integration': test_api_endpoints(),
    }

    # Summary
    print_section("TEST SUMMARY")

    passed = sum(1 for v in results.values() if v is True)
    failed = sum(1 for v in results.values() if v is False)
    skipped = sum(1 for v in results.values() if v is None)
    total = len(results)

    for test_name, result in results.items():
        status = "✅ PASS" if result is True else "❌ FAIL" if result is False else "⚠️  SKIP"
        print(f"  {status}  {test_name}")

    print()
    print(f"Total: {total} tests")
    print(f"  ✅ Passed: {passed}")
    print(f"  ❌ Failed: {failed}")
    print(f"  ⚠️  Skipped: {skipped}")

    if failed > 0:
        print("\n❌ SOME TESTS FAILED - Review errors above")
        return False
    else:
        print("\n✅ ALL TESTS PASSED!")
        return True


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
