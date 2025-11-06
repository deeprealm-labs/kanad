"""
Comprehensive Framework Testing Suite
Tests all core functionalities of the Kanad framework
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
from pyscf import gto
import numpy as np
from typing import Dict, List

# Import Kanad components
from kanad.core.hamiltonians.openfermion_jw import create_hamiltonian
from kanad.ansatze.governance_aware_ansatz import GovernanceAwareAnsatz
from kanad.ansatze.two_local_ansatz import TwoLocalAnsatz
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.utils.vqe_solver import VQESolver
from kanad.visualization.circuit_visualizer import create_preview_circuit


class IterationTracker:
    """Track iterations during VQE optimization"""
    def __init__(self):
        self.iterations = []
        self.energies = []
        self.parameters = []

    def callback(self, iteration, energy, parameters):
        self.iterations.append(iteration)
        self.energies.append(energy)
        self.parameters.append(parameters.copy() if hasattr(parameters, 'copy') else list(parameters))
        print(f"  Iteration {iteration}: Energy = {energy:.6f} Ha")

    def summary(self):
        return {
            'total_iterations': len(self.iterations),
            'initial_energy': self.energies[0] if self.energies else None,
            'final_energy': self.energies[-1] if self.energies else None,
            'converged': len(self.energies) > 1 and abs(self.energies[-1] - self.energies[-2]) < 1e-6
        }


def create_test_molecule(formula: str, basis: str = 'sto-3g'):
    """Create a test molecule"""
    geometries = {
        'H2': 'H 0 0 0; H 0 0 0.74',
        'LiH': 'Li 0 0 0; H 0 0 1.6',
        'H2O': 'O 0 0 0; H 0.757 0.586 0; H -0.757 0.586 0',
        'HeH+': 'He 0 0 0; H 0 0 0.9772'
    }

    charge = 1 if formula == 'HeH+' else 0
    mol = gto.M(atom=geometries[formula], basis=basis, charge=charge)
    return mol


def test_hamiltonian_generation():
    """Test 1: Hamiltonian Generation"""
    print("\n" + "="*60)
    print("TEST 1: Hamiltonian Generation")
    print("="*60)

    mol = create_test_molecule('H2')
    print(f"Molecule: H2, Basis: sto-3g")
    print(f"Number of electrons: {mol.nelectron}")
    print(f"Number of orbitals: {mol.nao_nr()}")

    hamiltonian, n_qubits = create_hamiltonian(mol)
    print(f"✓ Hamiltonian created successfully")
    print(f"  Number of qubits: {n_qubits}")
    print(f"  Number of terms: {len(hamiltonian.terms)}")

    assert n_qubits > 0, "Number of qubits should be positive"
    assert len(hamiltonian.terms) > 0, "Hamiltonian should have terms"

    return True


def test_vqe_with_governance_ansatz():
    """Test 2: VQE with Governance Ansatz"""
    print("\n" + "="*60)
    print("TEST 2: VQE with Governance Aware Ansatz")
    print("="*60)

    mol = create_test_molecule('H2')
    print(f"Molecule: H2")

    # Create Hamiltonian
    hamiltonian, n_qubits = create_hamiltonian(mol)
    print(f"Hamiltonian: {n_qubits} qubits")

    # Create ansatz
    ansatz = GovernanceAwareAnsatz(n_qubits=n_qubits, n_electrons=mol.nelectron)
    print(f"Ansatz: GovernanceAwareAnsatz")
    print(f"  Parameters: {ansatz.n_parameters}")

    # Run VQE with iteration tracking
    tracker = IterationTracker()
    solver = VQESolver(hamiltonian, ansatz)

    print("\nRunning VQE optimization...")
    result = solver.solve(
        maxiter=50,
        callback=tracker.callback
    )

    summary = tracker.summary()
    print(f"\n✓ VQE completed")
    print(f"  Total iterations: {summary['total_iterations']}")
    print(f"  Initial energy: {summary['initial_energy']:.6f} Ha")
    print(f"  Final energy: {summary['final_energy']:.6f} Ha")
    print(f"  Converged: {summary['converged']}")

    # Add nuclear repulsion
    nuclear_repulsion = mol.energy_nuc()
    total_energy = result['energy'] + nuclear_repulsion
    print(f"  Nuclear repulsion: {nuclear_repulsion:.6f} Ha")
    print(f"  Total energy: {total_energy:.6f} Ha")

    assert summary['total_iterations'] > 0, "Should have iterations"
    assert result['energy'] < 0, "Ground state energy should be negative"

    return True


def test_vqe_with_two_local_ansatz():
    """Test 3: VQE with TwoLocal Ansatz"""
    print("\n" + "="*60)
    print("TEST 3: VQE with TwoLocal Ansatz")
    print("="*60)

    mol = create_test_molecule('H2')
    print(f"Molecule: H2")

    hamiltonian, n_qubits = create_hamiltonian(mol)
    print(f"Hamiltonian: {n_qubits} qubits")

    ansatz = TwoLocalAnsatz(n_qubits=n_qubits, reps=2)
    print(f"Ansatz: TwoLocalAnsatz (reps=2)")
    print(f"  Parameters: {ansatz.n_parameters}")

    tracker = IterationTracker()
    solver = VQESolver(hamiltonian, ansatz)

    print("\nRunning VQE optimization...")
    result = solver.solve(
        maxiter=50,
        callback=tracker.callback
    )

    summary = tracker.summary()
    print(f"\n✓ VQE completed")
    print(f"  Total iterations: {summary['total_iterations']}")
    print(f"  Final energy: {summary['final_energy']:.6f} Ha")

    assert summary['total_iterations'] > 0, "Should have iterations"

    return True


def test_vqe_with_ucc_ansatz():
    """Test 4: VQE with UCC Ansatz"""
    print("\n" + "="*60)
    print("TEST 4: VQE with UCC Ansatz")
    print("="*60)

    mol = create_test_molecule('H2')
    print(f"Molecule: H2")

    hamiltonian, n_qubits = create_hamiltonian(mol)
    print(f"Hamiltonian: {n_qubits} qubits")

    ansatz = UCCAnsatz(n_qubits=n_qubits, n_electrons=mol.nelectron)
    print(f"Ansatz: UCCAnsatz")
    print(f"  Parameters: {ansatz.n_parameters}")

    tracker = IterationTracker()
    solver = VQESolver(hamiltonian, ansatz)

    print("\nRunning VQE optimization...")
    result = solver.solve(
        maxiter=50,
        callback=tracker.callback
    )

    summary = tracker.summary()
    print(f"\n✓ VQE completed")
    print(f"  Total iterations: {summary['total_iterations']}")
    print(f"  Final energy: {summary['final_energy']:.6f} Ha")

    assert summary['total_iterations'] > 0, "Should have iterations"

    return True


def test_different_basis_sets():
    """Test 5: Different Basis Sets"""
    print("\n" + "="*60)
    print("TEST 5: Different Basis Sets")
    print("="*60)

    basis_sets = ['sto-3g', '6-31g']
    results = {}

    for basis in basis_sets:
        print(f"\nTesting basis: {basis}")
        mol = create_test_molecule('H2', basis=basis)
        print(f"  Electrons: {mol.nelectron}")
        print(f"  Orbitals: {mol.nao_nr()}")

        hamiltonian, n_qubits = create_hamiltonian(mol)
        print(f"  Qubits: {n_qubits}")

        ansatz = GovernanceAwareAnsatz(n_qubits=n_qubits, n_electrons=mol.nelectron)
        solver = VQESolver(hamiltonian, ansatz)

        result = solver.solve(maxiter=30)
        nuclear_repulsion = mol.energy_nuc()
        total_energy = result['energy'] + nuclear_repulsion

        results[basis] = total_energy
        print(f"  ✓ Energy: {total_energy:.6f} Ha")

    print(f"\n✓ All basis sets tested successfully")

    return True


def test_ionic_molecule():
    """Test 6: Ionic Molecule (HeH+)"""
    print("\n" + "="*60)
    print("TEST 6: Ionic Molecule (HeH+)")
    print("="*60)

    mol = create_test_molecule('HeH+')
    print(f"Molecule: HeH+ (charge=+1)")
    print(f"Electrons: {mol.nelectron}")

    hamiltonian, n_qubits = create_hamiltonian(mol)
    print(f"Qubits: {n_qubits}")

    ansatz = GovernanceAwareAnsatz(n_qubits=n_qubits, n_electrons=mol.nelectron)
    tracker = IterationTracker()
    solver = VQESolver(hamiltonian, ansatz)

    print("\nRunning VQE optimization...")
    result = solver.solve(maxiter=50, callback=tracker.callback)

    summary = tracker.summary()
    nuclear_repulsion = mol.energy_nuc()
    total_energy = result['energy'] + nuclear_repulsion

    print(f"\n✓ VQE completed")
    print(f"  Iterations: {summary['total_iterations']}")
    print(f"  Total energy: {total_energy:.6f} Ha")

    assert result['energy'] < 0, "Energy should be negative"

    return True


def test_circuit_visualization():
    """Test 7: Circuit Visualization"""
    print("\n" + "="*60)
    print("TEST 7: Circuit Visualization")
    print("="*60)

    # Test circuit preview with H2
    molecule_config = {
        'atoms': [
            {'symbol': 'H', 'x': 0, 'y': 0, 'z': 0},
            {'symbol': 'H', 'x': 0, 'y': 0, 'z': 0.74}
        ],
        'basis': 'sto-3g',
        'charge': 0,
        'multiplicity': 1
    }

    backend_config = {
        'method': 'VQE',
        'ansatz': 'hardware_efficient',
        'mapper': 'jordan_wigner'
    }

    print("Creating circuit preview...")
    result = create_preview_circuit(molecule_config, backend_config)

    print(f"✓ Circuit preview created")
    print(f"  Qubits: {result['n_qubits']}")
    print(f"  Parameters: {result['n_parameters']}")
    print(f"  Depth: {result['depth']}")
    print(f"  ASCII diagram length: {len(result['ascii_diagram'])} chars")

    assert result['n_qubits'] > 0, "Should have qubits"
    assert result['ascii_diagram'], "Should have ASCII diagram"

    return True


def test_convergence_behavior():
    """Test 8: Convergence Behavior"""
    print("\n" + "="*60)
    print("TEST 8: Convergence Behavior Analysis")
    print("="*60)

    mol = create_test_molecule('LiH')
    print(f"Molecule: LiH")

    hamiltonian, n_qubits = create_hamiltonian(mol)
    ansatz = GovernanceAwareAnsatz(n_qubits=n_qubits, n_electrons=mol.nelectron)

    tracker = IterationTracker()
    solver = VQESolver(hamiltonian, ansatz)

    print("\nRunning VQE with convergence tracking...")
    result = solver.solve(maxiter=100, callback=tracker.callback)

    # Analyze convergence
    energies = tracker.energies
    if len(energies) > 2:
        improvements = [abs(energies[i] - energies[i-1]) for i in range(1, len(energies))]
        avg_improvement = np.mean(improvements)

        print(f"\n✓ Convergence analysis:")
        print(f"  Total iterations: {len(energies)}")
        print(f"  Energy range: {max(energies):.6f} to {min(energies):.6f} Ha")
        print(f"  Average improvement: {avg_improvement:.8f} Ha")
        print(f"  Monotonic decrease: {energies == sorted(energies, reverse=True)}")

    return True


def run_all_tests():
    """Run all framework tests"""
    print("\n" + "="*80)
    print(" KANAD FRAMEWORK COMPREHENSIVE TEST SUITE")
    print("="*80)

    tests = [
        ("Hamiltonian Generation", test_hamiltonian_generation),
        ("VQE with Governance Ansatz", test_vqe_with_governance_ansatz),
        ("VQE with TwoLocal Ansatz", test_vqe_with_two_local_ansatz),
        ("VQE with UCC Ansatz", test_vqe_with_ucc_ansatz),
        ("Different Basis Sets", test_different_basis_sets),
        ("Ionic Molecule", test_ionic_molecule),
        ("Circuit Visualization", test_circuit_visualization),
        ("Convergence Behavior", test_convergence_behavior),
    ]

    results = {}
    passed = 0
    failed = 0

    for name, test_func in tests:
        try:
            success = test_func()
            results[name] = "PASSED"
            passed += 1
        except Exception as e:
            results[name] = f"FAILED: {str(e)}"
            failed += 1
            print(f"\n❌ Test failed: {str(e)}")

    # Print summary
    print("\n" + "="*80)
    print(" TEST SUMMARY")
    print("="*80)

    for name, result in results.items():
        status = "✓" if result == "PASSED" else "✗"
        print(f"{status} {name}: {result}")

    print("\n" + "-"*80)
    print(f"Total: {len(tests)} | Passed: {passed} | Failed: {failed}")
    print("="*80)

    return passed, failed


if __name__ == "__main__":
    passed, failed = run_all_tests()
    sys.exit(0 if failed == 0 else 1)
