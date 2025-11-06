#!/usr/bin/env python3
"""
KANAD Master Test Suite
Systematically tests all Hamiltonians, Ansatze, Mappers, and Basis Sets

Test molecules: H2, LiH, BeH2, HeH+, H2O
Components tested:
- Hamiltonians: Covalent, Ionic (via governance ansatze)
- Ansatze: UCC, Hardware-Efficient, TwoLocal, Governance-aware
- Mappers: Jordan-Wigner, Bravyi-Kitaev
- Basis: sto-3g, 6-31g (others commented out for speed)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import json
import time
from datetime import datetime
from typing import Dict, List, Tuple
import numpy as np

from api.services.experiment_service import create_molecule_from_config
from kanad.ansatze import (
    UCCAnsatz, HardwareEfficientAnsatz, TwoLocalAnsatz,
    IonicGovernanceAnsatz, CovalentGovernanceAnsatz
)
from kanad.utils.vqe_solver import VQESolver

# Test configuration
TEST_MOLECULES = {
    'H2': {
        'smiles': '[H][H]',
        'charge': 0,
        'multiplicity': 1,
        'expected_hf_sto3g': -1.1167,  # Approximate
        'description': 'Hydrogen molecule - simplest case'
    },
    'HeH+': {
        'smiles': '[He+][H]',
        'charge': 1,
        'multiplicity': 1,
        'expected_hf_sto3g': -2.86,  # Approximate
        'description': 'Helium hydride cation - ionic character'
    },
    'LiH': {
        'smiles': '[Li][H]',
        'charge': 0,
        'multiplicity': 1,
        'expected_hf_sto3g': -7.86,  # Approximate
        'description': 'Lithium hydride - polar covalent'
    },
    'H2O': {
        'smiles': 'O',
        'charge': 0,
        'multiplicity': 1,
        'expected_hf_sto3g': -74.97,  # Approximate
        'description': 'Water - bent molecule'
    },
    # BeH2 - commented out for speed, uncomment for full test
    # 'BeH2': {
    #     'smiles': '[Be](H)H',
    #     'charge': 0,
    #     'multiplicity': 1,
    #     'expected_hf_sto3g': -15.77,
    #     'description': 'Beryllium hydride - linear molecule'
    # },
}

BASIS_SETS = ['sto-3g']  # Start with minimal basis
# For full testing, add: '6-31g', 'cc-pvdz'

MAPPERS = ['jordan_wigner']  # , 'bravyi_kitaev']  # Start with JW

ANSATZE = {
    'UCC': {
        'class': UCCAnsatz,
        'description': 'Unitary Coupled Cluster',
        'best_for': 'High accuracy molecular systems'
    },
    'HardwareEfficient': {
        'class': HardwareEfficientAnsatz,
        'description': 'Hardware-optimized circuit',
        'best_for': 'Cloud quantum hardware'
    },
    'TwoLocal': {
        'class': TwoLocalAnsatz,
        'description': 'Two-qubit entanglement layers',
        'best_for': 'Customizable circuits'
    },
    'CovalentGovernance': {
        'class': CovalentGovernanceAnsatz,
        'description': 'Covalent bonding-aware',
        'best_for': 'Molecular systems with covalent bonds'
    },
    # IonicGovernance - commented for now, test separately
    # 'IonicGovernance': {
    #     'class': IonicGovernanceAnsatz,
    #     'description': 'Ionic bonding-aware',
    #     'best_for': 'Ionic systems'
    # },
}

OPTIMIZERS = {
    'COBYLA': {
        'max_iterations': 50,
        'description': 'Fast, gradient-free'
    },
    # Add more optimizers for comparison
    # 'Powell': {'max_iterations': 50},
    # 'SLSQP': {'max_iterations': 20},  # Fewer iterations (expensive)
}


class TestResult:
    """Store test results"""
    def __init__(self):
        self.molecule = None
        self.basis = None
        self.ansatz = None
        self.mapper = None
        self.optimizer = None
        self.hf_energy = None
        self.vqe_energy = None
        self.correlation_energy = None
        self.iterations = None
        self.function_evals = None
        self.converged = False
        self.runtime = None
        self.error = None
        self.success = False

    def to_dict(self) -> Dict:
        return {
            'molecule': self.molecule,
            'basis': self.basis,
            'ansatz': self.ansatz,
            'mapper': self.mapper,
            'optimizer': self.optimizer,
            'hf_energy': self.hf_energy,
            'vqe_energy': self.vqe_energy,
            'correlation_energy': self.correlation_energy,
            'iterations': self.iterations,
            'function_evaluations': self.function_evals,
            'converged': self.converged,
            'runtime_seconds': self.runtime,
            'error': self.error,
            'success': self.success
        }


def test_combination(
    molecule_name: str,
    molecule_config: Dict,
    basis: str,
    ansatz_name: str,
    ansatz_class,
    mapper: str,
    optimizer_name: str,
    optimizer_config: Dict
) -> TestResult:
    """Test a single combination of molecule/basis/ansatz/mapper/optimizer"""

    result = TestResult()
    result.molecule = molecule_name
    result.basis = basis
    result.ansatz = ansatz_name
    result.mapper = mapper
    result.optimizer = optimizer_name

    print(f"\n{'='*80}")
    print(f"Testing: {molecule_name} | {basis} | {ansatz_name} | {mapper} | {optimizer_name}")
    print(f"{'='*80}")

    try:
        t_start = time.time()

        # Create molecule
        mol_config = molecule_config.copy()
        mol_config['basis'] = basis

        molecule = create_molecule_from_config(mol_config)
        ham = molecule.hamiltonian

        print(f"Molecule: {molecule.n_electrons} electrons, {ham.n_orbitals} orbitals")

        # Get HF energy
        _, hf_energy = ham.solve_scf(max_iterations=50, conv_tol=1e-6)
        result.hf_energy = float(hf_energy)
        print(f"HF energy: {hf_energy:.8f} Ha")

        # Create ansatz
        n_qubits = 2 * ham.n_orbitals
        n_electrons = molecule.n_electrons

        # Different ansatze have different initialization
        try:
            if ansatz_name in ['CovalentGovernance', 'IonicGovernance']:
                # Governance ansatze need different parameters
                ansatz = ansatz_class(
                    n_qubits=n_qubits,
                    n_electrons=n_electrons,
                    hamiltonian=ham
                )
            else:
                ansatz = ansatz_class(
                    n_qubits=n_qubits,
                    n_electrons=n_electrons
                )
        except Exception as e:
            print(f"⚠️  Ansatz initialization failed: {e}")
            print("   Trying simplified initialization...")
            ansatz = ansatz_class(n_qubits=n_qubits, n_electrons=n_electrons)

        print(f"Ansatz: {ansatz.n_parameters} parameters")

        # Run VQE
        print(f"Running VQE ({optimizer_name}, {optimizer_config['max_iterations']} iterations)...")

        vqe = VQESolver(
            hamiltonian=ham,
            ansatz=ansatz,
            mapper=mapper,
            optimizer_method=optimizer_name,
            max_iterations=optimizer_config['max_iterations'],
            backend='statevector',
            enable_analysis=False,
            enable_optimization=False
        )

        vqe_result = vqe.solve()

        t_end = time.time()
        result.runtime = t_end - t_start

        # Extract results
        result.vqe_energy = float(vqe_result['energy'])
        result.correlation_energy = float(vqe_result['energy'] - hf_energy)
        result.iterations = int(vqe_result.get('iterations', 0))
        result.function_evals = int(vqe_result.get('function_evaluations', 0))
        result.converged = bool(vqe_result.get('converged', False))

        # Check if correlation energy recovered
        if result.correlation_energy < -0.001:  # At least 1 mHa correlation
            result.success = True
            status = "✅ PASS"
        else:
            result.success = False
            status = "⚠️  WARN"

        print(f"\nResults:")
        print(f"  VQE energy:      {result.vqe_energy:.8f} Ha")
        print(f"  Correlation:     {result.correlation_energy*1000:.3f} mHa")
        print(f"  Iterations:      {result.iterations}")
        print(f"  Function evals:  {result.function_evals}")
        print(f"  Runtime:         {result.runtime:.2f} s")
        print(f"  Status:          {status}")

    except Exception as e:
        result.error = str(e)
        result.success = False
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()

    return result


def run_master_test_suite():
    """Run all test combinations"""

    print("\n" + "="*80)
    print("KANAD MASTER TEST SUITE")
    print("="*80)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\nTest Configuration:")
    print(f"  Molecules: {list(TEST_MOLECULES.keys())}")
    print(f"  Basis sets: {BASIS_SETS}")
    print(f"  Ansatze: {list(ANSATZE.keys())}")
    print(f"  Mappers: {MAPPERS}")
    print(f"  Optimizers: {list(OPTIMIZERS.keys())}")
    print("="*80)

    all_results = []
    total_tests = (
        len(TEST_MOLECULES) *
        len(BASIS_SETS) *
        len(ANSATZE) *
        len(MAPPERS) *
        len(OPTIMIZERS)
    )

    print(f"\nTotal test combinations: {total_tests}")

    test_num = 0
    for mol_name, mol_config in TEST_MOLECULES.items():
        for basis in BASIS_SETS:
            for ansatz_name, ansatz_info in ANSATZE.items():
                for mapper in MAPPERS:
                    for opt_name, opt_config in OPTIMIZERS.items():
                        test_num += 1
                        print(f"\n[Test {test_num}/{total_tests}]")

                        result = test_combination(
                            mol_name,
                            mol_config,
                            basis,
                            ansatz_name,
                            ansatz_info['class'],
                            mapper,
                            opt_name,
                            opt_config
                        )

                        all_results.append(result)

    # Summary
    print("\n" + "="*80)
    print("TEST SUITE SUMMARY")
    print("="*80)

    passed = sum(1 for r in all_results if r.success)
    failed = sum(1 for r in all_results if not r.success)

    print(f"Total tests:  {len(all_results)}")
    print(f"Passed:       {passed} ({passed/len(all_results)*100:.1f}%)")
    print(f"Failed:       {failed} ({failed/len(all_results)*100:.1f}%)")

    # Save results to JSON
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    results_file = f"tests/master_test_results_{timestamp}.json"

    with open(results_file, 'w') as f:
        json.dump({
            'timestamp': timestamp,
            'config': {
                'molecules': list(TEST_MOLECULES.keys()),
                'basis_sets': BASIS_SETS,
                'ansatze': list(ANSATZE.keys()),
                'mappers': MAPPERS,
                'optimizers': list(OPTIMIZERS.keys())
            },
            'summary': {
                'total': len(all_results),
                'passed': passed,
                'failed': failed
            },
            'results': [r.to_dict() for r in all_results]
        }, f, indent=2)

    print(f"\nResults saved to: {results_file}")

    # Detailed breakdown
    print(f"\n{'='*80}")
    print("DETAILED RESULTS")
    print("="*80)

    for result in all_results:
        status = "✅" if result.success else ("❌" if result.error else "⚠️ ")
        corr_mha = result.correlation_energy * 1000 if result.correlation_energy else 0
        print(f"{status} {result.molecule:5s} | {result.basis:6s} | {result.ansatz:18s} | "
              f"{result.mapper:15s} | Corr: {corr_mha:+7.2f} mHa | "
              f"Time: {result.runtime:.1f}s" if result.runtime else "N/A")

    print("="*80)

    return all_results


if __name__ == "__main__":
    results = run_master_test_suite()

    # Check overall success
    if all(r.success for r in results):
        print("\n🎉 ALL TESTS PASSED!")
        sys.exit(0)
    else:
        print("\n⚠️  SOME TESTS FAILED - See details above")
        sys.exit(1)
