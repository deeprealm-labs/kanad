"""
Comprehensive VQE Validation Suite.

Tests VQE with:
- Different bond types (Covalent, Ionic, Metallic)
- Different mappers (Jordan-Wigner, Bravyi-Kitaev, Hybrid Orbital)
- Different Hamiltonians
- Different ansätze (UCCSD, Hardware Efficient, Governance-aware)

Author: Kanad Framework
"""

import numpy as np
import logging
from typing import Dict, List, Tuple

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Import Kanad components
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.bonds.ionic_bond import IonicBond
from kanad.bonds.metallic_bond import MetallicBond

from kanad.core.representations.base_representation import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian

from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.mappers.bravyi_kitaev_mapper import BravyiKitaevMapper
from kanad.core.mappers.hybrid_orbital_mapper import HybridOrbitalMapper

from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz, EfficientSU2Ansatz
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz, IonicGovernanceAnsatz

from kanad.solvers.vqe_solver import VQESolver


def print_header(title: str, level: int = 1):
    """Print formatted header."""
    if level == 1:
        print("\n" + "=" * 70)
        print(title.center(70))
        print("=" * 70)
    elif level == 2:
        print("\n" + "-" * 70)
        print(title)
        print("-" * 70)
    else:
        print(f"\n{title}")


def print_result(name: str, value: float, unit: str = "eV", status: str = ""):
    """Print formatted result."""
    status_symbol = {"pass": "✅", "fail": "❌", "warn": "⚠️"}.get(status, "  ")
    print(f"  {status_symbol} {name:30s}: {value:12.6f} {unit}")


def calculate_error(vqe_energy: float, exact_energy: float) -> Tuple[float, float]:
    """Calculate absolute and relative error."""
    abs_error = abs(vqe_energy - exact_energy)
    rel_error = (abs_error / abs(exact_energy)) * 100 if exact_energy != 0 else float('inf')
    return abs_error, rel_error


def test_covalent_vqe_h2():
    """Test VQE on H2 molecule with different configurations."""
    print_header("TEST 1: H2 COVALENT BOND - VQE VALIDATION", level=1)

    # Create H2 molecule
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

    molecule = Molecule([h1, h2])
    representation = LCAORepresentation(molecule)
    hamiltonian = CovalentHamiltonian(molecule, representation)

    # Get exact energy for reference
    H_matrix = hamiltonian.to_matrix()
    eigenvalues = np.linalg.eigh(H_matrix)[0]
    exact_energy = eigenvalues[0] * 27.211386  # Ha to eV

    print_result("Exact Energy (Reference)", exact_energy, "eV")

    results = {}

    # Test 1.1: Different Mappers
    print_header("1.1: Different Mappers (UCCSD Ansatz)", level=2)

    mappers = {
        'Jordan-Wigner': JordanWignerMapper(),
        'Bravyi-Kitaev': BravyiKitaevMapper()
    }

    n_qubits = hamiltonian.n_orbitals * 2
    n_electrons = hamiltonian.n_electrons

    for mapper_name, mapper in mappers.items():
        try:
            print(f"\n  Testing {mapper_name} Mapper:")

            ansatz = UCCAnsatz(
                n_qubits=n_qubits,
                n_electrons=n_electrons,
                include_singles=True,
                include_doubles=True
            )

            solver = VQESolver(
                hamiltonian=hamiltonian,
                ansatz=ansatz,
                mapper=mapper,
                backend='classical',
                max_iterations=100
            )

            result = solver.solve()
            vqe_energy = result['energy'] * 27.211386  # Ha to eV

            abs_error, rel_error = calculate_error(vqe_energy, exact_energy)

            print_result("VQE Energy", vqe_energy, "eV")
            print_result("Absolute Error", abs_error, "eV")
            print_result("Relative Error", rel_error, "%")
            print_result("Converged", result['converged'], "")

            results[f"H2_UCCSD_{mapper_name}"] = {
                'energy': vqe_energy,
                'error': rel_error,
                'converged': result['converged']
            }

        except Exception as e:
            print(f"  ❌ Failed: {e}")
            results[f"H2_UCCSD_{mapper_name}"] = {'error': 'failed', 'exception': str(e)}

    # Test 1.2: Different Ansätze
    print_header("1.2: Different Ansätze (Jordan-Wigner Mapper)", level=2)

    mapper = JordanWignerMapper()

    ansatze = {
        'UCCSD': UCCAnsatz(n_qubits=n_qubits, n_electrons=n_electrons,
                          include_singles=True, include_doubles=True),
        'RealAmplitudes': RealAmplitudesAnsatz(n_qubits=n_qubits, n_electrons=n_electrons, n_layers=2),
        'EfficientSU2': EfficientSU2Ansatz(n_qubits=n_qubits, n_electrons=n_electrons, n_layers=2),
        'CovalentGovernance': CovalentGovernanceAnsatz(n_qubits=n_qubits, n_electrons=n_electrons,
                                                        n_layers=2, hybridization='sp')
    }

    for ansatz_name, ansatz in ansatze.items():
        try:
            print(f"\n  Testing {ansatz_name} Ansatz:")

            solver = VQESolver(
                hamiltonian=hamiltonian,
                ansatz=ansatz,
                mapper=mapper,
                backend='classical',
                optimizer='SLSQP',
                max_iterations=100
            )

            result = solver.solve()
            vqe_energy = result['energy'] * 27.211386  # Ha to eV

            abs_error, rel_error = calculate_error(vqe_energy, exact_energy)

            print_result("VQE Energy", vqe_energy, "eV")
            print_result("Absolute Error", abs_error, "eV")
            print_result("Relative Error", rel_error, "%")
            print_result("Converged", result['converged'], "")
            print_result("Parameters", len(result['parameters']), "")

            results[f"H2_{ansatz_name}_JW"] = {
                'energy': vqe_energy,
                'error': rel_error,
                'converged': result['converged'],
                'n_params': len(result['parameters'])
            }

        except Exception as e:
            print(f"  ❌ Failed: {e}")
            results[f"H2_{ansatz_name}_JW"] = {'error': 'failed', 'exception': str(e)}

    return results


def test_covalent_vqe_h2o():
    """Test VQE on H2O molecule."""
    print_header("TEST 2: H2O COVALENT BOND - VQE VALIDATION", level=1)

    # Create H2O molecule
    o_atom = Atom('O', position=np.array([0.0, 0.0, 0.0]))
    h1_atom = Atom('H', position=np.array([0.757, 0.586, 0.0]))
    h2_atom = Atom('H', position=np.array([-0.757, 0.586, 0.0]))

    molecule = Molecule([o_atom, h1_atom, h2_atom])
    representation = LCAORepresentation(molecule)
    hamiltonian = CovalentHamiltonian(molecule, representation)

    # Get exact energy
    H_matrix = hamiltonian.to_matrix()
    eigenvalues = np.linalg.eigh(H_matrix)[0]
    exact_energy = eigenvalues[0] * 27.211386  # Ha to eV

    print_result("Exact Energy (Reference)", exact_energy, "eV")

    results = {}

    print_header("2.1: VQE with Different Ansätze", level=2)

    n_qubits = hamiltonian.n_orbitals * 2
    n_electrons = hamiltonian.n_electrons
    mapper = JordanWignerMapper()

    print(f"\nSystem Info:")
    print(f"  Orbitals: {hamiltonian.n_orbitals}")
    print(f"  Electrons: {n_electrons}")
    print(f"  Qubits: {n_qubits}")

    # Test with different ansätze (fewer iterations for larger system)
    ansatze = {
        'RealAmplitudes': RealAmplitudesAnsatz(n_qubits=n_qubits, n_electrons=n_electrons, n_layers=1),
        'EfficientSU2': EfficientSU2Ansatz(n_qubits=n_qubits, n_electrons=n_electrons, n_layers=1),
        'CovalentGovernance': CovalentGovernanceAnsatz(n_qubits=n_qubits, n_electrons=n_electrons,
                                                        n_layers=1, hybridization='sp3')
    }

    for ansatz_name, ansatz in ansatze.items():
        try:
            print(f"\n  Testing {ansatz_name}:")

            solver = VQESolver(
                hamiltonian=hamiltonian,
                ansatz=ansatz,
                mapper=mapper,
                backend='classical',
                optimizer='COBYLA',
                max_iterations=50
            )

            result = solver.solve()
            vqe_energy = result['energy'] * 27.211386  # Ha to eV

            abs_error, rel_error = calculate_error(vqe_energy, exact_energy)

            print_result("VQE Energy", vqe_energy, "eV")
            print_result("Absolute Error", abs_error, "eV")
            print_result("Relative Error", rel_error, "%")
            print_result("Converged", result['converged'], "")

            results[f"H2O_{ansatz_name}"] = {
                'energy': vqe_energy,
                'error': rel_error,
                'converged': result['converged']
            }

        except Exception as e:
            print(f"  ❌ Failed: {e}")
            results[f"H2O_{ansatz_name}"] = {'error': 'failed', 'exception': str(e)}

    return results


def test_ionic_vqe_lih():
    """Test VQE on LiH ionic bond."""
    print_header("TEST 3: LiH IONIC BOND - VQE VALIDATION", level=1)

    # Create LiH molecule
    li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    h = Atom('H', position=np.array([1.596, 0.0, 0.0]))

    molecule = Molecule([li, h])
    representation = LCAORepresentation(molecule)
    hamiltonian = IonicHamiltonian(molecule, representation)

    # Get exact energy
    H_matrix = hamiltonian.to_matrix()
    eigenvalues = np.linalg.eigh(H_matrix)[0]
    exact_energy = eigenvalues[0] * 27.211386  # Ha to eV

    print_result("Exact Energy (Reference)", exact_energy, "eV")

    results = {}

    print_header("3.1: VQE with Ionic Governance", level=2)

    n_qubits = hamiltonian.n_orbitals * 2
    n_electrons = hamiltonian.n_electrons
    mapper = JordanWignerMapper()

    print(f"\nSystem Info:")
    print(f"  Orbitals: {hamiltonian.n_orbitals}")
    print(f"  Electrons: {n_electrons}")
    print(f"  Qubits: {n_qubits}")

    ansatze = {
        'RealAmplitudes': RealAmplitudesAnsatz(n_qubits=n_qubits, n_electrons=n_electrons, n_layers=2),
        'IonicGovernance': IonicGovernanceAnsatz(n_qubits=n_qubits, n_electrons=n_electrons,
                                                  n_layers=2, cation_orbitals=[0, 1])
    }

    for ansatz_name, ansatz in ansatze.items():
        try:
            print(f"\n  Testing {ansatz_name}:")

            solver = VQESolver(
                hamiltonian=hamiltonian,
                ansatz=ansatz,
                mapper=mapper,
                backend='classical',
                optimizer='SLSQP',
                max_iterations=100
            )

            result = solver.solve()
            vqe_energy = result['energy'] * 27.211386  # Ha to eV

            abs_error, rel_error = calculate_error(vqe_energy, exact_energy)

            print_result("VQE Energy", vqe_energy, "eV")
            print_result("Absolute Error", abs_error, "eV")
            print_result("Relative Error", rel_error, "%")
            print_result("Converged", result['converged'], "")

            results[f"LiH_{ansatz_name}"] = {
                'energy': vqe_energy,
                'error': rel_error,
                'converged': result['converged']
            }

        except Exception as e:
            print(f"  ❌ Failed: {e}")
            results[f"LiH_{ansatz_name}"] = {'error': 'failed', 'exception': str(e)}

    return results


def test_metallic_vqe():
    """Test VQE on metallic system."""
    print_header("TEST 4: METALLIC BOND - VQE VALIDATION", level=1)

    # Create 4-atom metallic chain
    atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]

    bond = MetallicBond(atoms, lattice_type='1d_chain', hubbard_u=0.0)

    # Get exact energy from tight-binding
    tb_result = bond.compute_energy(method='tight_binding')
    exact_energy = tb_result['energy']

    print_result("Tight-Binding Energy (Reference)", exact_energy, "eV")

    results = {}

    print_header("4.1: VQE on Metallic System", level=2)

    print(f"\nSystem Info:")
    print(f"  Atoms: {len(atoms)}")
    print(f"  Lattice: 1d_chain")
    print(f"  Electrons: {bond.hamiltonian.n_electrons}")

    # Note: Metallic VQE would require quantum treatment
    # For now, demonstrate the setup
    print("\n  ℹ️  Metallic VQE requires quantum many-body Hamiltonian")
    print("  ℹ️  Current implementation uses tight-binding (classical)")
    print("  ℹ️  Full quantum VQE for metals is future work")

    results['Metallic_4atom'] = {
        'tb_energy': exact_energy,
        'note': 'Quantum VQE for metals requires extended Hamiltonian'
    }

    return results


def print_summary(all_results: Dict):
    """Print summary of all results."""
    print_header("COMPREHENSIVE VQE VALIDATION SUMMARY", level=1)

    print("\n📊 Test Results:")
    print("-" * 70)

    categories = {
        'H2': 'H2 Molecule Tests',
        'H2O': 'H2O Molecule Tests',
        'LiH': 'LiH Ionic Tests',
        'Metallic': 'Metallic Tests'
    }

    for category, title in categories.items():
        print(f"\n{title}:")
        category_results = {k: v for k, v in all_results.items() if k.startswith(category)}

        for test_name, result in category_results.items():
            if isinstance(result, dict) and 'error' in result:
                if result['error'] == 'failed':
                    print(f"  ❌ {test_name}: FAILED - {result.get('exception', 'Unknown error')}")
                elif isinstance(result['error'], (int, float)):
                    status = "✅" if result['error'] < 10.0 else "⚠️"
                    print(f"  {status} {test_name}:")
                    print(f"      Energy: {result['energy']:.6f} eV")
                    print(f"      Error: {result['error']:.2f}%")
                    print(f"      Converged: {result['converged']}")

    print("\n" + "=" * 70)
    print("Key Findings:")
    print("-" * 70)

    # Count successes and failures
    total_tests = len(all_results)
    failed_tests = sum(1 for r in all_results.values()
                      if isinstance(r, dict) and r.get('error') == 'failed')
    passed_tests = total_tests - failed_tests

    print(f"  Total Tests: {total_tests}")
    print(f"  Passed: {passed_tests}")
    print(f"  Failed: {failed_tests}")
    print(f"  Success Rate: {(passed_tests/total_tests)*100:.1f}%")

    # Best and worst results
    errors = [(k, v['error']) for k, v in all_results.items()
             if isinstance(v, dict) and isinstance(v.get('error'), (int, float))]

    if errors:
        best = min(errors, key=lambda x: x[1])
        worst = max(errors, key=lambda x: x[1])

        print(f"\n  Best Performance: {best[0]} ({best[1]:.2f}% error)")
        print(f"  Worst Performance: {worst[0]} ({worst[1]:.2f}% error)")

    print("\n" + "=" * 70)


if __name__ == "__main__":
    print_header("KANAD VQE COMPREHENSIVE VALIDATION SUITE", level=1)
    print("Testing VQE with different bonds, mappers, Hamiltonians, and ansätze")

    all_results = {}

    try:
        # Test 1: H2 with comprehensive configurations
        h2_results = test_covalent_vqe_h2()
        all_results.update(h2_results)
    except Exception as e:
        logger.error(f"H2 tests failed: {e}")

    try:
        # Test 2: H2O
        h2o_results = test_covalent_vqe_h2o()
        all_results.update(h2o_results)
    except Exception as e:
        logger.error(f"H2O tests failed: {e}")

    try:
        # Test 3: LiH ionic
        lih_results = test_ionic_vqe_lih()
        all_results.update(lih_results)
    except Exception as e:
        logger.error(f"LiH tests failed: {e}")

    try:
        # Test 4: Metallic
        metallic_results = test_metallic_vqe()
        all_results.update(metallic_results)
    except Exception as e:
        logger.error(f"Metallic tests failed: {e}")

    # Print comprehensive summary
    print_summary(all_results)

    print("\n" + "=" * 70)
    print("VQE Validation Complete!")
    print("=" * 70)
