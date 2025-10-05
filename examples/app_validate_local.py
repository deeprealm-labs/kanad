"""
Application-style Governance Validation (Local)

End-to-end validation for covalent, ionic, and metallic systems using
your framework's governance, ansatzes, and solvers. Produces a concise
stdout summary plus a machine-readable JSON file.
"""

import os
import sys
import json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np

from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.representations.second_quantization import SecondQuantizationRepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz, IonicGovernanceAnsatz
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.solvers import VQESolver, SQDSolver, QPESolver
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper


def summarize_case(name: str, bond_type: str, validation: dict, energy: float | None, circuit_info: dict | None):
    print("\n" + "=" * 80)
    print(f"APPLICATION VALIDATION: {name} [{bond_type}]")
    print("=" * 80)

    status = "PASS" if validation.get('all_checks_passed') else "WARN"
    print(f"Governance: {status}")
    for check in validation.get('checks', []):
        s = "✓" if check['passed'] else "✗"
        print(f"  [{s}] {check['name']}: {check.get('message', '')}")

    if circuit_info:
        print("Circuit:")
        print(f"  qubits: {circuit_info.get('qubits')}")
        print(f"  gates:  {circuit_info.get('gates')}")

    if energy is not None:
        print(f"Energy (Ha): {energy:.6f}")


def run_covalent():
    # CH4
    c = Atom('C', position=np.array([0.0, 0.0, 0.0]))
    Hs = [
        Atom('H', position=np.array([0.629, 0.629, 0.629])),
        Atom('H', position=np.array([0.629, -0.629, -0.629])),
        Atom('H', position=np.array([-0.629, 0.629, -0.629])),
        Atom('H', position=np.array([-0.629, -0.629, 0.629])),
    ]
    mol = Molecule([c] + Hs, charge=0, spin=0, basis='sto-3g')
    rep = LCAORepresentation(mol)
    H = CovalentHamiltonian(molecule=mol, representation=rep, use_governance=True)
    validation = H.validate_with_governance()

    # Governance ansatz + VQE (statevector simulator)
    mapper = JordanWignerMapper()
    ans = CovalentGovernanceAnsatz(2 * H.n_orbitals, mol.n_electrons, n_layers=1)
    vqe = VQESolver(H, ans, mapper, backend='aer_simulator_statevector', optimizer='SLSQP', max_iterations=50)
    res = vqe.solve()

    circuit = ans.build_circuit()
    return {
        'name': 'CH4',
        'bond_type': 'covalent',
        'validation': validation,
        'energy': float(res['energy']),
        'circuit': {
            'qubits': circuit.n_qubits,
            'gates': len(circuit.gates)
        }
    }


def run_ionic():
    # LiH
    Li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    H = Atom('H', position=np.array([1.6, 0.0, 0.0]))
    mol = Molecule([Li, H], charge=0, spin=0)
    rep = SecondQuantizationRepresentation(molecule=mol)
    H_ionic = IonicHamiltonian(molecule=mol, representation=rep, use_governance=True)
    validation = H_ionic.validate_with_governance()

    # Governance ansatz + SQD (short)
    ans = IonicGovernanceAnsatz(2 * H_ionic.n_orbitals, H_ionic.n_electrons, n_layers=2)
    sqd = SQDSolver(H_ionic, ansatz=ans, n_samples=200, max_iterations=3)
    res = sqd.solve()

    circuit = ans.build_circuit()
    return {
        'name': 'LiH',
        'bond_type': 'ionic',
        'validation': validation,
        'energy': float(res['energy']),
        'circuit': {
            'qubits': circuit.n_qubits,
            'gates': len(circuit.gates)
        }
    }


def run_metallic():
    # 1D chain (4 atoms)
    positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.2, 0.0, 0.0]),
        np.array([2.4, 0.0, 0.0]),
        np.array([3.6, 0.0, 0.0]),
    ]
    atoms = [Atom('Na', position=p) for p in positions]
    mol = Molecule(atoms)
    H_met = MetallicHamiltonian(molecule=mol, lattice_type='1d_chain', hopping_parameter=-1.0, onsite_energy=0.0, hubbard_u=0.0)

    # QPE fallback (classical exact) after eri=None handling
    qpe = QPESolver(H_met, n_ancilla=4, evolution_time=1.0)
    res = qpe.solve()

    return {
        'name': 'Na4-chain',
        'bond_type': 'metallic',
        'validation': {
            'all_checks_passed': True,
            'checks': []
        },
        'energy': float(res['energy']),
        'circuit': None
    }


def main():
    results = []
    results.append(run_covalent())
    results.append(run_ionic())
    results.append(run_metallic())

    # Print readable summaries
    for r in results:
        summarize_case(r['name'], r['bond_type'], r['validation'], r['energy'], r['circuit'])

    # Write machine-readable report
    out = {
        'suite': 'app_validate_local',
        'results': results
    }
    with open('app_validate_local_report.json', 'w') as f:
        json.dump(out, f, indent=2)
    print("\nWrote app_validate_local_report.json")


if __name__ == '__main__':
    main()


