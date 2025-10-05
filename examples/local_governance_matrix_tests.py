"""
Local Governance Matrix Tests

Runs fast, local validations across Hamiltonians and solvers without cloud.
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np

from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.representations.second_quantization import SecondQuantizationRepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian

from kanad.ansatze.governance_aware_ansatz import (
    CovalentGovernanceAnsatz,
    IonicGovernanceAnsatz,
    AdaptiveGovernanceAnsatz,
)
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.solvers import VQESolver, SQDSolver, QPESolver
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper


def header(title: str):
    print("\n" + "=" * 80)
    print(f"{title}")
    print("=" * 80)


def test_covalent_h2_vqe_local():
    header("LOCAL: H2 Covalent + VQE (governance ansatz)")

    H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
    mol = Molecule([H1, H2], charge=0, spin=0, basis='sto-3g')

    rep = LCAORepresentation(mol)
    H = CovalentHamiltonian(molecule=mol, representation=rep, use_governance=True)

    ansatz = CovalentGovernanceAnsatz(n_qubits=2 * H.n_orbitals, n_electrons=mol.n_electrons, n_layers=1)
    mapper = JordanWignerMapper()

    vqe = VQESolver(H, ansatz, mapper, backend='classical', optimizer='SLSQP', max_iterations=200)
    res = vqe.solve()
    print(f"Energy (Ha): {res['energy']:.6f}")
    return res


def test_ionic_lih_sqd_local():
    header("LOCAL: LiH Ionic + SQD (governance ansatz)")

    Li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    H = Atom('H', position=np.array([1.6, 0.0, 0.0]))
    mol = Molecule([Li, H], charge=0, spin=0)

    rep = SecondQuantizationRepresentation(molecule=mol, include_spin=True)
    H_ionic = IonicHamiltonian(molecule=mol, representation=rep, use_governance=True)

    ansatz = IonicGovernanceAnsatz(n_qubits=2 * H_ionic.n_orbitals, n_electrons=H_ionic.n_electrons, n_layers=2)
    sqd = SQDSolver(H_ionic, ansatz=ansatz, n_samples=200, max_iterations=3)
    res = sqd.solve()
    print(f"Energy (Ha): {res['energy']:.6f}")
    return res


def test_metallic_chain_qpe_local():
    header("LOCAL: Metallic chain + QPE (tight-binding)")

    # Simple 1D chain of 4 sites
    positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([2.0, 0.0, 0.0]),
        np.array([3.0, 0.0, 0.0]),
    ]
    atoms = [Atom('Na', position=p) for p in positions]
    mol = Molecule(atoms)

    H_met = MetallicHamiltonian(
        molecule=mol,
        lattice_type='1d_chain',
        hopping_parameter=-1.0,
        onsite_energy=0.0,
        hubbard_u=0.0,
        periodic=True,
    )

    # Use classical fallback QPE (works with eri=None after fix)
    qpe = QPESolver(H_met, n_ancilla=4, evolution_time=1.0, initial_state=None)
    res = qpe.solve()
    # Report electronic-only and total if available
    e_total = res['energy']
    print(f"Energy (Ha): {e_total:.6f}")
    return res


def main():
    test_covalent_h2_vqe_local()
    test_ionic_lih_sqd_local()
    test_metallic_chain_qpe_local()


if __name__ == "__main__":
    main()


