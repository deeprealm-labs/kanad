"""
Local short optimization runs for CH4 (covalent) and LiH (ionic)
Compares governance-aware ansatz vs UCC with VQE and SQD.
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
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz, IonicGovernanceAnsatz
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.solvers import VQESolver, SQDSolver
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper


def header(t: str):
    print("\n" + "=" * 80)
    print(t)
    print("=" * 80)


def run_ch4():
    header("LOCAL OPT: CH4 covalent (VQE)")
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

    mapper = JordanWignerMapper()

    cov_ans = CovalentGovernanceAnsatz(2 * H.n_orbitals, mol.n_electrons, n_layers=1)
    vqe_cov = VQESolver(H, cov_ans, mapper, backend='aer_simulator_statevector', optimizer='COBYLA', max_iterations=30)
    r_cov = vqe_cov.solve()
    print(f"CovalentGovernanceAnsatz  Energy (Ha): {r_cov['energy']:.6f}")

    ucc = UCCAnsatz(2 * H.n_orbitals, mol.n_electrons, include_singles=True, include_doubles=True)
    vqe_ucc = VQESolver(H, ucc, mapper, backend='aer_simulator_statevector', optimizer='COBYLA', max_iterations=30)
    r_ucc = vqe_ucc.solve()
    print(f"UCCSD                   Energy (Ha): {r_ucc['energy']:.6f}")


def run_lih():
    header("LOCAL OPT: LiH ionic (SQD)")
    Li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    H = Atom('H', position=np.array([1.6, 0.0, 0.0]))
    mol = Molecule([Li, H], charge=0, spin=0)
    rep = SecondQuantizationRepresentation(molecule=mol)
    H_ionic = IonicHamiltonian(molecule=mol, representation=rep, use_governance=True)

    ion_ans = IonicGovernanceAnsatz(2 * H_ionic.n_orbitals, H_ionic.n_electrons, n_layers=2)
    sqd_ion = SQDSolver(H_ionic, ansatz=ion_ans, n_samples=200, max_iterations=3)
    r_ion = sqd_ion.solve()
    print(f"IonicGovernanceAnsatz   Energy (Ha): {r_ion['energy']:.6f}")

    ucc = UCCAnsatz(2 * H_ionic.n_orbitals, H_ionic.n_electrons, include_singles=True, include_doubles=True)
    sqd_ucc = SQDSolver(H_ionic, ansatz=ucc, n_samples=200, max_iterations=3)
    r_ucc = sqd_ucc.solve()
    print(f"UCCSD                   Energy (Ha): {r_ucc['energy']:.6f}")


def main():
    run_ch4()
    run_lih()


if __name__ == "__main__":
    main()


