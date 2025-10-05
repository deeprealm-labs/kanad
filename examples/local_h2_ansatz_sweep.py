"""
Local H2 ansatz sweep across governance-aware and standard ansatzes
Runs VQE locally with different ansatz choices and compares energies.
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np

from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.ansatze.governance_aware_ansatz import (
    CovalentGovernanceAnsatz,
    AdaptiveGovernanceAnsatz,
)
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.solvers import VQESolver
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper


def header(t: str):
    print("\n" + "=" * 80)
    print(t)
    print("=" * 80)


def run_vqe(hamiltonian, ansatz, label: str):
    mapper = JordanWignerMapper()
    vqe = VQESolver(hamiltonian, ansatz, mapper, backend='classical', optimizer='SLSQP', max_iterations=200)
    res = vqe.solve()
    print(f"{label:28s}  Energy (Ha): {res['energy']:.6f}")
    return res['energy']


def main():
    header("LOCAL: H2 ansatz sweep (governance vs UCC)")

    H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
    mol = Molecule([H1, H2], charge=0, spin=0, basis='sto-3g')

    rep = LCAORepresentation(mol)
    H_cov = CovalentHamiltonian(molecule=mol, representation=rep, use_governance=True)

    energies = {}

    # Governance-aware covalent ansatz
    cov_ansatz = CovalentGovernanceAnsatz(n_qubits=2 * H_cov.n_orbitals, n_electrons=mol.n_electrons, n_layers=1)
    energies['covalent_governance'] = run_vqe(H_cov, cov_ansatz, 'CovalentGovernanceAnsatz')

    # Adaptive governance (auto; covalent for H2)
    adaptive_ansatz = AdaptiveGovernanceAnsatz(n_qubits=2 * H_cov.n_orbitals, n_electrons=mol.n_electrons, n_layers=1, bonding_type='auto')
    energies['adaptive_governance'] = run_vqe(H_cov, adaptive_ansatz, 'AdaptiveGovernanceAnsatz')

    # UCCSD
    ucc_ansatz = UCCAnsatz(n_qubits=2 * H_cov.n_orbitals, n_electrons=mol.n_electrons, include_singles=True, include_doubles=True)
    energies['uccsd'] = run_vqe(H_cov, ucc_ansatz, 'UCCSD')

    # Summary
    header("SUMMARY: H2 ansatz sweep")
    for name, e in energies.items():
        print(f"{name:24s}  {e:.6f} Ha")


if __name__ == "__main__":
    main()


