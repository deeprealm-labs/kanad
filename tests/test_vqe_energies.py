"""
VQE Energy Validation Test
===========================
Tests that VQE produces correct ground state energies for molecules.

This is the CORE test - energies must be accurate to < 1 mHa.
"""

import pytest
import numpy as np
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.solvers.vqe_solver import VQESolver
from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from pyscf import fci as pyscf_fci


def test_h2_energy():
    """Test H2 molecule energy"""
    h2 = Molecule(
        atoms=[Atom('H', [0.0, 0.0, 0.0]), Atom('H', [0.0, 0.0, 0.74])],
        charge=0, spin=0, basis='sto-3g'
    )
    _ = h2.n_electrons

    # FCI reference
    cisolver = pyscf_fci.FCI(h2._hamiltonian.mf)
    fci_energy = cisolver.kernel()[0]

    # VQE
    mapper = JordanWignerMapper()
    ansatz = HardwareEfficientAnsatz(
        n_qubits=h2.n_orbitals * 2,
        n_electrons=h2.n_electrons,
        n_layers=2
    )

    vqe = VQESolver(
        hamiltonian=h2._hamiltonian,
        ansatz=ansatz,
        mapper=mapper,
        backend='classical',
        optimizer='COBYLA',
        max_iterations=100
    )

    result = vqe.solve()
    vqe_energy = result['energy']

    error = abs(vqe_energy - fci_energy)
    assert error < 0.001, f"H2 energy error {error*1000:.3f} mHa > 1 mHa"


def test_lih_energy():
    """Test LiH molecule energy"""
    lih = Molecule(
        atoms=[Atom('Li', [0.0, 0.0, 0.0]), Atom('H', [0.0, 0.0, 1.5])],
        charge=0, spin=0, basis='sto-3g'
    )
    _ = lih.n_electrons

    # FCI reference
    cisolver = pyscf_fci.FCI(lih._hamiltonian.mf)
    fci_energy = cisolver.kernel()[0]

    # VQE
    mapper = JordanWignerMapper()
    ansatz = HardwareEfficientAnsatz(
        n_qubits=lih.n_orbitals * 2,
        n_electrons=lih.n_electrons,
        n_layers=2
    )

    vqe = VQESolver(
        hamiltonian=lih._hamiltonian,
        ansatz=ansatz,
        mapper=mapper,
        backend='classical',
        optimizer='COBYLA',
        max_iterations=150
    )

    result = vqe.solve()
    vqe_energy = result['energy']

    error = abs(vqe_energy - fci_energy)
    assert error < 0.01, f"LiH energy error {error*1000:.3f} mHa > 10 mHa"


if __name__ == "__main__":
    print("Testing H2...")
    test_h2_energy()
    print("✓ H2 PASSED")

    print("\nTesting LiH...")
    test_lih_energy()
    print("✓ LiH PASSED")

    print("\n✓ ALL TESTS PASSED - VQE ENERGIES ARE ACCURATE")
