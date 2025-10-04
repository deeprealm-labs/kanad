#!/usr/bin/env python3
"""
Solver tests - VQE, QPE, SQD, FCI, etc.
"""

import pytest
import numpy as np
import sys
sys.path.insert(0, '/home/mk/deeprealm/kanad')

from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.solvers.vqe_solver import VQESolver
from kanad.solvers.qpe_solver import QPESolver
from kanad.solvers.sqd_solver import SQDSolver
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper


class TestVQE:
    """Test VQE Solver."""

    def test_vqe_h2(self):
        """Test VQE on H2."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        bond = CovalentBond(h1, h2)
        hamiltonian = bond.hamiltonian

        n_qubits = 2 * hamiltonian.n_orbitals
        ansatz = UCCAnsatz(n_qubits, hamiltonian.n_electrons)
        mapper = JordanWignerMapper()

        vqe = VQESolver(hamiltonian, ansatz, mapper)
        result = vqe.solve()

        assert 'energy' in result
        assert result['energy'] < 0  # Bound state


class TestQPE:
    """Test QPE Solver."""

    def test_qpe_h2(self):
        """Test QPE on H2."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        bond = CovalentBond(h1, h2)
        hamiltonian = bond.hamiltonian

        qpe = QPESolver(hamiltonian, n_ancilla=4)
        result = qpe.solve()

        assert 'energy' in result
        assert isinstance(result['energy'], (int, float))


class TestSQD:
    """Test SQD Solver."""

    def test_sqd_h2(self):
        """Test SQD on H2."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

        bond = CovalentBond(h1, h2)
        hamiltonian = bond.hamiltonian

        sqd = SQDSolver(hamiltonian, n_samples=100)
        result = sqd.solve()

        assert 'energy' in result
        assert isinstance(result['energy'], (int, float))


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
