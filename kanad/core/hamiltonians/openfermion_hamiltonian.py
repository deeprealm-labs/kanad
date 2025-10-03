"""
OpenFermion-based Hamiltonian builder for Kanad framework.

This module uses OpenFermion's validated fermionic-to-qubit transformations
to build molecularHamiltonians correctly, avoiding the bugs in custom
pauli_converter.py.

Requires: pip install openfermion
"""

import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian

try:
    from openfermion import FermionOperator, jordan_wigner, bravyi_kitaev
    from qiskit.quantum_info import SparsePauliOp
    OPENFERMION_AVAILABLE = True
except ImportError:
    OPENFERMION_AVAILABLE = False


class OpenFermionHamiltonian:
    """
    Build quantum Hamiltonians using OpenFermion.

    This class provides a validated alternative to the custom pauli_converter,
    using OpenFermion's battle-tested implementations of Jordan-Wigner and
    other fermion-to-qubit transformations.

    Usage:
        >>> from kanad.bonds import CovalentBond
        >>> from kanad.core.atom import Atom
        >>>
        >>> h1 = Atom('H', [0.0, 0.0, 0.0])
        >>> h2 = Atom('H', [0.0, 0.0, 0.74])
        >>> bond = CovalentBond(h1, h2)
        >>>
        >>> pauli_op = OpenFermionHamiltonian.to_pauli_operator(
        >>>     bond.hamiltonian, mapper='jordan_wigner'
        >>> )
    """

    @staticmethod
    def to_pauli_operator(hamiltonian: 'MolecularHamiltonian', mapper: str = 'jordan_wigner') -> SparsePauliOp:
        """
        Convert molecular Hamiltonian to Qiskit Pauli operator using OpenFermion.

        Args:
            hamiltonian: Molecular Hamiltonian with h_core, eri, nuclear_repulsion
            mapper: Transformation type ('jordan_wigner' or 'bravyi_kitaev')

        Returns:
            SparsePauliOp representing the molecular Hamiltonian

        Raises:
            ImportError: If OpenFermion is not installed
        """
        if not OPENFERMION_AVAILABLE:
            raise ImportError(
                "OpenFermion not installed. Install with: pip install openfermion"
            )

        # Build fermionic Hamiltonian
        fermion_hamiltonian = OpenFermionHamiltonian._build_fermion_hamiltonian(hamiltonian)

        # Apply transformation
        if mapper == 'jordan_wigner':
            qubit_hamiltonian = jordan_wigner(fermion_hamiltonian)
        elif mapper == 'bravyi_kitaev':
            qubit_hamiltonian = bravyi_kitaev(fermion_hamiltonian)
        else:
            raise ValueError(f"Unknown mapper: {mapper}")

        # Convert to Qiskit
        n_spin_orbitals = hamiltonian.n_orbitals * 2
        pauli_op = OpenFermionHamiltonian._openfermion_to_qiskit(
            qubit_hamiltonian, n_qubits=n_spin_orbitals
        )

        return pauli_op

    @staticmethod
    def _build_fermion_hamiltonian(hamiltonian: 'MolecularHamiltonian') -> FermionOperator:
        """
        Build OpenFermion FermionOperator from molecular integrals.

        This correctly handles the conversion from spatial orbital integrals
        to spin-orbital fermionic operators.

        Args:
            hamiltonian: Molecular Hamiltonian with integrals

        Returns:
            FermionOperator representing the molecular Hamiltonian
        """
        h_core = hamiltonian.h_core
        eri = hamiltonian.eri if hasattr(hamiltonian, 'eri') else None
        nuclear_repulsion = hamiltonian.nuclear_repulsion

        n_spatial = hamiltonian.n_orbitals

        # Initialize Hamiltonian
        ferm_ham = FermionOperator()

        # 1. Nuclear repulsion (constant term)
        ferm_ham += FermionOperator('', nuclear_repulsion)

        # 2. One-electron terms: h_ij a†_i a_j (for both spins)
        # Spin orbital convention: 0↑, 0↓, 1↑, 1↓, ... (interleaved)
        for i in range(n_spatial):
            for j in range(n_spatial):
                if abs(h_core[i, j]) > 1e-12:
                    # Spin-up: 2*i
                    ferm_ham += FermionOperator(f'{2*i}^ {2*j}', h_core[i, j])
                    # Spin-down: 2*i+1
                    ferm_ham += FermionOperator(f'{2*i+1}^ {2*j+1}', h_core[i, j])

        # 3. Two-electron terms: (1/2) Σ_{pqrs} ⟨pq|rs⟩ a†_p a†_q a_s a_r
        if eri is not None:
            for p in range(n_spatial):
                for q in range(n_spatial):
                    for r in range(n_spatial):
                        for s in range(n_spatial):
                            # Kanad ERI matches PySCF convention: eri[i,j,k,l] = (ij|kl)
                            # So ⟨pq|rs⟩ = (pq|rs) = eri[p,q,r,s]
                            v_pqrs = eri[p, q, r, s]  # ⟨pq|rs⟩ in chemist notation

                            if abs(v_pqrs) > 1e-12:
                                # Hamiltonian: (1/2) Σ ⟨pq|rs⟩ a†_p a†_q a_s a_r
                                # Add all spin combinations (same spin for p,r and for q,s)
                                # ↑↑
                                ferm_ham += FermionOperator(
                                    f'{2*p}^ {2*q}^ {2*s} {2*r}',
                                    0.5 * v_pqrs
                                )
                                # ↑↓
                                ferm_ham += FermionOperator(
                                    f'{2*p}^ {2*q+1}^ {2*s+1} {2*r}',
                                    0.5 * v_pqrs
                                )
                                # ↓↑
                                ferm_ham += FermionOperator(
                                    f'{2*p+1}^ {2*q}^ {2*s} {2*r+1}',
                                    0.5 * v_pqrs
                                )
                                # ↓↓
                                ferm_ham += FermionOperator(
                                    f'{2*p+1}^ {2*q+1}^ {2*s+1} {2*r+1}',
                                    0.5 * v_pqrs
                                )

        return ferm_ham

    @staticmethod
    def _openfermion_to_qiskit(qubit_op, n_qubits: int) -> SparsePauliOp:
        """
        Convert OpenFermion QubitOperator to Qiskit SparsePauliOp.

        Args:
            qubit_op: OpenFermion QubitOperator
            n_qubits: Number of qubits

        Returns:
            Qiskit SparsePauliOp
        """
        pauli_dict = {}

        for term, coeff in qubit_op.terms.items():
            if len(term) == 0:
                # Identity term
                pauli_str = 'I' * n_qubits
            else:
                # Build Pauli string
                pauli_list = ['I'] * n_qubits
                for qubit_idx, pauli_op in term:
                    pauli_list[qubit_idx] = pauli_op
                pauli_str = ''.join(pauli_list)

            # Accumulate coefficients
            if pauli_str in pauli_dict:
                pauli_dict[pauli_str] += coeff
            else:
                pauli_dict[pauli_str] = coeff

        # Filter negligible terms
        pauli_dict = {k: v for k, v in pauli_dict.items() if abs(v) > 1e-12}

        # Convert to Qiskit
        if len(pauli_dict) == 0:
            return SparsePauliOp(['I' * n_qubits], [0.0])

        pauli_strings = list(pauli_dict.keys())
        coefficients = list(pauli_dict.values())

        return SparsePauliOp(pauli_strings, coefficients)
