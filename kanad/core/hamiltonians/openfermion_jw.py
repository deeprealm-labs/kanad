"""
Jordan-Wigner transformation using OpenFermion (validated reference).

OpenFermion provides a well-tested implementation of the Jordan-Wigner
transformation that we can use as our ground truth.
"""

import numpy as np
from typing import Dict
from qiskit.quantum_info import SparsePauliOp
import logging

logger = logging.getLogger(__name__)


def openfermion_jordan_wigner(
    h_mo: np.ndarray,
    eri_mo: np.ndarray,
    nuclear_repulsion: float,
    n_electrons: int
) -> SparsePauliOp:
    """
    Build Hamiltonian using OpenFermion's Jordan-Wigner transformation.

    Args:
        h_mo: One-electron integrals in MO basis
        eri_mo: Two-electron integrals in MO basis (chemist notation ⟨pq|rs⟩)
        nuclear_repulsion: Nuclear repulsion energy
        n_electrons: Number of electrons

    Returns:
        Qiskit SparsePauliOp
    """
    try:
        from openfermion import FermionOperator, jordan_wigner
        from openfermion.transforms import get_fermion_operator
        from openfermion.ops import InteractionOperator
    except ImportError:
        raise ImportError("OpenFermion not installed. Install with: pip install openfermion")

    n_orbitals = h_mo.shape[0]
    n_qubits = 2 * n_orbitals

    logger.info(f"Using OpenFermion Jordan-Wigner: {n_orbitals} orbitals → {n_qubits} qubits")

    # Build FermionOperator directly
    fermion_ham = FermionOperator()

    # Add constant term (nuclear repulsion)
    print(f"➕ Adding nuclear repulsion constant: {nuclear_repulsion:.8f} Ha")
    fermion_ham += FermionOperator((), nuclear_repulsion)
    print(f"   Fermion Hamiltonian now has {len(fermion_ham.terms)} terms")

    # One-body terms: Σ_{pq} h_{pq} a†_p a_q
    # For both spin-up and spin-down
    for p in range(n_orbitals):
        for q in range(n_orbitals):
            if abs(h_mo[p, q]) < 1e-12:
                continue

            # Spin-up: indices 2*p, 2*q
            fermion_ham += FermionOperator(
                ((2*p, 1), (2*q, 0)),  # a†_{2p} a_{2q}
                h_mo[p, q]
            )

            # Spin-down: indices 2*p+1, 2*q+1
            fermion_ham += FermionOperator(
                ((2*p+1, 1), (2*q+1, 0)),  # a†_{2p+1} a_{2q+1}
                h_mo[p, q]
            )

    # Two-body terms: (1/2) Σ_{pqrs} ⟨pq|rs⟩ a†_p a†_r a_s a_q
    # For chemist notation ⟨pq|rs⟩, the operator is a†_p a†_r a_s a_q
    for p in range(n_orbitals):
        for q in range(n_orbitals):
            for r in range(n_orbitals):
                for s in range(n_orbitals):
                    if abs(eri_mo[p, q, r, s]) < 1e-12:
                        continue

                    coeff = 0.5 * eri_mo[p, q, r, s]

                    # ↑↑: a†_{p↑} a†_{r↑} a_{s↑} a_{q↑}
                    fermion_ham += FermionOperator(
                        ((2*p, 1), (2*r, 1), (2*s, 0), (2*q, 0)),
                        coeff
                    )

                    # ↑↓: a†_{p↑} a†_{r↓} a_{s↓} a_{q↑}
                    fermion_ham += FermionOperator(
                        ((2*p, 1), (2*r+1, 1), (2*s+1, 0), (2*q, 0)),
                        coeff
                    )

                    # ↓↑: a†_{p↓} a†_{r↑} a_{s↑} a_{q↓}
                    fermion_ham += FermionOperator(
                        ((2*p+1, 1), (2*r, 1), (2*s, 0), (2*q+1, 0)),
                        coeff
                    )

                    # ↓↓: a†_{p↓} a†_{r↓} a_{s↓} a_{q↓}
                    fermion_ham += FermionOperator(
                        ((2*p+1, 1), (2*r+1, 1), (2*s+1, 0), (2*q+1, 0)),
                        coeff
                    )

    # Apply Jordan-Wigner transformation
    print(f"🔄 Applying Jordan-Wigner transformation...")
    qubit_ham = jordan_wigner(fermion_ham)
    print(f"   Qubit Hamiltonian has {len(qubit_ham.terms)} terms")

    # Convert to Qiskit SparsePauliOp
    pauli_op = _openfermion_to_qiskit(qubit_ham, n_qubits)

    print(f"✅ Built Hamiltonian with {len(pauli_op)} Pauli terms")

    # DEBUG: Check if nuclear repulsion is preserved
    identity_coeff = 0.0
    for i, pauli_str in enumerate(pauli_op.paulis):
        if set(str(pauli_str).replace('I', '')) == set():
            print(f"   Identity term: {str(pauli_str)} = {pauli_op.coeffs[i]:.8f}")
            identity_coeff += pauli_op.coeffs[i]

    print(f"")
    print(f"🔍 NUCLEAR REPULSION CHECK:")
    print(f"   Input nuclear repulsion: {nuclear_repulsion:.8f} Ha")
    print(f"   Total identity coefficient: {identity_coeff:.8f} Ha")
    print(f"   Note: Identity includes nuclear repulsion + electronic contributions")
    print(f"")

    return pauli_op


def _openfermion_to_qiskit(qubit_operator, n_qubits: int) -> SparsePauliOp:
    """
    Convert OpenFermion QubitOperator to Qiskit SparsePauliOp.

    Args:
        qubit_operator: OpenFermion QubitOperator
        n_qubits: Number of qubits

    Returns:
        Qiskit SparsePauliOp
    """
    pauli_dict = {}

    for term, coeff in qubit_operator.terms.items():
        # term is a tuple of (qubit_index, pauli_char) pairs
        # e.g., ((0, 'X'), (1, 'Y'))

        # Build Pauli string
        pauli_str = ['I'] * n_qubits
        for qubit_idx, pauli_char in term:
            pauli_str[qubit_idx] = pauli_char

        # OpenFermion and Qiskit use the SAME qubit ordering
        # (both are little-endian: rightmost qubit is index 0)
        pauli_str_qiskit = ''.join(pauli_str)

        if pauli_str_qiskit in pauli_dict:
            pauli_dict[pauli_str_qiskit] += coeff
        else:
            pauli_dict[pauli_str_qiskit] = coeff

    # Filter negligible terms
    pauli_dict = {k: v for k, v in pauli_dict.items() if abs(v) > 1e-12}

    if len(pauli_dict) == 0:
        return SparsePauliOp(['I' * n_qubits], [0.0])

    pauli_strings = list(pauli_dict.keys())
    coefficients = [pauli_dict[s].real if isinstance(pauli_dict[s], complex) else pauli_dict[s]
                    for s in pauli_strings]

    return SparsePauliOp(pauli_strings, coefficients)
