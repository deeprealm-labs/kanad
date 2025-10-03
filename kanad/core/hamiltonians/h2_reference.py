"""
Reference H2 Hamiltonian implementation.

This is a minimal, correct implementation for H2 molecule that we can use
to validate and understand the correct formula.
"""

import numpy as np
from qiskit.quantum_info import SparsePauliOp


def build_h2_hamiltonian(h_core: np.ndarray, eri: np.ndarray, nuc_rep: float) -> SparsePauliOp:
    """
    Build H2 Hamiltonian using correct second quantization formula.

    For H2 with 2 spatial orbitals mapped to 4 spin orbitals (qubits):
    - Qubit 0: orbital 0, spin-up
    - Qubit 1: orbital 1, spin-up
    - Qubit 2: orbital 0, spin-down
    - Qubit 3: orbital 1, spin-down

    Args:
        h_core: One-electron integrals (2x2)
        eri: Two-electron integrals [i,j,k,l] = (ij|kl) (2x2x2x2)
        nuc_rep: Nuclear repulsion energy

    Returns:
        SparsePauliOp

    Reference:
        Szabo & Ostlund, "Modern Quantum Chemistry"
        Second quantization Hamiltonian in spatial orbital basis
    """
    pauli_dict = {}

    # Nuclear repulsion
    pauli_dict['IIII'] = nuc_rep

    # One-electron terms: Σ h_pq a†_p a_q
    # For each spatial orbital (p,q), we have spin-up and spin-down
    for p in range(2):
        for q in range(2):
            h_pq = h_core[p, q]
            if abs(h_pq) < 1e-12:
                continue

            # Number operator: a†_p a_p = n_p = (I - Z_p)/2
            if p == q:
                # Spin-up (qubit p)
                _add_to_dict(pauli_dict, 'IIII', h_pq / 2)
                z_str = ['I'] * 4
                z_str[p] = 'Z'
                _add_to_dict(pauli_dict, ''.join(z_str), -h_pq / 2)

                # Spin-down (qubit p+2)
                _add_to_dict(pauli_dict, 'IIII', h_pq / 2)
                z_str = ['I'] * 4
                z_str[p + 2] = 'Z'
                _add_to_dict(pauli_dict, ''.join(z_str), -h_pq / 2)

            # Hopping: a†_p a_q (p ≠ q)
            else:
                # Jordan-Wigner: a†_p a_q = (1/4)(XX + YY) + (i/4)(XY - YX) with Z string
                # For qubits 0,1 (spin-up) or 2,3 (spin-down)

                # Spin-up: qubits 0, 1
                if p < q:
                    z_qubits = list(range(p + 1, q))
                else:
                    z_qubits = list(range(q + 1, p)  )

                # XX term
                xx_str = _build_pauli_exc(p, q, 'X', 'X', z_qubits)
                _add_to_dict(pauli_dict, xx_str, h_pq / 4)

                # YY term
                yy_str = _build_pauli_exc(p, q, 'Y', 'Y', z_qubits)
                _add_to_dict(pauli_dict, yy_str, h_pq / 4)

                # XY term (imaginary)
                xy_str = _build_pauli_exc(p, q, 'X', 'Y', z_qubits)
                _add_to_dict(pauli_dict, xy_str, 1j * h_pq / 4)

                # YX term (imaginary)
                yx_str = _build_pauli_exc(p, q, 'Y', 'X', z_qubits)
                _add_to_dict(pauli_dict, yx_str, -1j * h_pq / 4)

                # Spin-down: qubits 2, 3
                p_down = p + 2
                q_down = q + 2

                if p_down < q_down:
                    z_qubits_down = list(range(p_down + 1, q_down))
                else:
                    z_qubits_down = list(range(q_down + 1, p_down))

                xx_str = _build_pauli_exc(p_down, q_down, 'X', 'X', z_qubits_down)
                _add_to_dict(pauli_dict, xx_str, h_pq / 4)

                yy_str = _build_pauli_exc(p_down, q_down, 'Y', 'Y', z_qubits_down)
                _add_to_dict(pauli_dict, yy_str, h_pq / 4)

                xy_str = _build_pauli_exc(p_down, q_down, 'X', 'Y', z_qubits_down)
                _add_to_dict(pauli_dict, xy_str, 1j * h_pq / 4)

                yx_str = _build_pauli_exc(p_down, q_down, 'Y', 'X', z_qubits_down)
                _add_to_dict(pauli_dict, yx_str, -1j * h_pq / 4)

    # Two-electron terms: (1/2) Σ (pq|rs) a†_p a†_q a_s a_r
    # Only implement number-number interactions for now (no excitations)
    # This gives Coulomb and exchange repulsion

    for p in range(2):
        for q in range(2):
            for r in range(2):
                for s in range(2):
                    v_pqrs = eri[p, q, r, s]  # (pq|rs)
                    if abs(v_pqrs) < 1e-12:
                        continue

                    # Only number operators: when p=r and q=s
                    # a†_p a†_q a_q a_p = n_p n_q
                    if p == r and q == s and p != q:
                        # n_p n_q for all spin combinations
                        # Spin-up, spin-up
                        _add_number_number(pauli_dict, p, q, 0.5 * v_pqrs)

                        # Spin-up, spin-down
                        _add_number_number(pauli_dict, p, q + 2, 0.5 * v_pqrs)

                        # Spin-down, spin-up
                        _add_number_number(pauli_dict, p + 2, q, 0.5 * v_pqrs)

                        # Spin-down, spin-down
                        _add_number_number(pauli_dict, p + 2, q + 2, 0.5 * v_pqrs)

    # Convert to SparsePauliOp
    pauli_strings = list(pauli_dict.keys())
    coefficients = list(pauli_dict.values())

    return SparsePauliOp(pauli_strings, coefficients)


def _add_to_dict(pauli_dict, pauli_str, coeff):
    """Add coefficient to Pauli dictionary."""
    if pauli_str in pauli_dict:
        pauli_dict[pauli_str] += coeff
    else:
        pauli_dict[pauli_str] = coeff


def _build_pauli_exc(q_from, q_to, pauli_from, pauli_to, z_qubits):
    """Build Pauli string for excitation with Z string."""
    pauli_list = ['I'] * 4
    pauli_list[q_from] = pauli_from
    pauli_list[q_to] = pauli_to

    for z_q in z_qubits:
        pauli_list[z_q] = 'Z'

    return ''.join(pauli_list)


def _add_number_number(pauli_dict, q_i, q_j, coeff):
    """
    Add n_i n_j term to Pauli dictionary.

    n_i n_j = [(I - Z_i)/2] [(I - Z_j)/2]
            = (1/4)[I - Z_i - Z_j + Z_i Z_j]
    """
    # I term
    _add_to_dict(pauli_dict, 'IIII', coeff / 4)

    # -Z_i term
    z_i_str = ['I'] * 4
    z_i_str[q_i] = 'Z'
    _add_to_dict(pauli_dict, ''.join(z_i_str), -coeff / 4)

    # -Z_j term
    z_j_str = ['I'] * 4
    z_j_str[q_j] = 'Z'
    _add_to_dict(pauli_dict, ''.join(z_j_str), -coeff / 4)

    # +Z_i Z_j term
    z_ij_str = ['I'] * 4
    z_ij_str[q_i] = 'Z'
    z_ij_str[q_j] = 'Z'
    _add_to_dict(pauli_dict, ''.join(z_ij_str), coeff / 4)
