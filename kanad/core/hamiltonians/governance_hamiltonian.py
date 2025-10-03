"""
Governance-Based Quantum Hamiltonian Builder.

Instead of using generic fermionic operators, build Hamiltonians based on
the specific physics and governance protocol of each bond type:
- Covalent: Molecular orbital (bonding/antibonding) picture
- Ionic: Ionic site occupation picture
- Metallic: Band structure picture
"""

import numpy as np
from typing import Dict
from qiskit.quantum_info import SparsePauliOp


def build_covalent_mo_hamiltonian(
    molecular_orbitals: np.ndarray,
    mo_energies: np.ndarray,
    eri_mo: np.ndarray,
    nuclear_repulsion: float,
    n_electrons: int
) -> SparsePauliOp:
    """
    Build Hamiltonian in molecular orbital basis for covalent bonds.

    Uses direct MO picture: H = Σ ε_i n_i + Σ V_ijkl n_i n_j + E_nuc

    For H2 minimal basis:
    - MO 0 (σ bonding): ε_0 < 0
    - MO 1 (σ* antibonding): ε_1 > 0
    - 2 electrons → both in MO 0 (ground state)

    Args:
        molecular_orbitals: MO coefficient matrix [n_orbitals, n_ao]
        mo_energies: MO energies [n_orbitals]
        eri_mo: Two-electron integrals in MO basis [i,j,k,l]
        nuclear_repulsion: Nuclear repulsion energy
        n_electrons: Number of electrons

    Returns:
        Qiskit SparsePauliOp for the Hamiltonian
    """
    n_orbitals = len(mo_energies)
    n_qubits = n_orbitals * 2  # spin orbitals: 2 per spatial MO

    pauli_dict = {}

    # 1. One-electron terms: Σ ε_i (n_i↑ + n_i↓)
    for i in range(n_orbitals):
        eps_i = mo_energies[i]

        if abs(eps_i) > 1e-10:
            # Spin-up occupation: n_i↑ = (I - Z_i)/2
            i_up = 2 * i
            pauli_I = 'I' * n_qubits
            pauli_Z = list('I' * n_qubits)
            pauli_Z[i_up] = 'Z'
            pauli_Z = ''.join(pauli_Z)

            if pauli_I in pauli_dict:
                pauli_dict[pauli_I] += 0.5 * eps_i
            else:
                pauli_dict[pauli_I] = 0.5 * eps_i

            if pauli_Z in pauli_dict:
                pauli_dict[pauli_Z] += -0.5 * eps_i
            else:
                pauli_dict[pauli_Z] = -0.5 * eps_i

            # Spin-down occupation: n_i↓ = (I - Z_{i+1})/2
            i_down = 2 * i + 1
            pauli_Z = list('I' * n_qubits)
            pauli_Z[i_down] = 'Z'
            pauli_Z = ''.join(pauli_Z)

            if pauli_I in pauli_dict:
                pauli_dict[pauli_I] += 0.5 * eps_i
            else:
                pauli_dict[pauli_I] = 0.5 * eps_i

            if pauli_Z in pauli_dict:
                pauli_dict[pauli_Z] += -0.5 * eps_i
            else:
                pauli_dict[pauli_Z] = -0.5 * eps_i

    # 2. Two-electron terms: (1/2) Σ ⟨ij|kl⟩_MO n_i n_j
    # Using MO basis ERIs (already transformed)
    # For same-spin pairs: ⟨ij|ij⟩ n_i↑ n_j↑ (coulomb - exchange)
    # For opposite-spin: ⟨ij|ji⟩ n_i↑ n_j↓ (coulomb only)

    for i in range(n_orbitals):
        for j in range(n_orbitals):
            # Coulomb integral in MO basis: ⟨ii|jj⟩
            # This is eri_mo[i,i,j,j] for chemist notation
            J_ij = eri_mo[i, i, j, j]

            if abs(J_ij) > 1e-10:
                # n_i↑ n_j↑ term
                pauli_ii_jj_upup = _build_number_number_operator(
                    2 * i, 2 * j, n_qubits
                )
                for pauli_str, coeff in pauli_ii_jj_upup.items():
                    full_coeff = 0.5 * J_ij * coeff
                    if pauli_str in pauli_dict:
                        pauli_dict[pauli_str] += full_coeff
                    else:
                        pauli_dict[pauli_str] = full_coeff

                # n_i↓ n_j↓ term
                pauli_ii_jj_downdown = _build_number_number_operator(
                    2 * i + 1, 2 * j + 1, n_qubits
                )
                for pauli_str, coeff in pauli_ii_jj_downdown.items():
                    full_coeff = 0.5 * J_ij * coeff
                    if pauli_str in pauli_dict:
                        pauli_dict[pauli_str] += full_coeff
                    else:
                        pauli_dict[pauli_str] = full_coeff

                # n_i↑ n_j↓ term
                pauli_ii_jj_updown = _build_number_number_operator(
                    2 * i, 2 * j + 1, n_qubits
                )
                for pauli_str, coeff in pauli_ii_jj_updown.items():
                    full_coeff = 0.5 * J_ij * coeff
                    if pauli_str in pauli_dict:
                        pauli_dict[pauli_str] += full_coeff
                    else:
                        pauli_dict[pauli_str] = full_coeff

                # n_i↓ n_j↑ term
                pauli_ii_jj_downup = _build_number_number_operator(
                    2 * i + 1, 2 * j, n_qubits
                )
                for pauli_str, coeff in pauli_ii_jj_downup.items():
                    full_coeff = 0.5 * J_ij * coeff
                    if pauli_str in pauli_dict:
                        pauli_dict[pauli_str] += full_coeff
                    else:
                        pauli_dict[pauli_str] = full_coeff

    # 3. Nuclear repulsion (constant)
    identity = 'I' * n_qubits
    if identity in pauli_dict:
        pauli_dict[identity] += nuclear_repulsion
    else:
        pauli_dict[identity] = nuclear_repulsion

    # 4. Clean up and convert to SparsePauliOp
    pauli_dict = {k: v for k, v in pauli_dict.items() if abs(v) > 1e-12}

    if len(pauli_dict) == 0:
        return SparsePauliOp(['I' * n_qubits], [0.0])

    pauli_strings = list(pauli_dict.keys())
    coefficients = list(pauli_dict.values())

    return SparsePauliOp(pauli_strings, coefficients)


def _build_number_number_operator(i: int, j: int, n_qubits: int) -> Dict[str, float]:
    """
    Build n_i n_j = (I - Z_i)/2 × (I - Z_j)/2 = (II - IZ - ZI + ZZ)/4.

    Args:
        i: First qubit index
        j: Second qubit index
        n_qubits: Total number of qubits

    Returns:
        Dictionary of Pauli strings to coefficients
    """
    result = {}

    # II term: 1/4
    pauli_II = 'I' * n_qubits
    result[pauli_II] = 0.25

    # -IZ term: -1/4
    pauli_IZ = list('I' * n_qubits)
    pauli_IZ[j] = 'Z'
    pauli_IZ = ''.join(pauli_IZ)
    result[pauli_IZ] = -0.25

    # -ZI term: -1/4
    pauli_ZI = list('I' * n_qubits)
    pauli_ZI[i] = 'Z'
    pauli_ZI = ''.join(pauli_ZI)
    result[pauli_ZI] = -0.25

    # ZZ term: 1/4
    pauli_ZZ = list('I' * n_qubits)
    pauli_ZZ[i] = 'Z'
    pauli_ZZ[j] = 'Z'
    pauli_ZZ = ''.join(pauli_ZZ)
    result[pauli_ZZ] = 0.25

    return result
