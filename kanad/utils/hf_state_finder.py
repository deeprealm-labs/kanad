"""
Utility to find the Hartree-Fock state representation for different qubit mappings.

Different fermion-to-qubit mappings (Jordan-Wigner, Bravyi-Kitaev, Parity) use
different qubit encodings. This means the HF state has a different bit-string
representation in each encoding.

This module provides a function to find the correct HF state for a given mapper.
"""

import numpy as np
from typing import List, Tuple
import logging

logger = logging.getLogger(__name__)


def find_hf_state(
    hamiltonian_pauli,
    n_electrons: int,
    hf_energy: float,
    tolerance: float = 1e-5
) -> Tuple[int, List[int]]:
    """
    Find the Hartree-Fock state representation for a given Pauli Hamiltonian.

    Args:
        hamiltonian_pauli: SparsePauliOp Hamiltonian
        n_electrons: Number of electrons
        hf_energy: Expected HF energy (from SCF)
        tolerance: Energy tolerance for matching HF energy

    Returns:
        (state_index, occupation_list) where:
            state_index: Integer index in state vector
            occupation_list: Binary occupation [q0, q1, q2, ...] where 1 = occupied

    Raises:
        ValueError: If HF state not found in computational basis
    """
    n_qubits = hamiltonian_pauli.num_qubits
    H_matrix = hamiltonian_pauli.to_matrix()

    # Search all computational basis states with correct electron count
    for i in range(2**n_qubits):
        # Check if this state has the right number of electrons
        bits = format(i, f'0{n_qubits}b')
        if bits.count('1') != n_electrons:
            continue

        # Compute energy of this state
        state = np.zeros(2**n_qubits, dtype=complex)
        state[i] = 1.0
        energy = np.real(state.conj() @ H_matrix @ state)

        # Check if it matches HF energy
        if abs(energy - hf_energy) < tolerance:
            # Convert to occupation list (little-endian: rightmost = qubit 0)
            occupation = [int(bits[n_qubits - 1 - j]) for j in range(n_qubits)]
            logger.info(f"Found HF state: index {i} = |{bits}⟩, occupation {occupation}")
            return i, occupation

    # HF state not found in computational basis - it's a superposition
    raise ValueError(
        f"HF state (E={hf_energy:.6f} Ha) not found in computational basis!\n"
        f"This means the HF state is a superposition in this qubit encoding.\n"
        f"Mapper may not support simple initial state preparation."
    )


def get_hf_initial_state(
    hamiltonian,
    molecule,
    mapper: str = 'jordan_wigner'
) -> List[int]:
    """
    Get the correct HF initial state occupation for a given mapper.

    This function automatically determines which qubits should be occupied
    to prepare the HF state for the specified mapper.

    Args:
        hamiltonian: Molecular Hamiltonian object
        molecule: Molecule object
        mapper: Mapper type ('jordan_wigner' or 'bravyi_kitaev')

    Returns:
        occupation_list: [q0, q1, q2, ...] where 1 = occupied, 0 = empty

    Example:
        >>> occupation = get_hf_initial_state(ham, mol, 'bravyi_kitaev')
        >>> # occupation might be [0, 0, 0, 1] for BK instead of [1, 1, 0, 0] for JW
    """
    # Get HF energy from SCF
    hf_energy = hamiltonian.solve_scf()[1]

    # Build Pauli Hamiltonian with specified mapper
    if hasattr(hamiltonian, 'to_sparse_hamiltonian'):
        pauli_ham = hamiltonian.to_sparse_hamiltonian(mapper=mapper)
    else:
        raise AttributeError(f"Hamiltonian {type(hamiltonian)} doesn't support mapper parameter")

    # Find HF state
    _, occupation = find_hf_state(pauli_ham, molecule.n_electrons, hf_energy)

    logger.info(f"HF state for {mapper}: {occupation}")
    return occupation
