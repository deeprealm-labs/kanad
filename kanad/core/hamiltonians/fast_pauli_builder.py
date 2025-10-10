"""
Fast direct construction of Pauli operators from molecular integrals.

Uses Qiskit's FermionicOp for CORRECT Jordan-Wigner transformation.
Works for ALL bonding types with ZERO accuracy loss.
"""

import numpy as np
from qiskit.quantum_info import SparsePauliOp
import logging

logger = logging.getLogger(__name__)


def build_molecular_hamiltonian_pauli(
    h_core: np.ndarray,
    eri: np.ndarray,
    nuclear_repulsion: float,
    n_orbitals: int
) -> SparsePauliOp:
    """
    Build molecular Hamiltonian as Pauli operators using Qiskit's FermionicOp.

    This uses Qiskit's built-in Jordan-Wigner transformation which is:
    - Correct (handles all fermion operator cases)
    - Fast (optimized implementation)
    - Works for all bonding types

    Args:
        h_core: One-electron integrals
        eri: Two-electron integrals
        nuclear_repulsion: Nuclear repulsion energy
        n_orbitals: Number of spatial orbitals

    Returns:
        SparsePauliOp representing the Hamiltonian
    """
    # Use OpenFermion for Jordan-Wigner transformation (our validated implementation)
    # This avoids the qiskit-nature dependency which has been removed
    try:
        from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner
        logger.info("Using OpenFermion Jordan-Wigner transformation (validated reference)")

        # Use our OpenFermion-based implementation directly
        return openfermion_jordan_wigner(
            h_mo=h_core,
            eri_mo=eri,
            nuclear_repulsion=nuclear_repulsion,
            n_electrons=0  # Not needed for Hamiltonian construction
        )
    except ImportError as e:
        raise ImportError(
            f"OpenFermion not available: {e}\n"
            "Install with: pip install openfermion"
        )
