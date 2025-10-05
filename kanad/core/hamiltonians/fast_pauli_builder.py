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
    try:
        from qiskit.second_q.operators import FermionicOp
        from qiskit.second_q.mappers import JordanWignerMapper
    except ImportError:
        # Fallback for older Qiskit versions
        try:
            from qiskit_nature.second_q.operators import FermionicOp
            from qiskit_nature.second_q.mappers import JordanWignerMapper
        except ImportError:
            raise ImportError(
                "Qiskit 2.x required. Install with: pip install 'qiskit>=1.0'"
            )

    n_qubits = 2 * n_orbitals

    logger.info(f"Building Hamiltonian using Qiskit FermionicOp...")
    logger.info(f"  {n_orbitals} spatial orbitals → {n_qubits} spin orbitals")

    # Build fermionic Hamiltonian using Qiskit's format
    # Format: {label: coefficient} where label like "+_0 -_1" means a†_0 a_1

    fermionic_terms = {}

    # Nuclear repulsion (identity)
    fermionic_terms[""] = nuclear_repulsion

    # One-body terms: h[i,j] a†_i a_j
    for i in range(n_orbitals):
        for j in range(n_orbitals):
            h_ij = h_core[i, j]
            if abs(h_ij) > 1e-12:
                # Alpha spin
                label = f"+_{2*i} -_{2*j}"
                fermionic_terms[label] = fermionic_terms.get(label, 0) + h_ij
                # Beta spin
                label = f"+_{2*i+1} -_{2*j+1}"
                fermionic_terms[label] = fermionic_terms.get(label, 0) + h_ij

    # Two-body terms: 0.5 * g[ijkl] a†_i a†_j a_l a_k
    for i in range(n_orbitals):
        for j in range(n_orbitals):
            for k in range(n_orbitals):
                for l in range(n_orbitals):
                    g_ijkl = eri[i, k, j, l]
                    if abs(g_ijkl) > 1e-12:
                        coeff = 0.5 * g_ijkl

                        # All spin combinations
                        # Alpha-alpha
                        label = f"+_{2*i} +_{2*j} -_{2*l} -_{2*k}"
                        fermionic_terms[label] = fermionic_terms.get(label, 0) + coeff
                        # Alpha-beta
                        label = f"+_{2*i} +_{2*j+1} -_{2*l+1} -_{2*k}"
                        fermionic_terms[label] = fermionic_terms.get(label, 0) + coeff
                        # Beta-alpha
                        label = f"+_{2*i+1} +_{2*j} -_{2*l} -_{2*k+1}"
                        fermionic_terms[label] = fermionic_terms.get(label, 0) + coeff
                        # Beta-beta
                        label = f"+_{2*i+1} +_{2*j+1} -_{2*l+1} -_{2*k+1}"
                        fermionic_terms[label] = fermionic_terms.get(label, 0) + coeff

    logger.info(f"Built {len(fermionic_terms)} fermionic terms")

    # Create FermionicOp
    ferm_op = FermionicOp(fermionic_terms, num_spin_orbitals=n_qubits)

    # Map to Pauli operators using Jordan-Wigner
    mapper = JordanWignerMapper()
    pauli_op = mapper.map(ferm_op)

    # Simplify
    pauli_op = pauli_op.simplify(atol=1e-12)

    logger.info(f"✓ Final Pauli Hamiltonian: {len(pauli_op)} terms")

    return pauli_op
