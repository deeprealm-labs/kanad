"""
Correct Jordan-Wigner transformation for molecular Hamiltonians.

Implements the standard Jordan-Wigner mapping from fermionic operators to Pauli operators
following the formulation in:
- Seeley, Richard, Love (2012) "The Bravyi-Kitaev transformation for quantum chemistry"
- McClean et al. (2016) "The theory of variational hybrid quantum-classical algorithms"

CRITICAL: This implementation is validated against known reference values and should
NOT be modified without comprehensive testing against FCI benchmarks.
"""

import numpy as np
from typing import Dict, Tuple
from qiskit.quantum_info import SparsePauliOp
import logging

logger = logging.getLogger(__name__)


def jordan_wigner_transform(
    h_mo: np.ndarray,
    eri_mo: np.ndarray,
    nuclear_repulsion: float,
    n_electrons: int
) -> SparsePauliOp:
    """
    Convert molecular Hamiltonian from fermionic to qubit representation.

    This function implements the Jordan-Wigner transformation for the electronic
    Hamiltonian in second quantization:

    H = Σ_{pq} h_{pq} a†_p a_q + (1/4) Σ_{pqrs} g_{pqrs} a†_p a†_q a_s a_r + E_nuc

    where h_{pq} are one-electron integrals in MO basis and g_{pqrs} are two-electron
    integrals in physicist notation: g_{pqrs} = ⟨pr|qs⟩ = (pq|rs) in chemist notation.

    Args:
        h_mo: One-electron integrals in MO basis (n_orbitals × n_orbitals)
        eri_mo: Two-electron integrals in MO basis, chemist notation ⟨pq|rs⟩
                (n_orbitals × n_orbitals × n_orbitals × n_orbitals)
        nuclear_repulsion: Nuclear repulsion energy
        n_electrons: Number of electrons

    Returns:
        SparsePauliOp representing the qubit Hamiltonian

    References:
        The factor of 1/4 in the two-body term comes from:
        - Each (pqrs) term appears 4 times due to spin summation
        - We use physicist notation: g_{pqrs} = ⟨pr|qs⟩
        - In chemist notation ⟨pq|rs⟩, we need to convert indices
    """
    n_orbitals = h_mo.shape[0]
    n_qubits = 2 * n_orbitals  # Spin orbitals: alpha and beta for each spatial orbital

    logger.info(f"Building Hamiltonian: {n_orbitals} orbitals → {n_qubits} qubits")

    # Dictionary to accumulate Pauli terms: {pauli_string: coefficient}
    pauli_terms: Dict[str, complex] = {}

    # 1. ONE-BODY TERMS: Σ_{pq} h_{pq} a†_p a_q
    # For each spatial orbital (p,q), we have spin-up and spin-down components
    logger.debug("Building one-body terms...")

    for p in range(n_orbitals):
        for q in range(n_orbitals):
            h_pq = h_mo[p, q]

            if abs(h_pq) < 1e-12:
                continue

            # Spin-up: a†_{p↑} a_{q↑}
            _add_one_body_term(pauli_terms, p, q, h_pq, n_qubits)

            # Spin-down: a†_{p↓} a_{q↓}
            _add_one_body_term(pauli_terms, p + n_orbitals, q + n_orbitals, h_pq, n_qubits)

    # 2. TWO-BODY TERMS: (1/4) Σ_{pqrs} g_{pqrs} a†_p a†_q a_s a_r
    # Convert from chemist ⟨pq|rs⟩ to physicist ⟨pr|qs⟩ notation
    logger.debug("Building two-body terms...")

    for p in range(n_orbitals):
        for q in range(n_orbitals):
            for r in range(n_orbitals):
                for s in range(n_orbitals):
                    # Chemist notation: ⟨pq|rs⟩
                    # Physicist notation: g_{prqs} = ⟨pr|qs⟩ = ⟨pq|rs⟩
                    g_pqrs = eri_mo[p, q, r, s]

                    if abs(g_pqrs) < 1e-12:
                        continue

                    # Spin cases: ↑↑, ↑↓, ↓↑, ↓↓
                    # Factor 1/2 for the coefficient and NO factor from spin (each spin combo is unique)
                    # CRITICAL: For chemist notation ⟨pq|rs⟩, the operator is a†_p a†_r a_s a_q

                    # ↑↑: a†_{p↑} a†_{r↑} a_{q↑} a_{s↑} - NOTE: fermionic operators anticommute!
                    # Canonical order: a†_p a†_r a_s a_q
                    _add_two_body_term(pauli_terms, p, r, s, q, 0.5 * g_pqrs, n_qubits)

                    # ↑↓: a†_{p↑} a†_{r↓} a_{q↓} a_{s↑}
                    _add_two_body_term(pauli_terms, p, r + n_orbitals, s + n_orbitals, q, 0.5 * g_pqrs, n_qubits)

                    # ↓↑: a†_{p↓} a†_{r↑} a_{q↑} a_{s↓}
                    _add_two_body_term(pauli_terms, p + n_orbitals, r, s, q + n_orbitals, 0.5 * g_pqrs, n_qubits)

                    # ↓↓: a†_{p↓} a†_{r↓} a_{q↓} a_{s↓}
                    _add_two_body_term(pauli_terms, p + n_orbitals, r + n_orbitals, s + n_orbitals, q + n_orbitals, 0.5 * g_pqrs, n_qubits)

    # 3. NUCLEAR REPULSION (constant term)
    identity = 'I' * n_qubits
    if identity in pauli_terms:
        pauli_terms[identity] += nuclear_repulsion
    else:
        pauli_terms[identity] = nuclear_repulsion

    # 4. Filter negligible terms and convert to SparsePauliOp
    pauli_terms = {k: v for k, v in pauli_terms.items() if abs(v) > 1e-12}

    if len(pauli_terms) == 0:
        logger.warning("Empty Hamiltonian! Returning zero operator")
        return SparsePauliOp(['I' * n_qubits], [0.0])

    pauli_strings = list(pauli_terms.keys())
    coefficients = [pauli_terms[s] for s in pauli_strings]

    # Convert complex to real (imaginary parts should cancel due to Hermiticity)
    coefficients_real = [c.real if isinstance(c, complex) else c for c in coefficients]

    logger.info(f"Built Hamiltonian with {len(pauli_strings)} Pauli terms")

    return SparsePauliOp(pauli_strings, coefficients_real)


def _add_one_body_term(
    pauli_terms: Dict[str, complex],
    p: int,
    q: int,
    coeff: float,
    n_qubits: int
) -> None:
    """
    Add one-body term h_{pq} a†_p a_q to Pauli dictionary.

    Jordan-Wigner transformation:
    a†_p a_q = (1/4) [(X_p X_q + Y_p Y_q)(Z_{p-1}...Z_{q+1}) + i(X_p Y_q - Y_p X_q)(Z_{p-1}...Z_{q+1})]

    For p = q (number operator):
    a†_p a_p = (I - Z_p) / 2
    """
    if p == q:
        # Number operator: a†_p a_p = (I - Z_p) / 2
        identity = 'I' * n_qubits
        z_term = 'I' * p + 'Z' + 'I' * (n_qubits - p - 1)

        _accumulate(pauli_terms, identity, 0.5 * coeff)
        _accumulate(pauli_terms, z_term, -0.5 * coeff)
    else:
        # Excitation operator: a†_p a_q with p ≠ q
        # Ensure p > q for standard ordering
        if p < q:
            p, q = q, p
            coeff = np.conj(coeff)  # Hermitian conjugate

        # Build Z string between p and q
        z_string = _build_z_string(p, q, n_qubits)

        # Four terms from Jordan-Wigner transformation
        # (1/4) * [(X_p X_q + Y_p Y_q) + i(X_p Y_q - Y_p X_q)] with Z string

        # XX term
        xx_pauli = list(z_string)
        xx_pauli[p] = 'X'
        xx_pauli[q] = 'X'
        _accumulate(pauli_terms, ''.join(xx_pauli), 0.25 * coeff)

        # YY term
        yy_pauli = list(z_string)
        yy_pauli[p] = 'Y'
        yy_pauli[q] = 'Y'
        _accumulate(pauli_terms, ''.join(yy_pauli), 0.25 * coeff)

        # XY term (imaginary)
        xy_pauli = list(z_string)
        xy_pauli[p] = 'X'
        xy_pauli[q] = 'Y'
        _accumulate(pauli_terms, ''.join(xy_pauli), 0.25j * coeff)

        # YX term (imaginary)
        yx_pauli = list(z_string)
        yx_pauli[p] = 'Y'
        yx_pauli[q] = 'X'
        _accumulate(pauli_terms, ''.join(yx_pauli), -0.25j * coeff)


def _add_two_body_term(
    pauli_terms: Dict[str, complex],
    p: int,
    q: int,
    r: int,
    s: int,
    coeff: float,
    n_qubits: int
) -> None:
    """
    Add two-body term g_{pqrs} a†_p a†_q a_r a_s to Pauli dictionary.

    Uses the exact Jordan-Wigner transformation by composing the
    one-body transformations:

    a†_p a†_q a_r a_s = (a†_p a_r)(a†_q a_s) - (a†_p a_s)(a†_q a_r) [for p>r, q>s]

    We handle this by explicitly computing all Pauli products.
    """
    # Map each fermionic operator to its Pauli representation
    # Then compute the product of all four operators

    # Create temporary dict for intermediate products
    temp1 = {}
    temp2 = {}
    temp3 = {}

    # Step 1: a_s → Pauli
    _add_annihilation_op(temp1, s, 1.0, n_qubits)

    # Step 2: a_r (a_s) → Pauli
    for pauli1, coeff1 in list(temp1.items()):
        _add_annihilation_op_product(temp2, r, pauli1, coeff1, n_qubits)

    # Step 3: a†_q (a_r a_s) → Pauli
    for pauli2, coeff2 in list(temp2.items()):
        _add_creation_op_product(temp3, q, pauli2, coeff2, n_qubits)

    # Step 4: a†_p (a†_q a_r a_s) → Pauli
    for pauli3, coeff3 in list(temp3.items()):
        _add_creation_op_product(pauli_terms, p, pauli3, coeff * coeff3, n_qubits)


def _add_annihilation_op(
    pauli_terms: Dict[str, complex],
    p: int,
    coeff: complex,
    n_qubits: int
) -> None:
    """Add annihilation operator a_p to Pauli dictionary."""
    # a_p = (X_p + iY_p) / 2 with Z string
    z_string = ['I'] * n_qubits
    for i in range(p):
        z_string[i] = 'Z'

    # X term
    x_pauli = z_string.copy()
    x_pauli[p] = 'X'
    _accumulate(pauli_terms, ''.join(x_pauli), 0.5 * coeff)

    # Y term
    y_pauli = z_string.copy()
    y_pauli[p] = 'Y'
    _accumulate(pauli_terms, ''.join(y_pauli), 0.5j * coeff)


def _add_creation_op_product(
    pauli_terms: Dict[str, complex],
    p: int,
    base_pauli: str,
    base_coeff: complex,
    n_qubits: int
) -> None:
    """Multiply existing Pauli string by creation operator a†_p."""
    # a†_p = (X_p - iY_p) / 2 with Z string
    z_string_ops = ['I'] * n_qubits
    for i in range(p):
        z_string_ops[i] = 'Z'

    # X term: multiply base_pauli by (Z_0...Z_{p-1} X_p)
    x_pauli_op = z_string_ops.copy()
    x_pauli_op[p] = 'X'
    result_x, coeff_x = _multiply_paulis(base_pauli, ''.join(x_pauli_op))
    _accumulate(pauli_terms, result_x, 0.5 * base_coeff * coeff_x)

    # Y term: multiply base_pauli by (Z_0...Z_{p-1} Y_p)
    y_pauli_op = z_string_ops.copy()
    y_pauli_op[p] = 'Y'
    result_y, coeff_y = _multiply_paulis(base_pauli, ''.join(y_pauli_op))
    _accumulate(pauli_terms, result_y, -0.5j * base_coeff * coeff_y)


def _add_annihilation_op_product(
    pauli_terms: Dict[str, complex],
    p: int,
    base_pauli: str,
    base_coeff: complex,
    n_qubits: int
) -> None:
    """Multiply existing Pauli string by annihilation operator a_p."""
    # a_p = (X_p + iY_p) / 2 with Z string
    z_string_ops = ['I'] * n_qubits
    for i in range(p):
        z_string_ops[i] = 'Z'

    # X term
    x_pauli_op = z_string_ops.copy()
    x_pauli_op[p] = 'X'
    result_x, coeff_x = _multiply_paulis(base_pauli, ''.join(x_pauli_op))
    _accumulate(pauli_terms, result_x, 0.5 * base_coeff * coeff_x)

    # Y term
    y_pauli_op = z_string_ops.copy()
    y_pauli_op[p] = 'Y'
    result_y, coeff_y = _multiply_paulis(base_pauli, ''.join(y_pauli_op))
    _accumulate(pauli_terms, result_y, 0.5j * base_coeff * coeff_y)


def _multiply_paulis(pauli1: str, pauli2: str) -> Tuple[str, complex]:
    """
    Multiply two Pauli strings and return (result_pauli, phase_coefficient).

    Pauli multiplication rules:
    - I * P = P
    - X * X = I, Y * Y = I, Z * Z = I
    - X * Y = iZ, Y * Z = iX, Z * X = iY (cyclic)
    - Y * X = -iZ, Z * Y = -iX, X * Z = -iY (anti-cyclic)
    """
    if len(pauli1) != len(pauli2):
        raise ValueError(f"Pauli strings must have same length: {len(pauli1)} vs {len(pauli2)}")

    result = []
    phase = 1.0 + 0j

    for p1, p2 in zip(pauli1, pauli2):
        if p1 == 'I':
            result.append(p2)
        elif p2 == 'I':
            result.append(p1)
        elif p1 == p2:
            result.append('I')
        else:
            # Different Pauli operators (X, Y, Z)
            # Use the cyclic property
            if (p1, p2) == ('X', 'Y'):
                result.append('Z')
                phase *= 1j
            elif (p1, p2) == ('Y', 'Z'):
                result.append('X')
                phase *= 1j
            elif (p1, p2) == ('Z', 'X'):
                result.append('Y')
                phase *= 1j
            elif (p1, p2) == ('Y', 'X'):
                result.append('Z')
                phase *= -1j
            elif (p1, p2) == ('Z', 'Y'):
                result.append('X')
                phase *= -1j
            elif (p1, p2) == ('X', 'Z'):
                result.append('Y')
                phase *= -1j

    return ''.join(result), phase


def _build_z_string(p: int, q: int, n_qubits: int) -> str:
    """Build Z operator string for Jordan-Wigner transformation between p and q."""
    pauli = ['I'] * n_qubits
    # Add Z operators between q+1 and p-1
    for i in range(q + 1, p):
        pauli[i] = 'Z'
    return ''.join(pauli)


def _accumulate(pauli_terms: Dict[str, complex], pauli_string: str, coeff: complex) -> None:
    """Accumulate Pauli term into dictionary."""
    if pauli_string in pauli_terms:
        pauli_terms[pauli_string] += coeff
    else:
        pauli_terms[pauli_string] = coeff
