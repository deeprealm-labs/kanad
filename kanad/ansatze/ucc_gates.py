"""
Correct UCC excitation gate implementations using Pauli evolution.

This module implements fermionic single and double excitations by:
1. Using OpenFermion to get the Jordan-Wigner Pauli decomposition
2. Using Qiskit's PauliEvolutionGate to correctly exponentiate these operators

This ensures the gates create the correct quantum interference needed for VQE
to lower energy below Hartree-Fock.

References:
- OpenFermion: https://github.com/quantumlib/OpenFermion
- Qiskit Evolution: https://docs.quantum.ibm.com/api/qiskit/qiskit.circuit.library.PauliEvolutionGate
"""

import numpy as np
from typing import List
from kanad.ansatze.base_ansatz import QuantumCircuit, Parameter

# Import Qiskit evolution for correct Pauli exponentiation
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qiskit.synthesis import LieTrotter

# Import OpenFermion for correct Jordan-Wigner transformation
try:
    from openfermion.ops import FermionOperator
    from openfermion.transforms import jordan_wigner
    HAS_OPENFERMION = True
except ImportError:
    HAS_OPENFERMION = False


def _openfermion_to_qiskit_pauli(of_op, n_qubits: int) -> SparsePauliOp:
    """
    Convert OpenFermion QubitOperator to Qiskit SparsePauliOp.

    Args:
        of_op: OpenFermion QubitOperator
        n_qubits: Number of qubits

    Returns:
        Qiskit SparsePauliOp
    """
    from qiskit.quantum_info import SparsePauliOp

    pauli_list = []
    coeffs = []

    for term, coeff in of_op.terms.items():
        # Build Pauli string (Qiskit uses big-endian for Pauli strings)
        pauli_str = ['I'] * n_qubits

        for qubit_idx, pauli_op in term:
            # term is tuple of (qubit_index, 'X'/'Y'/'Z')
            pauli_str[n_qubits - 1 - qubit_idx] = pauli_op

        pauli_list.append(''.join(pauli_str))
        coeffs.append(coeff)

    if not pauli_list:
        # Identity operator
        pauli_list = ['I' * n_qubits]
        coeffs = [0.0]

    return SparsePauliOp(pauli_list, coeffs)


def apply_single_excitation_jw(
    circuit: QuantumCircuit,
    p: int,
    q: int,
    theta: Parameter,
    n_qubits: int
):
    """
    Apply single excitation operator exp(θ(a†_p a_q - a†_q a_p)) using Jordan-Wigner.

    This uses OpenFermion to get the correct Pauli decomposition and Qiskit's
    PauliEvolutionGate to exponentiate it correctly.

    Args:
        circuit: Quantum circuit to modify
        p: Creation orbital index (qubit index in Qiskit little-endian)
        q: Annihilation orbital index (qubit index in Qiskit little-endian)
        theta: Variational parameter (can be numeric or Parameter)
        n_qubits: Total number of qubits
    """
    if not HAS_OPENFERMION:
        raise ImportError("OpenFermion is required for UCC gates. Install with: pip install openfermion")

    if p == q:
        return  # No excitation

    # Convert from Qiskit little-endian qubit indices to spin-orbital indices
    # Qiskit: qubit_idx = n_qubits - 1 - spin_orbital_idx
    # So: spin_orbital_idx = n_qubits - 1 - qubit_idx
    p_spin = n_qubits - 1 - p
    q_spin = n_qubits - 1 - q

    # Create fermionic single excitation: a†_p a_q - a†_q a_p
    fermion_op = (
        FermionOperator(f'{p_spin}^ {q_spin}', 1.0) -
        FermionOperator(f'{q_spin}^ {p_spin}', 1.0)
    )

    # Apply Jordan-Wigner transformation
    qubit_op = jordan_wigner(fermion_op)

    # Convert to Qiskit Pauli operator
    pauli_op = _openfermion_to_qiskit_pauli(qubit_op, n_qubits)

    # Handle complex coefficients from Jordan-Wigner
    # The JW transformation gives: T = -i/2 (X_p Z... Y_q - Y_p Z... X_q)
    # We want to implement exp(theta * T_real) where T_real is the Hermitian form

    # Remove imaginary unit from coefficients (JW gives imaginary coefficients)
    # The excitation operator is: i*T where T is Hermitian
    # So we multiply by -i to get the Hermitian operator
    real_coeffs = []
    for coeff in pauli_op.coeffs:
        if abs(coeff.imag) > abs(coeff.real):
            # Coefficient is primarily imaginary
            real_coeffs.append(-1j * coeff)
        else:
            real_coeffs.append(coeff)

    # Create real Pauli operator
    pauli_op_real = SparsePauliOp(pauli_op.paulis, np.array(real_coeffs))

    # Handle parameter type - convert custom Parameter to value or Qiskit Parameter
    from kanad.ansatze.base_ansatz import Parameter as KanadParameter
    from qiskit.circuit import Parameter as QiskitParameter

    # Track mapping between Qiskit parameters and Kanad parameters
    qiskit_to_kanad_param = {}

    if isinstance(theta, KanadParameter):
        if theta._is_symbolic:
            # Create Qiskit Parameter
            qiskit_param = QiskitParameter(theta.name)
            # Map it back to the Kanad parameter for tracking
            qiskit_to_kanad_param[qiskit_param] = theta
            # CRITICAL FIX: PauliEvolutionGate uses exp(-i*t*H), but we want exp(+i*θ*H)
            # So negate the parameter
            time_param = -qiskit_param
        else:
            # Use bound value, negated
            time_param = -theta.value
    else:
        # Already numeric or Qiskit Parameter, negate it
        if isinstance(theta, (int, float)):
            time_param = -theta
        else:
            # Qiskit Parameter
            time_param = -theta

    # Now create evolution gate with the Hermitian operator and decompose it
    from qiskit import QuantumCircuit as QiskitCircuit

    evolution_gate = PauliEvolutionGate(pauli_op_real, time=time_param, synthesis=LieTrotter())

    # Create a temporary Qiskit circuit to decompose the evolution gate
    temp_circuit = QiskitCircuit(n_qubits)
    temp_circuit.append(evolution_gate, range(n_qubits))

    # Decompose to basic gates
    decomposed = temp_circuit.decompose().decompose()  # Double decompose for full expansion

    # Add the decomposed gates to our circuit
    for instruction in decomposed.data:
        gate = instruction.operation
        qubits = [decomposed.qubits.index(q) for q in instruction.qubits]

        gate_name = gate.name.lower()

        # Convert Qiskit parameters back to Kanad parameters for proper tracking
        converted_params = []
        for p in gate.params:
            if isinstance(p, (int, float, complex)):
                converted_params.append(p)
            elif isinstance(p, QiskitParameter):
                # Direct Qiskit Parameter - check if it maps to a Kanad parameter
                if p in qiskit_to_kanad_param:
                    converted_params.append(qiskit_to_kanad_param[p])
                else:
                    # Check by name
                    found = False
                    for qparam, kparam in qiskit_to_kanad_param.items():
                        if hasattr(p, 'name') and hasattr(qparam, 'name') and p.name == qparam.name:
                            converted_params.append(kparam)
                            found = True
                            break
                    if not found:
                        # Create new Kanad parameter
                        new_kanad_param = KanadParameter(p.name if hasattr(p, 'name') else f'param_{len(circuit.parameters)}')
                        qiskit_to_kanad_param[p] = new_kanad_param
                        converted_params.append(new_kanad_param)
            elif hasattr(p, 'parameters'):
                # This is a ParameterExpression (e.g., -theta)
                # Extract the underlying Qiskit parameters
                underlying_qparams = list(p.parameters)
                if len(underlying_qparams) == 1:
                    # Single parameter expression - find the corresponding Kanad parameter
                    underlying_qparam = underlying_qparams[0]
                    if underlying_qparam in qiskit_to_kanad_param:
                        converted_params.append(qiskit_to_kanad_param[underlying_qparam])
                    else:
                        # Check by name
                        found = False
                        for qparam, kparam in qiskit_to_kanad_param.items():
                            if qparam.name == underlying_qparam.name:
                                converted_params.append(kparam)
                                qiskit_to_kanad_param[underlying_qparam] = kparam
                                found = True
                                break
                        if not found:
                            # Create new Kanad parameter
                            new_kanad_param = KanadParameter(underlying_qparam.name)
                            qiskit_to_kanad_param[underlying_qparam] = new_kanad_param
                            converted_params.append(new_kanad_param)
                else:
                    # Multiple parameters in expression - not supported, just append the expression
                    converted_params.append(p)
            else:
                converted_params.append(p)

        # Map to our circuit methods
        if gate_name == 'h':
            circuit.h(qubits[0])
        elif gate_name == 'x':
            circuit.x(qubits[0])
        elif gate_name == 'y':
            circuit.y(qubits[0])
        elif gate_name == 'z':
            circuit.z(qubits[0])
        elif gate_name == 'rx':
            circuit.rx(converted_params[0], qubits[0])
        elif gate_name == 'ry':
            circuit.ry(converted_params[0], qubits[0])
        elif gate_name == 'rz':
            circuit.rz(converted_params[0], qubits[0])
        elif gate_name == 'cx' or gate_name == 'cnot':
            circuit.cx(qubits[0], qubits[1])
        elif gate_name == 'barrier':
            pass  # Skip barriers
        else:
            # For any unknown gate, add as generic gate
            circuit.add_gate(gate_name, qubits, converted_params if converted_params else None)


def apply_double_excitation_jw(
    circuit: QuantumCircuit,
    p: int,
    q: int,
    r: int,
    s: int,
    theta: Parameter,
    n_qubits: int
):
    """
    Apply double excitation operator exp(θ(a†_p a†_q a_r a_s - a†_s a†_r a_q a_p))
    using Jordan-Wigner.

    This uses OpenFermion to get the correct Pauli decomposition (8 terms for JW)
    and Qiskit's PauliEvolutionGate to exponentiate it correctly.

    Args:
        circuit: Quantum circuit to modify
        p, q: Creation orbital indices (qubits in Qiskit little-endian)
        r, s: Annihilation orbital indices (qubits in Qiskit little-endian)
        theta: Variational parameter
        n_qubits: Total number of qubits

    Note:
        For H2 with HF state |1100⟩:
        - Main double excitation: qubits (3,2,1,0) in little-endian
        - This corresponds to spin-orbitals (0,1,2,3) in physical ordering
        - Excites both electrons from orbital 0 to orbital 1
    """
    if not HAS_OPENFERMION:
        raise ImportError("OpenFermion is required for UCC gates. Install with: pip install openfermion")

    # Convert from Qiskit little-endian qubit indices to spin-orbital indices
    p_spin = n_qubits - 1 - p
    q_spin = n_qubits - 1 - q
    r_spin = n_qubits - 1 - r
    s_spin = n_qubits - 1 - s

    # Create fermionic double excitation: a†_p a†_q a_r a_s - a†_s a†_r a_q a_p
    # Note: OpenFermion uses normal ordering automatically
    fermion_op = (
        FermionOperator(f'{p_spin}^ {q_spin}^ {r_spin} {s_spin}', 1.0) -
        FermionOperator(f'{s_spin}^ {r_spin}^ {q_spin} {p_spin}', 1.0)
    )

    # Apply Jordan-Wigner transformation
    qubit_op = jordan_wigner(fermion_op)

    # Convert to Qiskit Pauli operator
    pauli_op = _openfermion_to_qiskit_pauli(qubit_op, n_qubits)

    # Handle complex coefficients from Jordan-Wigner (same as single excitation)
    real_coeffs = []
    for coeff in pauli_op.coeffs:
        if abs(coeff.imag) > abs(coeff.real):
            # Coefficient is primarily imaginary
            real_coeffs.append(-1j * coeff)
        else:
            real_coeffs.append(coeff)

    # Create real Pauli operator
    pauli_op_real = SparsePauliOp(pauli_op.paulis, np.array(real_coeffs))

    # Handle parameter type - convert custom Parameter to value or Qiskit Parameter
    from kanad.ansatze.base_ansatz import Parameter as KanadParameter
    from qiskit.circuit import Parameter as QiskitParameter

    # Track mapping between Qiskit parameters and Kanad parameters
    qiskit_to_kanad_param = {}

    if isinstance(theta, KanadParameter):
        if theta._is_symbolic:
            # Create Qiskit Parameter
            qiskit_param = QiskitParameter(theta.name)
            # Map it back to the Kanad parameter for tracking
            qiskit_to_kanad_param[qiskit_param] = theta
            # CRITICAL FIX: negate for correct evolution direction
            time_param = -qiskit_param
        else:
            # Use bound value, negated
            time_param = -theta.value
    else:
        # Already numeric or Qiskit Parameter, negate it
        if isinstance(theta, (int, float)):
            time_param = -theta
        else:
            # Qiskit Parameter
            time_param = -theta

    # Create evolution gate: exp(theta * T) where T is the double excitation operator
    from qiskit import QuantumCircuit as QiskitCircuit

    evolution_gate = PauliEvolutionGate(pauli_op_real, time=time_param, synthesis=LieTrotter())

    # Create a temporary Qiskit circuit to decompose the evolution gate
    temp_circuit = QiskitCircuit(n_qubits)
    temp_circuit.append(evolution_gate, range(n_qubits))

    # Decompose to basic gates
    decomposed = temp_circuit.decompose().decompose()  # Double decompose for full expansion

    # Add the decomposed gates to our circuit
    for instruction in decomposed.data:
        gate = instruction.operation
        qubits = [decomposed.qubits.index(q) for q in instruction.qubits]

        gate_name = gate.name.lower()

        # Convert Qiskit parameters back to Kanad parameters for proper tracking
        converted_params = []
        for p in gate.params:
            if isinstance(p, (int, float, complex)):
                converted_params.append(p)
            elif isinstance(p, QiskitParameter):
                # Direct Qiskit Parameter - check if it maps to a Kanad parameter
                if p in qiskit_to_kanad_param:
                    converted_params.append(qiskit_to_kanad_param[p])
                else:
                    # Check by name
                    found = False
                    for qparam, kparam in qiskit_to_kanad_param.items():
                        if hasattr(p, 'name') and hasattr(qparam, 'name') and p.name == qparam.name:
                            converted_params.append(kparam)
                            found = True
                            break
                    if not found:
                        # Create new Kanad parameter
                        new_kanad_param = KanadParameter(p.name if hasattr(p, 'name') else f'param_{len(circuit.parameters)}')
                        qiskit_to_kanad_param[p] = new_kanad_param
                        converted_params.append(new_kanad_param)
            elif hasattr(p, 'parameters'):
                # This is a ParameterExpression (e.g., -theta)
                # Extract the underlying Qiskit parameters
                underlying_qparams = list(p.parameters)
                if len(underlying_qparams) == 1:
                    # Single parameter expression - find the corresponding Kanad parameter
                    underlying_qparam = underlying_qparams[0]
                    if underlying_qparam in qiskit_to_kanad_param:
                        converted_params.append(qiskit_to_kanad_param[underlying_qparam])
                    else:
                        # Check by name
                        found = False
                        for qparam, kparam in qiskit_to_kanad_param.items():
                            if qparam.name == underlying_qparam.name:
                                converted_params.append(kparam)
                                qiskit_to_kanad_param[underlying_qparam] = kparam
                                found = True
                                break
                        if not found:
                            # Create new Kanad parameter
                            new_kanad_param = KanadParameter(underlying_qparam.name)
                            qiskit_to_kanad_param[underlying_qparam] = new_kanad_param
                            converted_params.append(new_kanad_param)
                else:
                    # Multiple parameters in expression - not supported, just append the expression
                    converted_params.append(p)
            else:
                converted_params.append(p)

        # Map to our circuit methods
        if gate_name == 'h':
            circuit.h(qubits[0])
        elif gate_name == 'x':
            circuit.x(qubits[0])
        elif gate_name == 'y':
            circuit.y(qubits[0])
        elif gate_name == 'z':
            circuit.z(qubits[0])
        elif gate_name == 'rx':
            circuit.rx(converted_params[0], qubits[0])
        elif gate_name == 'ry':
            circuit.ry(converted_params[0], qubits[0])
        elif gate_name == 'rz':
            circuit.rz(converted_params[0], qubits[0])
        elif gate_name == 'cx' or gate_name == 'cnot':
            circuit.cx(qubits[0], qubits[1])
        elif gate_name == 'barrier':
            pass  # Skip barriers
        else:
            # For any unknown gate, add as generic gate
            circuit.add_gate(gate_name, qubits, converted_params if converted_params else None)


def apply_fermionic_excitation_efficient(
    circuit: QuantumCircuit,
    occupied: List[int],
    virtual: List[int],
    theta: Parameter,
    n_qubits: int
):
    """
    Apply fermionic excitation with correct Pauli evolution.

    This dispatches to the appropriate single or double excitation function.

    Args:
        circuit: Quantum circuit
        occupied: List of occupied orbital indices (Qiskit qubit indices)
        virtual: List of virtual orbital indices (Qiskit qubit indices)
        theta: Variational parameter
        n_qubits: Total qubits
    """
    if len(occupied) == 1 and len(virtual) == 1:
        # Single excitation
        apply_single_excitation_jw(
            circuit, virtual[0], occupied[0], theta, n_qubits
        )
    elif len(occupied) == 2 and len(virtual) == 2:
        # Double excitation
        # Standard convention: virtual[0] > virtual[1], occupied[0] > occupied[1]
        apply_double_excitation_jw(
            circuit, virtual[0], virtual[1], occupied[1], occupied[0], theta, n_qubits
        )
    else:
        raise ValueError(f"Unsupported excitation: {len(occupied)}->{len(virtual)}")
