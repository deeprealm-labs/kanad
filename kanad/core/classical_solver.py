"""
Classical diagonalization for Hi-VQE.

Instead of measuring all Pauli terms on quantum hardware, Hi-VQE:
1. Samples configurations from quantum state (Z measurement only)
2. Projects Hamiltonian into configuration subspace
3. Solves exactly using classical eigenvalue solver

This gives exact energy in the subspace with NO quantum measurement overhead!
"""

import numpy as np
from typing import Tuple, List, Optional
import logging
from qiskit.quantum_info import SparsePauliOp, Statevector

from kanad.core.configuration import Configuration, ConfigurationSubspace

logger = logging.getLogger(__name__)


class SubspaceHamiltonianBuilder:
    """
    Projects Hamiltonian into configuration subspace for classical diagonalization.

    Key insight: We only need matrix elements ⟨config_i|H|config_j⟩ for configs in subspace!
    """

    def __init__(self, hamiltonian: SparsePauliOp):
        """
        Initialize with Hamiltonian in Pauli form.

        Args:
            hamiltonian: Qiskit SparsePauliOp (sum of Pauli terms)
        """
        self.hamiltonian = hamiltonian
        self.n_qubits = hamiltonian.num_qubits

        logger.info(f"SubspaceHamiltonianBuilder: {len(hamiltonian)} Pauli terms, {self.n_qubits} qubits")

    def project(self, subspace: ConfigurationSubspace) -> np.ndarray:
        """
        Project Hamiltonian into configuration subspace.

        Returns:
            H_sub: n×n matrix where n = len(subspace)
                   H_sub[i,j] = ⟨config_i|H|config_j⟩
        """
        n = len(subspace)
        H_sub = np.zeros((n, n), dtype=complex)

        logger.info(f"Projecting Hamiltonian into subspace ({n} configurations)...")

        # Compute matrix elements
        for i in range(n):
            config_i = subspace[i]

            for j in range(i, n):  # Use symmetry: H is Hermitian
                config_j = subspace[j]

                # Compute ⟨config_i|H|config_j⟩
                matrix_element = self._compute_matrix_element(config_i, config_j)

                H_sub[i, j] = matrix_element
                if i != j:
                    H_sub[j, i] = np.conj(matrix_element)  # Hermitian symmetry

        logger.info(f"✓ Projected Hamiltonian: {n}×{n} matrix")

        return H_sub

    def _compute_matrix_element(self, config_i: Configuration, config_j: Configuration) -> complex:
        """
        Compute ⟨config_i|H|config_j⟩.

        For Hamiltonian H = Σ_k α_k P_k (sum of Pauli terms):
            ⟨i|H|j⟩ = Σ_k α_k ⟨i|P_k|j⟩

        We evaluate each Pauli term on the computational basis states.
        """
        # Create statevectors for configurations
        state_i = self._config_to_statevector(config_i)
        state_j = self._config_to_statevector(config_j)

        # Compute matrix element: ⟨i|H|j⟩ = ⟨i| H |j⟩
        # Using Qiskit: expectation = state_i.dagger() @ H @ state_j
        result = state_i.conjugate().T @ (self.hamiltonian.to_matrix() @ state_j)

        return complex(result)

    def _config_to_statevector(self, config: Configuration) -> np.ndarray:
        """
        Convert configuration to statevector.

        Configuration |110⟩ → statevector with 1.0 at index 6 (binary 110 = 6)
        """
        n_states = 2 ** self.n_qubits
        statevector = np.zeros(n_states, dtype=complex)
        statevector[config.to_int()] = 1.0
        return statevector

    def project_fast(self, subspace: ConfigurationSubspace) -> np.ndarray:
        """
        Fast projection using Pauli term evaluation (avoids full matrix construction).

        This is much faster for large systems!
        """
        n = len(subspace)
        H_sub = np.zeros((n, n), dtype=complex)

        # Validate qubit count matches
        if len(subspace) > 0:
            config_qubits = subspace[0].n_qubits
            if config_qubits != self.n_qubits:
                raise ValueError(
                    f"Qubit count mismatch: Hamiltonian has {self.n_qubits} qubits "
                    f"but configurations have {config_qubits} qubits. "
                    f"Ensure active space matches Hamiltonian construction."
                )

        logger.info(f"Fast projection into subspace ({n} configs, {len(self.hamiltonian)} Pauli terms)...")

        # For each Pauli term, compute contribution to all matrix elements
        for pauli_term in self.hamiltonian:
            pauli_str = pauli_term.paulis[0]  # Get Pauli string (e.g., "XXYZI")
            coeff = pauli_term.coeffs[0]  # Get coefficient

            # Compute ⟨config_i|P|config_j⟩ for this Pauli term
            for i in range(n):
                for j in range(i, n):
                    contrib = self._pauli_matrix_element(
                        subspace[i], subspace[j], pauli_str
                    )

                    H_sub[i, j] += coeff * contrib
                    if i != j:
                        H_sub[j, i] += np.conj(coeff * contrib)

        logger.info(f"✓ Fast projected Hamiltonian: {n}×{n} matrix")

        return H_sub

    def _pauli_matrix_element(self, config_i: Configuration, config_j: Configuration, pauli_str) -> complex:
        """
        Compute ⟨config_i|P|config_j⟩ for a Pauli string P.

        Pauli matrices in computational basis:
        - I: ⟨0|I|0⟩=1, ⟨1|I|1⟩=1, off-diagonal=0
        - Z: ⟨0|Z|0⟩=1, ⟨1|Z|1⟩=-1, off-diagonal=0
        - X: ⟨0|X|1⟩=1, ⟨1|X|0⟩=1, diagonal=0
        - Y: ⟨0|Y|1⟩=-i, ⟨1|Y|0⟩=i, diagonal=0

        For product P = P_0 ⊗ P_1 ⊗ ... ⊗ P_n:
            ⟨i|P|j⟩ = ∏_k ⟨i_k|P_k|j_k⟩

        IMPORTANT: Qiskit uses little-endian qubit ordering:
            - Bitstring "1100" means q3=1, q2=1, q1=0, q0=0 (rightmost is q0)
            - Pauli "ZIII" acts as Z₃⊗I₂⊗I₁⊗I₀ (leftmost acts on highest qubit)
        """
        bits_i = config_i.bitstring
        bits_j = config_j.bitstring

        result = 1.0 + 0j

        # Convert Pauli string to list (leftmost character is highest qubit index)
        pauli_list = list(str(pauli_str))

        # Qiskit convention: both Pauli and bitstring use left-to-right = high-to-low qubit index
        # Pauli "ZIII" position [0,1,2,3] acts on qubits [n-1, n-2, ..., 1, 0]
        # Bitstring "1100" position [0,1,2,3] represents qubits [n-1, n-2, ..., 1, 0]
        # Therefore: pauli_list[k] acts on the same qubit as bits[k]
        for k, pauli in enumerate(pauli_list):
            bit_i = bits_i[k]
            bit_j = bits_j[k]

            if pauli == 'I':
                # Identity: ⟨i|I|j⟩ = δ_ij (Kronecker delta)
                if bit_i != bit_j:
                    return 0.0  # Orthogonal states
                # else: contributes 1.0 (continue to next qubit)
                continue

            elif pauli == 'Z':
                if bit_i != bit_j:
                    return 0.0  # Off-diagonal of Z is zero
                elif bit_i == '1':
                    result *= -1.0  # ⟨1|Z|1⟩ = -1

            elif pauli == 'X':
                if bit_i == bit_j:
                    return 0.0  # Diagonal of X is zero
                # ⟨0|X|1⟩ = ⟨1|X|0⟩ = 1
                # result *= 1.0 (no change)

            elif pauli == 'Y':
                if bit_i == bit_j:
                    return 0.0  # Diagonal of Y is zero
                # ⟨0|Y|1⟩ = -i, ⟨1|Y|0⟩ = i
                if bit_i == '0' and bit_j == '1':
                    result *= -1j
                else:  # bit_i == '1' and bit_j == '0'
                    result *= 1j

        return result


def diagonalize_subspace(
    H_sub: np.ndarray,
    n_states: int = 1
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Diagonalize subspace Hamiltonian.

    Args:
        H_sub: Projected Hamiltonian matrix (n×n)
        n_states: Number of lowest eigenstates to return

    Returns:
        Tuple of (energies, eigenvectors)
        - energies: Array of n_states lowest eigenvalues
        - eigenvectors: Matrix of eigenvectors (n × n_states)
    """
    logger.info(f"Diagonalizing {H_sub.shape[0]}×{H_sub.shape[0]} subspace Hamiltonian...")

    # Hermitian eigenvalue problem
    eigenvalues, eigenvectors = np.linalg.eigh(H_sub)

    # Sort by energy (should already be sorted, but ensure it)
    sorted_indices = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[sorted_indices]
    eigenvectors = eigenvectors[:, sorted_indices]

    # Return requested number of states
    energies = eigenvalues[:n_states]
    vectors = eigenvectors[:, :n_states]

    logger.info(f"✓ Ground state energy: {energies[0]:.8f} Ha")

    if n_states > 1:
        logger.info(f"  Excited state energies: {energies[1:n_states]}")

    return energies, vectors


def compute_subspace_energy(
    hamiltonian: SparsePauliOp,
    subspace: ConfigurationSubspace,
    use_fast: bool = True
) -> Tuple[float, np.ndarray]:
    """
    Compute exact ground state energy in configuration subspace.

    This is the core of Hi-VQE: instead of measuring Pauli terms,
    we project H into subspace and solve classically!

    Args:
        hamiltonian: Qiskit SparsePauliOp
        subspace: ConfigurationSubspace with sampled configs
        use_fast: Use fast Pauli-term evaluation (recommended)

    Returns:
        Tuple of (ground_energy, ground_state_amplitudes)
    """
    # Project Hamiltonian
    builder = SubspaceHamiltonianBuilder(hamiltonian)

    if use_fast:
        H_sub = builder.project_fast(subspace)
    else:
        H_sub = builder.project(subspace)

    # Diagonalize
    energies, eigenvectors = diagonalize_subspace(H_sub, n_states=1)

    ground_energy = energies[0]
    ground_amplitudes = eigenvectors[:, 0]

    return ground_energy, ground_amplitudes


def get_important_configurations(
    subspace: ConfigurationSubspace,
    amplitudes: np.ndarray,
    threshold: float = 0.1
) -> List[Tuple[Configuration, float]]:
    """
    Get configurations with large amplitudes.

    These are the important configurations to generate excitations from.

    Args:
        subspace: Configuration subspace
        amplitudes: Amplitudes of configurations in ground state
        threshold: Return configs with |amplitude| > threshold

    Returns:
        List of (config, amplitude) tuples sorted by |amplitude|
    """
    important = []

    for i, config in enumerate(subspace):
        amp = amplitudes[i]
        if abs(amp) > threshold:
            important.append((config, amp))

    # Sort by amplitude magnitude
    important.sort(key=lambda x: abs(x[1]), reverse=True)

    logger.info(f"Found {len(important)} important configurations (threshold={threshold})")

    return important


def select_configurations_by_gradient(
    hamiltonian: SparsePauliOp,
    subspace: ConfigurationSubspace,
    ground_state: np.ndarray,
    candidate_pool: List[Configuration],
    k: int = 2
) -> List[Tuple[Configuration, float]]:
    """
    Select top-k configurations from pool based on energy gradient.

    This is the KEY to matching literature Hi-VQE accuracy!

    Gradient for configuration |φ⟩:
        ∇E_φ = 2 * |⟨ψ_current|H|φ⟩|

    Where ψ_current is current ground state approximation in subspace.

    Physical interpretation:
    - Large gradient → strong coupling with current ground state
    - Adding this config will significantly lower the energy
    - This is how literature Hi-VQE achieves <1 mHa error!

    Args:
        hamiltonian: Molecular Hamiltonian
        subspace: Current configuration subspace
        ground_state: Current ground state amplitudes in subspace
        candidate_pool: Pool of candidate configurations (not yet in subspace)
        k: Number of configurations to select

    Returns:
        List of (config, gradient) tuples for top-k candidates
        Sorted by gradient magnitude (largest first)

    Example:
        For H2, after HF reference:
        - |0011⟩ (double exc) has gradient ≈ 0.18 Ha
        - |1010⟩ (single exc) has gradient ≈ 0.00 Ha (no coupling!)
        → Select |0011⟩, skip singles, get exact energy!
    """
    logger.info(f"Computing gradients for {len(candidate_pool)} candidate configurations...")

    builder = SubspaceHamiltonianBuilder(hamiltonian)

    # Current ground state in subspace: |ψ⟩ = Σ_i c_i |config_i⟩
    # We have amplitudes c_i in `ground_state`

    gradients = []

    for candidate in candidate_pool:
        # Skip if already in subspace
        if candidate in subspace:
            continue

        # Compute ⟨ψ|H|candidate⟩ = Σ_i c_i * ⟨config_i|H|candidate⟩
        coupling = 0.0

        for i, config_i in enumerate(subspace):
            # Compute matrix element ⟨config_i|H|candidate⟩
            matrix_elem = 0.0

            for pauli_term in hamiltonian:
                pauli_str = pauli_term.paulis[0]
                coeff = pauli_term.coeffs[0]

                contrib = builder._pauli_matrix_element(config_i, candidate, pauli_str)
                matrix_elem += coeff * contrib

            # Weight by ground state amplitude
            coupling += ground_state[i] * matrix_elem

        # Gradient magnitude (factor of 2 from perturbation theory)
        gradient = 2 * abs(coupling)

        gradients.append((candidate, float(gradient)))

    # Sort by gradient magnitude (largest first)
    gradients.sort(key=lambda x: x[1], reverse=True)

    # Select top-k
    selected = gradients[:k]

    logger.info(f"Selected top-{k} configurations by gradient:")
    for config, grad in selected:
        logger.info(f"  {config}: gradient = {grad:.8f} Ha")

    return selected
