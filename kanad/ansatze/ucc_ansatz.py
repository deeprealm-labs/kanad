"""
Unitary Coupled Cluster (UCC) ansatz.

⚠️ DEPRECATED: This ansatz is deprecated and will be removed in v2.0.
Use CovalentGovernanceAnsatz, IonicGovernanceAnsatz, or HardwareEfficientAnsatz instead.

Reason: UCC family shows poor performance (0 mHa correlation) compared to governance
ansätze which achieve 49x better results on ionic molecules.
"""

from typing import List, Tuple, Optional
import numpy as np
import warnings
from kanad.ansatze.base_ansatz import BaseAnsatz, QuantumCircuit, Parameter


class UCCAnsatz(BaseAnsatz):
    """
    Unitary Coupled Cluster Singles and Doubles (UCCSD) ansatz.

    Based on traditional coupled cluster theory:
        |ψ⟩ = exp(T - T†) |HF⟩

    where T = T₁ + T₂ (singles + doubles excitations)

    ADVANTAGES:
    - Chemically inspired
    - Systematically improvable
    - Captures electron correlation

    DISADVANTAGES:
    - Many parameters (scales as O(N⁴) for UCCSD)
    - Deep circuits
    """

    def __init__(
        self,
        n_qubits: int,
        n_electrons: int,
        excitations: Optional[List[Tuple[List[int], List[int]]]] = None,
        include_singles: bool = True,
        include_doubles: bool = True
    ):
        """
        Initialize UCC ansatz.

        ⚠️ DEPRECATED: UCCAnsatz is deprecated and will be removed in v2.0.

        Recommended alternatives:
        - For ionic molecules: CovalentGovernanceAnsatz (49x better performance)
        - For covalent molecules: HardwareEfficientAnsatz (80% success rate)
        - For mixed character: AdaptiveGovernanceAnsatz

        Args:
            n_qubits: Number of qubits
            n_electrons: Number of electrons
            excitations: List of (occupied, virtual) excitation tuples
            include_singles: Include single excitations
            include_doubles: Include double excitations
        """
        warnings.warn(
            "UCCAnsatz is deprecated and will be removed in v2.0. "
            "Use CovalentGovernanceAnsatz (49x better on ionic molecules), "
            "IonicGovernanceAnsatz, or HardwareEfficientAnsatz instead. "
            "UCC family shows 0 mHa correlation in our implementation.",
            DeprecationWarning,
            stacklevel=2
        )
        super().__init__(n_qubits, n_electrons)
        self.include_singles = include_singles
        self.include_doubles = include_doubles

        if excitations is None:
            self.excitations = self._generate_excitations()
        else:
            self.excitations = excitations

    @property
    def n_parameters(self) -> int:
        """Number of variational parameters (one per excitation)."""
        return len(self.excitations)

    def _generate_excitations(self) -> List[Tuple[List[int], List[int]]]:
        """
        Generate all singles and doubles excitations.

        CRITICAL: Must determine occupied/virtual qubits from HF state,
        not assume they're the first n_electrons/2 qubits!

        Returns:
            List of (occupied_orbitals, virtual_orbitals) tuples
        """
        # Get HF state to determine which qubits are occupied
        hf_state = self._hartree_fock_state()

        # Find occupied and virtual qubits
        occupied = [i for i, occ in enumerate(hf_state) if occ == 1]
        virtual = [i for i, occ in enumerate(hf_state) if occ == 0]

        excitations = []

        # Singles: |occ⟩ → |virt⟩
        if self.include_singles:
            for i in occupied:
                for a in virtual:
                    excitations.append(([i], [a]))

        # Doubles: |occ₁ occ₂⟩ → |virt₁ virt₂⟩
        if self.include_doubles:
            for idx_i, i in enumerate(occupied):
                for j in occupied[idx_i + 1:]:
                    for idx_a, a in enumerate(virtual):
                        for b in virtual[idx_a + 1:]:
                            excitations.append(([i, j], [a, b]))

        return excitations

    def build_circuit(self, initial_state: Optional[List[int]] = None) -> QuantumCircuit:
        """
        Build UCCSD circuit.

        Circuit structure:
        1. Prepare Hartree-Fock reference state
        2. Apply single excitations
        3. Apply double excitations

        Args:
            initial_state: Initial state (default: HF state)

        Returns:
            UCC circuit
        """
        circuit = QuantumCircuit(self.n_qubits)

        # 1. Prepare reference state (Hartree-Fock)
        if initial_state is None:
            initial_state = self._hartree_fock_state()

        # CRITICAL: Qiskit uses little-endian qubit ordering (qubit 0 is rightmost)
        # So we need to reverse the occupation list when applying X gates
        # Example: HF state [1,1,0,0] means occupy spin-orbitals 0,1
        # In Qiskit little-endian, this is |1100⟩ which means qubits 2,3 are |1⟩
        # So we apply X to qubits at indices: n_qubits - 1 - i for each occupied i
        for i, occupation in enumerate(initial_state):
            if occupation == 1:
                # Map spin-orbital i to qubit (n_qubits - 1 - i) for little-endian
                qubit_index = self.n_qubits - 1 - i
                circuit.x(qubit_index)

        circuit.barrier()

        # 2. Apply excitation operators
        for exc_idx, (occupied, virtual) in enumerate(self.excitations):
            if len(occupied) == 1:
                # Single excitation
                self._apply_single_excitation(circuit, occupied[0], virtual[0], exc_idx)
            elif len(occupied) == 2:
                # Double excitation
                self._apply_double_excitation(
                    circuit,
                    occupied[0], occupied[1],
                    virtual[0], virtual[1],
                    exc_idx
                )

        self.circuit = circuit
        return circuit

    def _hartree_fock_state(self) -> List[int]:
        """
        Generate Hartree-Fock reference state for spin orbitals.

        CRITICAL: Uses OpenFermion/Qiskit spin-orbital ordering:
        - Spin orbitals are paired: [0↑, 0↓, 1↑, 1↓, 2↑, 2↓, ...]
        - Qubit index: 2*orbital_idx for spin-up, 2*orbital_idx+1 for spin-down

        For H2 (2 electrons, 2 spatial orbitals → 4 qubits):
        - HF state: both electrons in lowest orbital (orbital 0)
        - Qubits: [0↑, 0↓, 1↑, 1↓] = [1, 1, 0, 0]
        - This is |1100⟩ in qubit basis

        For closed-shell systems:
        - Fill lowest orbitals with electron pairs (↑↓)

        Returns:
            Occupation list for spin orbitals
        """
        state = [0] * self.n_qubits
        n_orbitals = self.n_qubits // 2

        # For closed-shell: assume even number of electrons, paired spins
        if self.n_electrons % 2 == 0:
            n_pairs = self.n_electrons // 2
            # Fill lowest n_pairs orbitals with both spin-up and spin-down
            for i in range(min(n_pairs, n_orbitals)):
                state[2*i] = 1      # Spin-up: qubit 2*i
                state[2*i + 1] = 1  # Spin-down: qubit 2*i+1
        else:
            # Open-shell: fill pairs first, then add unpaired electron
            n_pairs = self.n_electrons // 2
            for i in range(min(n_pairs, n_orbitals)):
                state[2*i] = 1      # Spin-up
                state[2*i + 1] = 1  # Spin-down

            # Add unpaired electron (spin-up in next orbital)
            if n_pairs < n_orbitals:
                state[2*n_pairs] = 1

        return state

    def _spin_orbital_to_qubit(self, spin_orbital_idx: int) -> int:
        """
        Map spin-orbital index to Qiskit qubit index.

        Accounts for Qiskit's little-endian qubit ordering.

        Args:
            spin_orbital_idx: Spin-orbital index (0 to n_qubits-1)

        Returns:
            Qubit index in Qiskit's little-endian convention
        """
        return self.n_qubits - 1 - spin_orbital_idx

    def _apply_single_excitation(
        self,
        circuit: QuantumCircuit,
        occ: int,
        virt: int,
        exc_idx: int
    ):
        """
        Apply single excitation operator using correct Jordan-Wigner gates.

        Implements: exp(θ (a†_virt a_occ - a†_occ a_virt))

        Args:
            circuit: Circuit to modify
            occ: Occupied spin-orbital index
            virt: Virtual spin-orbital index
            exc_idx: Excitation index
        """
        from kanad.ansatze.ucc_gates import apply_single_excitation_jw

        theta = Parameter(f'θ_s_{exc_idx}')

        # Map spin-orbital indices to qubit indices (little-endian)
        occ_qubit = self._spin_orbital_to_qubit(occ)
        virt_qubit = self._spin_orbital_to_qubit(virt)

        # Apply correct fermionic excitation using Jordan-Wigner
        apply_single_excitation_jw(
            circuit, virt_qubit, occ_qubit, theta, self.n_qubits
        )

    def _apply_double_excitation(
        self,
        circuit: QuantumCircuit,
        occ1: int,
        occ2: int,
        virt1: int,
        virt2: int,
        exc_idx: int
    ):
        """
        Apply double excitation operator using correct Jordan-Wigner gates.

        Implements: exp(θ (a†_virt1 a†_virt2 a_occ2 a_occ1 - a†_occ1 a†_occ2 a_virt2 a_virt1))

        Args:
            circuit: Circuit to modify
            occ1: First occupied spin-orbital index
            occ2: Second occupied spin-orbital index
            virt1: First virtual spin-orbital index
            virt2: Second virtual spin-orbital index
            exc_idx: Excitation index
        """
        from kanad.ansatze.ucc_gates import apply_double_excitation_jw

        theta = Parameter(f'θ_d_{exc_idx}')

        # Map spin-orbital indices to qubit indices (little-endian)
        occ1_q = self._spin_orbital_to_qubit(occ1)
        occ2_q = self._spin_orbital_to_qubit(occ2)
        virt1_q = self._spin_orbital_to_qubit(virt1)
        virt2_q = self._spin_orbital_to_qubit(virt2)

        # Apply correct fermionic double excitation using Jordan-Wigner
        apply_double_excitation_jw(
            circuit, virt1_q, virt2_q, occ2_q, occ1_q, theta, self.n_qubits
        )

    def get_excitation_list(self) -> List[Tuple[List[int], List[int]]]:
        """Get list of excitations."""
        return self.excitations

    def count_singles(self) -> int:
        """Count number of single excitations."""
        return sum(1 for exc in self.excitations if len(exc[0]) == 1)

    def count_doubles(self) -> int:
        """Count number of double excitations."""
        return sum(1 for exc in self.excitations if len(exc[0]) == 2)

    def __repr__(self) -> str:
        return (f"UCCAnsatz(n_qubits={self.n_qubits}, n_electrons={self.n_electrons}, "
                f"singles={self.count_singles()}, doubles={self.count_doubles()})")


class UCC_S_Ansatz(UCCAnsatz):
    """UCC with singles only (faster, less accurate)."""

    def __init__(self, n_qubits: int, n_electrons: int):
        super().__init__(n_qubits, n_electrons, include_singles=True, include_doubles=False)


class UCC_D_Ansatz(UCCAnsatz):
    """UCC with doubles only."""

    def __init__(self, n_qubits: int, n_electrons: int):
        super().__init__(n_qubits, n_electrons, include_singles=False, include_doubles=True)
