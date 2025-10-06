"""
Unitary Coupled Cluster (UCC) ansatz.

Gold standard ansatz for quantum chemistry VQE.
"""

from typing import List, Tuple, Optional
import numpy as np
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

        Args:
            n_qubits: Number of qubits
            n_electrons: Number of electrons
            excitations: List of (occupied, virtual) excitation tuples
            include_singles: Include single excitations
            include_doubles: Include double excitations
        """
        super().__init__(n_qubits, n_electrons)
        self.include_singles = include_singles
        self.include_doubles = include_doubles

        if excitations is None:
            self.excitations = self._generate_excitations()
        else:
            self.excitations = excitations

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

        for qubit, occupation in enumerate(initial_state):
            if occupation == 1:
                circuit.x(qubit)

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

        Uses interleaved spin ordering (Jordan-Wigner convention):
        - Qubits 0, 1, 2, ..., N-1: spin-up orbitals (orbital 0↑, 1↑, 2↑, ...)
        - Qubits N, N+1, ..., 2N-1: spin-down orbitals (orbital 0↓, 1↓, 2↓, ...)

        For H2 (2 electrons, 2 orbitals → 4 qubits):
        - Correct HF state: [1,0,1,0] (qubit 0=orb0↑, qubit 2=orb0↓)

        Returns:
            Occupation list for spin orbitals
        """
        state = [0] * self.n_qubits
        n_orbitals = self.n_qubits // 2

        # Calculate spin-up and spin-down electron counts
        n_up = (self.n_electrons + 1) // 2  # Ceiling division
        n_down = self.n_electrons // 2       # Floor division

        # Fill spin-up orbitals (qubits 0, 1, 2, ...)
        for i in range(min(n_up, n_orbitals)):
            state[i] = 1

        # Fill spin-down orbitals (qubits n_orbitals, n_orbitals+1, ...)
        for i in range(min(n_down, n_orbitals)):
            state[n_orbitals + i] = 1

        return state

    def _apply_single_excitation(
        self,
        circuit: QuantumCircuit,
        occ: int,
        virt: int,
        exc_idx: int
    ):
        """
        Apply single excitation operator.

        Implements: exp(θ (a†_virt a_occ - a†_occ a_virt))

        Using Givens rotation decomposition.

        Args:
            circuit: Circuit to modify
            occ: Occupied orbital
            virt: Virtual orbital
            exc_idx: Excitation index
        """
        theta = Parameter(f'θ_s_{exc_idx}')

        # Givens rotation: mixing occupied and virtual orbitals
        # Implemented as RY rotation between qubits

        # Entangle the orbitals
        circuit.cx(occ, virt)
        circuit.ry(theta, virt)
        circuit.cx(occ, virt)

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
        Apply double excitation operator.

        Implements: exp(θ (a†_virt1 a†_virt2 a_occ2 a_occ1 - h.c.))

        Args:
            circuit: Circuit to modify
            occ1: First occupied orbital
            occ2: Second occupied orbital
            virt1: First virtual orbital
            virt2: Second virtual orbital
            exc_idx: Excitation index
        """
        theta = Parameter(f'θ_d_{exc_idx}')

        # Simplified double excitation using RY rotations
        # Full implementation would use fermionic simulation

        # Create entanglement pattern
        circuit.cx(occ1, occ2)
        circuit.cx(virt1, virt2)
        circuit.cx(occ2, virt1)

        # Parametrized rotation
        circuit.ry(theta, virt1)

        # Uncompute entanglement
        circuit.cx(occ2, virt1)
        circuit.cx(virt1, virt2)
        circuit.cx(occ1, occ2)

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
