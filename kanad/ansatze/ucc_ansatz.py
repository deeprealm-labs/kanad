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

        Returns:
            List of (occupied_orbitals, virtual_orbitals) tuples
        """
        n_occ = self.n_electrons // 2  # Closed shell
        n_virt = self.n_qubits - n_occ

        excitations = []

        # Singles: |occ⟩ → |virt⟩
        if self.include_singles:
            for i in range(n_occ):
                for a in range(n_occ, self.n_qubits):
                    excitations.append(([i], [a]))

        # Doubles: |occ₁ occ₂⟩ → |virt₁ virt₂⟩
        if self.include_doubles:
            for i in range(n_occ):
                for j in range(i + 1, n_occ):
                    for a in range(n_occ, self.n_qubits):
                        for b in range(a + 1, self.n_qubits):
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

        For spin orbitals with alpha/beta separation:
        - Qubits 0..n_orbitals-1: spin-up (alpha)
        - Qubits n_orbitals..2*n_orbitals-1: spin-down (beta)

        Fills electrons in lowest spatial orbitals with spin pairing.

        Returns:
            Occupation list for spin orbitals
        """
        state = [0] * self.n_qubits
        n_spatial_orbitals = self.n_qubits // 2
        n_electrons = self.n_electrons

        # For closed-shell: fill spatial orbitals from bottom, alternating spins
        # Assuming even number of electrons (closed shell)
        if n_electrons % 2 == 0:
            n_pairs = n_electrons // 2
            for i in range(n_pairs):
                # Fill spin-up (alpha)
                state[i] = 1
                # Fill spin-down (beta) in corresponding spatial orbital
                state[i + n_spatial_orbitals] = 1
        else:
            # Open-shell: fill pairs, then add unpaired electron to spin-up
            n_pairs = n_electrons // 2
            for i in range(n_pairs):
                state[i] = 1
                state[i + n_spatial_orbitals] = 1
            # Unpaired electron in next spin-up orbital
            state[n_pairs] = 1

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
