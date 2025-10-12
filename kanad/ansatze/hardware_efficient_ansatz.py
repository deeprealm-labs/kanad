"""
Hardware-Efficient ansatz for near-term quantum devices.

Optimized for shallow circuits with native gate sets.
"""

from typing import List, Optional
import numpy as np
from kanad.ansatze.base_ansatz import BaseAnsatz, QuantumCircuit, Parameter


def get_hf_state_qubits(n_qubits: int, n_electrons: int, mapper: str = 'jordan_wigner') -> List[int]:
    """
    Get which qubits to set for HF state preparation in different mappers.

    Args:
        n_qubits: Total number of qubits
        n_electrons: Number of electrons
        mapper: Mapper type ('jordan_wigner', 'bravyi_kitaev', 'parity')

    Returns:
        List of qubit indices to apply X gates to

    Examples:
        H2 (2e, 4 qubits):
        - JW: |1100⟩ → qubits [2, 3]
        - BK: |1000⟩ → qubits [3]
    """
    mapper = mapper.lower()

    if mapper == 'bravyi_kitaev':
        # Bravyi-Kitaev uses binary tree encoding
        # For H2 (2e, 4 qubits): HF state is |1000⟩ (only highest qubit)
        # General pattern: Set the highest ceil(log2(n_electrons + 1)) qubits
        #
        # Actually for BK, the HF state depends on the specific transformation
        # For H2: |1000⟩ means qubit 3 only
        # For now, use empirical rule: highest qubit only for 2 electrons
        if n_electrons == 2 and n_qubits == 4:
            return [3]  # |1000⟩ in little-endian
        elif n_electrons == 2:
            return [n_qubits - 1]  # Highest qubit only
        else:
            # For other electron counts, use JW as fallback
            # TODO: Implement general BK HF state finder
            start_qubit = n_qubits - n_electrons
            return list(range(start_qubit, n_qubits))

    elif mapper in ['jordan_wigner', 'parity']:
        # Jordan-Wigner and Parity use blocked spin ordering: [0↑, 1↑, ..., 0↓, 1↓, ...]
        # Fill first n_electrons spin-orbitals
        # In Qiskit little-endian: |1100⟩ for 2e = qubits [2, 3]
        start_qubit = n_qubits - n_electrons
        return list(range(start_qubit, n_qubits))

    else:
        # Default to JW
        start_qubit = n_qubits - n_electrons
        return list(range(start_qubit, n_qubits))


class HardwareEfficientAnsatz(BaseAnsatz):
    """
    Hardware-Efficient ansatz.

    Uses repeating layers of:
    1. Single-qubit rotations (RY gates)
    2. Entangling gates (CNOT, CZ, or custom)

    ADVANTAGES:
    - Shallow circuits
    - Native gate set
    - Fast execution on hardware
    - Trainable on NISQ devices

    DISADVANTAGES:
    - Not chemically motivated
    - May require many layers
    - Barren plateau issues
    """

    def __init__(
        self,
        n_qubits: int,
        n_electrons: int,
        n_layers: int = 2,
        entanglement: str = 'linear',
        rotation_gates: List[str] = None,
        entangling_gate: str = 'cx',
        mapper: str = 'jordan_wigner'
    ):
        """
        Initialize hardware-efficient ansatz.

        Args:
            n_qubits: Number of qubits
            n_electrons: Number of electrons
            n_layers: Number of ansatz layers
            entanglement: Entanglement pattern ('linear', 'circular', 'full')
            rotation_gates: List of rotation gates ['ry', 'rz', 'rx']
            entangling_gate: Entangling gate type ('cx', 'cz', 'rxx')
            mapper: Fermion-to-qubit mapper type ('jordan_wigner', 'bravyi_kitaev', 'parity')
        """
        super().__init__(n_qubits, n_electrons)
        self.n_layers = n_layers
        self.entanglement = entanglement
        self.rotation_gates = rotation_gates or ['ry']
        self.entangling_gate = entangling_gate
        self.mapper = mapper.lower()

    @property
    def n_parameters(self) -> int:
        """Number of variational parameters."""
        # Each layer has: n_qubits * len(rotation_gates) rotations
        # (Note: no second rotation layer after entangling gates)
        return self.n_layers * self.n_qubits * len(self.rotation_gates)

    def build_circuit(self, initial_state: Optional[List[int]] = None) -> QuantumCircuit:
        """
        Build hardware-efficient circuit.

        Structure:
        [Initial state] - [Layer 1] - [Layer 2] - ... - [Layer n]

        Each layer:
        [Rotations on all qubits] - [Entangling gates] - [Rotations]

        Args:
            initial_state: Initial state preparation

        Returns:
            Hardware-efficient circuit
        """
        circuit = QuantumCircuit(self.n_qubits)

        # 1. Prepare initial state (optional)
        if initial_state is not None:
            for qubit, occupation in enumerate(initial_state):
                if occupation == 1:
                    circuit.x(qubit)
            circuit.barrier()
        elif self.n_electrons > 0:
            # Default: Hartree-Fock state with MAPPER-AWARE preparation
            # Different fermion-to-qubit mappings use different HF state encodings:
            #   - Jordan-Wigner: |1100⟩ for H2 (qubits [2, 3])
            #   - Bravyi-Kitaev: |1000⟩ for H2 (qubit [3] only)
            #   - Parity: Similar to JW

            # CRITICAL: Use mapper-aware HF state preparation
            hf_qubits = get_hf_state_qubits(self.n_qubits, self.n_electrons, self.mapper)
            for qubit in hf_qubits:
                circuit.x(qubit)

            circuit.barrier()

        # 2. Apply layers
        for layer_idx in range(self.n_layers):
            self._apply_layer(circuit, layer_idx)

        self.circuit = circuit
        return circuit

    def _apply_layer(self, circuit: QuantumCircuit, layer_idx: int):
        """
        Apply one layer of the ansatz.

        Args:
            circuit: Circuit to modify
            layer_idx: Layer index
        """
        # Single-qubit rotations
        for gate_type in self.rotation_gates:
            for qubit in range(self.n_qubits):
                param = Parameter(f'θ_{layer_idx}_{gate_type}_{qubit}')

                if gate_type == 'ry':
                    circuit.ry(param, qubit)
                elif gate_type == 'rz':
                    circuit.rz(param, qubit)
                elif gate_type == 'rx':
                    circuit.rx(param, qubit)

        # Entangling gates
        self._apply_entanglement(circuit, layer_idx)

    def _apply_entanglement(self, circuit: QuantumCircuit, layer_idx: int):
        """
        Apply entangling layer.

        Args:
            circuit: Circuit to modify
            layer_idx: Layer index
        """
        if self.entanglement == 'linear':
            # Linear chain: 0-1, 1-2, 2-3, ...
            for i in range(self.n_qubits - 1):
                self._apply_entangling_gate(circuit, i, i + 1, layer_idx)

        elif self.entanglement == 'circular':
            # Circular: linear + connect last to first
            for i in range(self.n_qubits - 1):
                self._apply_entangling_gate(circuit, i, i + 1, layer_idx)
            # Close the loop
            self._apply_entangling_gate(circuit, self.n_qubits - 1, 0, layer_idx)

        elif self.entanglement == 'full':
            # All-to-all entanglement
            for i in range(self.n_qubits):
                for j in range(i + 1, self.n_qubits):
                    self._apply_entangling_gate(circuit, i, j, layer_idx)

        elif self.entanglement == 'pairwise':
            # Pairs: (0,1), (2,3), (4,5), ...
            for i in range(0, self.n_qubits - 1, 2):
                self._apply_entangling_gate(circuit, i, i + 1, layer_idx)

    def _apply_entangling_gate(
        self,
        circuit: QuantumCircuit,
        qubit1: int,
        qubit2: int,
        layer_idx: int
    ):
        """Apply entangling gate between two qubits."""
        if self.entangling_gate == 'cx' or self.entangling_gate == 'cnot':
            circuit.cx(qubit1, qubit2)

        elif self.entangling_gate == 'cz':
            circuit.cz(qubit1, qubit2)

        elif self.entangling_gate == 'rxx':
            param = Parameter(f'θ_{layer_idx}_rxx_{qubit1}_{qubit2}')
            circuit.rxx(param, qubit1, qubit2)

        elif self.entangling_gate == 'ryy':
            param = Parameter(f'θ_{layer_idx}_ryy_{qubit1}_{qubit2}')
            circuit.ryy(param, qubit1, qubit2)

        elif self.entangling_gate == 'rzz':
            param = Parameter(f'θ_{layer_idx}_rzz_{qubit1}_{qubit2}')
            circuit.rzz(param, qubit1, qubit2)

    def __repr__(self) -> str:
        return (f"HardwareEfficientAnsatz(n_qubits={self.n_qubits}, "
                f"layers={self.n_layers}, entanglement='{self.entanglement}')")


class RealAmplitudesAnsatz(HardwareEfficientAnsatz):
    """
    Real-amplitudes ansatz (RY rotations only).

    Uses only RY gates + entanglers → real amplitudes.
    Commonly used in VQE applications.
    """

    def __init__(
        self,
        n_qubits: int,
        n_electrons: int,
        n_layers: int = 2,
        entanglement: str = 'linear'
    ):
        super().__init__(
            n_qubits,
            n_electrons,
            n_layers=n_layers,
            entanglement=entanglement,
            rotation_gates=['ry'],
            entangling_gate='cx'
        )


class EfficientSU2Ansatz(HardwareEfficientAnsatz):
    """
    Efficient SU(2) ansatz.

    Uses RY and RZ rotations to span SU(2).
    More expressive than real amplitudes.
    """

    def __init__(
        self,
        n_qubits: int,
        n_electrons: int,
        n_layers: int = 2,
        entanglement: str = 'linear'
    ):
        super().__init__(
            n_qubits,
            n_electrons,
            n_layers=n_layers,
            entanglement=entanglement,
            rotation_gates=['ry', 'rz'],
            entangling_gate='cx'
        )
