"""
TwoLocal ansatz - hardware-efficient variational form.

Simpler than UCC, uses only rotation and entanglement gates.
Well-suited for NISQ devices and often achieves good results.
"""

from typing import List, Optional
import numpy as np
from kanad.ansatze.base_ansatz import BaseAnsatz, QuantumCircuit, Parameter


class TwoLocalAnsatz(BaseAnsatz):
    """
    TwoLocal ansatz with rotation and entanglement layers.

    Structure:
    1. Initial rotation layer (RY gates on all qubits)
    2. Entanglement layer (CX gates in pattern)
    3. Rotation layer
    4. Entanglement layer
    5. ... repeat for n_layers
    6. Final rotation layer

    This is the ansatz used by Qiskit Nature and often achieves
    near-FCI accuracy with fewer parameters than UCC.
    """

    def __init__(
        self,
        n_qubits: int,
        n_electrons: int,
        n_layers: int = 2,
        rotation_gates: str = 'ry',
        entanglement: str = 'linear',
        entanglement_gate: str = 'cx'
    ):
        """
        Initialize TwoLocal ansatz.

        Args:
            n_qubits: Number of qubits
            n_electrons: Number of electrons
            n_layers: Number of repetition layers (depth)
            rotation_gates: Single-qubit rotation ('ry', 'rx', 'rz', or 'full')
            entanglement: Entanglement pattern ('linear', 'full', 'circular')
            entanglement_gate: Two-qubit gate ('cx', 'cz', 'rxx', 'ryy', 'rzz')
        """
        super().__init__(n_qubits, n_electrons)
        self.n_layers = n_layers
        self.rotation_gates = rotation_gates
        self.entanglement = entanglement
        self.entanglement_gate = entanglement_gate

    def build_circuit(self, initial_state: Optional[List[int]] = None) -> QuantumCircuit:
        """
        Build TwoLocal circuit.

        Args:
            initial_state: Initial state (default: HF state)

        Returns:
            TwoLocal circuit
        """
        circuit = QuantumCircuit(self.n_qubits)

        # 1. Prepare initial state (HF reference)
        if initial_state is None:
            initial_state = self._hartree_fock_state()

        for qubit, occupation in enumerate(initial_state):
            if occupation == 1:
                circuit.x(qubit)

        circuit.barrier()

        # 2. Build layers
        param_idx = 0

        for layer in range(self.n_layers):
            # Rotation layer
            param_idx = self._add_rotation_layer(circuit, param_idx)

            # Entanglement layer
            self._add_entanglement_layer(circuit)

        # Final rotation layer
        param_idx = self._add_rotation_layer(circuit, param_idx)

        self.circuit = circuit
        return circuit

    def _add_rotation_layer(self, circuit: QuantumCircuit, param_idx: int) -> int:
        """
        Add rotation layer to circuit.

        Args:
            circuit: Circuit to modify
            param_idx: Current parameter index

        Returns:
            Updated parameter index
        """
        if self.rotation_gates == 'ry':
            # Single RY rotation per qubit
            for qubit in range(self.n_qubits):
                theta = Parameter(f'θ_{param_idx}')
                circuit.ry(theta, qubit)
                param_idx += 1

        elif self.rotation_gates == 'rx':
            for qubit in range(self.n_qubits):
                theta = Parameter(f'θ_{param_idx}')
                circuit.rx(theta, qubit)
                param_idx += 1

        elif self.rotation_gates == 'rz':
            for qubit in range(self.n_qubits):
                theta = Parameter(f'θ_{param_idx}')
                circuit.rz(theta, qubit)
                param_idx += 1

        elif self.rotation_gates == 'full':
            # Full rotation: RZ-RY-RZ per qubit
            for qubit in range(self.n_qubits):
                for gate_type in ['rz', 'ry', 'rz']:
                    theta = Parameter(f'θ_{param_idx}')
                    if gate_type == 'rz':
                        circuit.rz(theta, qubit)
                    elif gate_type == 'ry':
                        circuit.ry(theta, qubit)
                    param_idx += 1
        else:
            raise ValueError(f"Unknown rotation_gates: {self.rotation_gates}")

        return param_idx

    def _add_entanglement_layer(self, circuit: QuantumCircuit):
        """
        Add entanglement layer to circuit.

        Args:
            circuit: Circuit to modify
        """
        if self.entanglement == 'linear':
            # Linear chain: 0-1, 1-2, 2-3, ...
            for i in range(self.n_qubits - 1):
                self._add_entangling_gate(circuit, i, i + 1)

        elif self.entanglement == 'full':
            # All-to-all: every pair
            for i in range(self.n_qubits):
                for j in range(i + 1, self.n_qubits):
                    self._add_entangling_gate(circuit, i, j)

        elif self.entanglement == 'circular':
            # Circular: 0-1, 1-2, ..., (n-1)-0
            for i in range(self.n_qubits):
                j = (i + 1) % self.n_qubits
                self._add_entangling_gate(circuit, i, j)
        else:
            raise ValueError(f"Unknown entanglement: {self.entanglement}")

    def _add_entangling_gate(self, circuit: QuantumCircuit, qubit1: int, qubit2: int):
        """
        Add entangling gate between two qubits.

        Args:
            circuit: Circuit to modify
            qubit1: First qubit
            qubit2: Second qubit
        """
        if self.entanglement_gate == 'cx':
            circuit.cx(qubit1, qubit2)
        elif self.entanglement_gate == 'cz':
            circuit.cz(qubit1, qubit2)
        elif self.entanglement_gate == 'rxx':
            # Fixed angle for simplicity
            circuit.rxx(np.pi/4, qubit1, qubit2)
        elif self.entanglement_gate == 'ryy':
            circuit.ryy(np.pi/4, qubit1, qubit2)
        elif self.entanglement_gate == 'rzz':
            circuit.rzz(np.pi/4, qubit1, qubit2)
        else:
            raise ValueError(f"Unknown entanglement_gate: {self.entanglement_gate}")

    def _hartree_fock_state(self) -> List[int]:
        """
        Generate HF reference state.

        Qubit ordering: interleaved (q0=orb0↑, q1=orb1↑, q2=orb0↓, q3=orb1↓, ...)
        For H2 (2 electrons in orbital 0): |0101⟩ = [1,0,1,0]

        Returns:
            Occupation list
        """
        state = [0] * self.n_qubits
        n_orbitals = self.n_qubits // 2
        n_up = (self.n_electrons + 1) // 2  # Spin-up electrons
        n_down = self.n_electrons // 2      # Spin-down electrons

        # Fill spin-up orbitals (qubits 0, 1, 2, ...)
        for i in range(min(n_up, n_orbitals)):
            state[i] = 1

        # Fill spin-down orbitals (qubits n_orbitals, n_orbitals+1, ...)
        for i in range(min(n_down, n_orbitals)):
            state[n_orbitals + i] = 1

        return state

    def __repr__(self) -> str:
        return (f"TwoLocalAnsatz(n_qubits={self.n_qubits}, n_layers={self.n_layers}, "
                f"rotation={self.rotation_gates}, entanglement={self.entanglement})")
