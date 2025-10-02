"""
Base classes for variational quantum ansätze.

Ansätze define the parametrized quantum circuits used in VQE.
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Optional, Tuple
import numpy as np


class Parameter:
    """
    Variational parameter for quantum circuits.

    Lightweight parameter class for circuit construction.
    """

    def __init__(self, name: str, value: Optional[float] = None):
        """
        Initialize parameter.

        Args:
            name: Parameter name
            value: Initial value (optional)
        """
        self.name = name
        self.value = value if value is not None else 0.0

    def __repr__(self) -> str:
        return f"Parameter('{self.name}', value={self.value:.4f})"


class QuantumCircuit:
    """
    Lightweight quantum circuit representation.

    Stores gate sequence for ansatz construction.
    """

    def __init__(self, n_qubits: int):
        """
        Initialize circuit.

        Args:
            n_qubits: Number of qubits
        """
        self.n_qubits = n_qubits
        self.gates: List[Dict] = []
        self.parameters: List[Parameter] = []
        self.depth = 0

    def add_gate(self, gate_type: str, qubits: List[int], params: Optional[List] = None):
        """Add gate to circuit."""
        self.gates.append({
            'type': gate_type,
            'qubits': qubits,
            'params': params or []
        })
        self.depth += 1

        # Track parameters
        if params:
            for p in params:
                if isinstance(p, Parameter) and p not in self.parameters:
                    self.parameters.append(p)

    # Single-qubit gates
    def h(self, qubit: int):
        """Hadamard gate."""
        self.add_gate('h', [qubit])

    def x(self, qubit: int):
        """Pauli X gate."""
        self.add_gate('x', [qubit])

    def y(self, qubit: int):
        """Pauli Y gate."""
        self.add_gate('y', [qubit])

    def z(self, qubit: int):
        """Pauli Z gate."""
        self.add_gate('z', [qubit])

    def rx(self, theta, qubit: int):
        """Rotation around X axis."""
        self.add_gate('rx', [qubit], [theta])

    def ry(self, theta, qubit: int):
        """Rotation around Y axis."""
        self.add_gate('ry', [qubit], [theta])

    def rz(self, theta, qubit: int):
        """Rotation around Z axis."""
        self.add_gate('rz', [qubit], [theta])

    # Two-qubit gates
    def cx(self, control: int, target: int):
        """CNOT gate."""
        self.add_gate('cx', [control, target])

    def cnot(self, control: int, target: int):
        """CNOT gate (alias)."""
        self.cx(control, target)

    def cz(self, control: int, target: int):
        """Controlled-Z gate."""
        self.add_gate('cz', [control, target])

    def rxx(self, theta, qubit1: int, qubit2: int):
        """XX rotation gate."""
        self.add_gate('rxx', [qubit1, qubit2], [theta])

    def ryy(self, theta, qubit1: int, qubit2: int):
        """YY rotation gate."""
        self.add_gate('ryy', [qubit1, qubit2], [theta])

    def rzz(self, theta, qubit1: int, qubit2: int):
        """ZZ rotation gate."""
        self.add_gate('rzz', [qubit1, qubit2], [theta])

    def swap(self, qubit1: int, qubit2: int):
        """SWAP gate."""
        self.add_gate('swap', [qubit1, qubit2])

    def barrier(self):
        """Barrier (for visualization)."""
        self.add_gate('barrier', list(range(self.n_qubits)))

    def get_num_parameters(self) -> int:
        """Get number of variational parameters."""
        return len(self.parameters)

    def bind_parameters(self, values: np.ndarray):
        """
        Bind parameter values.

        Args:
            values: Array of parameter values
        """
        if len(values) != len(self.parameters):
            raise ValueError(f"Expected {len(self.parameters)} values, got {len(values)}")

        for param, value in zip(self.parameters, values):
            param.value = value

    def __repr__(self) -> str:
        return f"QuantumCircuit(n_qubits={self.n_qubits}, depth={self.depth}, parameters={len(self.parameters)})"


class BaseAnsatz(ABC):
    """
    Abstract base class for variational ansätze.

    Defines the interface for all ansatz implementations.
    """

    def __init__(self, n_qubits: int, n_electrons: int):
        """
        Initialize ansatz.

        Args:
            n_qubits: Number of qubits
            n_electrons: Number of electrons
        """
        self.n_qubits = n_qubits
        self.n_electrons = n_electrons
        self.circuit: Optional[QuantumCircuit] = None

    @abstractmethod
    def build_circuit(self, **kwargs) -> QuantumCircuit:
        """
        Build the ansatz circuit.

        Returns:
            Parametrized quantum circuit
        """
        pass

    def get_num_parameters(self) -> int:
        """Get number of variational parameters."""
        if self.circuit is None:
            self.circuit = self.build_circuit()
        return self.circuit.get_num_parameters()

    def get_circuit_depth(self) -> int:
        """Get circuit depth."""
        if self.circuit is None:
            self.circuit = self.build_circuit()
        return self.circuit.depth

    def initialize_parameters(self, strategy: str = 'random') -> np.ndarray:
        """
        Initialize parameter values.

        Args:
            strategy: Initialization strategy ('random', 'zeros', 'small_random')

        Returns:
            Initial parameter values
        """
        n_params = self.get_num_parameters()

        if strategy == 'random':
            return np.random.uniform(-np.pi, np.pi, n_params)
        elif strategy == 'zeros':
            return np.zeros(n_params)
        elif strategy == 'small_random':
            return np.random.uniform(-0.1, 0.1, n_params)
        else:
            raise ValueError(f"Unknown initialization strategy: {strategy}")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(n_qubits={self.n_qubits}, n_electrons={self.n_electrons})"
