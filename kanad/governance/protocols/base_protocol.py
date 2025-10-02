"""
Base governance protocol for bonding-specific quantum circuit construction.

This is THE CORE INNOVATION of Kanad: bonding type determines quantum representation,
circuit topology, and allowed operations through governance rules.
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional, Callable
from enum import Enum
from dataclasses import dataclass
import numpy as np


class BondingType(Enum):
    """Types of chemical bonding."""
    IONIC = "ionic"
    COVALENT = "covalent"
    METALLIC = "metallic"
    MIXED = "mixed"


@dataclass
class GovernanceRule:
    """
    A single governance rule for quantum circuit construction.

    Rules dictate:
    - Which operators are physically valid
    - What circuit topology is allowed
    - How entanglement should be structured
    - Which symmetries must be preserved
    """

    name: str
    description: str
    condition: Callable[[Any, Dict], bool]  # (state, context) -> bool
    action: Callable[[Any, Dict], Any]  # (state, context) -> modified_state
    priority: int = 0  # Higher priority rules execute first
    required: bool = True  # Whether rule must be satisfied

    def applies(self, state: Any, context: Dict) -> bool:
        """
        Check if this rule applies to the current state.

        Args:
            state: Current quantum state or circuit
            context: Additional context information

        Returns:
            True if rule should be applied
        """
        try:
            return self.condition(state, context)
        except Exception as e:
            print(f"Warning: Rule '{self.name}' condition failed: {e}")
            return False

    def execute(self, state: Any, context: Dict) -> Any:
        """
        Execute the rule's action on the state.

        Args:
            state: Current state
            context: Context information

        Returns:
            Modified state
        """
        try:
            return self.action(state, context)
        except Exception as e:
            print(f"Error: Rule '{self.name}' execution failed: {e}")
            return state


class BaseGovernanceProtocol(ABC):
    """
    Abstract base class for governance protocols.

    Each bonding type has its own governance protocol that enforces
    physically appropriate quantum circuit construction.

    Philosophy:
        "The physics of bonding should dictate the quantum representation,
         not the other way around."
    """

    def __init__(self, bond_type: BondingType):
        """
        Initialize governance protocol.

        Args:
            bond_type: Type of bonding this protocol governs
        """
        self.bond_type = bond_type
        self.rules: List[GovernanceRule] = []
        self._initialize_rules()

    @abstractmethod
    def _initialize_rules(self):
        """
        Initialize bonding-specific governance rules.

        Each subclass must define its own rules based on the
        physical requirements of that bonding type.
        """
        pass

    @abstractmethod
    def validate_operator(self, operator: 'QuantumOperator') -> bool:
        """
        Validate if an operator is physically appropriate for this bonding type.

        Args:
            operator: Quantum operator to validate

        Returns:
            True if operator is valid for this bonding type
        """
        pass

    @abstractmethod
    def construct_ansatz(self, representation: Any) -> 'QuantumCircuit':
        """
        Construct physically appropriate variational ansatz circuit.

        The ansatz structure is determined by the bonding physics:
        - Ionic: localized gates, minimal entanglement
        - Covalent: paired gates for bonding orbitals
        - Metallic: collective gates, maximal entanglement

        Args:
            representation: Quantum representation object

        Returns:
            QuantumCircuit optimized for this bonding type
        """
        pass

    @abstractmethod
    def enforce_constraints(self, circuit: 'QuantumCircuit') -> 'QuantumCircuit':
        """
        Apply physical constraints to quantum circuit.

        Args:
            circuit: Input quantum circuit

        Returns:
            Modified circuit satisfying all constraints
        """
        pass

    def apply_governance(
        self,
        initial_state: Any,
        context: Dict[str, Any]
    ) -> Any:
        """
        Apply all governance rules in priority order.

        Args:
            initial_state: Initial quantum state or circuit
            context: Context dictionary with additional information

        Returns:
            Governed state satisfying all rules
        """
        state = initial_state

        # Sort rules by priority (highest first)
        sorted_rules = sorted(self.rules, key=lambda r: r.priority, reverse=True)

        # Apply each applicable rule
        for rule in sorted_rules:
            if rule.applies(state, context):
                state = rule.execute(state, context)
                # Note: state may be modified in-place, so identity check is not reliable

        return state

    def get_allowed_operators(self) -> List[str]:
        """
        Get list of operator types allowed for this bonding type.

        Returns:
            List of allowed operator names
        """
        return []  # Override in subclasses

    def get_forbidden_operators(self) -> List[str]:
        """
        Get list of operator types forbidden for this bonding type.

        Returns:
            List of forbidden operator names
        """
        return []  # Override in subclasses

    def check_symmetry_preservation(self, circuit: 'QuantumCircuit') -> bool:
        """
        Check if circuit preserves required symmetries.

        Args:
            circuit: Quantum circuit to check

        Returns:
            True if all symmetries are preserved
        """
        # Placeholder - override in subclasses
        return True

    def count_violations(self, circuit: 'QuantumCircuit') -> int:
        """
        Count number of governance rule violations in circuit.

        Args:
            circuit: Circuit to analyze

        Returns:
            Number of violations
        """
        violations = 0
        context = {'circuit': circuit}

        for rule in self.rules:
            if rule.required and not rule.applies(None, context):
                violations += 1

        return violations

    def get_entanglement_strategy(self) -> str:
        """
        Get the entanglement strategy for this bonding type.

        Returns:
            Description of entanglement strategy
        """
        strategies = {
            BondingType.IONIC: "minimal (nearest-neighbor only)",
            BondingType.COVALENT: "paired (bonding orbital pairs)",
            BondingType.METALLIC: "collective (GHZ-like states)"
        }
        return strategies.get(self.bond_type, "unknown")

    def __repr__(self) -> str:
        """String representation."""
        return (f"{self.__class__.__name__}("
                f"bond_type={self.bond_type.value}, "
                f"n_rules={len(self.rules)})")


class QuantumCircuitState:
    """
    Represents the state of a quantum circuit during construction.

    Used by governance rules to track circuit properties.
    """

    def __init__(self, n_qubits: int):
        """
        Initialize circuit state.

        Args:
            n_qubits: Number of qubits
        """
        self.n_qubits = n_qubits
        self.gates: List[Dict] = []
        self.entanglement_graph: Dict[int, List[int]] = {i: [] for i in range(n_qubits)}
        self.is_hybridized: bool = False
        self.has_mo_pairs: bool = False
        self.is_paired: bool = False
        self.in_k_space: bool = False
        self.is_localized: bool = True
        self.is_collectively_entangled: bool = False
        self.depth: int = 0
        self.metadata: Dict[str, Any] = {}

    def add_gate(self, gate_type: str, qubits: List[int], params: Optional[List[float]] = None):
        """Add a gate to the circuit."""
        self.gates.append({
            'type': gate_type,
            'qubits': qubits,
            'params': params or []
        })
        self.depth += 1

        # Update entanglement graph for 2-qubit gates
        if len(qubits) == 2:
            q1, q2 = qubits
            if q2 not in self.entanglement_graph[q1]:
                self.entanglement_graph[q1].append(q2)
            if q1 not in self.entanglement_graph[q2]:
                self.entanglement_graph[q2].append(q1)

    def get_entanglement_degree(self, qubit: int) -> int:
        """Get number of qubits entangled with given qubit."""
        return len(self.entanglement_graph[qubit])

    def max_entanglement_degree(self) -> int:
        """Get maximum entanglement degree in circuit."""
        return max(self.get_entanglement_degree(q) for q in range(self.n_qubits))

    def is_sparse(self) -> bool:
        """Check if circuit has sparse connectivity."""
        avg_degree = sum(self.get_entanglement_degree(q) for q in range(self.n_qubits)) / self.n_qubits
        return avg_degree < 2.0

    def __repr__(self) -> str:
        """String representation."""
        return (f"CircuitState(n_qubits={self.n_qubits}, "
                f"gates={len(self.gates)}, depth={self.depth})")
