"""
Circuit Optimization for VQE

Reduces gate count and circuit depth for faster execution.
"""

import numpy as np
from typing import List, Tuple, Optional
from kanad.ansatze.base_ansatz import QuantumCircuit
import logging

logger = logging.getLogger(__name__)


class CircuitOptimizer:
    """
    Optimize quantum circuits for faster VQE execution.

    Techniques:
    1. Gate cancellation (X-X, CNOT-CNOT inverses)
    2. Rotation merging (RZ(θ₁)+RZ(θ₂) → RZ(θ₁+θ₂))
    3. Small angle removal (RZ(ε) ≈ I for small ε)
    4. Commutation rules (move gates to reduce depth)
    """

    def __init__(self, threshold: float = 1e-6):
        """
        Initialize circuit optimizer.

        Args:
            threshold: Remove rotations smaller than this (radians)
        """
        self.threshold = threshold

    def optimize(self, circuit: QuantumCircuit) -> QuantumCircuit:
        """
        Apply all optimization passes.

        Args:
            circuit: Input circuit

        Returns:
            Optimized circuit
        """
        optimized = circuit

        # Pass 1: Cancel inverse gates
        optimized = self._cancel_inverse_gates(optimized)

        # Pass 2: Merge consecutive rotations
        optimized = self._merge_rotations(optimized)

        # Pass 3: Remove small rotations
        optimized = self._remove_small_rotations(optimized)

        return optimized

    def _cancel_inverse_gates(self, circuit: QuantumCircuit) -> QuantumCircuit:
        """
        Cancel adjacent inverse gate pairs.

        X-X → I
        CNOT-CNOT → I
        H-H → I
        """
        gates = circuit.gates.copy()
        optimized_gates = []

        i = 0
        while i < len(gates):
            if i < len(gates) - 1:
                current = gates[i]
                next_gate = gates[i+1]

                # Check if gates are inverses on same qubits
                if self._are_inverses(current, next_gate):
                    # Skip both gates (they cancel)
                    i += 2
                    continue

            optimized_gates.append(gates[i])
            i += 1

        circuit.gates = optimized_gates
        return circuit

    def _are_inverses(self, gate1: dict, gate2: dict) -> bool:
        """Check if two gates are inverses."""
        # Same gate type and qubits
        if gate1['type'] != gate2['type']:
            return False
        if gate1['qubits'] != gate2['qubits']:
            return False

        # Self-inverse gates (X, CNOT, H, etc.)
        self_inverse = ['X', 'Y', 'Z', 'H', 'CX', 'CNOT', 'SWAP']
        if gate1['type'] in self_inverse:
            return True

        # Rotation gates with opposite angles
        if gate1['type'] in ['RX', 'RY', 'RZ']:
            if 'parameter' in gate1 and 'parameter' in gate2:
                angle1 = gate1['parameter'].value if hasattr(gate1['parameter'], 'value') else gate1['parameter']
                angle2 = gate2['parameter'].value if hasattr(gate2['parameter'], 'value') else gate2['parameter']
                return np.isclose(angle1, -angle2)

        return False

    def _merge_rotations(self, circuit: QuantumCircuit) -> QuantumCircuit:
        """
        Merge consecutive rotations on same qubit.

        RZ(θ₁) + RZ(θ₂) → RZ(θ₁+θ₂)
        """
        gates = circuit.gates.copy()
        optimized_gates = []

        i = 0
        while i < len(gates):
            if i < len(gates) - 1:
                current = gates[i]
                next_gate = gates[i+1]

                # Check if can merge rotations
                merged = self._try_merge_rotations(current, next_gate)
                if merged is not None:
                    optimized_gates.append(merged)
                    i += 2
                    continue

            optimized_gates.append(gates[i])
            i += 1

        circuit.gates = optimized_gates
        return circuit

    def _try_merge_rotations(self, gate1: dict, gate2: dict) -> Optional[dict]:
        """Try to merge two rotation gates."""
        # Same type and qubit
        if gate1['type'] != gate2['type']:
            return None
        if gate1['qubits'] != gate2['qubits']:
            return None

        # Must be rotation gates
        if gate1['type'] not in ['RX', 'RY', 'RZ']:
            return None

        # Merge angles
        if 'parameter' in gate1 and 'parameter' in gate2:
            angle1 = gate1['parameter'].value if hasattr(gate1['parameter'], 'value') else gate1['parameter']
            angle2 = gate2['parameter'].value if hasattr(gate2['parameter'], 'value') else gate2['parameter']

            merged_angle = angle1 + angle2

            # Create merged gate
            merged = gate1.copy()
            if hasattr(gate1['parameter'], 'value'):
                merged['parameter'].value = merged_angle
            else:
                merged['parameter'] = merged_angle

            return merged

        return None

    def _remove_small_rotations(self, circuit: QuantumCircuit) -> QuantumCircuit:
        """
        Remove rotations with very small angles.

        RZ(ε) ≈ I for |ε| < threshold
        """
        gates = circuit.gates.copy()
        optimized_gates = []

        for gate in gates:
            # Check if rotation gate
            if gate['type'] in ['RX', 'RY', 'RZ']:
                if 'parameter' in gate:
                    angle = gate['parameter'].value if hasattr(gate['parameter'], 'value') else gate['parameter']

                    # Skip if angle is very small
                    if abs(angle) < self.threshold:
                        continue

            optimized_gates.append(gate)

        circuit.gates = optimized_gates
        return circuit

    def estimate_improvement(self, original: QuantumCircuit, optimized: QuantumCircuit) -> dict:
        """
        Estimate optimization improvement.

        Args:
            original: Original circuit
            optimized: Optimized circuit

        Returns:
            Dictionary with metrics
        """
        orig_gates = len(original.gates)
        opt_gates = len(optimized.gates)

        orig_depth = original.depth()
        opt_depth = optimized.depth()

        return {
            'original_gates': orig_gates,
            'optimized_gates': opt_gates,
            'gates_removed': orig_gates - opt_gates,
            'gate_reduction': (orig_gates - opt_gates) / orig_gates if orig_gates > 0 else 0,
            'original_depth': orig_depth,
            'optimized_depth': opt_depth,
            'depth_reduction': (orig_depth - opt_depth) / orig_depth if orig_depth > 0 else 0
        }
