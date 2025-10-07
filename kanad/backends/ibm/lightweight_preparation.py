"""
Lightweight IBM Quantum Job Preparation
========================================

Memory-efficient preparation that avoids creating large Hamiltonian matrices.
Uses sparse representations and lazy evaluation.
"""

import logging
from typing import Dict, Any, Optional, List, Tuple
import numpy as np

logger = logging.getLogger(__name__)


class LightweightIBMPreparation:
    """
    Memory-efficient preparation for IBM Quantum jobs.

    Avoids:
    - Full Hamiltonian matrix construction (2^n x 2^n)
    - Large ansatz circuits (UCC can be huge)
    - Unnecessary intermediate computations

    Uses:
    - Sparse Hamiltonian representation
    - Simple parametric ansatz
    - Direct observable construction
    """

    def __init__(
        self,
        bond: 'BaseBond',
        ansatz_type: str = 'hardware_efficient',
        n_layers: int = 2
    ):
        """
        Initialize lightweight preparation.

        Args:
            bond: Molecular bond
            ansatz_type: 'hardware_efficient' only (memory safe)
            n_layers: Number of ansatz layers
        """
        self.bond = bond
        self.hamiltonian = bond.hamiltonian
        self.ansatz_type = ansatz_type
        self.n_layers = n_layers

        # Basic info
        self.n_qubits = 2 * self.hamiltonian.n_orbitals
        self.n_electrons = self.hamiltonian.molecule.n_electrons

        # Get HF reference WITHOUT building full matrix
        logger.info(f"Computing HF reference energy...")
        self.density_matrix, self.hf_energy = self.hamiltonian.solve_scf()
        logger.info(f"  HF energy: {self.hf_energy:.8f} Ha")

        # Prepare lightweight ansatz
        logger.info(f"Building {ansatz_type} ansatz...")
        self._prepare_ansatz()
        logger.info(f"  Circuit depth: {self.circuit.depth()}")
        logger.info(f"  Parameters: {self.n_parameters}")

    def _prepare_ansatz(self):
        """Build simple hardware-efficient ansatz circuit."""
        from qiskit import QuantumCircuit
        from qiskit.circuit import Parameter

        qc = QuantumCircuit(self.n_qubits)

        params = []

        # Layer structure: RY rotations + CX entanglement
        for layer in range(self.n_layers):
            # Rotation layer
            for q in range(self.n_qubits):
                param = Parameter(f'θ_{layer}_{q}')
                params.append(param)
                qc.ry(param, q)

            # Entanglement layer (linear)
            for q in range(self.n_qubits - 1):
                qc.cx(q, q + 1)

        # Final rotation layer
        for q in range(self.n_qubits):
            param = Parameter(f'θ_final_{q}')
            params.append(param)
            qc.ry(param, q)

        self.circuit = qc
        self.parameters = params
        self.n_parameters = len(params)

    def _build_hamiltonian_observable(self):
        """
        Build Pauli observable from Hamiltonian.

        Uses sparse representation - doesn't build full matrix.
        """
        from qiskit.quantum_info import SparsePauliOp

        logger.info("Building Hamiltonian observable...")

        try:
            # Try to get Hamiltonian in Pauli basis
            if hasattr(self.hamiltonian, 'to_pauli_op'):
                pauli_op = self.hamiltonian.to_pauli_op()
                logger.info(f"  Using existing Pauli representation")
                return pauli_op

            # Otherwise, create simple Z-basis observable for testing
            # This is a placeholder - for real VQE, need proper Hamiltonian
            logger.warning("Full Hamiltonian not available, using Z-basis approximation")

            # Simple Z-observable (measures energy in computational basis)
            pauli_strings = []
            coeffs = []

            # Add Z terms for each qubit
            for q in range(self.n_qubits):
                pauli_str = 'I' * q + 'Z' + 'I' * (self.n_qubits - q - 1)
                pauli_strings.append(pauli_str)
                coeffs.append(-0.1)  # Small coefficient

            # Add identity (constant energy shift)
            pauli_strings.append('I' * self.n_qubits)
            coeffs.append(float(self.hf_energy))  # HF energy as constant

            observable = SparsePauliOp(pauli_strings, coeffs)
            logger.info(f"  Created Z-basis observable with {len(pauli_strings)} terms")

            return observable

        except Exception as e:
            logger.error(f"Failed to build observable: {e}")
            # Fallback: identity operator
            observable = SparsePauliOp(['I' * self.n_qubits], [float(self.hf_energy)])
            logger.warning("Using identity observable (HF energy)")
            return observable

    def prepare_vqe_circuits(
        self,
        initial_parameters: Optional[np.ndarray] = None
    ) -> Tuple[List, List]:
        """
        Prepare circuit and observable for VQE.

        Args:
            initial_parameters: Initial parameter values (uses small random if None)

        Returns:
            ([circuit], [observable])
        """
        # Use small random parameters if none provided
        if initial_parameters is None:
            initial_parameters = np.random.uniform(-0.1, 0.1, self.n_parameters)

        # Bind parameters
        bound_circuit = self.circuit.assign_parameters(initial_parameters)

        # Build observable
        observable = self._build_hamiltonian_observable()

        logger.info("VQE circuits prepared")
        logger.info(f"  Circuit: {self.n_qubits} qubits, depth {bound_circuit.depth()}")
        logger.info(f"  Observable: {observable.num_qubits} qubits")

        return [bound_circuit], [observable]

    def get_preparation_summary(self) -> Dict[str, Any]:
        """Get summary of preparation."""
        return {
            'ansatz_type': self.ansatz_type,
            'n_qubits': self.n_qubits,
            'n_electrons': self.n_electrons,
            'n_parameters': self.n_parameters,
            'n_layers': self.n_layers,
            'hf_energy': float(self.hf_energy),
            'circuit_depth': self.circuit.depth()
        }
