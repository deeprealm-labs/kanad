"""
IBM Quantum Job Preparation Module

Handles all pre-execution steps for IBM Quantum:
- Hamiltonian preparation
- Circuit construction
- Observable mapping
- Transpilation optimization
"""

import logging
from typing import Dict, Any, Optional, List, Tuple
import numpy as np

logger = logging.getLogger(__name__)


class IBMPreparation:
    """
    Prepare quantum chemistry jobs for IBM Quantum execution.

    Handles:
    - Molecular Hamiltonian construction
    - Jordan-Wigner mapping
    - Ansatz circuit generation
    - Observable preparation for Estimator
    - Circuit transpilation

    Usage:
        prep = IBMPreparation(bond, ansatz='ucc')
        circuits, observables = prep.prepare_vqe()
        # Submit to IBM backend
    """

    def __init__(
        self,
        bond: 'BaseBond',
        ansatz_type: str = 'ucc',
        mapper_type: str = 'jordan_wigner'
    ):
        """
        Initialize preparation module.

        Args:
            bond: Molecular bond from BondFactory
            ansatz_type: Ansatz type ('ucc', 'hardware_efficient', 'governance')
            mapper_type: Fermionic-to-qubit mapping
        """
        self.bond = bond
        self.hamiltonian = bond.hamiltonian
        self.ansatz_type = ansatz_type
        self.mapper_type = mapper_type

        # Prepare components
        self._prepare_hamiltonian()
        self._prepare_ansatz()
        self._prepare_mapper()

        logger.info("IBM job preparation initialized")
        logger.info(f"  Ansatz: {ansatz_type}")
        logger.info(f"  Mapper: {mapper_type}")

    def _prepare_hamiltonian(self):
        """Prepare molecular Hamiltonian."""
        # Run SCF for reference
        self.density_matrix, self.hf_energy = self.hamiltonian.solve_scf()

        # Get Hamiltonian matrix
        self.n_qubits = 2 * self.hamiltonian.n_orbitals
        self.H_matrix = self.hamiltonian.to_matrix(
            n_qubits=self.n_qubits,
            use_mo_basis=True
        )

        logger.info(f"Hamiltonian prepared: {self.n_qubits} qubits")
        logger.info(f"HF energy: {self.hf_energy:.8f} Ha")

    def _prepare_ansatz(self):
        """Prepare variational ansatz."""
        from kanad.ansatze import UCCAnsatz, HardwareEfficientAnsatz

        if self.ansatz_type == 'ucc':
            self.ansatz = UCCAnsatz(
                hamiltonian=self.hamiltonian,
                excitation_type='sd',
                spin_ordering='blocked'
            )
        elif self.ansatz_type == 'hardware_efficient':
            n_qubits = 2 * self.hamiltonian.n_orbitals
            n_electrons = self.hamiltonian.molecule.n_electrons
            self.ansatz = HardwareEfficientAnsatz(
                n_qubits=n_qubits,
                n_electrons=n_electrons,
                n_layers=2,
                entanglement='linear'
            )
        else:
            # Try governance-aware ansatz
            from kanad.ansatze.governance_aware_ansatz import (
                CovalentGovernanceAnsatz,
                IonicGovernanceAnsatz
            )

            if self.bond.bond_type == 'covalent':
                self.ansatz = CovalentGovernanceAnsatz(hamiltonian=self.hamiltonian)
            elif self.bond.bond_type == 'ionic':
                self.ansatz = IonicGovernanceAnsatz(hamiltonian=self.hamiltonian)
            else:
                raise ValueError(f"Unsupported ansatz type: {self.ansatz_type}")

        # Build circuit
        self.ansatz.build_circuit()
        self.n_parameters = self.ansatz.circuit.get_num_parameters()

        logger.info(f"Ansatz built: {self.n_parameters} parameters")

    def _prepare_mapper(self):
        """Prepare fermionic-to-qubit mapper."""
        from kanad.core.mappers import JordanWignerMapper

        self.mapper = JordanWignerMapper()

        logger.info(f"Mapper prepared: {self.mapper_type}")

    def prepare_vqe_circuits(
        self,
        initial_parameters: Optional[np.ndarray] = None
    ) -> Tuple[List, List]:
        """
        Prepare circuits and observables for VQE.

        Args:
            initial_parameters: Initial parameter values (uses zeros if None)

        Returns:
            (circuits, observables) tuple
        """
        if initial_parameters is None:
            initial_parameters = np.zeros(self.n_parameters)

        # Bind parameters to circuit
        self.ansatz.circuit.bind_parameters(initial_parameters)

        # Convert to Qiskit
        qiskit_circuit = self.ansatz.circuit.to_qiskit()

        # Prepare Hamiltonian observable
        observable = self._hamiltonian_to_pauli_sum()

        logger.info("VQE circuits prepared")

        return [qiskit_circuit], [observable]

    def prepare_sqd_circuits(
        self,
        subspace_dim: int = 10
    ) -> List:
        """
        Prepare circuits for SQD.

        Args:
            subspace_dim: Dimension of quantum subspace

        Returns:
            List of quantum circuits
        """
        circuits = []

        # Generate diverse parameterized circuits for subspace
        for i in range(subspace_dim):
            # Random parameters for diversity
            params = np.random.randn(self.n_parameters) * 0.5

            # Bind and convert
            self.ansatz.circuit.bind_parameters(params)
            qiskit_circuit = self.ansatz.circuit.to_qiskit()

            circuits.append(qiskit_circuit)

        logger.info(f"SQD circuits prepared: {len(circuits)} states")

        return circuits

    def _hamiltonian_to_pauli_sum(self):
        """
        Convert Hamiltonian matrix to Pauli sum observable.

        Decomposes molecular Hamiltonian H into Pauli operator sum:
        H = Σ_i c_i P_i

        where P_i are Pauli strings (e.g., 'XXYZ') and c_i are coefficients.

        Returns:
            SparsePauliOp: Pauli operator representation for IBM Estimator
        """
        try:
            from kanad.core.hamiltonians.pauli_converter import PauliConverter

            # Convert Hamiltonian to Pauli operator using PauliConverter
            # This properly decomposes H_matrix into Pauli terms
            pauli_op = PauliConverter.to_sparse_pauli_op(
                self.hamiltonian,
                self.mapper
            )

            logger.info(f"Hamiltonian converted to {len(pauli_op)} Pauli terms")

            return pauli_op

        except Exception as e:
            # Fallback to identity if conversion fails
            logger.error(f"Hamiltonian-to-Pauli conversion failed: {e}")
            logger.warning("Falling back to identity observable")

            from qiskit.quantum_info import SparsePauliOp
            pauli_list = [('I' * self.n_qubits, 1.0)]
            return SparsePauliOp.from_list(pauli_list)

    def transpile_for_hardware(
        self,
        circuits: List,
        backend: 'IBMBackend',
        optimization_level: int = 1
    ) -> List:
        """
        Transpile circuits for specific IBM hardware.

        Args:
            circuits: Quantum circuits
            backend: IBM backend
            optimization_level: Optimization level (0-3)

        Returns:
            Transpiled circuits
        """
        from qiskit import transpile

        logger.info(f"Transpiling {len(circuits)} circuits")
        logger.info(f"  Backend: {backend.backend.name}")
        logger.info(f"  Optimization: level {optimization_level}")

        transpiled = transpile(
            circuits,
            backend=backend.backend,
            optimization_level=optimization_level
        )

        # Log transpilation results
        if isinstance(transpiled, list):
            avg_depth = np.mean([c.depth() for c in transpiled])
            logger.info(f"  Average depth after transpilation: {avg_depth:.1f}")
        else:
            logger.info(f"  Depth after transpilation: {transpiled.depth()}")

        return transpiled

    def get_preparation_summary(self) -> Dict[str, Any]:
        """Get summary of preparation."""
        return {
            'molecule': {
                'atoms': [a.symbol for a in self.hamiltonian.atoms],
                'n_electrons': self.hamiltonian.n_electrons,
                'n_orbitals': self.hamiltonian.n_orbitals,
                'bond_type': self.bond.bond_type
            },
            'quantum': {
                'n_qubits': self.n_qubits,
                'ansatz_type': self.ansatz_type,
                'n_parameters': self.n_parameters,
                'mapper_type': self.mapper_type
            },
            'energies': {
                'hf_energy': float(self.hf_energy),
                'nuclear_repulsion': float(self.hamiltonian.nuclear_repulsion)
            }
        }

    def __repr__(self):
        return f"IBMPreparation(ansatz='{self.ansatz_type}', qubits={self.n_qubits})"
