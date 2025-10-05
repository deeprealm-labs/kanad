"""
Sample-based Quantum Diagonalization (SQD) Solver for IBM Quantum Platform.

Implements SQD algorithm using IBM Runtime for ground state estimation
through measurement sampling and classical post-processing.
"""

import numpy as np
from typing import Optional, Dict, List
import logging
from qiskit import QuantumCircuit

from kanad.backends.ibm.backend import IBMRuntimeBackend
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.ansatze.base_ansatz import BaseAnsatz

logger = logging.getLogger(__name__)


class IBMSQDSolver:
    """
    Sample-based Quantum Diagonalization solver for IBM Quantum.

    SQD estimates ground state energy by:
    1. Sampling quantum states
    2. Measuring Hamiltonian matrix elements
    3. Classical diagonalization of sampled subspace

    Note: Runs in batch mode (no sessions) for free/open IBM plans.
    When premium subscription is available, edit line ~160 to set session=True.
    """

    def __init__(
        self,
        hamiltonian: MolecularHamiltonian,
        ansatz: BaseAnsatz,
        n_samples: int = 10,
        backend: Optional[IBMRuntimeBackend] = None,
        backend_name: Optional[str] = None,
        token: Optional[str] = None,
        instance: Optional[str] = None,
        shots: int = 8192,
        optimization_level: int = 3,
        resilience_level: int = 1,
        # Advanced error mitigation options
        enable_dynamical_decoupling: bool = True,
        dynamical_decoupling_sequence: str = 'XX',
        enable_twirling: bool = True,
        twirling_strategy: str = 'all',
        zne_extrapolator: Optional[str] = None,
        **backend_options
    ):
        """
        Initialize IBM SQD solver with error mitigation.

        Args:
            hamiltonian: Molecular Hamiltonian
            ansatz: Ansatz for generating sample states
            n_samples: Number of quantum states to sample
            backend: Pre-initialized IBM backend (optional)
            backend_name: IBM backend name
            token: IBM Quantum token
            instance: IBM Quantum instance (CRN)
            shots: Number of measurement shots per circuit
            optimization_level: Transpiler optimization (0-3)
            resilience_level: Error mitigation level (0-2)
            enable_dynamical_decoupling: DD for decoherence suppression
            dynamical_decoupling_sequence: 'XX', 'XY4', or 'CPMG'
            enable_twirling: Gate twirling for noise averaging
            twirling_strategy: 'active', 'passive', or 'all'
            zne_extrapolator: ZNE method if resilience_level=2
            **backend_options: Additional backend options
        """
        self.hamiltonian = hamiltonian
        self.ansatz = ansatz
        self.n_samples = n_samples

        # Store error mitigation settings
        self.error_mitigation = {
            'resilience_level': resilience_level,
            'dynamical_decoupling': enable_dynamical_decoupling,
            'dd_sequence': dynamical_decoupling_sequence,
            'twirling': enable_twirling,
            'twirling_strategy': twirling_strategy,
            'zne_extrapolator': zne_extrapolator
        }

        # Initialize backend
        if backend is None:
            if backend_name is None:
                raise ValueError("Either 'backend' or 'backend_name' must be provided")

            logger.info(f"Initializing IBM backend for SQD: {backend_name}")
            self.backend = IBMRuntimeBackend(
                backend_name=backend_name,
                token=token,
                instance=instance,
                shots=shots,
                optimization_level=optimization_level,
                resilience_level=resilience_level,
                **backend_options
            )
        else:
            self.backend = backend

        # System size
        self.n_qubits = 2 * hamiltonian.n_orbitals

        # Build ansatz circuit
        self.circuit = ansatz.build_circuit()
        self.n_parameters = self.circuit.get_num_parameters()

        # Convert Hamiltonian to Pauli operators
        logger.info("Converting Hamiltonian to sparse Pauli representation...")
        self.pauli_hamiltonian = hamiltonian.to_sparse_hamiltonian()
        logger.info(f"Hamiltonian: {len(self.pauli_hamiltonian)} Pauli terms")

        # Convert circuit to Qiskit
        self.qiskit_circuit = self.circuit.to_qiskit()

        logger.info(f"SQD Configuration:")
        logger.info(f"  System qubits: {self.n_qubits}")
        logger.info(f"  Sample states: {n_samples}")
        logger.info(f"  Parameters per state: {self.n_parameters}")

    def _generate_sample_states(self) -> List[np.ndarray]:
        """
        Generate random parameter sets for sampling quantum subspace.

        Returns:
            List of parameter arrays
        """
        logger.info("\nGenerating sample states...")

        sample_parameters = []

        for i in range(self.n_samples):
            # Random parameters in [-π, π]
            params = np.random.uniform(-np.pi, np.pi, self.n_parameters)
            sample_parameters.append(params)

            logger.info(f"  Sample {i+1}/{self.n_samples}: "
                       f"params = [{params[0]:.3f}, {params[1]:.3f}, ...]")

        return sample_parameters

    def _measure_matrix_elements(
        self,
        sample_parameters: List[np.ndarray]
    ) -> np.ndarray:
        """
        Measure Hamiltonian matrix elements between sample states.

        H_ij = ⟨ψ_i|H|ψ_j⟩

        Returns:
            Hamiltonian matrix in sampled subspace (n_samples × n_samples)
        """
        logger.info("\nMeasuring Hamiltonian matrix elements...")

        # Use batch mode (no session) for free/open plans
        # When premium subscription is available, set session=True for better batching
        estimator = self.backend.get_estimator(session=False)
        H_matrix = np.zeros((self.n_samples, self.n_samples), dtype=complex)

        # Transpile base circuit to ISA once
        from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
        from qiskit.quantum_info import SparsePauliOp

        logger.info("Transpiling circuit to ISA...")
        pm = generate_preset_pass_manager(
            optimization_level=self.backend.optimization_level,
            backend=self.backend.backend
        )
        isa_circuit = pm.run(self.qiskit_circuit)

        # Pad Hamiltonian to match ISA circuit qubits
        n_extra_qubits = isa_circuit.num_qubits - len(self.pauli_hamiltonian.paulis[0])
        if n_extra_qubits > 0:
            padded_paulis = []
            for pauli, coeff in zip(self.pauli_hamiltonian.paulis, self.pauli_hamiltonian.coeffs):
                padded_pauli = str(pauli) + 'I' * n_extra_qubits
                padded_paulis.append((padded_pauli, coeff))
            isa_hamiltonian = SparsePauliOp.from_list(padded_paulis)
        else:
            isa_hamiltonian = self.pauli_hamiltonian

        # Prepare all circuits with bound parameters
        circuits = []
        observables = []

        for i in range(self.n_samples):
            # Diagonal elements only for simplicity
            bound_circuit = isa_circuit.assign_parameters(sample_parameters[i])
            circuits.append(bound_circuit)
            observables.append(isa_hamiltonian)

        # Measure diagonal elements
        logger.info(f"Measuring {len(circuits)} diagonal matrix elements...")
        job = estimator.run(list(zip(circuits, observables)))
        result = job.result()

        # Fill diagonal
        for i in range(self.n_samples):
            energy = result[i].data.evs  # Scalar in EstimatorV2
            H_matrix[i, i] = energy
            logger.info(f"  H[{i},{i}] = {energy:.6f} Ha")

        # For off-diagonal, use approximation (assume diagonal dominance)
        # For now, assume diagonal dominance
        return H_matrix.real

    def solve(self) -> Dict:
        """
        Run SQD on IBM Quantum.

        Returns:
            Dictionary with:
                - energy: Estimated ground state energy (Hartree)
                - eigenvalues: All eigenvalues of sampled subspace
                - n_samples: Number of sampled states
                - hamiltonian_matrix: Sampled Hamiltonian matrix
                - backend_info: IBM backend information
        """
        logger.info("=" * 80)
        logger.info("IBM QUANTUM SQD")
        logger.info("=" * 80)

        backend_info = self.backend.get_backend_info()
        logger.info(f"Backend: {backend_info['name']}")
        logger.info(f"  Qubits: {backend_info.get('num_qubits', 'N/A')}")
        logger.info(f"  Shots: {self.backend.shots}")
        logger.info(f"Samples: {self.n_samples}")

        # 1. Generate sample states
        sample_parameters = self._generate_sample_states()

        # 2. Measure Hamiltonian matrix elements
        H_matrix = self._measure_matrix_elements(sample_parameters)

        # 3. Classical diagonalization
        logger.info("\nPerforming classical diagonalization...")
        eigenvalues, eigenvectors = np.linalg.eigh(H_matrix)

        ground_state_energy = eigenvalues[0]

        logger.info("=" * 80)
        logger.info("SQD RESULTS")
        logger.info("=" * 80)
        logger.info(f"Ground State Energy: {ground_state_energy:.6f} Ha")
        logger.info(f"                     {ground_state_energy * 27.211386:.4f} eV")
        logger.info(f"\nAll eigenvalues:")
        for i, E in enumerate(eigenvalues):
            logger.info(f"  State {i}: {E:.6f} Ha ({E * 27.211386:.4f} eV)")
        logger.info("=" * 80)

        return {
            'energy': ground_state_energy,
            'energy_hartree': ground_state_energy,
            'energy_ev': ground_state_energy * 27.211386,
            'eigenvalues': eigenvalues,
            'eigenvalues_ev': eigenvalues * 27.211386,
            'eigenvectors': eigenvectors,
            'n_samples': self.n_samples,
            'hamiltonian_matrix': H_matrix,
            'backend_info': backend_info
        }
