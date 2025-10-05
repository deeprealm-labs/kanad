"""
Quantum Phase Estimation Solver for IBM Quantum Platform.

Implements QPE algorithm using IBM Runtime for eigenvalue estimation
with high precision and error mitigation.
"""

import numpy as np
from typing import Optional, Dict, List
import logging
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.circuit.library import QFT

from kanad.backends.ibm.backend import IBMRuntimeBackend
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian

logger = logging.getLogger(__name__)


class IBMQPESolver:
    """
    Quantum Phase Estimation solver for IBM Quantum.

    QPE provides exponential advantage for eigenvalue estimation
    compared to classical methods. Suitable for high-precision
    ground state energy calculations.
    """

    def __init__(
        self,
        hamiltonian: MolecularHamiltonian,
        n_ancilla: int = 6,
        backend: Optional[IBMRuntimeBackend] = None,
        backend_name: Optional[str] = None,
        token: Optional[str] = None,
        instance: Optional[str] = None,
        shots: int = 8192,
        optimization_level: int = 3,
        resilience_level: int = 2,
        **backend_options
    ):
        """
        Initialize IBM QPE solver.

        Args:
            hamiltonian: Molecular Hamiltonian
            n_ancilla: Number of ancilla qubits for phase precision
            backend: Pre-initialized IBM backend (optional)
            backend_name: IBM backend name
            token: IBM Quantum token
            instance: IBM Quantum instance (CRN)
            shots: Number of measurement shots
            optimization_level: Transpiler optimization (0-3)
            resilience_level: Error mitigation level (0-2, higher for QPE)
            **backend_options: Additional backend options
        """
        self.hamiltonian = hamiltonian
        self.n_ancilla = n_ancilla

        # Initialize backend
        if backend is None:
            if backend_name is None:
                raise ValueError("Either 'backend' or 'backend_name' must be provided")

            logger.info(f"Initializing IBM backend for QPE: {backend_name}")
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
        self.total_qubits = self.n_ancilla + self.n_qubits

        logger.info(f"QPE Configuration:")
        logger.info(f"  System qubits: {self.n_qubits}")
        logger.info(f"  Ancilla qubits: {self.n_ancilla}")
        logger.info(f"  Total qubits: {self.total_qubits}")
        logger.info(f"  Phase precision: {2**(-n_ancilla):.2e}")

    def _build_time_evolution_circuit(self, time: float) -> QuantumCircuit:
        """
        Build time evolution circuit: exp(-iHt)

        For now, uses Trotter decomposition. Can be optimized
        with advanced circuit synthesis methods.
        """
        # Convert Hamiltonian to Pauli operators
        pauli_hamiltonian = self.hamiltonian.to_sparse_hamiltonian()

        # Build time evolution circuit using Qiskit
        from qiskit.synthesis import SuzukiTrotter

        qc = QuantumCircuit(self.n_qubits)

        # Use PauliEvolutionGate for time evolution
        from qiskit.circuit.library import PauliEvolutionGate

        evolution_gate = PauliEvolutionGate(
            pauli_hamiltonian,
            time=time,
            synthesis=SuzukiTrotter()
        )

        qc.append(evolution_gate, range(self.n_qubits))

        return qc

    def _build_qpe_circuit(self, initial_state: Optional[QuantumCircuit] = None) -> QuantumCircuit:
        """
        Build full QPE circuit.

        Circuit structure:
        1. Initialize ancilla in |+⟩ state
        2. Initialize system in trial state (HF or other)
        3. Controlled time evolution operations
        4. Inverse QFT on ancilla
        5. Measure ancilla to extract phase
        """
        # Create quantum registers
        ancilla = QuantumRegister(self.n_ancilla, 'ancilla')
        system = QuantumRegister(self.n_qubits, 'system')
        classical = ClassicalRegister(self.n_ancilla, 'c')

        qc = QuantumCircuit(ancilla, system, classical)

        # 1. Initialize ancilla in superposition
        qc.h(ancilla)

        # 2. Initialize system state
        if initial_state is not None:
            qc.compose(initial_state, qubits=system, inplace=True)
        else:
            # Default: Hartree-Fock state (fill lowest orbitals)
            n_electrons = self.hamiltonian.molecule.n_electrons if hasattr(self.hamiltonian, 'molecule') else self.n_qubits // 2
            for i in range(min(n_electrons, self.n_qubits)):
                qc.x(system[i])

        # 3. Controlled time evolution
        for i in range(self.n_ancilla):
            # Time for this ancilla qubit: 2^i * base_time
            evolution_time = 2 ** i * (2 * np.pi)

            # Build controlled time evolution
            evolution_circuit = self._build_time_evolution_circuit(evolution_time)
            controlled_evolution = evolution_circuit.control(1)

            # Apply controlled evolution
            qc.append(controlled_evolution, [ancilla[i]] + list(system))

        # 4. Inverse QFT on ancilla
        iqft = QFT(self.n_ancilla, inverse=True, do_swaps=False)
        qc.append(iqft, ancilla)

        # 5. Measure ancilla
        qc.measure(ancilla, classical)

        return qc

    def solve(
        self,
        initial_state: Optional[QuantumCircuit] = None
    ) -> Dict:
        """
        Run QPE on IBM Quantum.

        Args:
            initial_state: Initial system state circuit (HF if None)

        Returns:
            Dictionary with:
                - energy: Estimated ground state energy (Hartree)
                - phase: Measured quantum phase
                - counts: Measurement counts
                - precision: Phase estimation precision
                - backend_info: IBM backend information
        """
        logger.info("=" * 80)
        logger.info("IBM QUANTUM QPE")
        logger.info("=" * 80)

        backend_info = self.backend.get_backend_info()
        logger.info(f"Backend: {backend_info['name']}")
        logger.info(f"  Qubits: {backend_info.get('num_qubits', 'N/A')}")
        logger.info(f"  Shots: {self.backend.shots}")

        # Build QPE circuit
        logger.info("\nBuilding QPE circuit...")
        qpe_circuit = self._build_qpe_circuit(initial_state)

        logger.info(f"Circuit: {qpe_circuit.num_qubits} qubits, {qpe_circuit.depth()} depth")

        # Execute circuit
        logger.info("\nSubmitting to IBM Quantum...")
        # Use batch mode (no session) for free/open plans
        sampler = self.backend.get_sampler(session=False)

        job = sampler.run([qpe_circuit])
        result = job.result()

        # Extract measurement counts
        counts = result[0].data.c.get_counts()

        logger.info(f"Job completed. Top measurements:")
        sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        for bitstring, count in sorted_counts[:5]:
            logger.info(f"  {bitstring}: {count}")

        # Extract most likely phase
        most_likely_bitstring = sorted_counts[0][0]
        phase_int = int(most_likely_bitstring, 2)
        phase = phase_int / (2 ** self.n_ancilla)

        # Convert phase to energy
        # E = phase * 2π (depends on time evolution encoding)
        energy = phase * 2 * np.pi

        logger.info("=" * 80)
        logger.info("QPE RESULTS")
        logger.info("=" * 80)
        logger.info(f"Measured Phase: {phase:.6f}")
        logger.info(f"Ground State Energy: {energy:.6f} Ha")
        logger.info(f"                     {energy * 27.211386:.4f} eV")
        logger.info(f"Precision: ±{2**(-self.n_ancilla):.2e}")
        logger.info("=" * 80)

        return {
            'energy': energy,
            'energy_hartree': energy,
            'energy_ev': energy * 27.211386,
            'phase': phase,
            'counts': counts,
            'most_likely_bitstring': most_likely_bitstring,
            'precision': 2 ** (-self.n_ancilla),
            'n_ancilla': self.n_ancilla,
            'backend_info': backend_info,
            'circuit_depth': qpe_circuit.depth(),
            'circuit_qubits': qpe_circuit.num_qubits
        }
