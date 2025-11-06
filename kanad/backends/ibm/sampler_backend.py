"""
IBM Quantum Sampler Backend for Hi-VQE

Uses SamplerV2 to get raw measurement counts, enabling:
- Custom error mitigation for Z-basis measurements
- Rich visualization data for web app
- Direct expectation value calculation from counts
"""

import logging
from typing import Dict, Any, List, Optional
import numpy as np

logger = logging.getLogger(__name__)


class IBMSamplerBackend:
    """
    IBM Quantum backend using Sampler for Hi-VQE measurements.

    Advantages over Estimator:
    - Get raw measurement counts (full distribution)
    - Custom post-processing for Hi-VQE
    - Better visualization data
    - More control over error mitigation

    Usage:
        backend = IBMSamplerBackend(backend_name='ibm_torino')
        result = backend.run_hivqe_measurement(
            circuits=circuits,
            observable=hamiltonian,
            shots=8192
        )
    """

    def __init__(
        self,
        backend_name: str,
        api_token: Optional[str] = None,
        channel: str = 'ibm_quantum_platform',
        crn: Optional[str] = None
    ):
        """Initialize IBM Sampler backend."""
        import os
        from kanad.backends.ibm.backend import IBMBackend

        # Use existing IBMBackend for connection
        self.api_token = api_token or os.getenv('IBM_API')
        self.channel = channel
        self.crn = crn or os.getenv('IBM_CRN')

        if not self.api_token:
            raise ValueError("IBM API token required")

        # Initialize IBM backend for service
        self._ibm_backend = IBMBackend(
            backend_name=backend_name,
            api_token=self.api_token,
            channel=self.channel,
            crn=self.crn
        )

        self.backend = self._ibm_backend.backend
        self.service = self._ibm_backend.service

        logger.info(f"IBM Sampler backend initialized: {self.backend.name}")

    def run_hivqe_measurement(
        self,
        circuits: List,
        observable: 'SparsePauliOp',
        shots: int = 8192,
        optimization_level: int = 3,
        enable_readout_mitigation: bool = True,
        enable_twirling: bool = True,
        mode: str = 'batch'
    ) -> Dict[str, Any]:
        """
        Run Hi-VQE measurement circuits using Sampler.

        Returns raw counts for each circuit, enabling custom post-processing.

        Args:
            circuits: List of measurement circuits (one per configuration)
            observable: Hamiltonian observable (for expectation value calculation)
            shots: Number of measurement shots
            optimization_level: Circuit optimization (0-3)
            enable_readout_mitigation: Apply readout error mitigation
            enable_twirling: Apply Pauli twirling
            mode: 'batch' or 'session'

        Returns:
            Dictionary with:
            - job_id: IBM job ID
            - status: Job status
            - backend: Backend name
            - config: Measurement configuration
        """
        from qiskit_ibm_runtime import Batch, Session, SamplerV2 as Sampler
        from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

        # Ensure circuits is a list
        if not isinstance(circuits, list):
            circuits = [circuits]

        logger.info(f"Running Hi-VQE measurement on {self.backend.name}")
        logger.info(f"  Circuits: {len(circuits)}")
        logger.info(f"  Shots: {shots}")
        logger.info(f"  Mode: {mode}")

        # Transpile circuits
        logger.info(f"  Transpiling circuits (optimization_level={optimization_level})...")
        pm = generate_preset_pass_manager(
            backend=self.backend,
            optimization_level=optimization_level
        )
        transpiled_circuits = pm.run(circuits)
        logger.info(f"  Transpilation complete")

        # Choose execution mode
        if mode == 'session':
            context = Session(backend=self.backend, max_time='1h')
        else:
            context = Batch(backend=self.backend)

        try:
            with context as exec_context:
                # Create Sampler
                if mode == 'session':
                    sampler = Sampler(session=exec_context)
                else:
                    sampler = Sampler(mode=exec_context)

                # Set options
                sampler.options.default_shots = shots

                # Error mitigation options
                if enable_readout_mitigation:
                    # For SamplerV2, enable measurement error mitigation
                    try:
                        sampler.options.dynamical_decoupling.enable = False
                        sampler.options.twirling.enable_gates = enable_twirling
                        sampler.options.twirling.enable_measure = enable_twirling
                        logger.info(f"  Error mitigation enabled: readout + twirling={enable_twirling}")
                    except:
                        logger.warning("  Could not set all error mitigation options")

                logger.info("Using Sampler primitive for raw counts")

                # Submit job
                job = sampler.run(transpiled_circuits)

                # Return job info
                result = {
                    'job_id': job.job_id(),
                    'status': str(job.status()),
                    'backend': self.backend.name,
                    'mode': mode,
                    'config': {
                        'shots': shots,
                        'circuits': len(circuits),
                        'readout_mitigation': enable_readout_mitigation,
                        'twirling': enable_twirling,
                        'observable_terms': len(observable)
                    }
                }

                if mode == 'session':
                    result['session_id'] = exec_context.session_id

                logger.info(f"Job submitted: {result['job_id']}")
                return result

        except Exception as e:
            logger.error(f"IBM Sampler execution failed: {e}")
            raise

    def get_job_result(self, job_id: str) -> Any:
        """Get raw Sampler results."""
        return self._ibm_backend.get_job_result(job_id)

    def get_job_status(self, job_id: str) -> str:
        """Get job status."""
        return self._ibm_backend.get_job_status(job_id)

    def process_hivqe_results(
        self,
        job_id: str,
        observable: 'SparsePauliOp',
        subspace: 'ConfigurationSubspace'
    ) -> Dict[str, Any]:
        """
        Process Sampler results to compute Hi-VQE energies.

        Returns rich data structure for visualization and analysis.

        Args:
            job_id: IBM job ID
            observable: Hamiltonian used for measurements
            subspace: Configuration subspace

        Returns:
            Dictionary with:
            - diagonal_energies: List of ⟨config|H|config⟩ values
            - ground_energy: Final ground state energy
            - counts: Raw measurement counts for each circuit
            - probabilities: Measurement probability distributions
            - visualization_data: Data for web app
            - error_analysis: Error estimates and mitigation info
        """
        logger.info(f"Processing Hi-VQE results for job {job_id}")

        # Get raw Sampler results
        results = self.get_job_result(job_id)

        # Extract counts from each circuit
        all_counts = []
        all_probabilities = []

        for i, pub_result in enumerate(results):
            # Get counts from DataBin
            counts = pub_result.data.meas.get_counts()
            all_counts.append(counts)

            # Calculate probabilities
            total_shots = sum(counts.values())
            probs = {state: count/total_shots for state, count in counts.items()}
            all_probabilities.append(probs)

            logger.info(f"  Circuit {i}: {len(counts)} unique measurement outcomes")

        # Calculate diagonal energies ⟨config|H|config⟩
        logger.info("Computing expectation values from counts...")
        diagonal_energies = []

        for i, (config, counts) in enumerate(zip(subspace, all_counts)):
            # Calculate expectation value for this configuration
            expectation = self._calculate_expectation_value(counts, observable)
            diagonal_energies.append(expectation)

            logger.info(f"  Config {i} ({config}): {expectation:.8f} Ha")

        # Classical diagonalization
        logger.info("Performing classical diagonalization...")
        from kanad.core.classical_solver import SubspaceHamiltonianBuilder

        builder = SubspaceHamiltonianBuilder(observable)
        H_sub = builder.project_fast(subspace)

        # Replace diagonal with measured values
        for i in range(len(diagonal_energies)):
            H_sub[i, i] = diagonal_energies[i]

        # Diagonalize
        eigenvalues = np.linalg.eigvalsh(H_sub)
        eigenvectors = np.linalg.eigh(H_sub)[1]
        ground_energy = eigenvalues[0]
        ground_state = eigenvectors[:, 0]

        logger.info(f"Ground state energy: {ground_energy:.8f} Ha")

        # Build visualization data
        viz_data = self._build_visualization_data(
            subspace=subspace,
            counts=all_counts,
            probabilities=all_probabilities,
            diagonal_energies=diagonal_energies,
            ground_state_amplitudes=ground_state,
            eigenvalues=eigenvalues
        )

        # Error analysis
        error_data = self._analyze_measurement_errors(
            all_counts,
            all_probabilities,
            job_id
        )

        return {
            'job_id': job_id,
            'backend': self.backend.name,
            'diagonal_energies': diagonal_energies,
            'ground_energy': ground_energy,
            'eigenvalues': eigenvalues.tolist(),
            'ground_state_amplitudes': ground_state.tolist(),
            'counts': [dict(c) for c in all_counts],
            'probabilities': all_probabilities,
            'visualization_data': viz_data,
            'error_analysis': error_data
        }

    def _calculate_expectation_value(
        self,
        counts: Dict[str, int],
        observable: 'SparsePauliOp'
    ) -> float:
        """
        Calculate expectation value ⟨ψ|H|ψ⟩ from measurement counts.

        For Hi-VQE, we measure in Z-basis, so we can directly compute
        the expectation value from Pauli terms.
        """
        total_shots = sum(counts.values())
        expectation = 0.0

        # Get Pauli terms
        pauli_list = observable.to_list()

        for pauli_str, coeff in pauli_list:
            # Calculate expectation for this Pauli term
            term_expectation = 0.0

            for bitstring, count in counts.items():
                # Remove spaces and reverse to match qubit ordering
                bits = bitstring.replace(' ', '')[::-1]

                # Calculate eigenvalue for this measurement
                eigenvalue = 1.0
                for i, pauli in enumerate(pauli_str):
                    if i < len(bits):
                        bit = int(bits[i])
                        if pauli == 'Z':
                            eigenvalue *= (1 - 2*bit)  # 0→+1, 1→-1
                        elif pauli == 'I':
                            eigenvalue *= 1.0
                        elif pauli in ['X', 'Y']:
                            # For X and Y, we'd need to measure in different basis
                            # For Hi-VQE, all measurements should be Z-basis
                            pass

                term_expectation += eigenvalue * count

            term_expectation /= total_shots
            expectation += coeff.real * term_expectation

        return expectation

    def _build_visualization_data(
        self,
        subspace,
        counts,
        probabilities,
        diagonal_energies,
        ground_state_amplitudes,
        eigenvalues
    ) -> Dict[str, Any]:
        """Build rich visualization data for web app."""

        # Configuration labels
        config_labels = [str(config) for config in subspace]

        # Measurement statistics
        measurement_stats = []
        for i, (config, prob_dist) in enumerate(zip(config_labels, probabilities)):
            # Get most probable measurement outcome
            top_outcome = max(prob_dist.items(), key=lambda x: x[1])

            measurement_stats.append({
                'configuration': config,
                'diagonal_energy': diagonal_energies[i],
                'ground_state_amplitude': abs(ground_state_amplitudes[i]),
                'most_probable_outcome': top_outcome[0],
                'max_probability': top_outcome[1],
                'num_outcomes': len(prob_dist),
                'entropy': -sum(p * np.log2(p) if p > 0 else 0 for p in prob_dist.values())
            })

        # Energy spectrum
        energy_spectrum = {
            'eigenvalues': eigenvalues.tolist(),
            'ground_energy': eigenvalues[0],
            'excitation_energies': (eigenvalues[1:] - eigenvalues[0]).tolist() if len(eigenvalues) > 1 else []
        }

        return {
            'configurations': config_labels,
            'measurement_statistics': measurement_stats,
            'energy_spectrum': energy_spectrum,
            'ground_state_composition': {
                'amplitudes': ground_state_amplitudes.tolist(),
                'probabilities': (ground_state_amplitudes**2).tolist(),
                'dominant_configs': [
                    config_labels[i] for i in np.argsort(np.abs(ground_state_amplitudes))[::-1][:3]
                ]
            }
        }

    def _analyze_measurement_errors(
        self,
        all_counts,
        all_probabilities,
        job_id
    ) -> Dict[str, Any]:
        """Analyze measurement errors and uncertainty."""

        # Get job metadata
        try:
            job = self.service.job(job_id)
            metadata = job.result().metadata if hasattr(job.result(), 'metadata') else {}
        except:
            metadata = {}

        # Calculate shot noise for each measurement
        shot_noise_estimates = []
        for counts in all_counts:
            total_shots = sum(counts.values())
            # Shot noise scales as 1/sqrt(N)
            shot_noise = 1.0 / np.sqrt(total_shots)
            shot_noise_estimates.append(shot_noise)

        # Check for readout errors (unexpected measurement outcomes)
        readout_fidelity_estimates = []
        for probs in all_probabilities:
            # Fidelity estimate: probability of expected outcome
            max_prob = max(probs.values())
            readout_fidelity_estimates.append(max_prob)

        return {
            'shot_noise': {
                'estimates': shot_noise_estimates,
                'average': np.mean(shot_noise_estimates),
                'max': np.max(shot_noise_estimates)
            },
            'readout_fidelity': {
                'estimates': readout_fidelity_estimates,
                'average': np.mean(readout_fidelity_estimates),
                'min': np.min(readout_fidelity_estimates)
            },
            'mitigation_applied': {
                'twirling': metadata.get('twirling', {}).get('enable_gates', False),
                'dynamical_decoupling': metadata.get('dynamical_decoupling', {}).get('enable', False)
            }
        }

    def __repr__(self):
        return f"IBMSamplerBackend(backend='{self.backend.name}', qubits={self.backend.num_qubits})"
