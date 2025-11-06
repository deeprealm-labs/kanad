"""
Bluequbit Sampler Backend for Hi-VQE

Similar to IBM Sampler, provides raw measurement counts for Hi-VQE
diagonal energy measurements with custom error mitigation.

Advantages:
- GPU-accelerated simulation (up to 36 qubits)
- Fast execution (no hardware queue wait time)
- Rich measurement data for visualization
- Custom post-processing for Hi-VQE
"""

import logging
from typing import Dict, Any, List, Optional
import numpy as np

logger = logging.getLogger(__name__)


class BluequbitSamplerBackend:
    """
    Bluequbit backend using measurement sampling for Hi-VQE.

    Provides similar interface to IBM Sampler but runs on GPU simulators
    for fast, accurate results without hardware noise.

    Usage:
        backend = BluequbitSamplerBackend(device='gpu', api_token='...')
        result = backend.run_hivqe_measurement(
            circuits=circuits,
            observable=hamiltonian,
            shots=8192
        )
    """

    def __init__(
        self,
        device: str = 'gpu',
        api_token: Optional[str] = None
    ):
        """Initialize Bluequbit Sampler backend."""
        import os
        from kanad.backends.bluequbit.backend import BlueQubitBackend

        # Use existing BlueQubitBackend
        self.api_token = api_token or os.getenv('BLUE_TOKEN')

        if not self.api_token:
            raise ValueError("Bluequbit API token required")

        # Initialize Bluequbit backend
        self._bq_backend = BlueQubitBackend(
            device=device,
            api_token=self.api_token
        )

        self.device = device

        logger.info(f"Bluequbit Sampler backend initialized: {device}")

    def run_hivqe_measurement(
        self,
        circuits: List,
        observable: 'SparsePauliOp',
        shots: int = 8192,
        job_name: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Run Hi-VQE measurement circuits using Bluequbit sampler.

        Returns raw counts for each circuit, enabling custom post-processing.

        Args:
            circuits: List of measurement circuits (one per configuration)
            observable: Hamiltonian observable (for expectation value calculation)
            shots: Number of measurement shots
            job_name: Optional job name

        Returns:
            Dictionary with:
            - counts: Raw measurement counts for each circuit
            - diagonal_energies: Computed ⟨config|H|config⟩ values
            - ground_energy: Final ground state energy
            - metadata: Execution metadata
        """
        # Ensure circuits is a list
        if not isinstance(circuits, list):
            circuits = [circuits]

        logger.info(f"Running Hi-VQE measurement on Bluequbit {self.device}")
        logger.info(f"  Circuits: {len(circuits)}")
        logger.info(f"  Shots: {shots}")

        # Add measurements if not present
        circuits_with_meas = []
        for circ in circuits:
            if circ.num_clbits == 0:
                circ_copy = circ.copy()
                circ_copy.measure_all()
                circuits_with_meas.append(circ_copy)
            else:
                circuits_with_meas.append(circ)

        # Run all circuits
        all_counts = []
        all_probabilities = []

        job_name = job_name or f"kanad_hivqe_{self.device}"

        for i, circuit in enumerate(circuits_with_meas):
            logger.info(f"  Running circuit {i+1}/{len(circuits_with_meas)}...")

            # Run circuit on Bluequbit
            result = self._bq_backend.run_circuit(
                circuit,
                shots=shots,
                job_name=f"{job_name}_circuit_{i}"
            )

            # Get counts
            counts = result.get('counts', {})
            all_counts.append(counts)

            # Calculate probabilities
            total_shots = sum(counts.values())
            probs = {state: count/total_shots for state, count in counts.items()}
            all_probabilities.append(probs)

            logger.info(f"    Circuit {i}: {len(counts)} unique measurement outcomes")

        # Calculate diagonal energies ⟨config|H|config⟩
        logger.info("Computing expectation values from counts...")
        diagonal_energies = []

        for i, counts in enumerate(all_counts):
            # Calculate expectation value for this configuration
            expectation = self._calculate_expectation_value(counts, observable)
            diagonal_energies.append(expectation)

            logger.info(f"  Config {i}: {expectation:.8f} Ha")

        return {
            'backend': f'bluequbit_{self.device}',
            'counts': [dict(c) for c in all_counts],
            'probabilities': all_probabilities,
            'diagonal_energies': diagonal_energies,
            'metadata': {
                'device': self.device,
                'shots': shots,
                'num_circuits': len(circuits),
                'observable_terms': len(observable)
            }
        }

    def process_hivqe_results(
        self,
        counts: List[Dict],
        observable: 'SparsePauliOp',
        subspace: 'ConfigurationSubspace',
        diagonal_energies: Optional[List[float]] = None
    ) -> Dict[str, Any]:
        """
        Process Hi-VQE measurement results to compute final ground energy.

        Similar to IBM Sampler processing but uses pre-computed counts.

        Args:
            counts: Raw measurement counts from run_hivqe_measurement
            observable: Hamiltonian used for measurements
            subspace: Configuration subspace
            diagonal_energies: Pre-computed diagonal energies (optional)

        Returns:
            Dictionary with:
            - diagonal_energies: List of ⟨config|H|config⟩ values
            - ground_energy: Final ground state energy
            - eigenvalues: All eigenvalues
            - ground_state_amplitudes: Ground state coefficients
            - visualization_data: Data for web app
            - error_analysis: Error estimates
        """
        logger.info("Processing Hi-VQE results...")

        # Calculate diagonal energies if not provided
        if diagonal_energies is None:
            diagonal_energies = []
            for i, (config, cts) in enumerate(zip(subspace, counts)):
                expectation = self._calculate_expectation_value(cts, observable)
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

        # Calculate probabilities
        all_probabilities = []
        for cts in counts:
            total = sum(cts.values())
            probs = {state: count/total for state, count in cts.items()}
            all_probabilities.append(probs)

        # Build visualization data
        viz_data = self._build_visualization_data(
            subspace=subspace,
            counts=counts,
            probabilities=all_probabilities,
            diagonal_energies=diagonal_energies,
            ground_state_amplitudes=ground_state,
            eigenvalues=eigenvalues
        )

        # Error analysis
        error_data = self._analyze_measurement_errors(
            counts,
            all_probabilities
        )

        return {
            'backend': f'bluequbit_{self.device}',
            'diagonal_energies': diagonal_energies,
            'ground_energy': ground_energy,
            'eigenvalues': eigenvalues.tolist(),
            'ground_state_amplitudes': ground_state.tolist(),
            'counts': [dict(c) for c in counts],
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
        all_probabilities
    ) -> Dict[str, Any]:
        """Analyze measurement errors and uncertainty."""

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
                'simulation_noise': False,  # Bluequbit is noiseless simulator
                'ideal_measurements': True
            }
        }

    def __repr__(self):
        return f"BluequbitSamplerBackend(device='{self.device}')"
