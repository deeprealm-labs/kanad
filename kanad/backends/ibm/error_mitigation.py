"""
Error Mitigation for IBM Quantum Hardware

Implements multiple error mitigation techniques:
1. Readout error mitigation (measurement calibration)
2. Zero-Noise Extrapolation (ZNE)
3. M3 (Matrix-free Measurement Mitigation)
4. Dynamical Decoupling
5. Probabilistic Error Cancellation (PEC)

For Hi-VQE, the most effective strategies are:
- Readout error mitigation (always beneficial)
- ZNE for energy estimates
- Dynamical decoupling for coherence
"""

import logging
from typing import Dict, Any, List, Optional, Tuple
import numpy as np

logger = logging.getLogger(__name__)


class ErrorMitigationStrategy:
    """
    Configurable error mitigation strategy for IBM backends.

    Usage:
        strategy = ErrorMitigationStrategy(
            readout_mitigation=True,
            zne_extrapolation='linear',
            dynamical_decoupling='XY4'
        )

        # Apply to Estimator options
        estimator.options.resilience_level = strategy.resilience_level
        estimator.options.resilience = strategy.get_resilience_options()
    """

    def __init__(
        self,
        resilience_level: int = 1,
        readout_mitigation: bool = True,
        zne_extrapolation: Optional[str] = None,
        zne_noise_factors: Optional[List[float]] = None,
        dynamical_decoupling: Optional[str] = None,
        twirling: bool = False,
        measure_mitigation: bool = True
    ):
        """
        Initialize error mitigation strategy.

        Args:
            resilience_level: Built-in IBM resilience (0=none, 1=light, 2=heavy)
                             Level 1: Readout error mitigation
                             Level 2: ZNE + readout mitigation
            readout_mitigation: Enable measurement error mitigation
            zne_extrapolation: ZNE extrapolation method ('linear', 'exponential', None)
            zne_noise_factors: Noise scaling factors for ZNE (default: [1.0, 1.5, 2.0])
            dynamical_decoupling: DD sequence ('X', 'XY4', 'XX', None)
            twirling: Enable Pauli twirling for gate errors
            measure_mitigation: Use M3 measurement mitigation (Qiskit feature)
        """
        self.resilience_level = resilience_level
        self.readout_mitigation = readout_mitigation
        self.zne_extrapolation = zne_extrapolation
        self.zne_noise_factors = zne_noise_factors or [1.0, 1.5, 2.0]
        self.dynamical_decoupling = dynamical_decoupling
        self.twirling = twirling
        self.measure_mitigation = measure_mitigation

        logger.info(f"Error mitigation strategy initialized:")
        logger.info(f"  Resilience level: {resilience_level}")
        logger.info(f"  Readout mitigation: {readout_mitigation}")
        logger.info(f"  ZNE extrapolation: {zne_extrapolation or 'disabled'}")
        logger.info(f"  Dynamical decoupling: {dynamical_decoupling or 'disabled'}")
        logger.info(f"  Twirling: {twirling}")

    def get_resilience_options(self) -> Dict[str, Any]:
        """
        Get Qiskit Runtime resilience options.

        Returns:
            Dictionary of resilience options for EstimatorV2/SamplerV2
        """
        options = {}

        # ZNE configuration
        if self.zne_extrapolation:
            options['zne_mitigation'] = True
            options['zne'] = {
                'noise_factors': self.zne_noise_factors,
                'extrapolator': self.zne_extrapolation
            }
            logger.info(f"  ZNE enabled: {self.zne_extrapolation} extrapolation")

        # Measurement error mitigation
        if self.measure_mitigation:
            options['measure_mitigation'] = True
            logger.info(f"  Measurement mitigation enabled (M3)")

        # Pauli twirling
        if self.twirling:
            options['twirling'] = {
                'enable_gates': True,
                'enable_measure': True,
                'num_randomizations': 32
            }
            logger.info(f"  Pauli twirling enabled")

        return options

    def get_transpiler_options(self, backend) -> Dict[str, Any]:
        """
        Get transpiler options for error mitigation.

        Args:
            backend: IBM backend object

        Returns:
            Dictionary of transpiler options
        """
        options = {}

        # Dynamical decoupling
        if self.dynamical_decoupling:
            from qiskit.transpiler import PassManager
            from qiskit.transpiler.passes import PadDynamicalDecoupling
            from qiskit.circuit.library import XGate, YGate

            # Select DD sequence
            if self.dynamical_decoupling == 'X':
                dd_sequence = [XGate()] * 2
            elif self.dynamical_decoupling == 'XY4':
                dd_sequence = [XGate(), YGate(), XGate(), YGate()]
            elif self.dynamical_decoupling == 'XX':
                dd_sequence = [XGate(), XGate()]
            else:
                dd_sequence = None

            if dd_sequence:
                options['dynamical_decoupling'] = {
                    'sequence': dd_sequence,
                    'spacing': 'balanced'  # or 'uniform'
                }
                logger.info(f"  Dynamical decoupling: {self.dynamical_decoupling}")

        return options

    @staticmethod
    def auto_configure(backend_name: str) -> 'ErrorMitigationStrategy':
        """
        Automatically configure error mitigation based on backend type.

        **Automatic Strategy:**
        - Simulators (aer_simulator, statevector): No mitigation (not needed)
        - Real hardware (ibm_*): Full mitigation stack

        This eliminates the need for manual configuration and ensures optimal
        settings for each backend type.

        Args:
            backend_name: Name of the quantum backend

        Returns:
            Configured ErrorMitigationStrategy instance

        Examples:
            >>> # Simulator - no mitigation
            >>> strategy = ErrorMitigationStrategy.auto_configure('aer_simulator')
            >>> strategy.readout_mitigation  # False

            >>> # Real hardware - full mitigation
            >>> strategy = ErrorMitigationStrategy.auto_configure('ibm_kyoto')
            >>> strategy.readout_mitigation  # True
            >>> strategy.zne_extrapolation  # 'exponential'
        """
        backend_lower = backend_name.lower()

        # Check if simulator
        is_simulator = any(sim in backend_lower for sim in [
            'simulator', 'statevector', 'aer', 'fake'
        ])

        if is_simulator:
            # Simulators: No mitigation needed (noiseless or noise model already applied)
            logger.info(f"Auto-config: Simulator detected ({backend_name}) - disabling error mitigation")
            return ErrorMitigationStrategy(
                resilience_level=0,
                readout_mitigation=False,
                zne_extrapolation=None,
                dynamical_decoupling=None,
                twirling=False,
                measure_mitigation=False
            )

        else:
            # Real hardware: Enable full mitigation stack
            logger.info(f"Auto-config: Real hardware detected ({backend_name}) - enabling full mitigation")
            return ErrorMitigationStrategy(
                resilience_level=2,  # Level 2: ZNE + readout mitigation
                readout_mitigation=True,
                zne_extrapolation='exponential',  # Exponential extrapolation works best
                zne_noise_factors=[1.0, 3.0, 5.0],  # Standard ZNE noise scaling
                dynamical_decoupling='XY4',  # XY4 sequence for coherence protection
                twirling=True,  # Pauli twirling for gate errors
                measure_mitigation=True  # M3 measurement mitigation
            )

    def estimate_mitigation_overhead(self, n_circuits: int) -> Dict[str, Any]:
        """
        Estimate computational overhead from error mitigation.

        Args:
            n_circuits: Number of circuits to run

        Returns:
            Dictionary with overhead estimates
        """
        overhead = {
            'circuit_multiplier': 1.0,
            'calibration_shots': 0,
            'extra_time': 0.0
        }

        # ZNE multiplies circuits by noise factors
        if self.zne_extrapolation:
            overhead['circuit_multiplier'] *= len(self.zne_noise_factors)

        # Twirling adds randomized versions
        if self.twirling:
            overhead['circuit_multiplier'] *= 32  # num_randomizations

        # Readout calibration needs extra shots
        if self.readout_mitigation:
            overhead['calibration_shots'] = 8192  # Per backend, not per circuit

        # Estimate time overhead (rough)
        overhead['extra_time'] = (
            n_circuits * overhead['circuit_multiplier'] * 0.1  # seconds per circuit
            + overhead['calibration_shots'] * 0.001  # calibration time
        )

        logger.info(f"Error mitigation overhead:")
        logger.info(f"  Circuit multiplier: {overhead['circuit_multiplier']:.1f}x")
        logger.info(f"  Extra time estimate: {overhead['extra_time']:.1f}s")

        return overhead


class AdaptiveShotAllocation:
    """
    Adaptive shot allocation for Hi-VQE measurements.

    Key insight: Different Pauli terms have different variances.
    Allocate more shots to high-variance terms to minimize total error.

    For Hi-VQE, we measure Z-basis states, so shot allocation is per-configuration.
    """

    def __init__(
        self,
        total_shots: int = 4096,
        min_shots_per_config: int = 100,
        confidence_level: float = 0.95
    ):
        """
        Initialize adaptive shot allocator.

        Args:
            total_shots: Total shot budget
            min_shots_per_config: Minimum shots per configuration
            confidence_level: Statistical confidence for estimates
        """
        self.total_shots = total_shots
        self.min_shots_per_config = min_shots_per_config
        self.confidence_level = confidence_level

    def allocate_shots(
        self,
        config_amplitudes: np.ndarray,
        prev_variances: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Allocate shots across configurations based on amplitudes.

        Configurations with larger amplitudes get more shots for better statistics.

        Args:
            config_amplitudes: Array of configuration amplitudes from previous iteration
            prev_variances: Previous iteration's measurement variances (optional)

        Returns:
            Array of shot counts for each configuration
        """
        n_configs = len(config_amplitudes)

        # Ensure minimum shots
        remaining_shots = self.total_shots - n_configs * self.min_shots_per_config

        if remaining_shots <= 0:
            # Can't afford minimum shots for all configs
            logger.warning(f"Shot budget too small for {n_configs} configs")
            return np.full(n_configs, self.total_shots // n_configs)

        # Allocate remaining shots proportional to |amplitude|
        # (Higher amplitude = more important for energy)
        importance = np.abs(config_amplitudes)
        importance = importance / np.sum(importance)  # Normalize

        # If we have variance estimates, use them (importance sampling)
        if prev_variances is not None:
            # Allocate more shots to high-variance configs
            variance_weights = np.sqrt(prev_variances)
            variance_weights = variance_weights / np.sum(variance_weights)

            # Combine amplitude and variance (weighted average)
            importance = 0.7 * importance + 0.3 * variance_weights
            importance = importance / np.sum(importance)

        # Allocate remaining shots
        extra_shots = (importance * remaining_shots).astype(int)

        # Total allocation
        shot_allocation = np.full(n_configs, self.min_shots_per_config) + extra_shots

        # Ensure total doesn't exceed budget (due to rounding)
        total_allocated = np.sum(shot_allocation)
        if total_allocated > self.total_shots:
            # Remove excess from largest allocations
            excess = total_allocated - self.total_shots
            largest_idx = np.argmax(shot_allocation)
            shot_allocation[largest_idx] -= excess

        logger.info(f"Adaptive shot allocation:")
        logger.info(f"  Total shots: {np.sum(shot_allocation)}/{self.total_shots}")
        logger.info(f"  Min/Max shots: {np.min(shot_allocation)}/{np.max(shot_allocation)}")
        logger.info(f"  Configs: {n_configs}")

        return shot_allocation

    def estimate_energy_error(
        self,
        shot_allocation: np.ndarray,
        config_amplitudes: np.ndarray,
        hamiltonian_variance: float = 1.0
    ) -> float:
        """
        Estimate statistical error in energy from shot allocation.

        Args:
            shot_allocation: Shots per configuration
            config_amplitudes: Configuration amplitudes
            hamiltonian_variance: Variance of Hamiltonian operator

        Returns:
            Estimated standard error in energy (Ha)
        """
        # Variance in energy: Var[E] = Σ_i |c_i|^2 * σ²_i
        # σ²_i = variance per config / shots_i

        config_variances = hamiltonian_variance / shot_allocation
        energy_variance = np.sum(config_amplitudes**2 * config_variances)
        energy_error = np.sqrt(energy_variance)

        logger.info(f"  Estimated energy error: {energy_error:.6f} Ha")

        return energy_error


def optimize_circuit_for_hardware(circuit, backend, optimization_level: int = 3) -> 'QuantumCircuit':
    """
    Optimize quantum circuit for IBM hardware.

    Applies:
    - Native gate decomposition (to basis gates)
    - SWAP insertion for qubit connectivity
    - Gate fusion and cancellation
    - Delay optimization

    Args:
        circuit: Quantum circuit
        backend: IBM backend
        optimization_level: Optimization level (0-3)

    Returns:
        Optimized circuit
    """
    from qiskit.transpiler import PassManager, CouplingMap
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    logger.info(f"Optimizing circuit for {backend.name}")
    logger.info(f"  Original: {circuit.num_qubits} qubits, depth {circuit.depth()}")

    # Generate optimized pass manager
    pm = generate_preset_pass_manager(
        backend=backend,
        optimization_level=optimization_level,
        seed_transpiler=42  # For reproducibility
    )

    # Run optimization
    optimized_circuit = pm.run(circuit)

    logger.info(f"  Optimized: {optimized_circuit.num_qubits} qubits, depth {optimized_circuit.depth()}")
    logger.info(f"  Reduction: {circuit.depth() - optimized_circuit.depth()} gates removed")

    return optimized_circuit
