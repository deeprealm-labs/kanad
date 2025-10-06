"""
Uncertainty Analysis for Quantum Measurements.

Provides statistical analysis of shot noise, confidence intervals,
and convergence analysis for VQE/QPE results.
"""

import numpy as np
from typing import Dict, Any, List, Tuple, Optional, Callable
from scipy import stats
import logging

logger = logging.getLogger(__name__)


class UncertaintyAnalyzer:
    """
    Statistical uncertainty analysis for quantum measurements.

    Supports:
    - Shot noise estimation from finite measurements
    - Bootstrap confidence intervals
    - Convergence analysis (error vs shots)
    - Multiple run comparison
    - Publication-ready error bars
    """

    def __init__(self, backend: str = 'statevector'):
        """
        Initialize uncertainty analyzer.

        Args:
            backend: 'statevector' (exact) or 'qasm' (shot noise)
        """
        self.backend = backend

    def estimate_shot_noise(self,
                           pauli_expectations: Dict[str, float],
                           n_shots: int = 1024,
                           method: str = 'theoretical') -> Dict[str, Any]:
        """
        Estimate energy uncertainty from shot noise.

        Theory:
            σ²_E = Σ_P var(⟨P⟩) / n_shots × c_P²

        where P are Pauli strings, c_P are coefficients.

        Args:
            pauli_expectations: {pauli_string: expectation_value}
                Example: {'II': 1.0, 'ZI': 0.5, 'IZ': 0.3, 'ZZ': 0.1}
            n_shots: Number of measurement shots
            method: 'theoretical' or 'empirical'

        Returns:
            result: Dictionary with:
                - energy_variance: σ²_E
                - energy_std: σ_E (standard deviation)
                - standard_error: σ_E / √n_shots
                - confidence_interval_95: (lower, upper)
        """
        # Variance of Pauli expectation value measurement
        # For Pauli P: ⟨P⟩ ∈ [-1, 1]
        # var(⟨P⟩) ≈ (1 - ⟨P⟩²) for single Pauli measurement

        total_variance = 0.0

        for pauli_string, expectation in pauli_expectations.items():
            # Variance of individual Pauli measurement
            if method == 'theoretical':
                # Theoretical maximum: var(⟨P⟩) = 1 - ⟨P⟩²
                var_pauli = 1.0 - expectation**2
            else:
                # Conservative estimate
                var_pauli = 1.0

            # Contribution to energy variance
            # Assuming equal coefficients (simplified)
            total_variance += var_pauli

        # Divide by number of shots
        energy_variance = total_variance / n_shots

        energy_std = np.sqrt(energy_variance)
        standard_error = energy_std  # Already normalized by √n_shots in formula

        # 95% confidence interval (1.96 sigma)
        ci_95 = (
            -1.96 * standard_error,  # Relative to mean
            +1.96 * standard_error
        )

        return {
            'energy_variance': energy_variance,
            'energy_std': energy_std,
            'standard_error': standard_error,
            'confidence_interval_95': ci_95,
            'n_shots': n_shots,
            'n_paulis': len(pauli_expectations)
        }

    def bootstrap_confidence_interval(self,
                                     measurements: np.ndarray,
                                     energy_fn: Callable,
                                     n_bootstrap: int = 1000,
                                     confidence: float = 0.95) -> Dict[str, Any]:
        """
        Bootstrap resampling for confidence intervals.

        Args:
            measurements: (n_shots, n_qubits) array of measurement bitstrings
                         OR (n_shots,) array of energies
            energy_fn: Function to compute energy from measurements
                      (optional if measurements are already energies)
            n_bootstrap: Number of bootstrap samples
            confidence: Confidence level (0.95 → 95%)

        Returns:
            result: Dictionary with:
                - mean: Bootstrap mean
                - std: Bootstrap standard deviation
                - ci_lower: Lower confidence bound
                - ci_upper: Upper confidence bound
                - samples: All bootstrap samples
        """
        if measurements.ndim == 1 and energy_fn is None:
            # Measurements are already energies
            energies = measurements
        else:
            # Compute initial energy
            energies = energy_fn(measurements)

        n_shots = len(measurements)
        bootstrap_samples = []

        for _ in range(n_bootstrap):
            # Resample with replacement
            indices = np.random.choice(n_shots, size=n_shots, replace=True)

            if measurements.ndim == 1:
                resample = measurements[indices]
                if energy_fn is not None:
                    E_bootstrap = energy_fn(resample)
                else:
                    E_bootstrap = np.mean(resample)
            else:
                resample = measurements[indices]
                E_bootstrap = energy_fn(resample)

            bootstrap_samples.append(E_bootstrap)

        bootstrap_samples = np.array(bootstrap_samples)

        # Compute confidence interval
        alpha = 1 - confidence
        ci_lower = np.percentile(bootstrap_samples, alpha/2 * 100)
        ci_upper = np.percentile(bootstrap_samples, (1 - alpha/2) * 100)

        return {
            'mean': np.mean(bootstrap_samples),
            'std': np.std(bootstrap_samples),
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'confidence': confidence,
            'n_bootstrap': n_bootstrap,
            'samples': bootstrap_samples
        }

    def convergence_analysis(self,
                            energy_fn: Callable[[int], float],
                            shot_counts: List[int] = [64, 128, 256, 512, 1024, 2048, 4096],
                            n_trials: int = 10,
                            exact_energy: Optional[float] = None) -> Dict[str, Any]:
        """
        Analyze energy convergence with number of shots.

        Args:
            energy_fn: Function that takes n_shots and returns energy
                      Example: lambda n: vqe.solve(n_shots=n)['energy']
            shot_counts: List of shot counts to test
            n_trials: Number of trials per shot count
            exact_energy: True energy (if known) for error calculation

        Returns:
            result: Dictionary with:
                - shot_counts: Input shot counts
                - mean_energies: Average energy at each count
                - std_energies: Standard deviation
                - errors: |E_measured - E_exact| (if exact_energy given)
                - theoretical_se: Expected SE = σ/√n_shots
        """
        results = {
            'shot_counts': shot_counts,
            'mean_energies': [],
            'std_energies': [],
            'min_energies': [],
            'max_energies': [],
        }

        if exact_energy is not None:
            results['errors'] = []

        for n_shots in shot_counts:
            energies_at_n = []

            for trial in range(n_trials):
                E = energy_fn(n_shots)
                energies_at_n.append(E)

            energies_at_n = np.array(energies_at_n)

            results['mean_energies'].append(np.mean(energies_at_n))
            results['std_energies'].append(np.std(energies_at_n))
            results['min_energies'].append(np.min(energies_at_n))
            results['max_energies'].append(np.max(energies_at_n))

            if exact_energy is not None:
                error = np.abs(np.mean(energies_at_n) - exact_energy)
                results['errors'].append(error)

        # Convert to arrays
        for key in ['mean_energies', 'std_energies', 'min_energies', 'max_energies']:
            results[key] = np.array(results[key])

        if exact_energy is not None:
            results['errors'] = np.array(results['errors'])

        # Theoretical standard error: SE = σ/√n
        # Estimate σ from highest shot count
        sigma_est = results['std_energies'][-1] * np.sqrt(shot_counts[-1])
        results['theoretical_se'] = sigma_est / np.sqrt(np.array(shot_counts))

        return results

    def compare_runs(self,
                    energy_results: List[Dict[str, float]],
                    labels: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Statistical comparison of multiple VQE/QPE runs.

        Args:
            energy_results: List of result dictionaries with 'energy' key
            labels: Optional labels for each run

        Returns:
            comparison: Dictionary with:
                - energies: List of energies
                - mean: Overall mean
                - std: Overall standard deviation
                - sem: Standard error of mean
                - ci_95: 95% confidence interval
                - consistency: Measure of consistency (0-1)
                - outliers: Indices of outlier runs
        """
        energies = np.array([r['energy'] for r in energy_results])
        n_runs = len(energies)

        if labels is None:
            labels = [f"Run {i+1}" for i in range(n_runs)]

        mean_energy = np.mean(energies)
        std_energy = np.std(energies)
        sem = std_energy / np.sqrt(n_runs)  # Standard error of mean

        # 95% confidence interval using t-distribution
        t_critical = stats.t.ppf(0.975, df=n_runs-1)
        ci_95 = (mean_energy - t_critical * sem, mean_energy + t_critical * sem)

        # Detect outliers (3-sigma rule)
        outliers = np.where(np.abs(energies - mean_energy) > 3 * std_energy)[0]

        # Consistency measure: 1 - CV (coefficient of variation)
        cv = std_energy / np.abs(mean_energy) if mean_energy != 0 else np.inf
        consistency = max(0, 1 - cv)

        return {
            'energies': energies,
            'labels': labels,
            'mean': mean_energy,
            'std': std_energy,
            'sem': sem,
            'ci_95': ci_95,
            'consistency': consistency,
            'outliers': outliers.tolist(),
            'n_runs': n_runs,
            'coefficient_of_variation': cv
        }

    def estimate_required_shots(self,
                               target_precision: float,
                               hamiltonian_variance: float = 1.0) -> int:
        """
        Estimate number of shots required for target precision.

        Args:
            target_precision: Desired standard error (e.g., 0.01 Ha)
            hamiltonian_variance: Variance of Hamiltonian expectation

        Returns:
            n_shots: Required number of shots

        Formula:
            SE = σ/√n → n = (σ/SE)²
        """
        n_shots = int(np.ceil((hamiltonian_variance / target_precision)**2))

        logger.info(f"To achieve SE = {target_precision:.4f}, "
                    f"need n_shots = {n_shots}")

        return n_shots

    def statistical_test(self,
                        energy1: float,
                        energy2: float,
                        std1: float,
                        std2: float,
                        n1: int = 1,
                        n2: int = 1) -> Dict[str, Any]:
        """
        Test if two energies are statistically different.

        Args:
            energy1, energy2: Mean energies
            std1, std2: Standard deviations
            n1, n2: Number of measurements

        Returns:
            result: Dictionary with:
                - t_statistic: t-test statistic
                - p_value: Probability of observing difference by chance
                - significant: Whether difference is significant (p < 0.05)
                - conclusion: Text interpretation
        """
        # Welch's t-test (unequal variances)
        se1 = std1 / np.sqrt(n1)
        se2 = std2 / np.sqrt(n2)
        se_diff = np.sqrt(se1**2 + se2**2)

        t_stat = (energy1 - energy2) / se_diff if se_diff > 0 else np.inf

        # Degrees of freedom (Welch-Satterthwaite equation)
        df = ((se1**2 + se2**2)**2 /
              (se1**4 / (n1-1) + se2**4 / (n2-1))) if n1 > 1 and n2 > 1 else 1

        p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df=df))

        significant = p_value < 0.05

        if significant:
            conclusion = f"Energies are significantly different (p={p_value:.4f})"
        else:
            conclusion = f"No significant difference (p={p_value:.4f})"

        return {
            't_statistic': t_stat,
            'p_value': p_value,
            'degrees_of_freedom': df,
            'significant': significant,
            'conclusion': conclusion,
            'difference': energy1 - energy2,
            'se_difference': se_diff
        }

    def plot_convergence(self,
                        convergence_result: Dict[str, Any],
                        save_path: Optional[str] = None):
        """
        Plot convergence analysis results.

        Args:
            convergence_result: Output from convergence_analysis()
            save_path: Path to save figure (None = show)
        """
        import matplotlib.pyplot as plt

        shot_counts = convergence_result['shot_counts']
        mean_energies = convergence_result['mean_energies']
        std_energies = convergence_result['std_energies']

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Plot 1: Energy vs shots
        ax1.errorbar(shot_counts, mean_energies, yerr=std_energies,
                    marker='o', capsize=5, capthick=2, linewidth=2,
                    label='Measured')

        if 'errors' in convergence_result:
            ax2_twin = ax1.twinx()
            ax2_twin.semilogy(shot_counts, convergence_result['errors'],
                             'r--', marker='s', label='Error')
            ax2_twin.set_ylabel('|E - E_exact| (Ha)', color='r', fontsize=12)
            ax2_twin.tick_params(axis='y', labelcolor='r')

        ax1.set_xlabel('Number of Shots', fontsize=12)
        ax1.set_ylabel('Energy (Ha)', fontsize=12)
        ax1.set_title('Energy Convergence', fontsize=14)
        ax1.set_xscale('log')
        ax1.grid(alpha=0.3)
        ax1.legend()

        # Plot 2: Standard deviation vs shots (log-log)
        ax2.loglog(shot_counts, std_energies, 'bo-', linewidth=2,
                  markersize=8, label='Measured σ')

        # Theoretical 1/√n scaling
        if 'theoretical_se' in convergence_result:
            ax2.loglog(shot_counts, convergence_result['theoretical_se'],
                      'r--', linewidth=2, label='Theoretical σ/√n')

        ax2.set_xlabel('Number of Shots', fontsize=12)
        ax2.set_ylabel('Standard Deviation (Ha)', fontsize=12)
        ax2.set_title('Shot Noise Scaling', fontsize=14)
        ax2.grid(alpha=0.3)
        ax2.legend()

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Convergence plot saved to {save_path}")
        else:
            plt.show()

        plt.close()

    def __repr__(self) -> str:
        """String representation."""
        return f"UncertaintyAnalyzer(backend='{self.backend}')"
