"""
Temperature Modulator

Modifies molecular Hamiltonians based on temperature, accounting for:
- Thermal population of excited states
- Temperature-dependent bond strengths
- Vibrational zero-point and thermal energy
- Boltzmann-weighted configuration mixing

Physical Basis:
    At temperature T, molecular properties are averaged over thermal
    populations:  ⟨A⟩_T = Σ_i p_i A_i  where p_i ∝ exp(-E_i/kT)

References:
    - McQuarrie "Statistical Mechanics" (2000)
    - Atkins "Physical Chemistry" (2018)
"""

import numpy as np
from typing import Dict, Any, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class TemperatureModulator:
    """
    Apply temperature effects to molecular Hamiltonians.

    Temperature affects:
    1. Bond dissociation: Bonds weaken at high T (Morse potential)
    2. Configuration entropy: More configurations accessible at high T
    3. Vibrational energy: ZPE + thermal contributions
    4. Electronic excitations: Population of excited states

    Example:
        >>> from kanad.bonds import BondFactory
        >>> from kanad.environment import TemperatureModulator
        >>>
        >>> bond = BondFactory.create_bond('H', 'H', distance=0.74)
        >>> temp_mod = TemperatureModulator()
        >>>
        >>> # Room temperature
        >>> result_298 = temp_mod.apply_temperature(bond, 298.15)
        >>> print(f"Energy at 298K: {result_298['energy']:.6f} Ha")
        >>>
        >>> # High temperature
        >>> result_1000 = temp_mod.apply_temperature(bond, 1000.0)
        >>> print(f"Energy at 1000K: {result_1000['energy']:.6f} Ha")
        >>> print(f"Bond weakening: {result_1000['bond_strength_factor']:.3f}")
    """

    # Physical constants (CODATA 2018)
    k_B_Ha = 3.1668115634556076e-6  # Boltzmann constant in Ha/K
    k_B_eV = 8.617333262e-5          # Boltzmann constant in eV/K
    Ha_to_kcal = 627.509474          # Hartree to kcal/mol

    def __init__(self):
        """Initialize temperature modulator."""
        self.reference_temp = 298.15  # K (standard conditions)
        logger.info("TemperatureModulator initialized")

    def apply_temperature(
        self,
        bond_or_molecule,
        temperature: float,
        include_vibrational: bool = True,
        include_electronic: bool = True,
        n_excited_states: int = 5
    ) -> Dict[str, Any]:
        """
        Apply temperature effects to molecular system.

        Args:
            bond_or_molecule: Bond or Molecule object
            temperature: Temperature in Kelvin
            include_vibrational: Include vibrational thermal energy
            include_electronic: Include electronic excited state populations
            n_excited_states: Number of excited states for Boltzmann averaging

        Returns:
            Dictionary with:
                energy: Temperature-corrected energy (Ha)
                free_energy: Helmholtz free energy A = E - TS (Ha)
                entropy: Entropy S (Ha/K)
                heat_capacity: Cv (Ha/K)
                bond_strength_factor: Multiplicative factor for bond strength
                thermal_population: Population distribution over states
                vibrational_energy: Thermal vibrational contribution (Ha)
        """
        logger.info(f"Applying temperature effects: T = {temperature:.2f} K")

        # Get base energy and Hamiltonian
        base_energy = self._get_base_energy(bond_or_molecule)

        # 1. Compute thermal bond weakening
        bond_strength_factor = self._compute_bond_strength_factor(temperature)

        # 2. Vibrational thermal energy
        if include_vibrational:
            E_vib = self._compute_vibrational_energy(
                bond_or_molecule, temperature
            )
        else:
            E_vib = 0.0

        # 3. Electronic thermal population
        if include_electronic and hasattr(bond_or_molecule, 'compute_excited_states'):
            electronic_correction, populations = self._compute_electronic_thermal_correction(
                bond_or_molecule, temperature, n_excited_states
            )
        else:
            electronic_correction = 0.0
            populations = [1.0]  # Only ground state

        # 4. Total thermal energy
        E_thermal = base_energy + E_vib + electronic_correction

        # 5. Entropy (approximation from partition function)
        S = self._compute_entropy(
            temperature, E_vib, electronic_correction
        )

        # 6. Free energy: A = E - TS
        A = E_thermal - temperature * S

        # 7. Heat capacity
        Cv = self._compute_heat_capacity(temperature, E_vib)

        return {
            'energy': E_thermal,
            'free_energy': A,
            'entropy': S,
            'heat_capacity': Cv,
            'bond_strength_factor': bond_strength_factor,
            'thermal_population': populations,
            'vibrational_energy': E_vib,
            'electronic_correction': electronic_correction,
            'temperature': temperature
        }

    def scan_temperature(
        self,
        bond_or_molecule,
        temp_range: Tuple[float, float] = (100, 1000),
        n_points: int = 20,
        **kwargs
    ) -> Dict[str, np.ndarray]:
        """
        Scan temperature and compute properties at each point.

        Args:
            bond_or_molecule: Bond or Molecule object
            temp_range: (T_min, T_max) in Kelvin
            n_points: Number of temperature points
            **kwargs: Additional arguments for apply_temperature

        Returns:
            Dictionary with arrays:
                temperatures: Temperature values (K)
                energies: Total energy vs T (Ha)
                free_energies: Helmholtz free energy vs T (Ha)
                entropies: Entropy vs T (Ha/K)
                heat_capacities: Cv vs T (Ha/K)
                bond_strengths: Bond strength factor vs T
        """
        T_min, T_max = temp_range
        temperatures = np.linspace(T_min, T_max, n_points)

        energies = []
        free_energies = []
        entropies = []
        heat_capacities = []
        bond_strengths = []

        for T in temperatures:
            result = self.apply_temperature(bond_or_molecule, T, **kwargs)
            energies.append(result['energy'])
            free_energies.append(result['free_energy'])
            entropies.append(result['entropy'])
            heat_capacities.append(result['heat_capacity'])
            bond_strengths.append(result['bond_strength_factor'])

        return {
            'temperatures': temperatures,
            'energies': np.array(energies),
            'free_energies': np.array(free_energies),
            'entropies': np.array(entropies),
            'heat_capacities': np.array(heat_capacities),
            'bond_strengths': np.array(bond_strengths)
        }

    def _get_base_energy(self, bond_or_molecule) -> float:
        """Get base ground state energy."""
        if hasattr(bond_or_molecule, 'energy'):
            return bond_or_molecule.energy
        elif hasattr(bond_or_molecule, 'hamiltonian'):
            # Compute ground state energy
            # For bonds, use cached energy if available
            if hasattr(bond_or_molecule, '_cached_energy'):
                return bond_or_molecule._cached_energy
            else:
                # Would compute VQE or FCI here
                logger.warning("Computing energy on the fly - consider caching")
                return 0.0  # Placeholder
        else:
            raise ValueError("Cannot extract energy from object")

    def _compute_bond_strength_factor(self, temperature: float) -> float:
        """
        Compute temperature-dependent bond strength factor.

        Physical basis: Morse potential depth decreases with T
        Empirical: D(T) ≈ D(0) × exp(-αT/T_ref)

        Args:
            temperature: Temperature in K

        Returns:
            Factor to multiply bond dissociation energy (0-1)
        """
        # Empirical temperature scaling (typical: α ~ 0.0003)
        alpha = 0.0003  # Adjustable parameter
        T_ref = self.reference_temp

        factor = np.exp(-alpha * (temperature - T_ref))

        # Constrain to reasonable range
        factor = np.clip(factor, 0.5, 1.2)

        return factor

    def _compute_vibrational_energy(
        self,
        bond_or_molecule,
        temperature: float
    ) -> float:
        """
        Compute vibrational thermal energy.

        For harmonic oscillator:
            E_vib = Σ_i [½ℏω_i + ℏω_i/(exp(ℏω_i/kT) - 1)]
                  = ZPE + thermal excitation

        Args:
            bond_or_molecule: Molecular system
            temperature: Temperature in K

        Returns:
            Vibrational thermal energy in Ha
        """
        # Get vibrational frequencies (if available)
        if hasattr(bond_or_molecule, 'frequencies'):
            frequencies = bond_or_molecule.frequencies  # cm^-1
        elif hasattr(bond_or_molecule, 'get_frequencies'):
            frequencies = bond_or_molecule.get_frequencies()
        else:
            # Estimate from bond type (very crude)
            logger.warning("No vibrational data - using estimate")
            # Typical C-H stretch ~ 3000 cm^-1, C-C ~ 1000 cm^-1
            frequencies = np.array([1000.0])  # Placeholder

        # Convert cm^-1 to Ha
        cm_to_Ha = 4.556335e-6  # conversion factor
        omega_Ha = frequencies * cm_to_Ha

        # Compute thermal energy for each mode
        E_vib = 0.0
        for omega in omega_Ha:
            # Zero-point energy
            E_ZPE = 0.5 * omega

            # Thermal excitation
            x = omega / (self.k_B_Ha * temperature)
            if x < 50:  # Avoid overflow
                E_thermal = omega / (np.exp(x) - 1)
            else:
                E_thermal = 0.0  # Frozen out

            E_vib += E_ZPE + E_thermal

        return E_vib

    def _compute_electronic_thermal_correction(
        self,
        molecule,
        temperature: float,
        n_states: int = 5
    ) -> Tuple[float, np.ndarray]:
        """
        Compute thermal population correction from excited states.

        Boltzmann-weighted energy:
            E = Σ_i p_i E_i  where p_i = exp(-E_i/kT) / Z

        Args:
            molecule: Molecule with excited state data
            temperature: Temperature in K
            n_states: Number of states to include

        Returns:
            (energy_correction, populations)
        """
        # Get excited state energies
        if hasattr(molecule, 'excited_energies'):
            energies = molecule.excited_energies[:n_states]
        else:
            logger.warning("No excited state data - skipping thermal correction")
            return 0.0, np.array([1.0])

        # Compute Boltzmann weights
        beta = 1.0 / (self.k_B_Ha * temperature)
        E_0 = energies[0]  # Ground state

        # Relative energies
        Delta_E = energies - E_0

        # Boltzmann factors
        weights = np.exp(-beta * Delta_E)
        Z = np.sum(weights)  # Partition function
        populations = weights / Z

        # Thermal average energy
        E_avg = np.sum(populations * energies)

        # Correction = thermal average - ground state
        correction = E_avg - E_0

        return correction, populations

    def _compute_entropy(
        self,
        temperature: float,
        E_vib: float,
        E_elec: float
    ) -> float:
        """
        Estimate entropy from thermal energies.

        Approximation: S ≈ (E_thermal - E_0) / T + k ln(Z)

        Args:
            temperature: Temperature in K
            E_vib: Vibrational thermal energy (Ha)
            E_elec: Electronic thermal energy (Ha)

        Returns:
            Entropy in Ha/K
        """
        # Simple approximation: S ~ E_thermal/T
        # More accurate would require partition functions

        if temperature > 0:
            S = (E_vib + E_elec) / temperature + self.k_B_Ha * np.log(2)
        else:
            S = 0.0

        return S

    def _compute_heat_capacity(
        self,
        temperature: float,
        E_vib: float
    ) -> float:
        """
        Compute heat capacity Cv from vibrational energy.

        Cv = dE/dT (at constant volume)

        Args:
            temperature: Temperature in K
            E_vib: Vibrational energy (Ha)

        Returns:
            Heat capacity in Ha/K
        """
        # Numerical derivative
        dT = 1.0  # K
        T1 = temperature - dT/2
        T2 = temperature + dT/2

        # Very crude approximation
        # Would need to recompute E_vib at T1 and T2
        # For now, use classical limit: Cv ~ 3Nk

        Cv = 3 * self.k_B_Ha  # Per atom

        return Cv

    def modify_hamiltonian_with_temperature(
        self,
        hamiltonian,
        temperature: float,
        bond_strength_factor: Optional[float] = None
    ):
        """
        Create temperature-modified Hamiltonian.

        Modifies:
        - Bond dissociation terms (weaker at high T)
        - Adds thermal fluctuations (optional)

        Args:
            hamiltonian: Original Hamiltonian
            temperature: Temperature in K
            bond_strength_factor: Pre-computed factor (or compute if None)

        Returns:
            Modified Hamiltonian object
        """
        if bond_strength_factor is None:
            bond_strength_factor = self._compute_bond_strength_factor(temperature)

        # For SparsePauliOp, scale coupling terms
        # This is a simplified approach - in practice would need
        # to identify bond-specific terms

        logger.info(f"Applying bond strength factor: {bond_strength_factor:.4f}")

        # Scale Hamiltonian (affects all interactions equally)
        # More sophisticated: identify and scale only bond terms
        H_thermal = hamiltonian * bond_strength_factor

        return H_thermal

    def plot_temperature_scan(
        self,
        scan_result: Dict[str, np.ndarray],
        properties: list = ['energies', 'free_energies', 'entropies'],
        save_path: Optional[str] = None
    ):
        """
        Plot temperature-dependent properties.

        Args:
            scan_result: Output from scan_temperature()
            properties: List of properties to plot
            save_path: Optional path to save figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed - cannot plot")
            return

        fig, axes = plt.subplots(len(properties), 1, figsize=(10, 4*len(properties)))
        if len(properties) == 1:
            axes = [axes]

        T = scan_result['temperatures']

        for ax, prop in zip(axes, properties):
            values = scan_result[prop]

            ax.plot(T, values, 'o-', linewidth=2, markersize=6)
            ax.set_xlabel('Temperature (K)', fontsize=12)
            ax.set_ylabel(prop.replace('_', ' ').title(), fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.set_xlim(T[0], T[-1])

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")

        plt.show()
