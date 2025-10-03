"""
Temperature effects in quantum chemistry and materials.

Implements:
- Fermi-Dirac distribution for electron occupation
- Thermal averaging of observables
- Alloy formation and mixing entropy
- Temperature-dependent properties
"""

from typing import Dict, List, Optional, Tuple, Any
import numpy as np
from scipy.constants import k as k_B_SI  # Boltzmann constant in J/K


class Temperature:
    """
    Temperature effects in quantum systems.

    Handles:
        - Fermi-Dirac statistics for electrons
        - Thermal occupation of states
        - Free energy and entropy
        - Temperature-dependent mixing (alloys)

    Units:
        - Temperature: Kelvin (K)
        - Energy: eV
        - k_B = 8.617333×10⁻⁵ eV/K
    """

    # Boltzmann constant in eV/K
    k_B = 8.617333e-5  # eV/K

    def __init__(self, T: float):
        """
        Initialize temperature.

        Args:
            T: Temperature in Kelvin
        """
        if T < 0:
            raise ValueError(f"Temperature must be non-negative, got {T} K")

        self.T = T  # Kelvin
        self.beta = 1.0 / (self.k_B * T) if T > 0 else np.inf  # 1/eV

    @classmethod
    def from_celsius(cls, T_celsius: float) -> 'Temperature':
        """Create from Celsius temperature."""
        return cls(T_celsius + 273.15)

    @classmethod
    def room_temperature(cls) -> 'Temperature':
        """Create at room temperature (298 K, ~25°C)."""
        return cls(298.0)

    @classmethod
    def zero(cls) -> 'Temperature':
        """Create at absolute zero (T = 0 K)."""
        return cls(0.0)

    def fermi_dirac(self, energy: float, mu: float) -> float:
        """
        Fermi-Dirac distribution.

        f(E) = 1 / (exp((E - μ) / k_B T) + 1)

        Args:
            energy: Energy level (eV)
            mu: Chemical potential / Fermi energy (eV)

        Returns:
            Occupation probability [0, 1]
        """
        if self.T == 0:
            # At T=0: step function
            return 1.0 if energy <= mu else 0.0

        x = self.beta * (energy - mu)

        # Avoid overflow for large x
        if x > 50:
            return 0.0
        elif x < -50:
            return 1.0
        else:
            return 1.0 / (np.exp(x) + 1.0)

    def thermal_occupation(
        self,
        energies: np.ndarray,
        fermi_energy: float
    ) -> np.ndarray:
        """
        Compute thermal occupation of energy levels.

        Args:
            energies: Array of energy levels (eV)
            fermi_energy: Fermi energy (eV)

        Returns:
            Occupation numbers for each level
        """
        if self.T == 0:
            # At T=0: fill up to Fermi level
            return np.where(energies <= fermi_energy, 1.0, 0.0)

        # Vectorized Fermi-Dirac
        x = self.beta * (energies - fermi_energy)

        # Clip to avoid overflow
        x_clipped = np.clip(x, -50, 50)
        occupation = 1.0 / (np.exp(x_clipped) + 1.0)

        return occupation

    def thermal_energy(
        self,
        energies: np.ndarray,
        fermi_energy: float,
        degeneracy: int = 2
    ) -> float:
        """
        Compute total thermal energy.

        E_thermal = Σ_n f(E_n) * E_n * g

        Args:
            energies: Energy levels (eV)
            fermi_energy: Fermi energy (eV)
            degeneracy: Spin degeneracy (2 for electrons)

        Returns:
            Thermal energy (eV)
        """
        occupation = self.thermal_occupation(energies, fermi_energy)
        return degeneracy * np.sum(occupation * energies)

    def entropy(
        self,
        energies: np.ndarray,
        fermi_energy: float,
        degeneracy: int = 2
    ) -> float:
        """
        Compute electronic entropy.

        S = -k_B Σ_n [f ln(f) + (1-f) ln(1-f)] * g

        Args:
            energies: Energy levels (eV)
            fermi_energy: Fermi energy (eV)
            degeneracy: Spin degeneracy

        Returns:
            Entropy in eV/K
        """
        if self.T == 0:
            return 0.0

        occupation = self.thermal_occupation(energies, fermi_energy)

        # Avoid log(0)
        f_safe = np.clip(occupation, 1e-50, 1.0 - 1e-50)
        one_minus_f = 1.0 - f_safe

        entropy_per_level = -(f_safe * np.log(f_safe) + one_minus_f * np.log(one_minus_f))
        total_entropy = degeneracy * self.k_B * np.sum(entropy_per_level)

        return total_entropy

    def free_energy(
        self,
        energies: np.ndarray,
        fermi_energy: float,
        degeneracy: int = 2
    ) -> float:
        """
        Compute Helmholtz free energy.

        F = E - TS

        Args:
            energies: Energy levels (eV)
            fermi_energy: Fermi energy (eV)
            degeneracy: Spin degeneracy

        Returns:
            Free energy (eV)
        """
        E = self.thermal_energy(energies, fermi_energy, degeneracy)
        S = self.entropy(energies, fermi_energy, degeneracy)
        F = E - self.T * S

        return F

    def boltzmann_factor(self, energy: float) -> float:
        """
        Boltzmann factor exp(-E / k_B T).

        Args:
            energy: Energy (eV)

        Returns:
            Boltzmann factor
        """
        if self.T == 0:
            return 0.0 if energy > 0 else 1.0

        x = -self.beta * energy
        x_clipped = np.clip(x, -50, 50)
        return np.exp(x_clipped)

    def maxwell_boltzmann(self, energy: float, partition_function: float) -> float:
        """
        Maxwell-Boltzmann distribution.

        P(E) = exp(-E / k_B T) / Z

        Args:
            energy: Energy (eV)
            partition_function: Z = Σ exp(-E_i / k_B T)

        Returns:
            Probability
        """
        if partition_function <= 0:
            raise ValueError("Partition function must be positive")

        return self.boltzmann_factor(energy) / partition_function

    def thermal_average(
        self,
        values: np.ndarray,
        energies: np.ndarray
    ) -> float:
        """
        Compute thermal average of observable.

        ⟨A⟩ = Σ A_i exp(-E_i / k_B T) / Z

        Args:
            values: Observable values
            energies: Corresponding energies (eV)

        Returns:
            Thermal average
        """
        if len(values) != len(energies):
            raise ValueError("values and energies must have same length")

        weights = np.array([self.boltzmann_factor(E) for E in energies])
        Z = np.sum(weights)

        if Z == 0:
            return 0.0

        return np.sum(values * weights) / Z

    def __repr__(self) -> str:
        """String representation."""
        return f"Temperature({self.T:.2f} K)"

    def __str__(self) -> str:
        """Human-readable string."""
        celsius = self.T - 273.15
        return f"{self.T:.2f} K ({celsius:.2f} °C)"


class AlloyFormation:
    """
    Model alloy formation with temperature-dependent mixing.

    Handles:
        - Binary and multi-component alloys
        - Mixing enthalpy and entropy
        - Phase diagrams (simplified)
        - Solid solution vs phase separation
    """

    def __init__(
        self,
        elements: List[str],
        compositions: List[float],
        temperature: Temperature
    ):
        """
        Initialize alloy.

        Args:
            elements: Element symbols (e.g., ['Cu', 'Zn'])
            compositions: Mole fractions (must sum to 1)
            temperature: Temperature object
        """
        if len(elements) != len(compositions):
            raise ValueError("elements and compositions must have same length")

        if abs(sum(compositions) - 1.0) > 1e-6:
            raise ValueError(f"Compositions must sum to 1, got {sum(compositions)}")

        self.elements = elements
        self.compositions = np.array(compositions)
        self.temperature = temperature

    def mixing_entropy(self) -> float:
        """
        Compute ideal mixing entropy.

        S_mix = -R Σ x_i ln(x_i)

        This favors mixing at high temperature.

        Returns:
            Mixing entropy (eV/K)
        """
        # Avoid log(0)
        x_safe = np.clip(self.compositions, 1e-50, 1.0)
        S_mix = -Temperature.k_B * np.sum(x_safe * np.log(x_safe))

        return S_mix

    def mixing_enthalpy(self, interaction_parameters: Dict[Tuple[str, str], float]) -> float:
        """
        Compute mixing enthalpy using regular solution model.

        H_mix = Σ_ij Ω_ij x_i x_j

        Ω_ij > 0: phase separation (unfavorable mixing)
        Ω_ij < 0: mixing favored

        Args:
            interaction_parameters: Dict of (element_i, element_j) -> Ω_ij (eV)

        Returns:
            Mixing enthalpy (eV)
        """
        H_mix = 0.0

        for i, elem_i in enumerate(self.elements):
            for j, elem_j in enumerate(self.elements):
                if i < j:  # Avoid double counting
                    key = (elem_i, elem_j)
                    key_rev = (elem_j, elem_i)

                    # Try both orderings
                    Omega = interaction_parameters.get(key, interaction_parameters.get(key_rev, 0.0))

                    H_mix += Omega * self.compositions[i] * self.compositions[j]

        return H_mix

    def mixing_free_energy(self, interaction_parameters: Dict[Tuple[str, str], float]) -> float:
        """
        Compute Gibbs free energy of mixing.

        ΔG_mix = H_mix - T S_mix

        ΔG_mix < 0: mixing favorable (solid solution)
        ΔG_mix > 0: phase separation

        Args:
            interaction_parameters: Interaction parameters

        Returns:
            Mixing free energy (eV/atom)
        """
        H_mix = self.mixing_enthalpy(interaction_parameters)
        S_mix = self.mixing_entropy()

        G_mix = H_mix - self.temperature.T * S_mix

        return G_mix

    def will_mix(self, interaction_parameters: Dict[Tuple[str, str], float]) -> bool:
        """
        Predict if alloy will form solid solution.

        Args:
            interaction_parameters: Interaction parameters

        Returns:
            True if mixing is thermodynamically favorable
        """
        G_mix = self.mixing_free_energy(interaction_parameters)
        return G_mix < 0.0

    # API compatibility aliases
    def get_mixing_entropy(self) -> float:
        """Alias for mixing_entropy() for API compatibility."""
        return self.mixing_entropy()

    def get_mixing_enthalpy(self, interaction_parameters: Dict[Tuple[str, str], float]) -> float:
        """Alias for mixing_enthalpy() for API compatibility."""
        return self.mixing_enthalpy(interaction_parameters)

    def get_free_energy_of_mixing(self, interaction_parameters: Dict[Tuple[str, str], float]) -> float:
        """Alias for mixing_free_energy() for API compatibility."""
        return self.mixing_free_energy(interaction_parameters)

    def is_alloy_stable(self, interaction_parameters: Dict[Tuple[str, str], float]) -> bool:
        """Alias for will_mix() for API compatibility."""
        return self.will_mix(interaction_parameters)

    def spinodal_temperature(
        self,
        composition: float,
        interaction_parameter: float
    ) -> float:
        """
        Estimate spinodal decomposition temperature (binary alloy).

        T_s = Ω / (2 k_B x (1-x))

        Above T_s: single phase
        Below T_s: phase separation

        Args:
            composition: Composition of first element (0-1)
            interaction_parameter: Ω (eV)

        Returns:
            Spinodal temperature (K)
        """
        if composition <= 0 or composition >= 1:
            return 0.0

        T_s = interaction_parameter / (2 * Temperature.k_B * composition * (1 - composition))

        return max(T_s, 0.0)

    def phase_diagram_point(self) -> Dict[str, Any]:
        """
        Get current point on phase diagram.

        Returns:
            Dictionary with phase information
        """
        return {
            'elements': self.elements,
            'compositions': self.compositions.tolist(),
            'temperature': self.temperature.T,
            'mixing_entropy': self.mixing_entropy(),
        }

    def __repr__(self) -> str:
        """String representation."""
        comp_str = '-'.join([f"{e}{x:.2f}" for e, x in zip(self.elements, self.compositions)])
        return f"Alloy({comp_str}, T={self.temperature.T:.1f} K)"
