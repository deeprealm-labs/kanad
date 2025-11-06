"""
Molecular Dynamics Thermostats

Temperature control methods for MD simulations to maintain canonical (NVT) ensemble.

Without a thermostat, MD samples the microcanonical (NVE) ensemble where total
energy is conserved. Thermostats couple the system to a heat bath to maintain
constant temperature.

Implemented Thermostats:
-----------------------
1. Berendsen: Simple velocity rescaling (weak coupling)
2. Velocity Rescaling: Improved Berendsen with correct statistics
3. Nose-Hoover: Deterministic canonical ensemble
4. Langevin: Stochastic dynamics with friction and random force

Choosing a Thermostat:
---------------------
- **Berendsen**: Fast equilibration, use for initial relaxation
- **Velocity Rescaling**: Better than Berendsen, correct temperature distribution
- **Nose-Hoover**: True canonical ensemble, use for production runs
- **Langevin**: Implicit solvent, includes friction/diffusion

References:
----------
- Berendsen: Berendsen et al. (1984) J. Chem. Phys. 81, 3684
- Velocity Rescaling: Bussi et al. (2007) J. Chem. Phys. 126, 014101
- Nose-Hoover: Hoover (1985) Phys. Rev. A 31, 1695
- Langevin: Lemons & Gythiel (1997) Am. J. Phys. 65, 1079
"""

import numpy as np
import logging
from typing import Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# Physical constants
K_BOLTZMANN = 3.1668105e-6  # Hartree/K


@dataclass
class ThermostatState:
    """
    Extended state for thermostats that require auxiliary variables.

    Attributes:
        xi: Nose-Hoover thermostat variable
        v_xi: Nose-Hoover thermostat velocity
    """
    xi: float = 0.0
    v_xi: float = 0.0


class BaseThermostat:
    """
    Base class for molecular dynamics thermostats.

    All thermostats modify velocities to maintain target temperature.
    """

    def __init__(self, target_temperature: float):
        """
        Initialize thermostat.

        Args:
            target_temperature: Target temperature in Kelvin
        """
        self.target_temperature = target_temperature  # K
        self.target_ke = None  # Will be computed based on degrees of freedom

        logger.debug(f"Initialized {self.__class__.__name__}")
        logger.debug(f"  Target temperature: {target_temperature:.2f} K")

    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        timestep: float,
        **kwargs
    ) -> np.ndarray:
        """
        Apply thermostat to modify velocities.

        Args:
            velocities: Atomic velocities (N_atoms, 3) in Bohr/fs
            masses: Atomic masses (N_atoms,) in amu
            timestep: Integration timestep in fs
            **kwargs: Additional parameters (e.g., current_ke, forces)

        Returns:
            Modified velocities (N_atoms, 3) in Bohr/fs
        """
        raise NotImplementedError("Subclasses must implement apply()")

    def compute_temperature(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        n_dof: Optional[int] = None
    ) -> float:
        """
        Compute instantaneous temperature from velocities.

        T = (2 * KE) / (k_B * N_dof)

        Args:
            velocities: (N_atoms, 3) in Bohr/fs
            masses: (N_atoms,) in amu
            n_dof: Number of degrees of freedom (default: 3*N_atoms - 6)

        Returns:
            Temperature in Kelvin
        """
        n_atoms = len(masses)
        if n_dof is None:
            # 3N - 6 for non-linear molecules (remove translation + rotation)
            n_dof = 3 * n_atoms - 6
            if n_atoms == 1:
                n_dof = 3  # Monoatomic
            elif n_atoms == 2:
                n_dof = 3 * n_atoms - 5  # Diatomic (1 rotation is degenerate)

        # Convert masses to electron masses
        masses_me = masses * 1822.888486

        # Convert velocities to a.u.
        velocities_au = velocities * 41.341

        # Kinetic energy
        ke = 0.5 * np.sum(masses_me[:, np.newaxis] * velocities_au**2)

        # Temperature: T = 2*KE / (k_B * N_dof)
        temperature = (2.0 * ke) / (K_BOLTZMANN * n_dof)

        return temperature


class BerendsenThermostat(BaseThermostat):
    """
    Berendsen weak-coupling thermostat.

    Algorithm:
    ---------
    Rescale velocities at each step:
        v' = λ * v
    where:
        λ = sqrt(1 + (dt/τ) * (T_target/T_current - 1))

    Properties:
    ----------
    - Simple and efficient
    - Fast equilibration
    - Does NOT sample canonical ensemble correctly
    - Use for initial equilibration, not production
    - Coupling time τ controls relaxation rate

    The Berendsen thermostat drives the system towards target temperature
    exponentially with time constant τ. Smaller τ = faster but less physical.

    Recommended τ: 0.1-1.0 ps for equilibration

    References:
    ----------
    - Berendsen, Postma, van Gunsteren, DiNola, Haak (1984)
      J. Chem. Phys. 81, 3684
    """

    def __init__(self, target_temperature: float, coupling_time: float = 0.5):
        """
        Initialize Berendsen thermostat.

        Args:
            target_temperature: Target temperature in K
            coupling_time: Coupling time constant τ in picoseconds
        """
        super().__init__(target_temperature)
        self.coupling_time = coupling_time * 1000.0  # Convert ps to fs
        logger.debug(f"  Coupling time: {coupling_time:.3f} ps")

    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        timestep: float,
        **kwargs
    ) -> np.ndarray:
        """
        Apply Berendsen weak coupling.

        Args:
            velocities: Current velocities
            masses: Atomic masses
            timestep: Integration timestep in fs

        Returns:
            Rescaled velocities
        """
        # Compute current temperature
        T_current = self.compute_temperature(velocities, masses)

        # Berendsen scaling factor
        # λ = sqrt(1 + (dt/τ) * (T_target/T_current - 1))
        ratio = self.target_temperature / T_current
        lambda_factor = np.sqrt(1.0 + (timestep / self.coupling_time) * (ratio - 1.0))

        # Rescale velocities
        v_new = lambda_factor * velocities

        logger.debug(f"Berendsen: T={T_current:.2f}K → {self.target_temperature:.2f}K, λ={lambda_factor:.4f}")

        return v_new


class VelocityRescaling(BaseThermostat):
    """
    Velocity rescaling thermostat with correct canonical distribution.

    Algorithm:
    ---------
    Improved version of Berendsen that samples the canonical ensemble correctly
    by adding stochastic noise to the rescaling factor.

    Properties:
    ----------
    - Samples correct canonical ensemble
    - Better than Berendsen for production runs
    - Fast equilibration like Berendsen
    - Minimal overhead vs. Berendsen

    This is essentially Berendsen + noise to fix the statistical ensemble.
    Recommended for most MD simulations.

    References:
    ----------
    - Bussi, Donadio, Parrinello (2007) J. Chem. Phys. 126, 014101
    - Bussi, Parrinello (2007) Phys. Rev. E 75, 056707
    """

    def __init__(self, target_temperature: float, coupling_time: float = 0.5):
        """
        Initialize velocity rescaling thermostat.

        Args:
            target_temperature: Target temperature in K
            coupling_time: Coupling time constant τ in ps
        """
        super().__init__(target_temperature)
        self.coupling_time = coupling_time * 1000.0  # ps to fs
        logger.debug(f"  Coupling time: {coupling_time:.3f} ps")

    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        timestep: float,
        **kwargs
    ) -> np.ndarray:
        """
        Apply stochastic velocity rescaling.

        Args:
            velocities: Current velocities
            masses: Atomic masses
            timestep: Integration timestep in fs

        Returns:
            Rescaled velocities with stochastic correction
        """
        n_atoms = len(masses)
        n_dof = 3 * n_atoms - 6

        # Current kinetic energy and temperature
        T_current = self.compute_temperature(velocities, masses, n_dof=n_dof)

        # Target kinetic energy
        KE_target = 0.5 * K_BOLTZMANN * self.target_temperature * n_dof

        # Current kinetic energy (in Hartree)
        masses_me = masses * 1822.888486
        velocities_au = velocities * 41.341
        KE_current = 0.5 * np.sum(masses_me[:, np.newaxis] * velocities_au**2)

        # Stochastic rescaling (Bussi et al. 2007)
        c = timestep / self.coupling_time
        R1 = np.random.normal(0, 1)
        R2_sum = np.sum(np.random.normal(0, 1, n_dof - 1)**2)

        # Stochastic target kinetic energy
        KE_stochastic = KE_target + (KE_current - KE_target) * np.exp(-c) + \
                        KE_target * (1 - np.exp(-c)) * (R1**2 / n_dof + R2_sum / n_dof) + \
                        2 * np.sqrt(KE_current * KE_target * (1 - np.exp(-c)) / n_dof) * R1

        # Scaling factor
        if KE_current > 0:
            alpha = np.sqrt(KE_stochastic / KE_current)
        else:
            alpha = 1.0

        # Rescale velocities
        v_new = alpha * velocities

        T_new = self.compute_temperature(v_new, masses, n_dof=n_dof)
        logger.debug(f"Velocity Rescaling: T={T_current:.2f}K → {T_new:.2f}K, α={alpha:.4f}")

        return v_new


class NoseHooverThermostat(BaseThermostat):
    """
    Nose-Hoover thermostat for canonical ensemble.

    Algorithm:
    ---------
    Extended Hamiltonian with thermostat variable ξ:
        dv/dt = F/m - ξ * v
        dξ/dt = (T_current - T_target) / Q

    where Q is the thermostat mass (controls oscillation frequency).

    Properties:
    ----------
    - Deterministic (no random forces)
    - Rigorously samples canonical ensemble
    - Requires careful choice of Q
    - Can show temperature oscillations
    - Use Nose-Hoover chains for better ergodicity

    The thermostat mass Q determines the coupling strength. Too small Q leads
    to fast oscillations, too large Q leads to weak coupling.

    Recommended Q ~ k_B * T * τ² where τ is desired relaxation time.

    References:
    ----------
    - Nose (1984) Mol. Phys. 52, 255
    - Hoover (1985) Phys. Rev. A 31, 1695
    - Martyna, Klein, Tuckerman (1992) J. Chem. Phys. 97, 2635 (chains)
    """

    def __init__(
        self,
        target_temperature: float,
        relaxation_time: float = 0.5,
        n_chains: int = 1
    ):
        """
        Initialize Nose-Hoover thermostat.

        Args:
            target_temperature: Target temperature in K
            relaxation_time: Characteristic time for temperature relaxation in ps
            n_chains: Number of Nose-Hoover chain thermostats (1=simple, 3-5=better)
        """
        super().__init__(target_temperature)
        self.relaxation_time = relaxation_time * 1000.0  # ps to fs
        self.n_chains = n_chains

        # Thermostat mass (in a.u. squared)
        # Q = k_B * T * τ²
        self.Q = K_BOLTZMANN * target_temperature * (self.relaxation_time * 41.341)**2

        # Initialize thermostat variables
        self.xi = np.zeros(n_chains)  # Thermostat positions
        self.v_xi = np.zeros(n_chains)  # Thermostat velocities

        logger.debug(f"  Relaxation time: {relaxation_time:.3f} ps")
        logger.debug(f"  Thermostat mass Q: {self.Q:.3e} a.u.²")
        logger.debug(f"  Number of chains: {n_chains}")

    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        timestep: float,
        **kwargs
    ) -> np.ndarray:
        """
        Apply Nose-Hoover thermostat.

        Args:
            velocities: Current velocities
            masses: Atomic masses
            timestep: Integration timestep in fs

        Returns:
            Modified velocities
        """
        n_atoms = len(masses)
        n_dof = 3 * n_atoms - 6

        # Current kinetic energy
        masses_me = masses * 1822.888486
        velocities_au = velocities * 41.341
        KE_current = 0.5 * np.sum(masses_me[:, np.newaxis] * velocities_au**2)

        # Target kinetic energy
        KE_target = 0.5 * K_BOLTZMANN * self.target_temperature * n_dof

        dt_au = timestep * 41.341

        # Update thermostat (first chain only for simplicity)
        # dξ/dt = (KE_current - KE_target) / Q
        self.v_xi[0] += 0.5 * dt_au * (KE_current - KE_target) / self.Q
        self.xi[0] += dt_au * self.v_xi[0]

        # Scale velocities
        # dv/dt = -ξ * v
        scale_factor = np.exp(-self.v_xi[0] * dt_au)
        v_new = scale_factor * velocities

        # Update thermostat again (leap-frog style)
        velocities_au_new = v_new * 41.341
        KE_new = 0.5 * np.sum(masses_me[:, np.newaxis] * velocities_au_new**2)
        self.v_xi[0] += 0.5 * dt_au * (KE_new - KE_target) / self.Q

        T_current = (2.0 * KE_current) / (K_BOLTZMANN * n_dof)
        T_new = (2.0 * KE_new) / (K_BOLTZMANN * n_dof)
        logger.debug(f"Nose-Hoover: T={T_current:.2f}K → {T_new:.2f}K, ξ={self.xi[0]:.4f}")

        return v_new


class LangevinThermostat(BaseThermostat):
    """
    Langevin thermostat - stochastic dynamics with friction.

    Algorithm:
    ---------
    Modified equation of motion:
        m * dv/dt = F - γ * m * v + R(t)

    where:
        γ = friction coefficient
        R(t) = random force (Gaussian white noise)
        ⟨R_i(t) R_j(t')⟩ = 2 * γ * m * k_B * T * δ_ij * δ(t-t')

    Properties:
    ----------
    - Stochastic (includes random forces)
    - Implicit solvent (friction mimics solvent)
    - Correct canonical ensemble
    - Can compute diffusion coefficients
    - Friction coefficient γ determines coupling

    The Langevin equation describes Brownian motion and is widely used for
    modeling implicit solvent effects. The friction γ should match the
    solvent viscosity.

    Recommended γ: 0.1-1.0 ps⁻¹ for water-like solvent

    References:
    ----------
    - Langevin (1908) C. R. Acad. Sci. Paris 146, 530
    - Lemons & Gythiel (1997) Am. J. Phys. 65, 1079
    - Bussi & Parrinello (2008) Comput. Phys. Commun. 179, 26
    """

    def __init__(
        self,
        target_temperature: float,
        friction_coefficient: float = 1.0
    ):
        """
        Initialize Langevin thermostat.

        Args:
            target_temperature: Target temperature in K
            friction_coefficient: Friction γ in ps⁻¹
        """
        super().__init__(target_temperature)
        self.gamma = friction_coefficient / 1000.0  # ps⁻¹ to fs⁻¹
        logger.debug(f"  Friction γ: {friction_coefficient:.3f} ps⁻¹")

    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        timestep: float,
        **kwargs
    ) -> np.ndarray:
        """
        Apply Langevin dynamics.

        Args:
            velocities: Current velocities
            masses: Atomic masses
            timestep: Integration timestep in fs

        Returns:
            Modified velocities with friction and random force
        """
        # Convert to atomic units
        masses_me = masses * 1822.888486
        dt_au = timestep * 41.341
        gamma_au = self.gamma / 41.341  # fs⁻¹ to a.u.⁻¹

        # Friction term: exp(-γ * dt)
        friction_factor = np.exp(-gamma_au * dt_au)

        # Random force term
        # σ² = k_B * T * (1 - exp(-2*γ*dt)) / m
        sigma_sq = K_BOLTZMANN * self.target_temperature * (1 - friction_factor**2) / masses_me

        # Convert velocities to a.u.
        v_au = velocities * 41.341

        # Apply Langevin update
        # v(t+dt) = v(t) * exp(-γ*dt) + sqrt(σ²) * N(0,1)
        random_force = np.random.normal(0, 1, v_au.shape)
        v_au_new = friction_factor * v_au + np.sqrt(sigma_sq)[:, np.newaxis] * random_force

        # Convert back to Bohr/fs
        v_new = v_au_new / 41.341

        T_old = self.compute_temperature(velocities, masses)
        T_new = self.compute_temperature(v_new, masses)
        logger.debug(f"Langevin: T={T_old:.2f}K → {T_new:.2f}K, γ={self.gamma*1000:.3f} ps⁻¹")

        return v_new


# Factory function
def create_thermostat(
    name: str,
    target_temperature: float,
    **kwargs
) -> BaseThermostat:
    """
    Factory function to create thermostats by name.

    Args:
        name: Thermostat name ('berendsen', 'velocity_rescaling', 'nose_hoover', 'langevin')
        target_temperature: Target temperature in K
        **kwargs: Additional parameters (coupling_time, friction_coefficient, etc.)

    Returns:
        Thermostat instance

    Raises:
        ValueError: If thermostat name not recognized
    """
    thermostats = {
        'berendsen': BerendsenThermostat,
        'velocity_rescaling': VelocityRescaling,
        'v_rescale': VelocityRescaling,  # Alias
        'nose_hoover': NoseHooverThermostat,
        'nh': NoseHooverThermostat,  # Alias
        'langevin': LangevinThermostat,
    }

    name_lower = name.lower()
    if name_lower not in thermostats:
        available = ', '.join(thermostats.keys())
        raise ValueError(f"Unknown thermostat '{name}'. Available: {available}")

    return thermostats[name_lower](target_temperature, **kwargs)
