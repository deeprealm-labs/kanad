"""
MD Initial Conditions and Equilibration

Provides methods for generating initial conditions for MD simulations:
- Maxwell-Boltzmann velocity distribution
- Center-of-mass motion removal
- System equilibration protocols

Proper initialization is critical for correct MD sampling. The velocity
distribution must match the target temperature and respect conservation laws
(total momentum = 0, total angular momentum = 0).

References:
----------
- Maxwell-Boltzmann: Maxwell (1860) Phil. Mag. 19, 19-32
- Equipartition theorem: k_B T = m ⟨v²⟩ / 2 per degree of freedom
- Allen & Tildesley (2017) Computer Simulation of Liquids, Chapter 3
"""

import numpy as np
import logging
from typing import Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# Physical constants
K_BOLTZMANN = 3.1668105e-6  # Hartree/K


@dataclass
class InitializationResult:
    """
    Result from initialization procedure.

    Attributes:
        positions: Initial positions (N_atoms, 3) in Bohr
        velocities: Initial velocities (N_atoms, 3) in Bohr/fs
        temperature: Actual temperature from velocities in K
        kinetic_energy: Kinetic energy in Hartree
        com_velocity: Center-of-mass velocity (should be ~0)
        angular_momentum: Total angular momentum (should be ~0)
    """
    positions: np.ndarray
    velocities: np.ndarray
    temperature: float
    kinetic_energy: float
    com_velocity: np.ndarray
    angular_momentum: np.ndarray


class MaxwellBoltzmannInitializer:
    """
    Generate velocities from Maxwell-Boltzmann distribution.

    For a system at temperature T, the velocity of each atom along each
    Cartesian direction is drawn from a Gaussian distribution:

        P(v_i) ∝ exp(-m_i v_i² / 2k_B T)

    The width of the distribution is:
        σ_i = sqrt(k_B T / m_i)

    This ensures the equipartition theorem: ⟨K.E.⟩ = (3/2) N k_B T
    """

    def __init__(self, temperature: float, remove_com: bool = True, seed: Optional[int] = None):
        """
        Initialize Maxwell-Boltzmann velocity generator.

        Args:
            temperature: Target temperature in Kelvin
            remove_com: Remove center-of-mass motion (recommended: True)
            seed: Random seed for reproducibility (None = random)
        """
        self.temperature = temperature
        self.remove_com = remove_com
        self.seed = seed

        if seed is not None:
            np.random.seed(seed)

        logger.debug(f"Initialized Maxwell-Boltzmann generator")
        logger.debug(f"  Target temperature: {temperature:.2f} K")
        logger.debug(f"  Remove COM motion: {remove_com}")

    def generate(
        self,
        masses: np.ndarray,
        positions: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Generate velocities from Maxwell-Boltzmann distribution.

        Args:
            masses: Atomic masses (N_atoms,) in amu
            positions: Atomic positions (N_atoms, 3) in Bohr (optional, for COM removal)

        Returns:
            Velocities (N_atoms, 3) in Bohr/fs
        """
        n_atoms = len(masses)

        # Convert masses to electron masses
        masses_me = masses * 1822.888486

        # Generate random velocities from Gaussian distribution
        # σ² = k_B T / m (in atomic units)
        sigma_au = np.sqrt(K_BOLTZMANN * self.temperature / masses_me)

        # Generate velocities in atomic units
        velocities_au = np.random.normal(0, 1, (n_atoms, 3)) * sigma_au[:, np.newaxis]

        # Convert from a.u. to Bohr/fs
        # 1 a.u. velocity = 1 Bohr / (1 a.u. time)
        # 1 a.u. time = 2.4189 × 10^-17 s = 0.024189 fs
        # So 1 a.u. velocity = 1 Bohr / 0.024189 fs = 41.341 Bohr/fs
        velocities = velocities_au / 41.341  # Bohr/fs

        # Remove center-of-mass motion
        if self.remove_com:
            velocities = remove_com_motion(velocities, masses, positions)

        # Rescale to exact target temperature
        velocities = self._rescale_to_temperature(velocities, masses, self.temperature)

        # Compute actual temperature
        T_actual = self._compute_temperature(velocities, masses)

        logger.debug(f"Generated MB velocities: T_target={self.temperature:.2f}K, T_actual={T_actual:.2f}K")

        return velocities

    def _compute_temperature(self, velocities: np.ndarray, masses: np.ndarray) -> float:
        """Compute temperature from velocities."""
        n_atoms = len(masses)
        n_dof = 3 * n_atoms - 6  # Remove translation + rotation

        # For small molecules (diatomic), n_dof can be 0
        # Use at least 1 DOF to avoid division by zero
        if n_dof <= 0:
            n_dof = 1

        masses_me = masses * 1822.888486
        velocities_au = velocities / 41.341  # v_au = v_fs / 41.341

        # Kinetic energy
        ke = 0.5 * np.sum(masses_me[:, np.newaxis] * velocities_au**2)

        # Temperature: T = 2*KE / (k_B * N_dof)
        temperature = (2.0 * ke) / (K_BOLTZMANN * n_dof)

        return temperature

    def _rescale_to_temperature(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        target_temperature: float
    ) -> np.ndarray:
        """
        Rescale velocities to exact target temperature.

        This ensures the initial temperature exactly matches the target,
        which is important for NVE simulations.
        """
        current_temperature = self._compute_temperature(velocities, masses)

        if current_temperature > 0:
            scale_factor = np.sqrt(target_temperature / current_temperature)
            velocities_scaled = scale_factor * velocities
        else:
            velocities_scaled = velocities

        return velocities_scaled


def remove_com_motion(
    velocities: np.ndarray,
    masses: np.ndarray,
    positions: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Remove center-of-mass translation and rotation.

    This ensures:
    1. Total momentum = 0 (no COM translation)
    2. Total angular momentum = 0 (no COM rotation)

    Conservation of momentum and angular momentum requires removing these
    unphysical motions that can arise from initial condition generation.

    Args:
        velocities: Atomic velocities (N_atoms, 3) in Bohr/fs
        masses: Atomic masses (N_atoms,) in amu
        positions: Atomic positions (N_atoms, 3) in Bohr (optional, for rotation removal)

    Returns:
        Velocities with COM motion removed (N_atoms, 3) in Bohr/fs
    """
    n_atoms = len(masses)

    # Step 1: Remove translational motion
    # V_COM = Σ m_i v_i / Σ m_i
    total_mass = np.sum(masses)
    com_velocity = np.sum(masses[:, np.newaxis] * velocities, axis=0) / total_mass

    # Subtract COM velocity from all atoms
    velocities_no_trans = velocities - com_velocity

    logger.debug(f"Removed COM translation: V_COM = {np.linalg.norm(com_velocity):.3e} Bohr/fs")

    # Step 2: Remove rotational motion (if positions provided)
    if positions is not None:
        velocities_no_rot = _remove_rotation(velocities_no_trans, masses, positions)

        # Check angular momentum
        L = _compute_angular_momentum(velocities_no_rot, masses, positions)
        logger.debug(f"Removed COM rotation: |L| = {np.linalg.norm(L):.3e} a.u.")

        return velocities_no_rot
    else:
        return velocities_no_trans


def _remove_rotation(
    velocities: np.ndarray,
    masses: np.ndarray,
    positions: np.ndarray
) -> np.ndarray:
    """
    Remove angular momentum (rotation about COM).

    The rotational velocity of atom i is:
        v_rot = ω × r_i
    where ω is the angular velocity and r_i is position relative to COM.

    We compute ω from the total angular momentum:
        L = I · ω
    where I is the moment of inertia tensor.

    Then subtract rotational component from each atom's velocity.
    """
    n_atoms = len(masses)
    total_mass = np.sum(masses)

    # Center of mass position
    com_position = np.sum(masses[:, np.newaxis] * positions, axis=0) / total_mass

    # Positions relative to COM
    r = positions - com_position

    # Compute moment of inertia tensor
    I = _compute_inertia_tensor(masses, r)

    # Compute angular momentum
    L = _compute_angular_momentum(velocities, masses, r)

    # Angular velocity: ω = I^(-1) · L
    try:
        I_inv = np.linalg.inv(I)
        omega = I_inv @ L
    except np.linalg.LinAlgError:
        # Singular matrix (e.g., linear molecule)
        logger.warning("Singular inertia tensor - cannot remove all rotational motion")
        return velocities

    # Subtract rotational velocity from each atom
    # v_rot = ω × r
    velocities_no_rot = velocities.copy()
    for i in range(n_atoms):
        v_rot = np.cross(omega, r[i])
        velocities_no_rot[i] -= v_rot

    return velocities_no_rot


def _compute_inertia_tensor(masses: np.ndarray, positions: np.ndarray) -> np.ndarray:
    """
    Compute moment of inertia tensor.

    I_ij = Σ_k m_k (r_k² δ_ij - r_k,i r_k,j)
    """
    I = np.zeros((3, 3))

    for m, r in zip(masses, positions):
        r_sq = np.dot(r, r)
        I += m * (r_sq * np.eye(3) - np.outer(r, r))

    return I


def _compute_angular_momentum(
    velocities: np.ndarray,
    masses: np.ndarray,
    positions: np.ndarray
) -> np.ndarray:
    """
    Compute total angular momentum.

    L = Σ_i m_i r_i × v_i
    """
    L = np.zeros(3)

    for m, r, v in zip(masses, positions, velocities):
        L += m * np.cross(r, v)

    return L


def equilibrate_system(
    positions: np.ndarray,
    velocities: np.ndarray,
    masses: np.ndarray,
    force_function,
    target_temperature: float,
    n_steps: int = 1000,
    timestep: float = 0.5,
    thermostat: str = 'berendsen',
    verbose: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Equilibrate system at target temperature.

    Runs a short NVT simulation to:
    1. Relax initial structure
    2. Equilibrate temperature
    3. Achieve thermal equilibrium

    Typical equilibration: 1-10 ps with strong thermostat coupling.

    Args:
        positions: Initial positions (N_atoms, 3) in Bohr
        velocities: Initial velocities (N_atoms, 3) in Bohr/fs
        masses: Atomic masses (N_atoms,) in amu
        force_function: Function to compute forces
        target_temperature: Target T in K
        n_steps: Number of equilibration steps
        timestep: Integration timestep in fs
        thermostat: Thermostat type ('berendsen', 'velocity_rescaling')
        verbose: Print progress

    Returns:
        (equilibrated_positions, equilibrated_velocities)
    """
    from kanad.dynamics.integrators import VelocityVerletIntegrator
    from kanad.dynamics.thermostats import create_thermostat

    if verbose:
        logger.info(f"Equilibrating system: {n_steps} steps, T={target_temperature:.1f}K")

    # Create integrator and thermostat
    integrator = VelocityVerletIntegrator(timestep)
    thermostat_obj = create_thermostat(thermostat, target_temperature, coupling_time=0.1)

    # Initial state
    r = positions.copy()
    v = velocities.copy()

    # Equilibration loop
    for step in range(n_steps):
        # Compute forces
        f, pot_energy = force_function(r)

        # Convert forces to accelerations
        conversion_factor = 974.878  # Bohr/fs² per (Ha/Bohr)/amu
        a = (f / masses[:, np.newaxis]) * conversion_factor

        # Velocity Verlet step
        v_half = v + 0.5 * timestep * a
        r = r + timestep * v_half

        # New forces
        f_new, pot_energy = force_function(r)
        a_new = (f_new / masses[:, np.newaxis]) * conversion_factor
        v = v_half + 0.5 * timestep * a_new

        # Apply thermostat
        v = thermostat_obj.apply(v, masses, timestep)

        # Progress
        if verbose and (step + 1) % (n_steps // 10) == 0:
            T = thermostat_obj.compute_temperature(v, masses)
            logger.info(f"  Step {step+1}/{n_steps}: T={T:.1f}K")

    if verbose:
        T_final = thermostat_obj.compute_temperature(v, masses)
        logger.info(f"Equilibration complete: T_final={T_final:.1f}K")

    return r, v


def generate_initial_conditions(
    positions: np.ndarray,
    masses: np.ndarray,
    temperature: float,
    force_function=None,
    equilibrate: bool = False,
    n_equil_steps: int = 1000,
    seed: Optional[int] = None
) -> InitializationResult:
    """
    Complete initialization procedure for MD simulation.

    Generates initial conditions with proper:
    - Maxwell-Boltzmann velocity distribution
    - COM motion removal
    - Optional equilibration

    Args:
        positions: Atomic positions (N_atoms, 3) in Bohr
        masses: Atomic masses (N_atoms,) in amu
        temperature: Target temperature in K
        force_function: Function to compute forces (for equilibration)
        equilibrate: Run equilibration (recommended: True)
        n_equil_steps: Number of equilibration steps
        seed: Random seed for reproducibility

    Returns:
        InitializationResult with initialized positions and velocities
    """
    logger.info(f"Generating initial conditions: T={temperature:.1f}K")

    # Generate Maxwell-Boltzmann velocities
    mb_init = MaxwellBoltzmannInitializer(temperature, remove_com=True, seed=seed)
    velocities = mb_init.generate(masses, positions)

    # Equilibrate if requested
    if equilibrate:
        if force_function is None:
            raise ValueError("force_function required for equilibration")

        positions, velocities = equilibrate_system(
            positions, velocities, masses, force_function,
            temperature, n_equil_steps
        )

    # Compute final properties
    T_final = mb_init._compute_temperature(velocities, masses)
    total_mass = np.sum(masses)
    com_velocity = np.sum(masses[:, np.newaxis] * velocities, axis=0) / total_mass

    # Kinetic energy
    masses_me = masses * 1822.888486
    velocities_au = velocities / 41.341  # v_au = v_fs / 41.341
    ke = 0.5 * np.sum(masses_me[:, np.newaxis] * velocities_au**2)

    # Angular momentum
    com_position = np.sum(masses[:, np.newaxis] * positions, axis=0) / total_mass
    r_rel = positions - com_position
    L = _compute_angular_momentum(velocities, masses, r_rel)

    result = InitializationResult(
        positions=positions,
        velocities=velocities,
        temperature=T_final,
        kinetic_energy=ke,
        com_velocity=com_velocity,
        angular_momentum=L
    )

    logger.info(f"Initialization complete:")
    logger.info(f"  Temperature: {T_final:.2f} K")
    logger.info(f"  KE: {ke:.6f} Ha")
    logger.info(f"  |V_COM|: {np.linalg.norm(com_velocity):.3e} Bohr/fs")
    logger.info(f"  |L|: {np.linalg.norm(L):.3e} a.u.")

    return result
