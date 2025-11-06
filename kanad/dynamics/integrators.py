"""
Molecular Dynamics Integrators

Time integration algorithms for propagating Newton's equations of motion:
    m_i * d²r_i/dt² = F_i = -∇_i V

All integrators are designed to:
- Conserve energy (for NVE ensemble)
- Be time-reversible (symplectic when possible)
- Handle varying forces (quantum or classical)
- Work with thermostats/barostats

Implemented Integrators:
-----------------------
1. Velocity Verlet: Most popular, 2nd order symplectic
2. Leapfrog: Alternative formulation, also symplectic
3. Runge-Kutta: 4th order accurate (not symplectic)

References:
----------
- Velocity Verlet: Swope et al. (1982) J. Chem. Phys. 76, 637
- Leapfrog: Hockney & Eastwood (1981) Computer Simulation Using Particles
- RK4: Press et al. (2007) Numerical Recipes, 3rd ed.
- Symplectic integrators: Hairer et al. (2006) Geometric Numerical Integration
"""

import numpy as np
import logging
from typing import Tuple, Callable, Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class IntegratorState:
    """
    State of the MD system at a given time.

    Attributes:
        positions: Atomic positions (N_atoms, 3) in Bohr
        velocities: Atomic velocities (N_atoms, 3) in Bohr/fs
        forces: Forces on atoms (N_atoms, 3) in Ha/Bohr
        masses: Atomic masses (N_atoms,) in amu
        time: Current simulation time in fs
        kinetic_energy: Kinetic energy in Hartree
        potential_energy: Potential energy in Hartree
    """
    positions: np.ndarray
    velocities: np.ndarray
    forces: np.ndarray
    masses: np.ndarray
    time: float
    kinetic_energy: float
    potential_energy: float


class BaseIntegrator:
    """
    Base class for molecular dynamics integrators.

    All integrators must implement the step() method which propagates
    the system by one timestep.
    """

    def __init__(self, timestep: float):
        """
        Initialize integrator.

        Args:
            timestep: Integration timestep in femtoseconds
        """
        self.timestep = timestep  # fs
        self.dt = timestep  # Alias

        # Convert timestep to atomic units (1 fs = 41.341 a.u.)
        self.dt_au = timestep * 41.341  # a.u. of time

        logger.debug(f"Initialized {self.__class__.__name__}")
        logger.debug(f"  Timestep: {timestep:.3f} fs ({self.dt_au:.2f} a.u.)")

    def step(
        self,
        state: IntegratorState,
        force_function: Callable
    ) -> IntegratorState:
        """
        Propagate system by one timestep.

        Args:
            state: Current system state
            force_function: Function to compute forces
                Signature: force_function(positions) -> (forces, potential_energy)

        Returns:
            New system state after one timestep
        """
        raise NotImplementedError("Subclasses must implement step()")

    def compute_kinetic_energy(
        self,
        velocities: np.ndarray,
        masses: np.ndarray
    ) -> float:
        """
        Compute kinetic energy: KE = 0.5 * Σ_i m_i v_i²

        Args:
            velocities: (N_atoms, 3) in Bohr/fs
            masses: (N_atoms,) in amu

        Returns:
            Kinetic energy in Hartree
        """
        # Convert masses from amu to electron masses
        # 1 amu = 1822.888486 m_e
        masses_me = masses * 1822.888486

        # Convert velocities from Bohr/fs to a.u. (Bohr / a.u. time)
        # 1 a.u. velocity = 41.341 Bohr/fs
        # Therefore: v_au = v_fs / 41.341
        velocities_au = velocities / 41.341

        # KE = 0.5 * Σ m v²
        ke = 0.5 * np.sum(masses_me[:, np.newaxis] * velocities_au**2)

        return ke  # Hartree


class VelocityVerletIntegrator(BaseIntegrator):
    """
    Velocity Verlet integrator - the gold standard for MD.

    Algorithm:
    ---------
    1. v(t + dt/2) = v(t) + (dt/2) * a(t)
    2. r(t + dt) = r(t) + dt * v(t + dt/2)
    3. Compute a(t + dt) from r(t + dt)
    4. v(t + dt) = v(t + dt/2) + (dt/2) * a(t + dt)

    Properties:
    ----------
    - Symplectic (conserves phase space volume)
    - Time-reversible
    - 2nd order accurate: O(dt³) local error, O(dt²) global
    - Minimal force evaluations (1 per step)
    - Excellent energy conservation

    This is the most widely used MD integrator and should be your default choice.

    References:
    ----------
    - Swope, Andersen, Berens, Wilson (1982) J. Chem. Phys. 76, 637
    - Frenkel & Smit (2002) Understanding Molecular Simulation
    """

    def step(
        self,
        state: IntegratorState,
        force_function: Callable
    ) -> IntegratorState:
        """
        Velocity Verlet step.

        Args:
            state: Current system state
            force_function: Function to compute forces

        Returns:
            New state after one timestep
        """
        # Extract current state
        r = state.positions.copy()
        v = state.velocities.copy()
        f = state.forces.copy()
        m = state.masses

        # Convert forces to accelerations: a = F/m
        # Forces in Ha/Bohr, masses in amu
        # Need to convert: 1 Ha/(Bohr * amu) = 974.878 Bohr/fs²
        conversion_factor = 974.878  # Bohr/fs² per (Ha/Bohr)/amu
        a = (f / m[:, np.newaxis]) * conversion_factor  # Bohr/fs²

        # Step 1: Half-step velocity update
        v_half = v + 0.5 * self.timestep * a

        # Step 2: Full-step position update
        r_new = r + self.timestep * v_half

        # Step 3: Compute forces at new positions
        f_new, pot_energy = force_function(r_new)

        # Convert new forces to accelerations
        a_new = (f_new / m[:, np.newaxis]) * conversion_factor

        # Step 4: Full-step velocity update
        v_new = v_half + 0.5 * self.timestep * a_new

        # Compute kinetic energy
        ke = self.compute_kinetic_energy(v_new, m)

        # Create new state
        new_state = IntegratorState(
            positions=r_new,
            velocities=v_new,
            forces=f_new,
            masses=m,
            time=state.time + self.timestep,
            kinetic_energy=ke,
            potential_energy=pot_energy
        )

        return new_state


class LeapfrogIntegrator(BaseIntegrator):
    """
    Leapfrog integrator - alternative symplectic method.

    Algorithm:
    ---------
    Velocities and positions are staggered by half a timestep:
    1. v(t + dt/2) = v(t - dt/2) + dt * a(t)
    2. r(t + dt) = r(t) + dt * v(t + dt/2)

    Properties:
    ----------
    - Symplectic
    - Time-reversible
    - 2nd order accurate
    - Velocities at half-steps (need synchronization for output)
    - Equivalent to Velocity Verlet but different formulation

    Note: For the first step, we approximate v(t - dt/2) using:
          v(t - dt/2) ≈ v(t) - (dt/2) * a(t)

    References:
    ----------
    - Hockney & Eastwood (1981) Computer Simulation Using Particles
    - Allen & Tildesley (2017) Computer Simulation of Liquids
    """

    def __init__(self, timestep: float):
        super().__init__(timestep)
        self.first_step = True
        self.v_half = None  # Store half-step velocity

    def step(
        self,
        state: IntegratorState,
        force_function: Callable
    ) -> IntegratorState:
        """
        Leapfrog step.

        Args:
            state: Current system state
            force_function: Function to compute forces

        Returns:
            New state after one timestep
        """
        r = state.positions.copy()
        v = state.velocities.copy()
        f = state.forces.copy()
        m = state.masses

        conversion_factor = 974.878  # Bohr/fs² per (Ha/Bohr)/amu
        a = (f / m[:, np.newaxis]) * conversion_factor

        # First step: Initialize v(t - dt/2)
        if self.first_step or self.v_half is None:
            self.v_half = v - 0.5 * self.timestep * a
            self.first_step = False

        # Step 1: Update velocity to v(t + dt/2)
        v_half_new = self.v_half + self.timestep * a

        # Step 2: Update position using half-step velocity
        r_new = r + self.timestep * v_half_new

        # Compute forces at new positions
        f_new, pot_energy = force_function(r_new)
        a_new = (f_new / m[:, np.newaxis]) * conversion_factor

        # Synchronize velocity to full step for output
        # v(t + dt) = v(t + dt/2) + (dt/2) * a(t + dt)
        v_new = v_half_new + 0.5 * self.timestep * a_new

        # Store half-step velocity for next iteration
        self.v_half = v_half_new

        # Compute kinetic energy
        ke = self.compute_kinetic_energy(v_new, m)

        new_state = IntegratorState(
            positions=r_new,
            velocities=v_new,
            forces=f_new,
            masses=m,
            time=state.time + self.timestep,
            kinetic_energy=ke,
            potential_energy=pot_energy
        )

        return new_state


class RungeKuttaIntegrator(BaseIntegrator):
    """
    4th-order Runge-Kutta integrator (RK4).

    Algorithm:
    ---------
    Classical RK4 for ordinary differential equations.
    More accurate than Velocity Verlet but NOT symplectic.

    Properties:
    ----------
    - 4th order accurate: O(dt⁴) local error, O(dt⁴) global
    - NOT symplectic (doesn't conserve phase space volume)
    - Requires 4 force evaluations per step (expensive!)
    - Better for high accuracy short trajectories
    - Use Velocity Verlet for long MD simulations

    Note: RK4 is included for comparison and high-accuracy applications,
          but Velocity Verlet is recommended for typical MD.

    References:
    ----------
    - Press et al. (2007) Numerical Recipes, 3rd ed.
    - Butcher (2008) Numerical Methods for ODEs
    """

    def step(
        self,
        state: IntegratorState,
        force_function: Callable
    ) -> IntegratorState:
        """
        RK4 step with 4 force evaluations.

        Args:
            state: Current system state
            force_function: Function to compute forces

        Returns:
            New state after one timestep
        """
        r = state.positions.copy()
        v = state.velocities.copy()
        f = state.forces.copy()
        m = state.masses
        dt = self.timestep

        conversion_factor = 974.878

        # Helper function to compute derivatives
        def derivatives(pos, vel):
            """Compute dr/dt and dv/dt"""
            forces, _ = force_function(pos)
            accel = (forces / m[:, np.newaxis]) * conversion_factor
            return vel, accel  # dr/dt = v, dv/dt = a

        # RK4 stages
        # k1
        k1_r, k1_v = derivatives(r, v)

        # k2
        r2 = r + 0.5 * dt * k1_r
        v2 = v + 0.5 * dt * k1_v
        k2_r, k2_v = derivatives(r2, v2)

        # k3
        r3 = r + 0.5 * dt * k2_r
        v3 = v + 0.5 * dt * k2_v
        k3_r, k3_v = derivatives(r3, v3)

        # k4
        r4 = r + dt * k3_r
        v4 = v + dt * k3_v
        k4_r, k4_v = derivatives(r4, v4)

        # Weighted average
        r_new = r + (dt / 6.0) * (k1_r + 2*k2_r + 2*k3_r + k4_r)
        v_new = v + (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)

        # Compute final forces and energy
        f_new, pot_energy = force_function(r_new)
        ke = self.compute_kinetic_energy(v_new, m)

        new_state = IntegratorState(
            positions=r_new,
            velocities=v_new,
            forces=f_new,
            masses=m,
            time=state.time + dt,
            kinetic_energy=ke,
            potential_energy=pot_energy
        )

        return new_state


# Factory function for creating integrators
def create_integrator(name: str, timestep: float) -> BaseIntegrator:
    """
    Factory function to create integrators by name.

    Args:
        name: Integrator name ('velocity_verlet', 'leapfrog', 'rk4')
        timestep: Integration timestep in fs

    Returns:
        Integrator instance

    Raises:
        ValueError: If integrator name not recognized
    """
    integrators = {
        'velocity_verlet': VelocityVerletIntegrator,
        'verlet': VelocityVerletIntegrator,  # Alias
        'leapfrog': LeapfrogIntegrator,
        'rk4': RungeKuttaIntegrator,
        'runge_kutta': RungeKuttaIntegrator,  # Alias
    }

    name_lower = name.lower()
    if name_lower not in integrators:
        available = ', '.join(integrators.keys())
        raise ValueError(f"Unknown integrator '{name}'. Available: {available}")

    return integrators[name_lower](timestep)
