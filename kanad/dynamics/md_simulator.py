"""
Main Molecular Dynamics Simulator

Orchestrates all MD components to run simulations with:
- Multiple force methods (HF, MP2, VQE, SQD)
- Various integrators (Velocity Verlet, Leapfrog, RK4)
- Temperature control (Berendsen, Nose-Hoover, Langevin)
- Trajectory storage and analysis
- Energy monitoring and conservation checks

This is the primary user interface for running MD simulations in Kanad.

Example Usage:
-------------
```python
from kanad.bonds import BondFactory
from kanad.dynamics import MDSimulator

# Create molecule
bond = BondFactory.create_bond('H', 'H', distance=0.74)

# Run classical MD with HF forces
md = MDSimulator(
    bond_or_molecule=bond,
    temperature=300.0,  # K
    timestep=0.5,  # fs
    integrator='velocity_verlet',
    thermostat='berendsen',
    force_method='hf'
)

result = md.run(n_steps=1000, save_trajectory=True)

# Access results
print(f"Final energy: {result.final_energy:.6f} Ha")
print(f"Average temperature: {result.avg_temperature:.2f} K")

# Or run quantum MD with VQE
md_quantum = MDSimulator(
    bond,
    temperature=300.0,
    timestep=0.5,
    force_method='vqe',
    use_governance=True
)

result = md_quantum.run(n_steps=100)
```

References:
----------
- Allen & Tildesley (2017) Computer Simulation of Liquids
- Frenkel & Smit (2002) Understanding Molecular Simulation
- Car & Parrinello (1985) Ab initio MD: Phys. Rev. Lett. 55, 2471
"""

import numpy as np
import logging
from typing import Optional, Callable, Dict, Any, Union
from dataclasses import dataclass, field
from pathlib import Path
import time as time_module

from kanad.dynamics.integrators import create_integrator, IntegratorState
from kanad.dynamics.thermostats import create_thermostat
from kanad.dynamics.trajectory import Trajectory, TrajectoryWriter
from kanad.dynamics.initialization import generate_initial_conditions

logger = logging.getLogger(__name__)


@dataclass
class MDResult:
    """
    Result from MD simulation.

    Attributes:
        trajectory: Full trajectory (if saved)
        final_positions: Final atomic positions (N_atoms, 3) Bohr
        final_velocities: Final velocities (N_atoms, 3) Bohr/fs
        final_energy: Final total energy (Hartree)
        avg_temperature: Average temperature (K)
        avg_kinetic_energy: Average KE (Hartree)
        avg_potential_energy: Average PE (Hartree)
        avg_total_energy: Average total energy (Hartree)
        energy_drift: Energy drift (Hartree)
        temperature_std: Temperature std deviation (K)
        n_steps_completed: Number of steps completed
        wall_time: Wall clock time (seconds)
        steps_per_second: Performance metric
        converged: Whether simulation completed successfully
        metadata: Additional simulation information
    """
    trajectory: Optional[Trajectory] = None
    final_positions: Optional[np.ndarray] = None
    final_velocities: Optional[np.ndarray] = None
    final_energy: float = 0.0
    avg_temperature: float = 0.0
    avg_kinetic_energy: float = 0.0
    avg_potential_energy: float = 0.0
    avg_total_energy: float = 0.0
    energy_drift: float = 0.0
    temperature_std: float = 0.0
    n_steps_completed: int = 0
    wall_time: float = 0.0
    steps_per_second: float = 0.0
    converged: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)


class MDSimulator:
    """
    Main molecular dynamics simulator.

    Coordinates integrator, thermostat, force computation, and trajectory storage
    to run MD simulations with various force methods (classical or quantum).
    """

    def __init__(
        self,
        bond_or_molecule,
        temperature: float,
        timestep: float = 0.5,
        integrator: str = 'velocity_verlet',
        thermostat: Optional[str] = 'berendsen',
        force_method: str = 'hf',
        use_governance: bool = False,
        backend: str = 'statevector',
        initial_velocities: Optional[np.ndarray] = None,
        random_seed: Optional[int] = None,
        **kwargs
    ):
        """
        Initialize MD simulator.

        Args:
            bond_or_molecule: Bond or Molecule object
            temperature: Target temperature in K
            timestep: Integration timestep in fs
            integrator: Integrator type ('velocity_verlet', 'leapfrog', 'rk4')
            thermostat: Thermostat type ('berendsen', 'nose_hoover', 'langevin', None for NVE)
            force_method: Force computation method ('hf', 'mp2', 'hivqe', 'vqe', 'sqd')
                         Note: 'hivqe' is recommended for quantum MD (more efficient than 'vqe')
            use_governance: Use governance protocols (for VQE/SQD)
            backend: Quantum backend ('statevector', 'aer', 'ibm')
            initial_velocities: Initial velocities (None = Maxwell-Boltzmann)
            random_seed: Random seed for initialization
            **kwargs: Additional parameters for integrator/thermostat/forces
        """
        self.bond_or_molecule = bond_or_molecule
        self.temperature = temperature

        # Initialize solver cache for quantum forces (CRITICAL for performance!)
        self.solver_cache = {}
        self.timestep = timestep
        self.force_method = force_method.lower()
        self.use_governance = use_governance
        self.backend = backend
        self.random_seed = random_seed

        # Extract atomic information
        self._setup_system()

        # Create integrator
        self.integrator = create_integrator(integrator, timestep)
        logger.info(f"Created integrator: {integrator}")

        # Create thermostat (if specified)
        if thermostat is not None:
            thermostat_kwargs = {}
            if 'coupling_time' in kwargs:
                thermostat_kwargs['coupling_time'] = kwargs['coupling_time']
            if 'friction_coefficient' in kwargs:
                thermostat_kwargs['friction_coefficient'] = kwargs['friction_coefficient']

            self.thermostat = create_thermostat(thermostat, temperature, **thermostat_kwargs)
            self.ensemble = 'NVT'
            logger.info(f"Created thermostat: {thermostat} (NVT ensemble)")
        else:
            self.thermostat = None
            self.ensemble = 'NVE'
            logger.info("No thermostat (NVE ensemble - energy conservation)")

        # Setup force computation
        self._setup_forces(**kwargs)

        # Initialize velocities
        if initial_velocities is not None:
            self.velocities = initial_velocities.copy()
            logger.info("Using provided initial velocities")
        else:
            # Will generate Maxwell-Boltzmann in run()
            self.velocities = None
            logger.info("Will generate Maxwell-Boltzmann velocities")

        # Statistics
        self.n_steps_run = 0

        logger.info(f"MDSimulator initialized:")
        logger.info(f"  System: {self.n_atoms} atoms")
        logger.info(f"  Temperature: {temperature:.1f} K")
        logger.info(f"  Timestep: {timestep:.3f} fs")
        logger.info(f"  Force method: {force_method}")
        logger.info(f"  Ensemble: {self.ensemble}")

    def _setup_system(self):
        """Extract atomic positions, masses, symbols from bond/molecule."""
        # Check if it's a bond or molecule
        if hasattr(self.bond_or_molecule, 'atom_1') and hasattr(self.bond_or_molecule, 'atom_2'):
            # It's a bond
            atoms = [self.bond_or_molecule.atom_1, self.bond_or_molecule.atom_2]
            self.is_bond = True
        elif hasattr(self.bond_or_molecule, 'atoms'):
            # It's a molecule
            atoms = self.bond_or_molecule.atoms
            self.is_bond = False
        else:
            raise ValueError("Input must be Bond or Molecule object")

        self.atoms = atoms
        self.n_atoms = len(atoms)

        # Extract positions, masses, symbols
        self.positions = np.array([atom.position for atom in atoms])
        self.masses = np.array([atom.atomic_mass for atom in atoms])
        self.symbols = [atom.symbol for atom in atoms]

        logger.debug(f"System setup: {self.n_atoms} atoms, total mass = {np.sum(self.masses):.2f} amu")

    def _setup_forces(self, **kwargs):
        """Setup force computation based on force_method."""
        if self.force_method in ['hf', 'mp2']:
            # Classical ab initio forces using PySCF gradients
            from kanad.core.gradients import GradientCalculator

            # Get hamiltonian
            if hasattr(self.bond_or_molecule, 'hamiltonian'):
                hamiltonian = self.bond_or_molecule.hamiltonian
            else:
                raise ValueError("Bond/Molecule must have hamiltonian for force calculation")

            # Create gradient calculator
            self.gradient_calc = GradientCalculator(
                self.bond_or_molecule,
                method=self.force_method.upper()
            )

            logger.info(f"Using {self.force_method.upper()} forces via PySCF gradients")

        elif self.force_method in ['hivqe', 'vqe', 'sqd']:
            # Quantum forces with solver caching for performance
            logger.info(f"Using quantum {self.force_method.upper()} forces")
            logger.info(f"  Backend: {self.backend}")
            logger.info(f"  Governance: {self.use_governance}")
            logger.info(f"  Solver caching: ENABLED (critical for performance!)")

            # Store parameters for quantum force computation
            self.quantum_params = {
                'backend': self.backend,
                'use_governance': self.use_governance,
                'solver_cache': self.solver_cache,  # Pass cache for reuse!
                **kwargs
            }

        else:
            raise ValueError(f"Unknown force_method: {self.force_method}")

    def compute_forces(self, positions: np.ndarray) -> tuple:
        """
        Compute forces on atoms at given positions.

        Args:
            positions: Atomic positions (N_atoms, 3) in Bohr

        Returns:
            (forces, potential_energy):
                forces: (N_atoms, 3) in Ha/Bohr
                potential_energy: Potential energy in Hartree
        """
        if self.force_method in ['hf', 'mp2']:
            # Update positions
            for i, atom in enumerate(self.atoms):
                atom.position = positions[i]

            # Compute gradient
            result = self.gradient_calc.compute_gradient()

            forces = result['forces']  # Already -gradient
            potential_energy = result['energy']

            return forces, potential_energy

        elif self.force_method in ['hivqe', 'vqe', 'sqd']:
            # Quantum forces - call quantum MD module with caching!
            from kanad.dynamics.quantum_md import compute_quantum_forces

            forces, potential_energy = compute_quantum_forces(
                positions,
                self.bond_or_molecule,
                method=self.force_method,
                **self.quantum_params
            )

            return forces, potential_energy

        else:
            raise ValueError(f"Unknown force_method: {self.force_method}")

    def run(
        self,
        n_steps: int,
        save_frequency: int = 10,
        save_trajectory: bool = True,
        output_file: Optional[str] = None,
        equilibrate: bool = False,
        n_equil_steps: int = 100,
        verbose: bool = True,
        check_energy: bool = True
    ) -> MDResult:
        """
        Run molecular dynamics simulation.

        Args:
            n_steps: Number of MD steps
            save_frequency: Save trajectory every N steps
            save_trajectory: Store trajectory in memory
            output_file: Save trajectory to file (HDF5 or XYZ)
            equilibrate: Run equilibration before production
            n_equil_steps: Equilibration steps
            verbose: Print progress
            check_energy: Monitor energy conservation (NVE only)

        Returns:
            MDResult with simulation results and trajectory
        """
        logger.info("=" * 70)
        logger.info("STARTING MOLECULAR DYNAMICS SIMULATION")
        logger.info("=" * 70)

        start_time = time_module.time()

        # Initialize velocities if not provided
        if self.velocities is None:
            if verbose:
                logger.info("Generating Maxwell-Boltzmann initial velocities...")

            init_result = generate_initial_conditions(
                self.positions,
                self.masses,
                self.temperature,
                force_function=self.compute_forces if equilibrate else None,
                equilibrate=equilibrate,
                n_equil_steps=n_equil_steps,
                seed=self.random_seed
            )

            self.positions = init_result.positions
            self.velocities = init_result.velocities

            if verbose:
                logger.info(f"Initial temperature: {init_result.temperature:.2f} K")
                logger.info(f"Initial KE: {init_result.kinetic_energy:.6f} Ha")

        # Initialize trajectory
        if save_trajectory:
            trajectory = Trajectory()
            trajectory.n_atoms = self.n_atoms
            trajectory.atom_symbols = self.symbols
            trajectory.atom_masses = self.masses
        else:
            trajectory = None

        # Compute initial forces
        forces, pot_energy = self.compute_forces(self.positions)

        # Compute initial kinetic energy
        ke = self.integrator.compute_kinetic_energy(self.velocities, self.masses)

        # Create initial state
        state = IntegratorState(
            positions=self.positions.copy(),
            velocities=self.velocities.copy(),
            forces=forces,
            masses=self.masses,
            time=0.0,
            kinetic_energy=ke,
            potential_energy=pot_energy
        )

        # Save initial frame
        if save_trajectory and (0 % save_frequency == 0):
            if self.thermostat:
                T = self.thermostat.compute_temperature(state.velocities, state.masses)
            else:
                T = 0.0

            trajectory.add_frame(
                state.positions,
                state.velocities,
                state.forces,
                state.kinetic_energy,
                state.potential_energy,
                T,
                state.time
            )

        # Statistics
        initial_energy = state.kinetic_energy + state.potential_energy
        energies = []
        temperatures = []

        # MD loop
        if verbose:
            logger.info(f"\nRunning {n_steps} MD steps...")
            logger.info(f"  Timestep: {self.timestep:.3f} fs")
            logger.info(f"  Total time: {n_steps * self.timestep:.2f} fs")
            logger.info(f"  Ensemble: {self.ensemble}")

        for step in range(1, n_steps + 1):
            # Integration step
            state = self.integrator.step(state, self.compute_forces)

            # Apply thermostat (if NVT)
            if self.thermostat is not None:
                state.velocities = self.thermostat.apply(
                    state.velocities,
                    state.masses,
                    self.timestep,
                    current_ke=state.kinetic_energy
                )

                # Recompute kinetic energy after thermostat
                state.kinetic_energy = self.integrator.compute_kinetic_energy(
                    state.velocities,
                    state.masses
                )

            # Compute temperature
            if self.thermostat:
                T = self.thermostat.compute_temperature(state.velocities, state.masses)
            else:
                n_dof = 3 * self.n_atoms - 6
                if n_dof <= 0:
                    n_dof = 1  # Avoid division by zero for small molecules
                T = (2.0 * state.kinetic_energy) / (3.1668105e-6 * n_dof)

            # Statistics
            total_energy = state.kinetic_energy + state.potential_energy
            energies.append(total_energy)
            temperatures.append(T)

            # Save frame
            if save_trajectory and (step % save_frequency == 0):
                trajectory.add_frame(
                    state.positions,
                    state.velocities,
                    state.forces,
                    state.kinetic_energy,
                    state.potential_energy,
                    T,
                    state.time
                )

            # Progress
            if verbose and (step % (n_steps // 10) == 0 or step == n_steps):
                energy_drift = total_energy - initial_energy
                logger.info(f"  Step {step}/{n_steps}: "
                          f"T={T:.1f}K, E_tot={total_energy:.6f}Ha, "
                          f"ΔE={energy_drift:.6e}Ha")

        # End timing
        wall_time = time_module.time() - start_time
        steps_per_second = n_steps / wall_time

        # Compute final statistics
        final_energy = energies[-1]
        avg_energy = np.mean(energies)
        energy_drift = final_energy - initial_energy

        avg_temp = np.mean(temperatures)
        temp_std = np.std(temperatures)

        # Create result
        result = MDResult(
            trajectory=trajectory,
            final_positions=state.positions,
            final_velocities=state.velocities,
            final_energy=final_energy,
            avg_temperature=avg_temp,
            avg_kinetic_energy=np.mean([e - p for e, p in zip(energies, [state.potential_energy]*len(energies))]),
            avg_potential_energy=state.potential_energy,  # Approximate
            avg_total_energy=avg_energy,
            energy_drift=energy_drift,
            temperature_std=temp_std,
            n_steps_completed=n_steps,
            wall_time=wall_time,
            steps_per_second=steps_per_second,
            converged=True,
            metadata={
                'force_method': self.force_method,
                'integrator': type(self.integrator).__name__,
                'thermostat': type(self.thermostat).__name__ if self.thermostat else 'None',
                'ensemble': self.ensemble,
                'timestep': self.timestep,
                'temperature': self.temperature,
            }
        )

        # Save trajectory to file
        if output_file is not None and trajectory is not None:
            output_path = Path(output_file)
            if output_path.suffix in ['.h5', '.hdf5']:
                writer = TrajectoryWriter('hdf5')
            elif output_path.suffix == '.xyz':
                writer = TrajectoryWriter('xyz')
            else:
                logger.warning(f"Unknown file extension: {output_path.suffix}, using HDF5")
                writer = TrajectoryWriter('hdf5')

            writer.write(trajectory, output_file)
            logger.info(f"Trajectory saved to {output_file}")

        # Final summary
        logger.info("=" * 70)
        logger.info("MD SIMULATION COMPLETE")
        logger.info("=" * 70)
        logger.info(f"Final results:")
        logger.info(f"  Steps completed: {n_steps}")
        logger.info(f"  Final energy: {final_energy:.6f} Ha")
        logger.info(f"  Average energy: {avg_energy:.6f} Ha")
        logger.info(f"  Energy drift: {energy_drift:.6e} Ha ({abs(energy_drift/initial_energy)*100:.3f}%)")
        logger.info(f"  Average temperature: {avg_temp:.2f} ± {temp_std:.2f} K")
        logger.info(f"  Wall time: {wall_time:.2f} s")
        logger.info(f"  Performance: {steps_per_second:.1f} steps/s")

        if check_energy and self.ensemble == 'NVE':
            energy_conservation = abs(energy_drift / initial_energy)
            if energy_conservation < 1e-4:
                logger.info(f"  ✅ Excellent energy conservation (<0.01%)")
            elif energy_conservation < 1e-3:
                logger.info(f"  ✅ Good energy conservation (<0.1%)")
            elif energy_conservation < 1e-2:
                logger.warning(f"  ⚠️  Moderate energy drift (>0.1%)")
            else:
                logger.warning(f"  ⚠️  Poor energy conservation (>1%) - reduce timestep!")

        return result
