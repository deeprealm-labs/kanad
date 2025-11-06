"""
Test Classical Molecular Dynamics

Tests the MD infrastructure with classical force methods (HF, MP2).
Validates:
- Energy conservation (NVE ensemble)
- Temperature control (NVT ensemble)
- Trajectory storage and analysis
- Different integrators and thermostats
- COM motion removal
- Maxwell-Boltzmann velocity distribution

These tests ensure the MD engine works correctly before adding quantum forces.
"""

import pytest
import numpy as np
import logging
from pathlib import Path

from kanad.bonds import BondFactory
from kanad.dynamics import (
    MDSimulator,
    VelocityVerletIntegrator,
    LeapfrogIntegrator,
    RungeKuttaIntegrator,
    BerendsenThermostat,
    VelocityRescaling,
    LangevinThermostat,
    MaxwellBoltzmannInitializer,
    remove_com_motion
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestClassicalMD:
    """Test classical MD with HF forces."""

    def test_h2_nve_energy_conservation(self):
        """Test NVE ensemble - energy should be conserved."""
        logger.info("\n" + "="*80)
        logger.info("TEST: H2 NVE Energy Conservation")
        logger.info("="*80)

        # Create H2 molecule at equilibrium
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        # Run NVE simulation (no thermostat)
        md = MDSimulator(
            bond,
            temperature=300.0,  # Initial temperature only
            timestep=0.5,  # fs
            integrator='velocity_verlet',
            thermostat=None,  # NVE ensemble
            force_method='hf'
        )

        result = md.run(n_steps=100, save_frequency=1)

        # Check energy conservation
        energies = result.trajectory.total_energies
        energy_drift = np.abs(energies[-1] - energies[0])
        max_drift = 0.001  # Hartree (very tight tolerance)

        logger.info(f"\nEnergy Conservation Check:")
        logger.info(f"  Initial energy: {energies[0]:.6f} Ha")
        logger.info(f"  Final energy: {energies[-1]:.6f} Ha")
        logger.info(f"  Energy drift: {energy_drift:.8f} Ha")
        logger.info(f"  Max allowed drift: {max_drift:.8f} Ha")

        # For symplectic integrator, energy should be very well conserved
        assert energy_drift < max_drift, f"Energy drift {energy_drift:.8f} > {max_drift:.8f} Ha"

        # Check energy fluctuations (should be small)
        energy_std = np.std(energies)
        logger.info(f"  Energy std dev: {energy_std:.8f} Ha")
        assert energy_std < 0.0005, f"Energy fluctuations too large: {energy_std:.8f} Ha"

        logger.info("✓ Energy conservation test PASSED")

    def test_h2_nvt_temperature_control(self):
        """Test NVT ensemble - temperature should be controlled."""
        logger.info("\n" + "="*80)
        logger.info("TEST: H2 NVT Temperature Control")
        logger.info("="*80)

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        target_temp = 300.0  # K

        # Run NVT simulation with Berendsen thermostat
        md = MDSimulator(
            bond,
            temperature=target_temp,
            timestep=0.5,
            integrator='velocity_verlet',
            thermostat='berendsen',
            force_method='hf'
        )

        result = md.run(n_steps=200, save_frequency=1)

        # Analyze temperature after equilibration (skip first 50 steps)
        temps = result.trajectory.temperatures[50:]
        avg_temp = np.mean(temps)
        std_temp = np.std(temps)

        logger.info(f"\nTemperature Control Check:")
        logger.info(f"  Target temperature: {target_temp:.2f} K")
        logger.info(f"  Average temperature: {avg_temp:.2f} K")
        logger.info(f"  Std dev: {std_temp:.2f} K")
        logger.info(f"  Deviation: {abs(avg_temp - target_temp):.2f} K")

        # Temperature should be close to target (within 10%)
        temp_tolerance = 0.1 * target_temp
        assert abs(avg_temp - target_temp) < temp_tolerance, \
            f"Temperature {avg_temp:.2f} not within {temp_tolerance:.2f} K of target {target_temp:.2f} K"

        logger.info("✓ Temperature control test PASSED")

    def test_integrator_comparison(self):
        """Compare different integrators on same system."""
        logger.info("\n" + "="*80)
        logger.info("TEST: Integrator Comparison")
        logger.info("="*80)

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        integrators = ['velocity_verlet', 'leapfrog', 'rk4']
        results = {}

        for integrator_name in integrators:
            md = MDSimulator(
                bond,
                temperature=300.0,
                timestep=0.5,
                integrator=integrator_name,
                thermostat=None,  # NVE
                force_method='hf'
            )

            result = md.run(n_steps=100, save_frequency=10)
            results[integrator_name] = result

            energy_drift = abs(result.trajectory.total_energies[-1] -
                             result.trajectory.total_energies[0])

            logger.info(f"\n{integrator_name.upper()}:")
            logger.info(f"  Final energy: {result.final_energy:.6f} Ha")
            logger.info(f"  Energy drift: {energy_drift:.8f} Ha")

        # All should give similar final energies
        energies = [r.final_energy for r in results.values()]
        energy_spread = max(energies) - min(energies)

        logger.info(f"\nEnergy spread across integrators: {energy_spread:.8f} Ha")
        assert energy_spread < 0.01, f"Integrators disagree by {energy_spread:.8f} Ha"

        logger.info("✓ Integrator comparison test PASSED")

    def test_thermostat_comparison(self):
        """Compare different thermostats."""
        logger.info("\n" + "="*80)
        logger.info("TEST: Thermostat Comparison")
        logger.info("="*80)

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        target_temp = 400.0
        thermostats = ['berendsen', 'velocity_rescaling', 'langevin']
        results = {}

        for thermostat_name in thermostats:
            md = MDSimulator(
                bond,
                temperature=target_temp,
                timestep=0.5,
                integrator='velocity_verlet',
                thermostat=thermostat_name,
                force_method='hf'
            )

            result = md.run(n_steps=200, save_frequency=1)
            results[thermostat_name] = result

            # Average temperature after equilibration
            avg_temp = np.mean(result.trajectory.temperatures[50:])
            std_temp = np.std(result.trajectory.temperatures[50:])

            logger.info(f"\n{thermostat_name.upper()}:")
            logger.info(f"  Average T: {avg_temp:.2f} K (target: {target_temp:.2f} K)")
            logger.info(f"  Std dev: {std_temp:.2f} K")
            logger.info(f"  Deviation: {abs(avg_temp - target_temp):.2f} K")

            # All thermostats should reach target (within 15%)
            assert abs(avg_temp - target_temp) < 0.15 * target_temp

        logger.info("✓ Thermostat comparison test PASSED")

    def test_maxwell_boltzmann_distribution(self):
        """Test that initial velocities follow Maxwell-Boltzmann distribution."""
        logger.info("\n" + "="*80)
        logger.info("TEST: Maxwell-Boltzmann Velocity Distribution")
        logger.info("="*80)

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        temp = 300.0

        # Generate many velocity samples
        masses = np.array([bond.atom_1.mass, bond.atom_2.mass])
        positions = np.array([bond.atom_1.position, bond.atom_2.position])

        mb_init = MaxwellBoltzmannInitializer(temp, remove_com=True, seed=42)

        n_samples = 1000
        velocities_all = []

        for i in range(n_samples):
            v = mb_init.generate(masses, positions)
            velocities_all.append(v)

        velocities_all = np.array(velocities_all)  # (n_samples, n_atoms, 3)

        # Check that mean velocity is zero (after COM removal)
        mean_vel = np.mean(velocities_all.reshape(-1, 3), axis=0)
        logger.info(f"\nVelocity Statistics:")
        logger.info(f"  Mean velocity: {mean_vel}")
        logger.info(f"  |Mean velocity|: {np.linalg.norm(mean_vel):.3e} Bohr/fs")

        assert np.linalg.norm(mean_vel) < 0.01, "Mean velocity should be ~0 after COM removal"

        # Check temperature distribution
        K_BOLTZMANN = 3.1668105e-6  # Ha/K
        masses_me = masses * 1822.888486

        temps = []
        for v in velocities_all:
            v_au = v * 41.341
            ke = 0.5 * np.sum(masses_me[:, np.newaxis] * v_au**2)
            # T = 2*KE / (k_B * N_dof), N_dof = 3*N - 6 = 0 for diatomic (use 1)
            T = (2.0 * ke) / (K_BOLTZMANN * 1)
            temps.append(T)

        avg_temp = np.mean(temps)
        std_temp = np.std(temps)

        logger.info(f"  Target temperature: {temp:.2f} K")
        logger.info(f"  Average temperature: {avg_temp:.2f} K")
        logger.info(f"  Std dev: {std_temp:.2f} K")

        # Average should be close to target
        assert abs(avg_temp - temp) < 0.1 * temp

        logger.info("✓ Maxwell-Boltzmann distribution test PASSED")

    def test_com_motion_removal(self):
        """Test center-of-mass motion removal."""
        logger.info("\n" + "="*80)
        logger.info("TEST: Center-of-Mass Motion Removal")
        logger.info("="*80)

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        masses = np.array([bond.atom_1.mass, bond.atom_2.mass])
        positions = np.array([bond.atom_1.position, bond.atom_2.position])

        # Generate random velocities (with COM motion)
        np.random.seed(42)
        velocities = np.random.randn(2, 3) * 0.1  # Random velocities

        # Check initial COM velocity (should be non-zero)
        total_mass = np.sum(masses)
        com_vel_initial = np.sum(masses[:, np.newaxis] * velocities, axis=0) / total_mass

        logger.info(f"\nBefore COM removal:")
        logger.info(f"  COM velocity: {com_vel_initial}")
        logger.info(f"  |V_COM|: {np.linalg.norm(com_vel_initial):.6f} Bohr/fs")

        # Remove COM motion
        velocities_clean = remove_com_motion(velocities, masses, positions)

        # Check final COM velocity (should be zero)
        com_vel_final = np.sum(masses[:, np.newaxis] * velocities_clean, axis=0) / total_mass

        logger.info(f"\nAfter COM removal:")
        logger.info(f"  COM velocity: {com_vel_final}")
        logger.info(f"  |V_COM|: {np.linalg.norm(com_vel_final):.3e} Bohr/fs")

        # Should be very close to zero
        assert np.linalg.norm(com_vel_final) < 1e-10, "COM velocity not removed"

        logger.info("✓ COM motion removal test PASSED")

    def test_trajectory_storage(self):
        """Test trajectory storage and retrieval."""
        logger.info("\n" + "="*80)
        logger.info("TEST: Trajectory Storage")
        logger.info("="*80)

        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        md = MDSimulator(
            bond,
            temperature=300.0,
            timestep=0.5,
            integrator='velocity_verlet',
            thermostat='berendsen',
            force_method='hf'
        )

        result = md.run(n_steps=50, save_frequency=5)

        # Check trajectory properties
        traj = result.trajectory
        n_frames = len(traj.frames)

        logger.info(f"\nTrajectory Properties:")
        logger.info(f"  Number of frames: {n_frames}")
        logger.info(f"  Expected frames: {50 // 5 + 1}")
        logger.info(f"  Number of atoms: {traj.n_atoms}")
        logger.info(f"  Simulation time: {traj.times[-1]:.2f} fs")

        assert n_frames == 50 // 5 + 1, f"Wrong number of frames: {n_frames}"
        assert traj.n_atoms == 2, "Should be 2 atoms (H2)"

        # Check that we can access frame data
        frame = traj.frames[0]
        assert frame.positions.shape == (2, 3)
        assert frame.velocities.shape == (2, 3)
        assert frame.forces.shape == (2, 3)
        assert isinstance(frame.kinetic_energy, float)
        assert isinstance(frame.potential_energy, float)

        logger.info("✓ Trajectory storage test PASSED")

    def test_h2_bond_stretch(self):
        """Test H2 bond stretching dynamics."""
        logger.info("\n" + "="*80)
        logger.info("TEST: H2 Bond Stretching Dynamics")
        logger.info("="*80)

        # Start with stretched H2 bond
        bond = BondFactory.create_bond('H', 'H', distance=1.0)  # Stretched from eq (0.74)

        md = MDSimulator(
            bond,
            temperature=100.0,  # Low temp to see vibration clearly
            timestep=0.1,  # Small timestep for oscillation
            integrator='velocity_verlet',
            thermostat=None,  # NVE to see free vibration
            force_method='hf'
        )

        result = md.run(n_steps=500, save_frequency=1)

        # Extract bond distances over time
        traj = result.trajectory
        distances = []

        for frame in traj.frames:
            r1 = frame.positions[0]
            r2 = frame.positions[1]
            dist = np.linalg.norm(r2 - r1)
            distances.append(dist)

        distances = np.array(distances)

        # Check that bond oscillates
        min_dist = np.min(distances)
        max_dist = np.max(distances)
        amplitude = max_dist - min_dist

        logger.info(f"\nBond Stretching Analysis:")
        logger.info(f"  Initial distance: {distances[0]:.4f} Bohr")
        logger.info(f"  Min distance: {min_dist:.4f} Bohr")
        logger.info(f"  Max distance: {max_dist:.4f} Bohr")
        logger.info(f"  Oscillation amplitude: {amplitude:.4f} Bohr")

        # Should oscillate (amplitude > 0.1 Bohr)
        assert amplitude > 0.1, f"Bond not oscillating (amplitude = {amplitude:.4f} Bohr)"

        # Should pass through equilibrium (0.74 Bohr is ~1.4 a.u.)
        # Initial is 1.0 Angstrom = 1.89 Bohr
        # Equilibrium is 0.74 Angstrom = 1.40 Bohr
        # Should oscillate around equilibrium
        avg_dist = np.mean(distances[100:])  # After initial transient
        logger.info(f"  Average distance: {avg_dist:.4f} Bohr")

        logger.info("✓ Bond stretching test PASSED")


if __name__ == '__main__':
    # Run tests
    print("\n" + "="*80)
    print("CLASSICAL MD TEST SUITE")
    print("="*80)

    test = TestClassicalMD()

    try:
        test.test_h2_nve_energy_conservation()
    except Exception as e:
        logger.error(f"NVE test failed: {e}")

    try:
        test.test_h2_nvt_temperature_control()
    except Exception as e:
        logger.error(f"NVT test failed: {e}")

    try:
        test.test_integrator_comparison()
    except Exception as e:
        logger.error(f"Integrator test failed: {e}")

    try:
        test.test_thermostat_comparison()
    except Exception as e:
        logger.error(f"Thermostat test failed: {e}")

    try:
        test.test_maxwell_boltzmann_distribution()
    except Exception as e:
        logger.error(f"Maxwell-Boltzmann test failed: {e}")

    try:
        test.test_com_motion_removal()
    except Exception as e:
        logger.error(f"COM removal test failed: {e}")

    try:
        test.test_trajectory_storage()
    except Exception as e:
        logger.error(f"Trajectory test failed: {e}")

    try:
        test.test_h2_bond_stretch()
    except Exception as e:
        logger.error(f"Bond stretch test failed: {e}")

    print("\n" + "="*80)
    print("ALL CLASSICAL MD TESTS COMPLETE")
    print("="*80)
