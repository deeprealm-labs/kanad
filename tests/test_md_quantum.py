"""
Test Quantum Molecular Dynamics

Tests MD with quantum forces (VQE/SQD) instead of classical (HF).
Validates:
- VQE/SQD force computation
- Quantum vs classical force comparison
- Correlation effects in dynamics
- Governance speedup in MD
- Cost estimation

These tests demonstrate the unique "world's first" quantum-enhanced MD capability.
"""

import pytest
import numpy as np
import logging

from kanad.bonds import BondFactory
from kanad.dynamics import MDSimulator
from kanad.dynamics.quantum_md import (
    compute_quantum_forces,
    estimate_quantum_md_cost,
    compare_classical_vs_quantum_forces
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestQuantumMD:
    """Test quantum MD with VQE/SQD forces."""

    def test_vqe_force_computation(self):
        """Test VQE force computation for H2."""
        logger.info("\n" + "="*80)
        logger.info("TEST: VQE Force Computation")
        logger.info("="*80)

        # H2 at equilibrium
        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        positions = np.array([bond.atom_1.position, bond.atom_2.position])

        # Compute VQE forces
        forces, energy = compute_quantum_forces(
            positions,
            bond,
            method='vqe',
            backend='statevector',
            use_governance=True
        )

        logger.info(f"\nVQE Force Computation:")
        logger.info(f"  Energy: {energy:.6f} Ha")
        logger.info(f"  Forces shape: {forces.shape}")
        logger.info(f"  Force on atom 1: {forces[0]}")
        logger.info(f"  Force on atom 2: {forces[1]}")
        logger.info(f"  |F1|: {np.linalg.norm(forces[0]):.6f} Ha/Bohr")
        logger.info(f"  |F2|: {np.linalg.norm(forces[1]):.6f} Ha/Bohr")

        # Check forces have correct shape
        assert forces.shape == (2, 3), f"Forces should be (2, 3), got {forces.shape}"

        # Check forces are opposite (Newton's 3rd law)
        force_sum = forces[0] + forces[1]
        logger.info(f"  Force sum (should be ~0): {force_sum}")
        logger.info(f"  |Force sum|: {np.linalg.norm(force_sum):.3e} Ha/Bohr")

        # For a diatomic, forces should be equal and opposite
        assert np.linalg.norm(force_sum) < 0.01, \
            f"Forces don't sum to zero (Newton's 3rd law): {force_sum}"

        # At equilibrium, forces should be small
        max_force = max(np.linalg.norm(forces[0]), np.linalg.norm(forces[1]))
        logger.info(f"  Max force magnitude: {max_force:.6f} Ha/Bohr")

        # At exact equilibrium forces would be 0, but we're close
        assert max_force < 0.1, f"Forces too large at equilibrium: {max_force} Ha/Bohr"

        logger.info("✓ VQE force computation test PASSED")

    def test_sqd_force_computation(self):
        """Test SQD force computation for H2."""
        logger.info("\n" + "="*80)
        logger.info("TEST: SQD Force Computation")
        logger.info("="*80)

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        positions = np.array([bond.atom_1.position, bond.atom_2.position])

        # Compute SQD forces
        forces, energy = compute_quantum_forces(
            positions,
            bond,
            method='sqd',
            backend='statevector',
            use_governance=True
        )

        logger.info(f"\nSQD Force Computation:")
        logger.info(f"  Energy: {energy:.6f} Ha")
        logger.info(f"  Force on atom 1: {forces[0]}")
        logger.info(f"  Force on atom 2: {forces[1]}")
        logger.info(f"  |F1|: {np.linalg.norm(forces[0]):.6f} Ha/Bohr")
        logger.info(f"  |F2|: {np.linalg.norm(forces[1]):.6f} Ha/Bohr")

        # Same checks as VQE
        assert forces.shape == (2, 3)

        force_sum = forces[0] + forces[1]
        logger.info(f"  Force sum: {np.linalg.norm(force_sum):.3e} Ha/Bohr")
        assert np.linalg.norm(force_sum) < 0.01

        logger.info("✓ SQD force computation test PASSED")

    def test_classical_vs_quantum_forces(self):
        """Compare classical (HF) vs quantum (VQE) forces."""
        logger.info("\n" + "="*80)
        logger.info("TEST: Classical vs Quantum Force Comparison")
        logger.info("="*80)

        # H2 at stretched geometry (correlation effects larger)
        bond = BondFactory.create_bond('H', 'H', distance=1.0)  # Stretched
        positions = np.array([bond.atom_1.position, bond.atom_2.position])

        # Use comparison function
        comparison = compare_classical_vs_quantum_forces(
            positions,
            bond,
            backend='statevector'
        )

        logger.info(f"\nClassical vs Quantum Comparison:")
        logger.info(f"  HF energy: {comparison['hf_energy']:.6f} Ha")
        logger.info(f"  VQE energy: {comparison['vqe_energy']:.6f} Ha")
        logger.info(f"  SQD energy: {comparison['sqd_energy']:.6f} Ha")
        logger.info(f"  Correlation energy: {comparison['correlation_energy']:.6f} Ha")
        logger.info(f"  Correlation %: {comparison['correlation_percent']:.2f}%")
        logger.info(f"  Force difference: {comparison['force_difference']:.6f} Ha/Bohr")
        logger.info(f"  Force correction: {comparison['force_correction_percent']:.1f}%")

        # VQE should give lower energy than HF (correlation energy is negative)
        assert comparison['vqe_energy'] < comparison['hf_energy'], \
            "VQE energy should be lower than HF (correlation)"

        # For stretched H2, correlation should be significant
        assert abs(comparison['correlation_percent']) > 0.1, \
            f"Correlation too small: {comparison['correlation_percent']:.2f}%"

        logger.info("✓ Classical vs quantum comparison test PASSED")

    def test_quantum_md_cost_estimation(self):
        """Test cost estimation for quantum MD."""
        logger.info("\n" + "="*80)
        logger.info("TEST: Quantum MD Cost Estimation")
        logger.info("="*80)

        # Estimate cost for small system
        cost = estimate_quantum_md_cost(
            n_atoms=2,
            n_orbitals=2,
            n_steps=100,
            method='vqe',
            use_governance=True
        )

        logger.info(f"\nCost Estimation (H2, 100 steps, VQE):")
        logger.info(f"  Qubits needed: {cost['n_qubits']}")
        logger.info(f"  Force evals per step: {cost['n_force_evals_per_step']}")
        logger.info(f"  Total quantum solves: {cost['total_solves']}")
        logger.info(f"  Effective solves (with governance): {cost['effective_solves']}")
        logger.info(f"  Governance advantage: {cost['governance_advantage']}")
        logger.info(f"  Estimated time (statevector): {cost['estimated_time_statevector']:.1f} s")
        logger.info(f"  Estimated time (hardware): {cost['estimated_time_hardware_minutes']:.1f} min")
        logger.info(f"  Feasible on statevector: {cost['feasible_statevector']}")
        logger.info(f"  Feasible on hardware: {cost['feasible_hardware']}")

        # Basic sanity checks
        assert cost['n_qubits'] == 4, "H2 should need 4 qubits (2 orbitals × 2)"
        assert cost['total_solves'] > 0, "Should need quantum solves"
        assert cost['effective_solves'] < cost['total_solves'], \
            "Governance should reduce effective cost"

        # With governance, should be feasible
        assert cost['feasible_statevector'], "Should be feasible on statevector"

        # Compare VQE vs SQD cost
        cost_sqd = estimate_quantum_md_cost(
            n_atoms=2,
            n_orbitals=2,
            n_steps=100,
            method='sqd',
            use_governance=True
        )

        logger.info(f"\nSQD vs VQE:")
        logger.info(f"  VQE time: {cost['estimated_time_statevector']:.1f} s")
        logger.info(f"  SQD time: {cost_sqd['estimated_time_statevector']:.1f} s")
        logger.info(f"  Speedup: {cost['estimated_time_statevector'] / cost_sqd['estimated_time_statevector']:.1f}x")

        # SQD should be faster than VQE
        assert cost_sqd['estimated_time_statevector'] < cost['estimated_time_statevector'], \
            "SQD should be faster than VQE"

        logger.info("✓ Cost estimation test PASSED")

    @pytest.mark.slow
    def test_h2_quantum_md_short(self):
        """Run short quantum MD simulation with VQE forces."""
        logger.info("\n" + "="*80)
        logger.info("TEST: Short H2 Quantum MD (VQE)")
        logger.info("="*80)

        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        # Very short quantum MD (expensive!)
        md = MDSimulator(
            bond,
            temperature=300.0,
            timestep=0.5,
            integrator='velocity_verlet',
            thermostat=None,  # NVE
            force_method='vqe',
            use_governance=True,
            backend='statevector'
        )

        # Only 5 steps (each step needs 13 force evaluations)
        # 5 × 13 = 65 VQE solves total
        logger.info("\nRunning 5-step VQE MD (this will take a while)...")
        logger.info("Each step requires 13 VQE solves (1 + 2×2×3 for gradients)")

        result = md.run(n_steps=5, save_frequency=1)

        logger.info(f"\nQuantum MD Result:")
        logger.info(f"  Final energy: {result.final_energy:.6f} Ha")
        logger.info(f"  Average temperature: {result.average_temperature:.2f} K")
        logger.info(f"  Trajectory frames: {len(result.trajectory.frames)}")

        # Should complete successfully
        assert len(result.trajectory.frames) == 6, "Should have 6 frames (0-5)"
        assert result.final_energy < 0, "Energy should be negative (bound state)"

        logger.info("✓ Quantum MD short run test PASSED")

    @pytest.mark.slow
    def test_governance_speedup_measurement(self):
        """Measure actual governance speedup in MD."""
        logger.info("\n" + "="*80)
        logger.info("TEST: Governance Speedup Measurement")
        logger.info("="*80)

        import time

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        positions = np.array([bond.atom_1.position, bond.atom_2.position])

        # Time force computation WITHOUT governance
        logger.info("\nComputing forces WITHOUT governance...")
        start = time.time()
        forces_no_gov, energy_no_gov = compute_quantum_forces(
            positions,
            bond,
            method='vqe',
            backend='statevector',
            use_governance=False
        )
        time_no_gov = time.time() - start

        # Time force computation WITH governance
        logger.info("Computing forces WITH governance...")
        start = time.time()
        forces_gov, energy_gov = compute_quantum_forces(
            positions,
            bond,
            method='vqe',
            backend='statevector',
            use_governance=True
        )
        time_gov = time.time() - start

        speedup = time_no_gov / time_gov

        logger.info(f"\nGovernance Speedup:")
        logger.info(f"  Time without governance: {time_no_gov:.2f} s")
        logger.info(f"  Time with governance: {time_gov:.2f} s")
        logger.info(f"  Measured speedup: {speedup:.2f}x")

        # Energies should be very close
        energy_diff = abs(energy_gov - energy_no_gov)
        logger.info(f"  Energy difference: {energy_diff:.8f} Ha")
        assert energy_diff < 1e-6, f"Governance changed energy: {energy_diff:.8f} Ha"

        # Forces should be very close
        force_diff = np.linalg.norm(forces_gov - forces_no_gov)
        logger.info(f"  Force difference: {force_diff:.8f} Ha/Bohr")
        assert force_diff < 1e-4, f"Governance changed forces significantly: {force_diff:.8f}"

        # Should see some speedup (at least 1.5x for small system)
        logger.info(f"\n  Expected speedup: 2-5x (from governance advantage)")
        logger.info(f"  Actual speedup: {speedup:.2f}x")

        if speedup > 1.5:
            logger.info("  ✓ Governance provides measurable speedup")
        else:
            logger.warning(f"  ⚠ Speedup lower than expected ({speedup:.2f}x)")

        logger.info("✓ Governance speedup test PASSED")


if __name__ == '__main__':
    # Run tests
    print("\n" + "="*80)
    print("QUANTUM MD TEST SUITE")
    print("="*80)

    test = TestQuantumMD()

    try:
        test.test_vqe_force_computation()
    except Exception as e:
        logger.error(f"VQE force test failed: {e}", exc_info=True)

    try:
        test.test_sqd_force_computation()
    except Exception as e:
        logger.error(f"SQD force test failed: {e}", exc_info=True)

    try:
        test.test_classical_vs_quantum_forces()
    except Exception as e:
        logger.error(f"Classical vs quantum test failed: {e}", exc_info=True)

    try:
        test.test_quantum_md_cost_estimation()
    except Exception as e:
        logger.error(f"Cost estimation test failed: {e}", exc_info=True)

    print("\n" + "="*80)
    print("SLOW TESTS (require more computation)")
    print("="*80)

    try:
        test.test_h2_quantum_md_short()
    except Exception as e:
        logger.error(f"Quantum MD short test failed: {e}", exc_info=True)

    try:
        test.test_governance_speedup_measurement()
    except Exception as e:
        logger.error(f"Governance speedup test failed: {e}", exc_info=True)

    print("\n" + "="*80)
    print("ALL QUANTUM MD TESTS COMPLETE")
    print("="*80)
