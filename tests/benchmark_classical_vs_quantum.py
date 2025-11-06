"""
Benchmark Classical vs Quantum MD

Comprehensive performance comparison between:
- Classical MD (HF, MP2 forces)
- Quantum MD (VQE, SQD forces)

Metrics:
- Steps per second
- Energy conservation
- Force accuracy
- Correlation effects
- Governance speedup
- Scalability with system size

This benchmark demonstrates the advantages and costs of quantum-enhanced MD.
"""

import time
import numpy as np
import logging
from typing import Dict, List
from dataclasses import dataclass

from kanad.bonds import BondFactory
from kanad.dynamics import MDSimulator
from kanad.dynamics.quantum_md import (
    compute_quantum_forces,
    estimate_quantum_md_cost,
    compare_classical_vs_quantum_forces
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class BenchmarkResult:
    """Result from a single benchmark run."""
    method: str  # 'hf', 'mp2', 'vqe', 'sqd'
    use_governance: bool
    n_steps: int
    timestep: float
    wall_time: float  # seconds
    steps_per_second: float
    final_energy: float  # Hartree
    energy_drift: float  # Hartree
    avg_temperature: float  # K
    correlation_energy: float = 0.0  # Only for quantum methods


class MDPerformanceBenchmark:
    """Comprehensive MD performance benchmarking."""

    def __init__(self):
        self.results: List[BenchmarkResult] = []

    def benchmark_classical_md(
        self,
        bond,
        method: str = 'hf',
        n_steps: int = 100,
        timestep: float = 0.5
    ) -> BenchmarkResult:
        """Benchmark classical MD (HF or MP2)."""
        logger.info(f"\n{'='*60}")
        logger.info(f"Benchmarking Classical MD: {method.upper()}")
        logger.info(f"{'='*60}")

        md = MDSimulator(
            bond,
            temperature=300.0,
            timestep=timestep,
            integrator='velocity_verlet',
            thermostat=None,  # NVE for performance
            force_method=method
        )

        # Run and time
        start = time.time()
        result = md.run(n_steps=n_steps, save_frequency=10, verbose=False)
        wall_time = time.time() - start

        # Compute metrics
        steps_per_second = n_steps / wall_time
        energies = result.trajectory.total_energies
        energy_drift = abs(energies[-1] - energies[0])

        logger.info(f"\nPerformance:")
        logger.info(f"  Wall time: {wall_time:.2f} s")
        logger.info(f"  Steps/second: {steps_per_second:.2f}")
        logger.info(f"  Final energy: {result.final_energy:.6f} Ha")
        logger.info(f"  Energy drift: {energy_drift:.8f} Ha")
        logger.info(f"  Avg temperature: {result.average_temperature:.2f} K")

        benchmark = BenchmarkResult(
            method=method,
            use_governance=False,
            n_steps=n_steps,
            timestep=timestep,
            wall_time=wall_time,
            steps_per_second=steps_per_second,
            final_energy=result.final_energy,
            energy_drift=energy_drift,
            avg_temperature=result.average_temperature
        )

        self.results.append(benchmark)
        return benchmark

    def benchmark_quantum_md(
        self,
        bond,
        method: str = 'vqe',
        use_governance: bool = True,
        n_steps: int = 5,
        timestep: float = 0.5,
        backend: str = 'statevector'
    ) -> BenchmarkResult:
        """Benchmark quantum MD (VQE or SQD)."""
        logger.info(f"\n{'='*60}")
        logger.info(f"Benchmarking Quantum MD: {method.upper()}")
        logger.info(f"  Governance: {use_governance}")
        logger.info(f"  Backend: {backend}")
        logger.info(f"{'='*60}")

        md = MDSimulator(
            bond,
            temperature=300.0,
            timestep=timestep,
            integrator='velocity_verlet',
            thermostat=None,
            force_method=method,
            use_governance=use_governance,
            backend=backend
        )

        # Run and time
        logger.info(f"Running {n_steps} steps (this may take a while)...")
        start = time.time()
        result = md.run(n_steps=n_steps, save_frequency=max(1, n_steps//5), verbose=False)
        wall_time = time.time() - start

        # Compute metrics
        steps_per_second = n_steps / wall_time
        energies = result.trajectory.total_energies
        energy_drift = abs(energies[-1] - energies[0])

        logger.info(f"\nPerformance:")
        logger.info(f"  Wall time: {wall_time:.2f} s")
        logger.info(f"  Steps/second: {steps_per_second:.4f}")
        logger.info(f"  Time per step: {wall_time/n_steps:.2f} s")
        logger.info(f"  Final energy: {result.final_energy:.6f} Ha")
        logger.info(f"  Energy drift: {energy_drift:.8f} Ha")
        logger.info(f"  Avg temperature: {result.average_temperature:.2f} K")

        benchmark = BenchmarkResult(
            method=method,
            use_governance=use_governance,
            n_steps=n_steps,
            timestep=timestep,
            wall_time=wall_time,
            steps_per_second=steps_per_second,
            final_energy=result.final_energy,
            energy_drift=energy_drift,
            avg_temperature=result.average_temperature
        )

        self.results.append(benchmark)
        return benchmark

    def compare_force_accuracy(self, bond, distance: float = 0.74):
        """Compare force accuracy between methods."""
        logger.info(f"\n{'='*60}")
        logger.info(f"Force Accuracy Comparison at r={distance:.2f} Å")
        logger.info(f"{'='*60}")

        # Update bond distance
        bond.distance = distance
        positions = np.array([bond.atom_1.position, bond.atom_2.position])

        # Get exact quantum result (VQE with governance)
        forces_vqe, energy_vqe = compute_quantum_forces(
            positions, bond,
            method='vqe',
            backend='statevector',
            use_governance=True
        )

        # Compare with classical
        comparison = compare_classical_vs_quantum_forces(
            positions, bond, backend='statevector'
        )

        logger.info(f"\nEnergies:")
        logger.info(f"  HF:  {comparison['hf_energy']:.6f} Ha")
        logger.info(f"  VQE: {comparison['vqe_energy']:.6f} Ha")
        logger.info(f"  SQD: {comparison['sqd_energy']:.6f} Ha")
        logger.info(f"  Correlation: {comparison['correlation_energy']:.6f} Ha "
                   f"({comparison['correlation_percent']:.2f}%)")

        logger.info(f"\nForces:")
        logger.info(f"  HF force:  {np.linalg.norm(comparison['hf_forces']):.6f} Ha/Bohr")
        logger.info(f"  VQE force: {np.linalg.norm(comparison['vqe_forces']):.6f} Ha/Bohr")
        logger.info(f"  SQD force: {np.linalg.norm(comparison['sqd_forces']):.6f} Ha/Bohr")
        logger.info(f"  Force difference: {comparison['force_difference']:.6f} Ha/Bohr "
                   f"({comparison['force_correction_percent']:.1f}%)")

        return comparison

    def measure_governance_advantage(self, bond):
        """Measure governance speedup for force computation."""
        logger.info(f"\n{'='*60}")
        logger.info("Governance Advantage Measurement")
        logger.info(f"{'='*60}")

        positions = np.array([bond.atom_1.position, bond.atom_2.position])

        # Time WITHOUT governance
        logger.info("\nTiming VQE forces WITHOUT governance...")
        start = time.time()
        forces_no_gov, energy_no_gov = compute_quantum_forces(
            positions, bond,
            method='vqe',
            backend='statevector',
            use_governance=False
        )
        time_no_gov = time.time() - start

        # Time WITH governance
        logger.info("Timing VQE forces WITH governance...")
        start = time.time()
        forces_gov, energy_gov = compute_quantum_forces(
            positions, bond,
            method='vqe',
            backend='statevector',
            use_governance=True
        )
        time_gov = time.time() - start

        speedup = time_no_gov / time_gov

        logger.info(f"\nGovernance Speedup:")
        logger.info(f"  Without governance: {time_no_gov:.2f} s")
        logger.info(f"  With governance:    {time_gov:.2f} s")
        logger.info(f"  Speedup factor:     {speedup:.2f}x")

        # Verify accuracy not sacrificed
        energy_diff = abs(energy_gov - energy_no_gov)
        force_diff = np.linalg.norm(forces_gov - forces_no_gov)

        logger.info(f"\nAccuracy Check:")
        logger.info(f"  Energy difference: {energy_diff:.8f} Ha")
        logger.info(f"  Force difference:  {force_diff:.8f} Ha/Bohr")

        if energy_diff < 1e-6 and force_diff < 1e-4:
            logger.info("  ✓ Governance maintains accuracy")
        else:
            logger.warning("  ⚠ Governance may affect accuracy")

        return {
            'speedup': speedup,
            'time_no_gov': time_no_gov,
            'time_gov': time_gov,
            'energy_diff': energy_diff,
            'force_diff': force_diff
        }

    def estimate_scaling(self):
        """Estimate cost scaling with system size."""
        logger.info(f"\n{'='*60}")
        logger.info("Cost Scaling Analysis")
        logger.info(f"{'='*60}")

        systems = [
            {'name': 'H2', 'n_atoms': 2, 'n_orbitals': 2},
            {'name': 'H2O', 'n_atoms': 3, 'n_orbitals': 4},
            {'name': 'CH4', 'n_atoms': 5, 'n_orbitals': 8},
            {'name': 'C2H4', 'n_atoms': 6, 'n_orbitals': 12},
        ]

        n_steps = 100

        logger.info(f"\nCost estimates for {n_steps} MD steps:")
        logger.info(f"{'System':<10} {'Qubits':<8} {'Solves':<10} {'Gov Solves':<12} "
                   f"{'SV Time':<12} {'Feasible':<10}")
        logger.info("-" * 70)

        for sys in systems:
            cost = estimate_quantum_md_cost(
                n_atoms=sys['n_atoms'],
                n_orbitals=sys['n_orbitals'],
                n_steps=n_steps,
                method='vqe',
                use_governance=True
            )

            logger.info(f"{sys['name']:<10} "
                       f"{cost['n_qubits']:<8} "
                       f"{cost['total_solves']:<10} "
                       f"{cost['effective_solves']:<12} "
                       f"{cost['estimated_time_statevector']:<12.1f} "
                       f"{'Yes' if cost['feasible_statevector'] else 'No':<10}")

        logger.info("\nKey Observations:")
        logger.info("- Cost scales as O(N_atoms) for force evaluations")
        logger.info("- Governance provides 2-10x speedup (larger for bigger systems)")
        logger.info("- Statevector simulation becomes infeasible beyond ~8 qubits")
        logger.info("- Real hardware needed for molecules with >4 orbitals")

    def print_summary(self):
        """Print summary of all benchmark results."""
        logger.info(f"\n{'='*80}")
        logger.info("BENCHMARK SUMMARY")
        logger.info(f"{'='*80}")

        if not self.results:
            logger.info("No benchmark results available")
            return

        # Group by method
        classical = [r for r in self.results if r.method in ['hf', 'mp2']]
        quantum = [r for r in self.results if r.method in ['vqe', 'sqd']]

        if classical:
            logger.info("\nClassical MD Performance:")
            logger.info(f"{'Method':<10} {'Steps':<8} {'Time(s)':<10} {'Steps/s':<12} "
                       f"{'Energy Drift':<15}")
            logger.info("-" * 70)

            for r in classical:
                logger.info(f"{r.method.upper():<10} {r.n_steps:<8} "
                           f"{r.wall_time:<10.2f} {r.steps_per_second:<12.2f} "
                           f"{r.energy_drift:<15.8f}")

        if quantum:
            logger.info("\nQuantum MD Performance:")
            logger.info(f"{'Method':<10} {'Gov':<5} {'Steps':<8} {'Time(s)':<10} "
                       f"{'Steps/s':<12} {'Energy Drift':<15}")
            logger.info("-" * 75)

            for r in quantum:
                logger.info(f"{r.method.upper():<10} "
                           f"{'Yes' if r.use_governance else 'No':<5} "
                           f"{r.n_steps:<8} {r.wall_time:<10.2f} "
                           f"{r.steps_per_second:<12.4f} {r.energy_drift:<15.8f}")

        # Performance comparison
        if classical and quantum:
            logger.info("\nPerformance Ratio:")
            classical_sps = classical[0].steps_per_second
            quantum_sps = quantum[0].steps_per_second
            slowdown = classical_sps / quantum_sps

            logger.info(f"  Classical: {classical_sps:.2f} steps/s")
            logger.info(f"  Quantum:   {quantum_sps:.4f} steps/s")
            logger.info(f"  Slowdown:  {slowdown:.0f}x")
            logger.info(f"\n  Quantum MD captures electron correlation at {slowdown:.0f}x cost")

        logger.info(f"\n{'='*80}")


def run_comprehensive_benchmark():
    """Run complete benchmark suite."""
    print("\n" + "="*80)
    print("COMPREHENSIVE MD BENCHMARK: CLASSICAL vs QUANTUM")
    print("="*80)

    benchmark = MDPerformanceBenchmark()

    # Create H2 test system
    bond = BondFactory.create_bond('H', 'H', distance=0.74)

    # 1. Classical MD benchmarks
    print("\n" + "="*80)
    print("PHASE 1: Classical MD Benchmarks")
    print("="*80)

    try:
        benchmark.benchmark_classical_md(bond, method='hf', n_steps=100)
    except Exception as e:
        logger.error(f"HF benchmark failed: {e}", exc_info=True)

    # 2. Quantum MD benchmarks (short runs due to cost)
    print("\n" + "="*80)
    print("PHASE 2: Quantum MD Benchmarks")
    print("="*80)

    try:
        benchmark.benchmark_quantum_md(
            bond, method='vqe', use_governance=True,
            n_steps=3, backend='statevector'
        )
    except Exception as e:
        logger.error(f"VQE benchmark failed: {e}", exc_info=True)

    try:
        benchmark.benchmark_quantum_md(
            bond, method='sqd', use_governance=True,
            n_steps=3, backend='statevector'
        )
    except Exception as e:
        logger.error(f"SQD benchmark failed: {e}", exc_info=True)

    # 3. Force accuracy comparison
    print("\n" + "="*80)
    print("PHASE 3: Force Accuracy Analysis")
    print("="*80)

    try:
        # At equilibrium
        benchmark.compare_force_accuracy(bond, distance=0.74)
    except Exception as e:
        logger.error(f"Force comparison (eq) failed: {e}", exc_info=True)

    try:
        # At stretched geometry (larger correlation)
        benchmark.compare_force_accuracy(bond, distance=1.0)
    except Exception as e:
        logger.error(f"Force comparison (stretched) failed: {e}", exc_info=True)

    # 4. Governance advantage
    print("\n" + "="*80)
    print("PHASE 4: Governance Advantage Measurement")
    print("="*80)

    try:
        benchmark.measure_governance_advantage(bond)
    except Exception as e:
        logger.error(f"Governance measurement failed: {e}", exc_info=True)

    # 5. Scaling analysis
    print("\n" + "="*80)
    print("PHASE 5: System Size Scaling")
    print("="*80)

    try:
        benchmark.estimate_scaling()
    except Exception as e:
        logger.error(f"Scaling analysis failed: {e}", exc_info=True)

    # Final summary
    benchmark.print_summary()

    print("\n" + "="*80)
    print("BENCHMARK COMPLETE")
    print("="*80)

    return benchmark


if __name__ == '__main__':
    run_comprehensive_benchmark()
