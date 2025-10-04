#!/usr/bin/env python3
"""
H2O VQE Benchmark - Comprehensive Performance Analysis

This script performs VQE simulation of water (H2O) molecule using:
- Qiskit Nature for molecular setup
- Multiple ansatz types (UCC, Hardware Efficient)
- Different optimization strategies
- Performance metrics (time, iterations, accuracy)
"""

import numpy as np
import time
from typing import Dict, Any
import matplotlib.pyplot as plt

# Qiskit imports
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper, ParityMapper
from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
try:
    from qiskit.primitives import StatevectorEstimator as Estimator
    from qiskit_algorithms.minimum_eigensolvers import VQE
    from qiskit_algorithms.optimizers import SLSQP, COBYLA, L_BFGS_B, SPSA
except ImportError:
    # Fallback for older Qiskit versions
    from qiskit.primitives import Estimator
    from qiskit_algorithms import VQE
    from qiskit_algorithms.optimizers import SLSQP, COBYLA, L_BFGS_B, SPSA

# For comparison
from pyscf import gto, scf as pyscf_scf


class H2OVQEBenchmark:
    """Benchmark VQE performance on H2O molecule."""

    def __init__(self, geometry='equilibrium'):
        """
        Initialize H2O molecule.

        Args:
            geometry: 'equilibrium' or 'stretched'
        """
        self.geometry = geometry
        self.results = {}

        # H2O geometries
        if geometry == 'equilibrium':
            # Equilibrium geometry (experimental)
            self.atom_string = "O 0.0000 0.0000 0.1173; H 0.0000 0.7572 -0.4692; H 0.0000 -0.7572 -0.4692"
            self.description = "H2O at equilibrium geometry"
        else:
            # Stretched O-H bonds (1.2x equilibrium)
            self.atom_string = "O 0.0000 0.0000 0.1173; H 0.0000 0.9086 -0.5630; H 0.0000 -0.9086 -0.5630"
            self.description = "H2O with stretched O-H bonds"

        print(f"\n{'='*70}")
        print(f"H2O VQE BENCHMARK")
        print(f"{'='*70}")
        print(f"Geometry: {self.description}")
        print(f"{'='*70}\n")

    def run_hartree_fock_reference(self):
        """Run classical Hartree-Fock for reference."""
        print("Running Hartree-Fock reference calculation...")
        start_time = time.time()

        # Setup molecule with PySCF
        mol = gto.Mole()
        mol.atom = self.atom_string
        mol.basis = 'sto-3g'
        mol.build()

        # Run HF
        mf = pyscf_scf.RHF(mol)
        hf_energy = mf.kernel()

        hf_time = time.time() - start_time

        # Convert to eV
        hartree_to_ev = 27.211386245988
        hf_energy_ev = hf_energy * hartree_to_ev

        self.results['HF'] = {
            'energy': hf_energy,
            'energy_ev': hf_energy_ev,
            'time': hf_time,
            'converged': mf.converged
        }

        print(f"  Energy: {hf_energy:.6f} Ha ({hf_energy_ev:.4f} eV)")
        print(f"  Time: {hf_time:.3f}s")
        print(f"  Converged: {mf.converged}")
        print()

        return hf_energy

    def setup_problem(self):
        """Setup the quantum chemistry problem with Qiskit Nature."""
        print("Setting up quantum problem with Qiskit Nature...")

        # Create driver
        driver = PySCFDriver(
            atom=self.atom_string,
            basis='sto-3g',
            charge=0,
            spin=0
        )

        # Get electronic structure problem
        self.problem = driver.run()

        # Get Hamiltonian (API changed in newer versions)
        try:
            # Newer API
            self.hamiltonian = self.problem.hamiltonian.second_q_op()
            self.num_particles = self.problem.num_particles
            self.num_spatial_orbitals = self.problem.num_spatial_orbitals
        except AttributeError:
            # Older API
            self.second_q_ops = self.problem.second_q_ops()
            self.hamiltonian = self.second_q_ops[0]  # First element
            self.num_particles = self.problem.num_particles
            self.num_spatial_orbitals = self.problem.num_spatial_orbitals

        print(f"  Spatial orbitals: {self.num_spatial_orbitals}")
        print(f"  Electrons: {sum(self.num_particles)}")
        print(f"  Alpha electrons: {self.num_particles[0]}")
        print(f"  Beta electrons: {self.num_particles[1]}")
        print()

    def run_vqe(self,
                mapper_name='jordan_wigner',
                optimizer_name='SLSQP',
                use_uccsd=True):
        """
        Run VQE with specified configuration.

        Args:
            mapper_name: 'jordan_wigner' or 'parity'
            optimizer_name: 'SLSQP', 'COBYLA', 'L_BFGS_B', 'SPSA'
            use_uccsd: Use UCCSD ansatz if True, else HF only
        """
        config_name = f"{mapper_name}_{optimizer_name}_{'UCCSD' if use_uccsd else 'HF'}"
        print(f"Running VQE: {config_name}")

        start_time = time.time()

        # Choose mapper
        if mapper_name == 'jordan_wigner':
            mapper = JordanWignerMapper()
        else:
            mapper = ParityMapper()

        # Map Hamiltonian to qubit operator
        qubit_op = mapper.map(self.hamiltonian)
        num_qubits = qubit_op.num_qubits

        # Create initial state (Hartree-Fock)
        init_state = HartreeFock(
            num_spatial_orbitals=self.num_spatial_orbitals,
            num_particles=self.num_particles,
            qubit_mapper=mapper
        )

        # Create ansatz
        if use_uccsd:
            ansatz = UCCSD(
                num_spatial_orbitals=self.num_spatial_orbitals,
                num_particles=self.num_particles,
                qubit_mapper=mapper,
                initial_state=init_state
            )
        else:
            ansatz = init_state  # Just Hartree-Fock state

        # Choose optimizer
        if optimizer_name == 'SLSQP':
            optimizer = SLSQP(maxiter=1000)
        elif optimizer_name == 'COBYLA':
            optimizer = COBYLA(maxiter=1000)
        elif optimizer_name == 'L_BFGS_B':
            optimizer = L_BFGS_B(maxiter=1000)
        else:  # SPSA
            optimizer = SPSA(maxiter=100)  # Fewer iterations for SPSA

        # Create VQE instance
        estimator = Estimator()
        vqe = VQE(estimator, ansatz, optimizer)

        # Run VQE
        vqe_result = vqe.compute_minimum_eigenvalue(qubit_op)

        vqe_time = time.time() - start_time

        # Extract results
        vqe_energy = vqe_result.eigenvalue.real
        hartree_to_ev = 27.211386245988
        vqe_energy_ev = vqe_energy * hartree_to_ev

        num_parameters = ansatz.num_parameters

        # Calculate error vs HF
        hf_energy = self.results['HF']['energy']
        error_ha = vqe_energy - hf_energy
        error_mha = error_ha * 1000  # milliHartree

        self.results[config_name] = {
            'energy': vqe_energy,
            'energy_ev': vqe_energy_ev,
            'time': vqe_time,
            'optimizer_evals': vqe_result.optimizer_evals,
            'num_qubits': num_qubits,
            'num_parameters': num_parameters,
            'mapper': mapper_name,
            'optimizer': optimizer_name,
            'ansatz': 'UCCSD' if use_uccsd else 'HF',
            'error_vs_hf_ha': error_ha,
            'error_vs_hf_mha': error_mha
        }

        print(f"  Energy: {vqe_energy:.6f} Ha ({vqe_energy_ev:.4f} eV)")
        print(f"  Error vs HF: {error_mha:.3f} mHa")
        print(f"  Qubits: {num_qubits}")
        print(f"  Parameters: {num_parameters}")
        print(f"  Optimizer evaluations: {vqe_result.optimizer_evals}")
        print(f"  Time: {vqe_time:.3f}s")
        print()

        return vqe_result

    def run_comprehensive_benchmark(self):
        """Run complete benchmark suite."""
        # 1. Reference HF
        self.run_hartree_fock_reference()

        # 2. Setup quantum problem
        self.setup_problem()

        # 3. VQE with different configurations
        configs = [
            # Jordan-Wigner mapper
            ('jordan_wigner', 'SLSQP', True),
            ('jordan_wigner', 'COBYLA', True),
            ('jordan_wigner', 'L_BFGS_B', True),

            # Parity mapper
            ('parity', 'SLSQP', True),

            # Quick comparison with HF-only (no UCCSD)
            ('jordan_wigner', 'SLSQP', False),
        ]

        for mapper, optimizer, use_uccsd in configs:
            try:
                self.run_vqe(mapper, optimizer, use_uccsd)
            except Exception as e:
                print(f"  ERROR: {e}\n")

    def print_summary(self):
        """Print summary comparison table."""
        print(f"\n{'='*70}")
        print("BENCHMARK SUMMARY")
        print(f"{'='*70}")
        print(f"{'Method':<30} {'Energy (eV)':<15} {'Time (s)':<12} {'Error (mHa)':<15}")
        print(f"{'-'*70}")

        # HF reference
        hf = self.results['HF']
        print(f"{'Hartree-Fock (Reference)':<30} {hf['energy_ev']:>14.4f} {hf['time']:>11.3f} {'0.000':>14}")

        # VQE results
        for name, result in self.results.items():
            if name == 'HF':
                continue

            method_name = f"VQE ({result['ansatz']}/{result['optimizer']})"
            print(f"{method_name:<30} {result['energy_ev']:>14.4f} {result['time']:>11.3f} {result['error_vs_hf_mha']:>14.3f}")

        print(f"{'='*70}\n")

    def plot_results(self, save_path='h2o_vqe_benchmark.png'):
        """Plot benchmark results."""
        # Extract data for VQE methods only
        vqe_methods = [name for name in self.results.keys() if name != 'HF']

        if not vqe_methods:
            print("No VQE results to plot")
            return

        energies = [self.results[m]['energy_ev'] for m in vqe_methods]
        times = [self.results[m]['time'] for m in vqe_methods]
        errors = [abs(self.results[m]['error_vs_hf_mha']) for m in vqe_methods]

        # Simplify method names for plotting
        method_labels = [
            f"{self.results[m]['ansatz']}\n{self.results[m]['optimizer']}"
            for m in vqe_methods
        ]

        # Create figure
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('H2O VQE Benchmark Results', fontsize=16, fontweight='bold')

        # Plot 1: Energies
        ax = axes[0, 0]
        hf_energy = self.results['HF']['energy_ev']
        ax.axhline(y=hf_energy, color='r', linestyle='--', label='HF Reference', linewidth=2)
        ax.bar(range(len(energies)), energies, color='steelblue', alpha=0.7)
        ax.set_ylabel('Energy (eV)', fontweight='bold')
        ax.set_title('Ground State Energy')
        ax.set_xticks(range(len(method_labels)))
        ax.set_xticklabels(method_labels, rotation=45, ha='right')
        ax.legend()
        ax.grid(axis='y', alpha=0.3)

        # Plot 2: Computation Time
        ax = axes[0, 1]
        ax.bar(range(len(times)), times, color='orange', alpha=0.7)
        ax.set_ylabel('Time (seconds)', fontweight='bold')
        ax.set_title('Computation Time')
        ax.set_xticks(range(len(method_labels)))
        ax.set_xticklabels(method_labels, rotation=45, ha='right')
        ax.grid(axis='y', alpha=0.3)

        # Plot 3: Error vs HF
        ax = axes[1, 0]
        ax.bar(range(len(errors)), errors, color='green', alpha=0.7)
        ax.set_ylabel('|Error| (mHa)', fontweight='bold')
        ax.set_title('Absolute Error vs Hartree-Fock')
        ax.set_xticks(range(len(method_labels)))
        ax.set_xticklabels(method_labels, rotation=45, ha='right')
        ax.grid(axis='y', alpha=0.3)

        # Plot 4: Efficiency (accuracy per time)
        ax = axes[1, 1]
        # Lower error and lower time = better, so plot 1/(error*time)
        efficiency = [1.0 / (e * t + 0.001) for e, t in zip(errors, times)]
        ax.bar(range(len(efficiency)), efficiency, color='purple', alpha=0.7)
        ax.set_ylabel('Efficiency (1/(error×time))', fontweight='bold')
        ax.set_title('Computational Efficiency')
        ax.set_xticks(range(len(method_labels)))
        ax.set_xticklabels(method_labels, rotation=45, ha='right')
        ax.grid(axis='y', alpha=0.3)

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Benchmark plot saved to: {save_path}")
        plt.close()

    def save_results(self, filename='h2o_vqe_results.txt'):
        """Save detailed results to file."""
        with open(filename, 'w') as f:
            f.write("="*70 + "\n")
            f.write("H2O VQE BENCHMARK - DETAILED RESULTS\n")
            f.write("="*70 + "\n")
            f.write(f"Geometry: {self.description}\n")
            f.write("="*70 + "\n\n")

            for name, result in self.results.items():
                f.write(f"{name}\n")
                f.write("-"*70 + "\n")
                for key, value in result.items():
                    f.write(f"  {key}: {value}\n")
                f.write("\n")

        print(f"Detailed results saved to: {filename}")


def main():
    """Run H2O VQE benchmark."""
    # Create benchmark
    benchmark = H2OVQEBenchmark(geometry='equilibrium')

    # Run comprehensive benchmark
    benchmark.run_comprehensive_benchmark()

    # Print summary
    benchmark.print_summary()

    # Save results
    benchmark.save_results('h2o_vqe_results.txt')

    # Plot results
    benchmark.plot_results('h2o_vqe_benchmark.png')

    print("\n" + "="*70)
    print("BENCHMARK COMPLETE!")
    print("="*70)
    print("\nKey Findings:")

    # Find best method by accuracy
    vqe_methods = {k: v for k, v in benchmark.results.items() if k != 'HF'}
    if vqe_methods:
        best_accuracy = min(vqe_methods.items(),
                           key=lambda x: abs(x[1]['error_vs_hf_mha']))
        print(f"  Most Accurate: {best_accuracy[0]}")
        print(f"    Error: {abs(best_accuracy[1]['error_vs_hf_mha']):.3f} mHa")

        # Find fastest method
        fastest = min(vqe_methods.items(), key=lambda x: x[1]['time'])
        print(f"  Fastest: {fastest[0]}")
        print(f"    Time: {fastest[1]['time']:.3f}s")

        # System info
        first_vqe = list(vqe_methods.values())[0]
        print(f"\nSystem Size:")
        print(f"  Qubits: {first_vqe['num_qubits']}")
        print(f"  Parameters (UCCSD): {first_vqe['num_parameters']}")

    print("="*70 + "\n")


if __name__ == "__main__":
    main()
