#!/usr/bin/env python3
"""
Real-Life Cloud Test: Custom Workflow with Multiple Backends
===========================================================

Demonstrates advanced usage:
- Custom Hamiltonian configurations
- Multiple mapper comparisons
- Custom ansatz design
- Parallel execution on multiple backends
- Result aggregation and analysis

Use Case: Comparing different quantum algorithms for the same molecule
"""

import os
import sys
import logging
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import json

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.mappers.bravyi_kitaev_mapper import BravyiKitaevMapper
from kanad.backends.ibm import IBMBackend, IBMPreparation, IBMRunner
from kanad.backends.bluequbit import BlueQubitBackend, BlueQubitRunner
from kanad.io import SMILESConverter


class CustomWorkflowRunner:
    """
    Advanced workflow for comparing quantum algorithms across backends.
    """

    def __init__(self):
        """Initialize workflow runner."""
        self.results = {}
        self.backends = {}

    def setup_backends(self):
        """Setup all available backends."""
        logger.info("Setting up quantum backends...")

        # BlueQubit
        if os.getenv('BLUE_TOKEN'):
            try:
                self.backends['bluequbit_gpu'] = BlueQubitBackend(
                    device='gpu',
                    api_token=os.getenv('BLUE_TOKEN')
                )
                logger.info("✓ BlueQubit GPU backend ready")
            except Exception as e:
                logger.warning(f"BlueQubit GPU failed: {e}")

            try:
                self.backends['bluequbit_cpu'] = BlueQubitBackend(
                    device='cpu',
                    api_token=os.getenv('BLUE_TOKEN')
                )
                logger.info("✓ BlueQubit CPU backend ready")
            except Exception as e:
                logger.warning(f"BlueQubit CPU failed: {e}")

        # IBM
        if os.getenv('IBM_API'):
            try:
                self.backends['ibm_simulator'] = IBMBackend(
                    backend_name='ibmq_qasm_simulator',
                    api_token=os.getenv('IBM_API')
                )
                logger.info("✓ IBM Simulator backend ready")
            except Exception as e:
                logger.warning(f"IBM Simulator failed: {e}")

        if not self.backends:
            logger.error("No backends available! Set IBM_API or BLUE_TOKEN")
            sys.exit(1)

        logger.info(f"\nTotal backends available: {len(self.backends)}")

    def create_molecule_from_smiles(self, smiles, name):
        """Create molecule from SMILES."""
        logger.info(f"\nCreating molecule: {name}")
        logger.info(f"  SMILES: {smiles}")

        converter = SMILESConverter(optimize_geometry=True)
        molecule = converter.smiles_to_molecule(smiles, charge=0, multiplicity=1, name=name)

        logger.info(f"  Atoms: {len(molecule.atoms)}")
        logger.info(f"  Electrons: {molecule.n_electrons}")

        return molecule

    def build_hamiltonians(self, molecule, basis_sets=['sto-3g', '6-31g']):
        """Build Hamiltonians with different basis sets."""
        hamiltonians = {}

        for basis in basis_sets:
            logger.info(f"\nBuilding Hamiltonian: {basis} basis")

            try:
                representation = LCAORepresentation(molecule)
                ham = CovalentHamiltonian(molecule, representation, basis_name=basis)
                hamiltonians[basis] = ham

                logger.info(f"  Orbitals: {ham.n_orbitals}")
                logger.info(f"  Qubits needed: {2 * ham.n_orbitals}")

            except Exception as e:
                logger.warning(f"  Failed to build {basis} Hamiltonian: {e}")

        return hamiltonians

    def build_mappers(self):
        """Build different qubit mappers."""
        mappers = {
            'jordan_wigner': JordanWignerMapper(),
            'bravyi_kitaev': BravyiKitaevMapper(),
        }

        logger.info("\nMappers configured:")
        for name in mappers:
            logger.info(f"  - {name}")

        return mappers

    def build_ansatze(self, n_qubits, n_electrons):
        """Build different ansatze."""
        ansatze = {}

        logger.info(f"\nBuilding ansatze for {n_qubits} qubits, {n_electrons} electrons:")

        # UCC ansatz
        try:
            ansatze['ucc'] = UCCAnsatz(
                n_qubits=n_qubits,
                n_electrons=n_electrons,
                include_singles=True,
                include_doubles=True
            )
            logger.info(f"  UCC: {ansatze['ucc'].n_parameters} parameters")
        except Exception as e:
            logger.warning(f"  UCC failed: {e}")

        # Hardware-efficient (shallow)
        try:
            ansatze['hw_eff_shallow'] = HardwareEfficientAnsatz(
                n_qubits=n_qubits,
                n_electrons=n_electrons,
                n_layers=2
            )
            logger.info(f"  HW-Efficient (2 layers): {ansatze['hw_eff_shallow'].n_parameters} parameters")
        except Exception as e:
            logger.warning(f"  HW-Efficient (shallow) failed: {e}")

        # Hardware-efficient (deep)
        try:
            ansatze['hw_eff_deep'] = HardwareEfficientAnsatz(
                n_qubits=n_qubits,
                n_electrons=n_electrons,
                n_layers=4
            )
            logger.info(f"  HW-Efficient (4 layers): {ansatze['hw_eff_deep'].n_parameters} parameters")
        except Exception as e:
            logger.warning(f"  HW-Efficient (deep) failed: {e}")

        return ansatze

    def run_configuration(
        self,
        config_name: str,
        hamiltonian,
        mapper,
        ansatz,
        backend_name: str,
        backend
    ):
        """
        Run a single configuration (Hamiltonian + Mapper + Ansatz + Backend).
        """
        logger.info(f"\n{'='*80}")
        logger.info(f"Running: {config_name}")
        logger.info(f"  Backend: {backend_name}")
        logger.info(f"{'='*80}")

        try:
            # Build circuit
            initial_params = np.zeros(ansatz.n_parameters)
            circuit = ansatz.build_circuit(initial_params)

            logger.info(f"Circuit built:")
            logger.info(f"  Depth: {circuit.depth()}")
            logger.info(f"  Gates: {len(circuit.data)}")

            # Map Hamiltonian
            qubit_ham = mapper.map_hamiltonian(hamiltonian)
            logger.info(f"  Pauli terms: {len(qubit_ham.terms)}")

            # Run on backend
            start_time = datetime.now()

            if 'bluequbit' in backend_name:
                result = backend.run_circuit(
                    circuit,
                    shots=1024,
                    job_name=config_name.replace(' ', '_')
                )

                # Extract energy
                if 'statevector' in result:
                    statevector = result['statevector']
                    pauli_sum = qubit_ham.to_pauli_sum()
                    H_matrix = pauli_sum.to_matrix()
                    energy = np.vdot(statevector, H_matrix @ statevector).real
                else:
                    energy = None

            elif 'ibm' in backend_name:
                observable = mapper.map_hamiltonian_to_observable(hamiltonian)
                result = backend.run_batch(
                    circuits=[circuit],
                    observables=[observable],
                    shots=1024,
                    optimization_level=1
                )
                energy = result['values'][0] if 'values' in result else None

            else:
                energy = None

            end_time = datetime.now()
            runtime = (end_time - start_time).total_seconds()

            # Get reference energy
            hf_energy = hamiltonian.compute_hf_energy()

            result_dict = {
                'config': config_name,
                'backend': backend_name,
                'n_qubits': circuit.num_qubits,
                'circuit_depth': circuit.depth(),
                'n_parameters': ansatz.n_parameters,
                'hf_energy': hf_energy,
                'vqe_energy': energy,
                'correlation': energy - hf_energy if energy else None,
                'runtime_seconds': runtime,
                'status': 'success',
                'timestamp': datetime.now().isoformat()
            }

            logger.info(f"\nResults:")
            logger.info(f"  HF Energy: {hf_energy:.6f} Ha")
            if energy:
                logger.info(f"  VQE Energy: {energy:.6f} Ha")
                logger.info(f"  Correlation: {energy - hf_energy:.6f} Ha")
            logger.info(f"  Runtime: {runtime:.2f}s")

            return result_dict

        except Exception as e:
            logger.error(f"Configuration failed: {e}", exc_info=True)
            return {
                'config': config_name,
                'backend': backend_name,
                'status': 'failed',
                'error': str(e)
            }

    def run_comprehensive_test(self, smiles, molecule_name):
        """
        Run comprehensive test on a molecule:
        - Multiple basis sets
        - Multiple mappers
        - Multiple ansatze
        - Multiple backends
        """
        logger.info(f"\n{'#'*80}")
        logger.info(f"# Comprehensive Test: {molecule_name}")
        logger.info(f"# SMILES: {smiles}")
        logger.info(f"{'#'*80}")

        # Create molecule
        molecule = self.create_molecule_from_smiles(smiles, molecule_name)

        # Build Hamiltonians
        hamiltonians = self.build_hamiltonians(molecule, basis_sets=['sto-3g'])

        # Build mappers
        mappers = self.build_mappers()

        # Results storage
        all_results = []

        # Run all combinations
        for basis_name, hamiltonian in hamiltonians.items():
            n_qubits = 2 * hamiltonian.n_orbitals

            # Skip if too many qubits
            if n_qubits > 30:
                logger.warning(f"Skipping {basis_name}: too many qubits ({n_qubits})")
                continue

            # Build ansatze
            ansatze = self.build_ansatze(n_qubits, molecule.n_electrons)

            # Run each configuration
            for mapper_name, mapper in mappers.items():
                for ansatz_name, ansatz in ansatze.items():
                    for backend_name, backend in self.backends.items():
                        config_name = f"{molecule_name}_{basis_name}_{mapper_name}_{ansatz_name}"

                        result = self.run_configuration(
                            config_name=config_name,
                            hamiltonian=hamiltonian,
                            mapper=mapper,
                            ansatz=ansatz,
                            backend_name=backend_name,
                            backend=backend
                        )

                        all_results.append(result)

        # Store results
        self.results[molecule_name] = all_results

        return all_results

    def save_results(self, filename='cloud_test_results.json'):
        """Save results to JSON file."""
        logger.info(f"\nSaving results to {filename}...")

        with open(filename, 'w') as f:
            json.dump(self.results, f, indent=2, default=str)

        logger.info(f"✓ Results saved")

    def print_comparative_analysis(self):
        """Print comparative analysis of all results."""
        logger.info(f"\n{'#'*80}")
        logger.info(f"# COMPARATIVE ANALYSIS")
        logger.info(f"{'#'*80}\n")

        for mol_name, results in self.results.items():
            logger.info(f"\n{mol_name}:")
            logger.info("=" * 80)

            # Group by backend
            by_backend = {}
            for r in results:
                if r['status'] == 'success':
                    backend = r['backend']
                    if backend not in by_backend:
                        by_backend[backend] = []
                    by_backend[backend].append(r)

            # Print backend comparison
            for backend, backend_results in by_backend.items():
                logger.info(f"\n  {backend}:")
                logger.info(f"  {'-' * 76}")

                for r in backend_results:
                    config = r['config'].replace(f"{mol_name}_", "")
                    logger.info(f"    {config}:")
                    logger.info(f"      HF: {r['hf_energy']:.6f} Ha")
                    if r.get('vqe_energy'):
                        logger.info(f"      VQE: {r['vqe_energy']:.6f} Ha")
                        logger.info(f"      Correlation: {r['correlation']:.6f} Ha")
                    logger.info(f"      Runtime: {r['runtime_seconds']:.2f}s")


def main():
    """Main entry point."""
    logger.info(f"\n{'#'*80}")
    logger.info(f"# Custom Workflow Cloud Testing")
    logger.info(f"# Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"{'#'*80}\n")

    # Initialize runner
    runner = CustomWorkflowRunner()

    # Setup backends
    runner.setup_backends()

    # Test molecules
    test_molecules = [
        ('O', 'Water'),
        ('[H][H]', 'Hydrogen'),
        ('c1ccccc1', 'Benzene'),  # More complex
    ]

    # Run tests on each molecule
    for smiles, name in test_molecules:
        try:
            runner.run_comprehensive_test(smiles, name)
        except Exception as e:
            logger.error(f"Failed to test {name}: {e}", exc_info=True)

    # Print analysis
    runner.print_comparative_analysis()

    # Save results
    runner.save_results('cloud_workflow_results.json')

    logger.info(f"\n{'#'*80}")
    logger.info(f"# Testing Complete!")
    logger.info(f"# Ended: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"{'#'*80}\n")


if __name__ == '__main__':
    main()
