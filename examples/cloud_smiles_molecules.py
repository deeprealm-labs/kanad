#!/usr/bin/env python3
"""
Real-Life Cloud Test: SMILES-based Molecule Library Testing
==========================================================

Demonstrates comprehensive SMILES molecule input and cloud execution.

Features:
- Common molecule library (water, benzene, drugs, etc.)
- Batch molecule processing
- Automatic geometry optimization
- Cloud backend comparison
- Result persistence

Use Cases:
- High-throughput virtual screening
- Drug discovery pipelines
- Chemical database quantum analysis
"""

import os
import sys
import logging
from datetime import datetime
import json
import numpy as np
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from kanad.io import SMILESConverter
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.backends.bluequbit import BlueQubitBackend
from kanad.backends.ibm import IBMBackend


class SMILESMoleculeLibraryTester:
    """
    Test suite for SMILES-based molecule input and cloud execution.
    """

    def __init__(self, backend_type='bluequbit', output_dir='results'):
        """
        Initialize tester.

        Args:
            backend_type: 'ibm' or 'bluequbit'
            output_dir: Directory for results
        """
        self.backend_type = backend_type
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        self.converter = SMILESConverter(optimize_geometry=True)
        self.results = []

        logger.info(f"SMILESMoleculeLibraryTester initialized")
        logger.info(f"  Backend: {backend_type}")
        logger.info(f"  Output: {output_dir}")

    def get_molecule_library(self):
        """
        Get comprehensive molecule library with SMILES.

        Organized by category for systematic testing.
        """
        library = {
            'simple_diatomics': {
                'H2': {'smiles': '[H][H]', 'description': 'Hydrogen molecule'},
                'N2': {'smiles': 'N#N', 'description': 'Nitrogen molecule'},
                'O2': {'smiles': 'O=O', 'description': 'Oxygen molecule'},
                'F2': {'smiles': 'F-F', 'description': 'Fluorine molecule'},
            },

            'small_molecules': {
                'H2O': {'smiles': 'O', 'description': 'Water'},
                'NH3': {'smiles': 'N', 'description': 'Ammonia'},
                'CH4': {'smiles': 'C', 'description': 'Methane'},
                'HF': {'smiles': 'F', 'description': 'Hydrogen fluoride'},
                'CO': {'smiles': '[C-]#[O+]', 'description': 'Carbon monoxide'},
                'CO2': {'smiles': 'O=C=O', 'description': 'Carbon dioxide'},
            },

            'organic_molecules': {
                'methanol': {'smiles': 'CO', 'description': 'Methanol'},
                'ethanol': {'smiles': 'CCO', 'description': 'Ethanol'},
                'acetone': {'smiles': 'CC(=O)C', 'description': 'Acetone'},
                'formaldehyde': {'smiles': 'C=O', 'description': 'Formaldehyde'},
                'formic_acid': {'smiles': 'C(=O)O', 'description': 'Formic acid'},
                'acetic_acid': {'smiles': 'CC(=O)O', 'description': 'Acetic acid'},
            },

            'aromatics': {
                'benzene': {'smiles': 'c1ccccc1', 'description': 'Benzene'},
                'phenol': {'smiles': 'c1ccc(cc1)O', 'description': 'Phenol'},
                'toluene': {'smiles': 'Cc1ccccc1', 'description': 'Toluene'},
                'aniline': {'smiles': 'c1ccc(cc1)N', 'description': 'Aniline'},
            },

            'drug_fragments': {
                'aspirin_core': {'smiles': 'c1ccc(cc1)O', 'description': 'Salicylic acid core'},
                'paracetamol': {'smiles': 'CC(=O)Nc1ccc(cc1)O', 'description': 'Paracetamol'},
                'caffeine_core': {'smiles': 'Cn1cnc2c1c(=O)[nH]c(=O)n2C', 'description': 'Caffeine core'},
            },

            'ionic_molecules': {
                'NaCl': {'smiles': '[Na+].[Cl-]', 'description': 'Sodium chloride'},
                'LiF': {'smiles': '[Li+].[F-]', 'description': 'Lithium fluoride'},
                'ammonium': {'smiles': '[NH4+]', 'description': 'Ammonium ion'},
            }
        }

        return library

    def initialize_backend(self):
        """Initialize cloud backend."""
        logger.info("\nInitializing backend...")

        if self.backend_type == 'bluequbit':
            self.backend = BlueQubitBackend(
                device='gpu',
                api_token=os.getenv('BLUE_TOKEN')
            )
            logger.info("✓ BlueQubit GPU backend ready")

        elif self.backend_type == 'ibm':
            self.backend = IBMBackend(
                backend_name='ibmq_qasm_simulator',
                api_token=os.getenv('IBM_API')
            )
            logger.info(f"✓ IBM backend ready: {self.backend.backend.name}")

    def test_molecule(
        self,
        smiles: str,
        name: str,
        description: str,
        basis: str = 'sto-3g',
        max_qubits: int = 20
    ):
        """
        Test a single molecule.

        Args:
            smiles: SMILES string
            name: Molecule name
            description: Description
            basis: Basis set
            max_qubits: Maximum qubits (skip if exceeded)

        Returns:
            Result dictionary
        """
        logger.info(f"\n{'='*80}")
        logger.info(f"Testing: {name}")
        logger.info(f"  Description: {description}")
        logger.info(f"  SMILES: {smiles}")
        logger.info(f"{'='*80}")

        try:
            # Convert SMILES to molecule
            logger.info("Converting SMILES to molecule...")
            molecule = self.converter.smiles_to_molecule(
                smiles,
                charge=0,
                multiplicity=1,
                name=name
            )

            logger.info(f"  Atoms: {len(molecule.atoms)}")
            logger.info(f"  Electrons: {molecule.n_electrons}")

            # Build Hamiltonian
            logger.info(f"Building Hamiltonian ({basis} basis)...")
            representation = LCAORepresentation(molecule)
            hamiltonian = CovalentHamiltonian(molecule, representation, basis_name=basis)
            n_qubits = 2 * hamiltonian.n_orbitals

            logger.info(f"  Orbitals: {hamiltonian.n_orbitals}")
            logger.info(f"  Qubits: {n_qubits}")

            # Check qubit limit
            if n_qubits > max_qubits:
                logger.warning(f"  Too many qubits ({n_qubits} > {max_qubits}), skipping...")
                return {
                    'name': name,
                    'smiles': smiles,
                    'status': 'skipped',
                    'reason': f'too_many_qubits ({n_qubits})'
                }

            # Build ansatz (use HW-efficient for larger systems)
            if n_qubits <= 10:
                logger.info("  Using UCC ansatz...")
                ansatz = UCCAnsatz(
                    n_qubits=n_qubits,
                    n_electrons=molecule.n_electrons,
                    include_singles=True,
                    include_doubles=True
                )
            else:
                logger.info("  Using hardware-efficient ansatz...")
                ansatz = HardwareEfficientAnsatz(
                    n_qubits=n_qubits,
                    n_electrons=molecule.n_electrons,
                    n_layers=2
                )

            logger.info(f"  Parameters: {ansatz.n_parameters}")

            # Build circuit
            mapper = JordanWignerMapper()
            qubit_ham = mapper.map_hamiltonian(hamiltonian)

            initial_params = np.zeros(ansatz.n_parameters)
            circuit = ansatz.build_circuit(initial_params)

            logger.info(f"  Circuit depth: {circuit.depth()}")
            logger.info(f"  Pauli terms: {len(qubit_ham.terms)}")

            # Run on cloud
            logger.info(f"\n  Submitting to {self.backend_type}...")

            if self.backend_type == 'bluequbit':
                result = self.backend.run_circuit(
                    circuit,
                    shots=1024,
                    job_name=f"{name.replace(' ', '_')}"
                )

                if 'statevector' in result:
                    statevector = result['statevector']
                    pauli_sum = qubit_ham.to_pauli_sum()
                    H_matrix = pauli_sum.to_matrix()
                    energy = np.vdot(statevector, H_matrix @ statevector).real
                else:
                    energy = None

            elif self.backend_type == 'ibm':
                observable = mapper.map_hamiltonian_to_observable(hamiltonian)
                result = self.backend.run_batch(
                    circuits=[circuit],
                    observables=[observable],
                    shots=1024
                )
                energy = result['values'][0] if 'values' in result else None

            # Reference energy
            hf_energy = hamiltonian.compute_hf_energy()

            result_dict = {
                'name': name,
                'smiles': smiles,
                'description': description,
                'n_atoms': len(molecule.atoms),
                'n_electrons': molecule.n_electrons,
                'n_qubits': n_qubits,
                'n_parameters': ansatz.n_parameters,
                'circuit_depth': circuit.depth(),
                'basis': basis,
                'hf_energy': float(hf_energy),
                'vqe_energy': float(energy) if energy else None,
                'correlation_energy': float(energy - hf_energy) if energy else None,
                'backend': self.backend_type,
                'status': 'success',
                'timestamp': datetime.now().isoformat()
            }

            logger.info(f"\n  ✓ Results:")
            logger.info(f"    HF Energy:  {hf_energy:.6f} Ha")
            if energy:
                logger.info(f"    VQE Energy: {energy:.6f} Ha")
                logger.info(f"    Correlation: {energy - hf_energy:.6f} Ha")

            return result_dict

        except Exception as e:
            logger.error(f"  ✗ Error: {e}", exc_info=True)
            return {
                'name': name,
                'smiles': smiles,
                'description': description,
                'status': 'failed',
                'error': str(e),
                'timestamp': datetime.now().isoformat()
            }

    def run_category(self, category_name, molecules, max_qubits=20):
        """Run all molecules in a category."""
        logger.info(f"\n{'#'*80}")
        logger.info(f"# Category: {category_name}")
        logger.info(f"{'#'*80}")

        category_results = []

        for mol_name, mol_data in molecules.items():
            result = self.test_molecule(
                smiles=mol_data['smiles'],
                name=mol_name,
                description=mol_data['description'],
                max_qubits=max_qubits
            )

            if result:
                category_results.append(result)
                self.results.append(result)

        return category_results

    def run_all_tests(self, categories=None, max_qubits=20):
        """
        Run tests on molecule library.

        Args:
            categories: List of categories to test (None = all)
            max_qubits: Maximum qubits allowed
        """
        logger.info(f"\n{'#'*80}")
        logger.info(f"# SMILES Molecule Library Testing")
        logger.info(f"# Backend: {self.backend_type.upper()}")
        logger.info(f"# Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"{'#'*80}")

        # Initialize backend
        self.initialize_backend()

        # Get library
        library = self.get_molecule_library()

        # Filter categories
        if categories:
            library = {k: v for k, v in library.items() if k in categories}

        # Run tests
        for category_name, molecules in library.items():
            self.run_category(category_name, molecules, max_qubits=max_qubits)

        # Save and summarize
        self.save_results()
        self.print_summary()

    def save_results(self):
        """Save results to JSON."""
        filename = self.output_dir / f'smiles_library_results_{self.backend_type}.json'

        logger.info(f"\nSaving results to {filename}...")

        with open(filename, 'w') as f:
            json.dump({
                'backend': self.backend_type,
                'timestamp': datetime.now().isoformat(),
                'total_molecules': len(self.results),
                'results': self.results
            }, f, indent=2)

        logger.info(f"✓ Results saved")

    def print_summary(self):
        """Print summary of results."""
        logger.info(f"\n{'#'*80}")
        logger.info(f"# SUMMARY")
        logger.info(f"{'#'*80}\n")

        # Statistics
        total = len(self.results)
        successful = len([r for r in self.results if r['status'] == 'success'])
        failed = len([r for r in self.results if r['status'] == 'failed'])
        skipped = len([r for r in self.results if r['status'] == 'skipped'])

        logger.info(f"Total molecules: {total}")
        logger.info(f"  Successful: {successful}")
        logger.info(f"  Failed: {failed}")
        logger.info(f"  Skipped: {skipped}")

        # Successful results
        if successful > 0:
            logger.info(f"\n{'='*80}")
            logger.info("Successful Calculations:")
            logger.info(f"{'='*80}")

            for r in [r for r in self.results if r['status'] == 'success']:
                logger.info(f"\n  {r['name']}:")
                logger.info(f"    SMILES: {r['smiles']}")
                logger.info(f"    Qubits: {r['n_qubits']}, Depth: {r['circuit_depth']}")
                logger.info(f"    HF: {r['hf_energy']:.6f} Ha")
                if r.get('vqe_energy'):
                    logger.info(f"    VQE: {r['vqe_energy']:.6f} Ha")
                    logger.info(f"    Correlation: {r['correlation_energy']:.6f} Ha")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description='SMILES Molecule Library Testing')
    parser.add_argument(
        '--backend',
        choices=['ibm', 'bluequbit'],
        default='bluequbit',
        help='Cloud backend'
    )
    parser.add_argument(
        '--category',
        nargs='+',
        help='Specific categories to test'
    )
    parser.add_argument(
        '--max-qubits',
        type=int,
        default=20,
        help='Maximum qubits allowed'
    )

    args = parser.parse_args()

    # Check environment
    if args.backend == 'ibm' and not os.getenv('IBM_API'):
        logger.error("IBM_API not set!")
        sys.exit(1)

    if args.backend == 'bluequbit' and not os.getenv('BLUE_TOKEN'):
        logger.error("BLUE_TOKEN not set!")
        sys.exit(1)

    # Run tests
    tester = SMILESMoleculeLibraryTester(backend_type=args.backend)

    if args.category:
        logger.info(f"Testing categories: {', '.join(args.category)}")
        tester.run_all_tests(categories=args.category, max_qubits=args.max_qubits)
    else:
        logger.info("Testing all categories")
        tester.run_all_tests(max_qubits=args.max_qubits)


if __name__ == '__main__':
    main()
