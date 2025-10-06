#!/usr/bin/env python3
"""
Real-Life Cloud Test: Drug Molecule VQE on IBM Quantum & BlueQubit
=================================================================

Tests VQE energy calculations for real drug molecules using cloud backends.

Molecules tested:
- Aspirin fragment (salicylic acid)
- Caffeine fragment (methylxanthine core)
- Ibuprofen fragment
- Paracetamol fragment

Backends:
- IBM Quantum (real hardware)
- BlueQubit GPU simulator
"""

import os
import sys
import logging
from datetime import datetime
import numpy as np

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Add kanad to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from kanad.io import SMILESConverter
from kanad.bonds import BondFactory
from kanad.solvers import VQESolver
from kanad.backends.ibm import IBMBackend
from kanad.backends.bluequbit import BlueQubitBackend
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper


class DrugMoleculeVQETester:
    """
    Test VQE on real drug molecule fragments using cloud backends.
    """

    def __init__(self, backend_type='bluequbit'):
        """
        Initialize tester.

        Args:
            backend_type: 'ibm' or 'bluequbit'
        """
        self.backend_type = backend_type
        self.results = []

        # Initialize SMILES converter
        self.converter = SMILESConverter(optimize_geometry=True)

        # Drug molecule fragments (small enough for current hardware)
        self.drug_molecules = {
            'aspirin_fragment': {
                'smiles': 'c1ccc(cc1)O',  # Phenol (aspirin core)
                'name': 'Aspirin Fragment (Phenol)',
                'description': 'Benzene ring with hydroxyl group'
            },
            'caffeine_fragment': {
                'smiles': 'Cn1cnc2c1c(=O)[nH]c(=O)n2C',  # Methylxanthine core
                'name': 'Caffeine Fragment',
                'description': 'Xanthine core structure'
            },
            'ibuprofen_fragment': {
                'smiles': 'CC(C)Cc1ccc(cc1)C',  # Isobutylbenzene
                'name': 'Ibuprofen Fragment',
                'description': 'Aromatic core with isobutyl group'
            },
            'paracetamol_fragment': {
                'smiles': 'CC(=O)Nc1ccc(cc1)O',  # Para-acetaminophen
                'name': 'Paracetamol',
                'description': 'Full paracetamol structure'
            },
            'simple_aromatic': {
                'smiles': 'c1ccccc1',  # Benzene
                'name': 'Benzene',
                'description': 'Simple aromatic ring'
            }
        }

        logger.info(f"DrugMoleculeVQETester initialized with {backend_type} backend")

    def initialize_backend(self):
        """Initialize cloud backend."""
        if self.backend_type == 'ibm':
            logger.info("Initializing IBM Quantum backend...")
            # Use simulator for testing (switch to real hardware for production)
            self.backend = IBMBackend(
                backend_name='ibmq_qasm_simulator',  # or 'ibm_brisbane' for real hardware
                api_token=os.getenv('IBM_API')
            )
            logger.info(f"IBM backend ready: {self.backend.backend.name}")

        elif self.backend_type == 'bluequbit':
            logger.info("Initializing BlueQubit backend...")
            self.backend = BlueQubitBackend(
                device='gpu',  # Use GPU for faster simulation
                api_token=os.getenv('BLUE_TOKEN')
            )
            logger.info("BlueQubit GPU backend ready")

        else:
            raise ValueError(f"Unknown backend: {self.backend_type}")

    def run_vqe_on_molecule(
        self,
        smiles: str,
        name: str,
        basis: str = 'sto-3g',
        ansatz_type: str = 'ucc',
        max_iterations: int = 100
    ):
        """
        Run VQE on a molecule.

        Args:
            smiles: SMILES string
            name: Molecule name
            basis: Basis set
            ansatz_type: Ansatz type ('ucc' or 'hardware_efficient')
            max_iterations: Max optimization iterations

        Returns:
            Result dictionary
        """
        logger.info(f"\n{'='*80}")
        logger.info(f"Running VQE on: {name}")
        logger.info(f"SMILES: {smiles}")
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
            logger.info(f"Building Hamiltonian with {basis} basis...")
            representation = LCAORepresentation(molecule)
            hamiltonian = CovalentHamiltonian(
                molecule=molecule,
                representation=representation,
                basis_name=basis
            )
            n_orbitals = hamiltonian.n_orbitals
            logger.info(f"  Orbitals: {n_orbitals}")
            logger.info(f"  Qubits needed: {2 * n_orbitals}")

            # Check qubit requirements
            n_qubits = 2 * n_orbitals
            if n_qubits > 30:
                logger.warning(f"Too many qubits ({n_qubits}), skipping...")
                return None

            # Build ansatz
            logger.info(f"Building {ansatz_type} ansatz...")
            if ansatz_type == 'ucc':
                ansatz = UCCAnsatz(
                    n_qubits=n_qubits,
                    n_electrons=molecule.n_electrons,
                    include_singles=True,
                    include_doubles=True
                )
            else:
                from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz
                ansatz = HardwareEfficientAnsatz(
                    n_qubits=n_qubits,
                    n_electrons=molecule.n_electrons,
                    n_layers=2
                )

            logger.info(f"  Parameters: {ansatz.n_parameters}")

            # Build mapper
            logger.info("Building Jordan-Wigner mapper...")
            mapper = JordanWignerMapper()

            # Map Hamiltonian to qubits
            logger.info("Mapping Hamiltonian to qubit operators...")
            qubit_hamiltonian = mapper.map_hamiltonian(hamiltonian)
            logger.info(f"  Pauli terms: {len(qubit_hamiltonian.terms)}")

            # Build circuit
            logger.info("Building quantum circuit...")
            initial_params = np.zeros(ansatz.n_parameters)
            circuit = ansatz.build_circuit(initial_params)
            logger.info(f"  Circuit depth: {circuit.depth()}")
            logger.info(f"  Circuit gates: {len(circuit.data)}")

            # Run on cloud backend
            logger.info(f"\nSubmitting to {self.backend_type} backend...")

            if self.backend_type == 'bluequbit':
                # BlueQubit execution
                result = self.backend.run_circuit(
                    circuit,
                    shots=1024,
                    job_name=f"vqe_{name.replace(' ', '_')}"
                )
                logger.info("Job completed successfully!")

                # Extract statevector if available
                if 'statevector' in result:
                    statevector = result['statevector']
                    # Compute energy expectation
                    energy = self._compute_energy_from_statevector(
                        statevector,
                        qubit_hamiltonian
                    )
                else:
                    energy = None

            elif self.backend_type == 'ibm':
                # IBM execution
                observable = mapper.map_hamiltonian_to_observable(hamiltonian)
                result = self.backend.run_batch(
                    circuits=[circuit],
                    observables=[observable],
                    shots=1024,
                    optimization_level=1,
                    resilience_level=1
                )
                logger.info(f"Job completed! Job ID: {result['job_id']}")

                # Extract energy
                energy = result['values'][0] if 'values' in result else None

            # Get HF reference for comparison
            hf_energy = hamiltonian.compute_hf_energy()

            result_dict = {
                'name': name,
                'smiles': smiles,
                'n_atoms': len(molecule.atoms),
                'n_electrons': molecule.n_electrons,
                'n_qubits': n_qubits,
                'n_parameters': ansatz.n_parameters,
                'circuit_depth': circuit.depth(),
                'hf_energy': hf_energy,
                'vqe_energy': energy,
                'correlation_energy': energy - hf_energy if energy else None,
                'backend': self.backend_type,
                'status': 'success'
            }

            logger.info(f"\n{'='*80}")
            logger.info(f"RESULTS for {name}:")
            logger.info(f"  HF Energy:          {hf_energy:.6f} Hartree")
            if energy:
                logger.info(f"  VQE Energy:         {energy:.6f} Hartree")
                logger.info(f"  Correlation:        {energy - hf_energy:.6f} Hartree")
            logger.info(f"{'='*80}\n")

            return result_dict

        except Exception as e:
            logger.error(f"Error running VQE on {name}: {e}", exc_info=True)
            return {
                'name': name,
                'smiles': smiles,
                'status': 'failed',
                'error': str(e)
            }

    def _compute_energy_from_statevector(self, statevector, hamiltonian):
        """Compute energy expectation from statevector."""
        # Convert Hamiltonian to matrix
        pauli_sum = hamiltonian.to_pauli_sum()
        H_matrix = pauli_sum.to_matrix()

        # Compute <ψ|H|ψ>
        energy = np.vdot(statevector, H_matrix @ statevector).real

        return energy

    def run_all_tests(self):
        """Run VQE on all drug molecules."""
        logger.info(f"\n{'#'*80}")
        logger.info(f"# Drug Molecule VQE Testing on {self.backend_type.upper()} Backend")
        logger.info(f"# Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"{'#'*80}\n")

        # Initialize backend
        self.initialize_backend()

        # Test each molecule
        for mol_id, mol_data in self.drug_molecules.items():
            result = self.run_vqe_on_molecule(
                smiles=mol_data['smiles'],
                name=mol_data['name'],
                basis='sto-3g',
                ansatz_type='ucc',
                max_iterations=100
            )

            if result:
                self.results.append(result)

        # Print summary
        self.print_summary()

    def print_summary(self):
        """Print summary of all results."""
        logger.info(f"\n{'#'*80}")
        logger.info(f"# SUMMARY OF RESULTS")
        logger.info(f"{'#'*80}\n")

        successful = [r for r in self.results if r['status'] == 'success']
        failed = [r for r in self.results if r['status'] == 'failed']

        logger.info(f"Total molecules tested: {len(self.results)}")
        logger.info(f"Successful: {len(successful)}")
        logger.info(f"Failed: {len(failed)}\n")

        if successful:
            logger.info("Successful Results:")
            logger.info("-" * 80)
            for r in successful:
                logger.info(f"\n{r['name']}:")
                logger.info(f"  Qubits: {r['n_qubits']}")
                logger.info(f"  Circuit depth: {r['circuit_depth']}")
                logger.info(f"  HF Energy: {r['hf_energy']:.6f} Ha")
                if r.get('vqe_energy'):
                    logger.info(f"  VQE Energy: {r['vqe_energy']:.6f} Ha")
                    logger.info(f"  Correlation: {r['correlation_energy']:.6f} Ha")

        if failed:
            logger.info("\n\nFailed Molecules:")
            logger.info("-" * 80)
            for r in failed:
                logger.info(f"\n{r['name']}: {r.get('error', 'Unknown error')}")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description='Drug Molecule VQE Cloud Testing')
    parser.add_argument(
        '--backend',
        choices=['ibm', 'bluequbit'],
        default='bluequbit',
        help='Cloud backend to use'
    )

    args = parser.parse_args()

    # Check environment variables
    if args.backend == 'ibm' and not os.getenv('IBM_API'):
        logger.error("IBM_API environment variable not set!")
        logger.error("Set it with: export IBM_API=your_token")
        sys.exit(1)

    if args.backend == 'bluequbit' and not os.getenv('BLUE_TOKEN'):
        logger.error("BLUE_TOKEN environment variable not set!")
        logger.error("Set it with: export BLUE_TOKEN=your_token")
        sys.exit(1)

    # Run tests
    tester = DrugMoleculeVQETester(backend_type=args.backend)
    tester.run_all_tests()


if __name__ == '__main__':
    main()
