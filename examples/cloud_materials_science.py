#!/usr/bin/env python3
"""
Real-Life Cloud Test: Materials Science - Battery & Catalyst Materials
=====================================================================

Tests quantum simulations of materials used in:
- Battery electrodes (Li-ion batteries)
- Catalysts (transition metal complexes)
- Semiconductors (band gap calculations)

Uses custom mapper and Hamiltonian configurations for materials.

Backends:
- IBM Quantum
- BlueQubit
"""

import os
import sys
import logging
from datetime import datetime
import numpy as np

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
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.mappers.bravyi_kitaev_mapper import BravyiKitaevMapper
from kanad.backends.ibm import IBMBackend
from kanad.backends.bluequbit import BlueQubitBackend
from kanad.solvers.vqe_solver import VQESolver


class MaterialsScienceTester:
    """
    Test quantum simulations of real materials for energy applications.
    """

    def __init__(self, backend_type='bluequbit'):
        """
        Initialize materials tester.

        Args:
            backend_type: 'ibm' or 'bluequbit'
        """
        self.backend_type = backend_type
        self.results = []

        logger.info(f"MaterialsScienceTester initialized with {backend_type} backend")

    def initialize_backend(self):
        """Initialize cloud backend."""
        if self.backend_type == 'ibm':
            logger.info("Initializing IBM Quantum backend...")
            self.backend = IBMBackend(
                backend_name='ibmq_qasm_simulator',
                api_token=os.getenv('IBM_API')
            )
            logger.info(f"IBM backend ready: {self.backend.backend.name}")

        elif self.backend_type == 'bluequbit':
            logger.info("Initializing BlueQubit backend...")
            self.backend = BlueQubitBackend(
                device='gpu',
                api_token=os.getenv('BLUE_TOKEN')
            )
            logger.info("BlueQubit GPU backend ready")

    def test_lithium_hydride_battery_material(self):
        """
        Test LiH - a model system for Li-ion battery electrode materials.

        LiH represents the lithium-hydrogen interaction important in
        solid-state batteries.
        """
        logger.info(f"\n{'='*80}")
        logger.info("TEST 1: Lithium Hydride (Battery Material)")
        logger.info(f"{'='*80}")

        try:
            # Create LiH molecule
            li = Atom('Li', [0.0, 0.0, 0.0])
            h = Atom('H', [0.0, 0.0, 1.5])  # Equilibrium distance
            lih = Molecule([li, h], charge=0, spin=0)

            logger.info(f"Molecule: LiH")
            logger.info(f"  Distance: 1.5 Angstrom")
            logger.info(f"  Electrons: {lih.n_electrons}")

            # Build Hamiltonian with different basis sets
            for basis in ['sto-3g', '6-31g']:
                logger.info(f"\nTesting with {basis} basis...")

                representation = LCAORepresentation(lih)
                hamiltonian = CovalentHamiltonian(lih, representation, basis_name=basis)
                n_qubits = 2 * hamiltonian.n_orbitals

                logger.info(f"  Orbitals: {hamiltonian.n_orbitals}")
                logger.info(f"  Qubits: {n_qubits}")

                # Test with different mappers
                for mapper_type in ['jordan_wigner', 'bravyi_kitaev']:
                    logger.info(f"\n  Testing {mapper_type} mapper...")

                    if mapper_type == 'jordan_wigner':
                        mapper = JordanWignerMapper()
                    else:
                        mapper = BravyiKitaevMapper()

                    # Build ansatz
                    ansatz = UCCAnsatz(
                        n_qubits=n_qubits,
                        n_electrons=lih.n_electrons,
                        include_singles=True,
                        include_doubles=True
                    )

                    logger.info(f"    Parameters: {ansatz.n_parameters}")

                    # Build and map
                    qubit_ham = mapper.map_hamiltonian(hamiltonian)
                    circuit = ansatz.build_circuit(np.zeros(ansatz.n_parameters))

                    logger.info(f"    Circuit depth: {circuit.depth()}")

                    # Run on cloud
                    if self.backend_type == 'bluequbit':
                        result = self.backend.run_circuit(circuit, shots=1024)
                        energy = self._compute_energy_from_result(result, qubit_ham)
                    else:
                        observable = mapper.map_hamiltonian_to_observable(hamiltonian)
                        result = self.backend.run_batch([circuit], [observable], shots=1024)
                        energy = result['values'][0]

                    hf_energy = hamiltonian.compute_hf_energy()

                    self.results.append({
                        'material': 'LiH (Battery)',
                        'basis': basis,
                        'mapper': mapper_type,
                        'n_qubits': n_qubits,
                        'hf_energy': hf_energy,
                        'vqe_energy': energy,
                        'correlation': energy - hf_energy if energy else None
                    })

                    logger.info(f"    HF Energy: {hf_energy:.6f} Ha")
                    if energy:
                        logger.info(f"    VQE Energy: {energy:.6f} Ha")
                        logger.info(f"    Correlation: {energy - hf_energy:.6f} Ha")

            return True

        except Exception as e:
            logger.error(f"Error testing LiH: {e}", exc_info=True)
            return False

    def test_transition_metal_catalyst(self):
        """
        Test FeH - model for Fe-based catalysts in Fischer-Tropsch synthesis.

        Iron hydrides are intermediates in catalytic processes.
        """
        logger.info(f"\n{'='*80}")
        logger.info("TEST 2: Iron Hydride (Catalyst Material)")
        logger.info(f"{'='*80}")

        try:
            # Create FeH molecule
            fe = Atom('Fe', [0.0, 0.0, 0.0])
            h = Atom('H', [0.0, 0.0, 1.6])
            feh = Molecule([fe, h], charge=0, spin=3)  # Quartet state

            logger.info(f"Molecule: FeH")
            logger.info(f"  Spin state: Quartet (S=3/2)")
            logger.info(f"  Electrons: {feh.n_electrons}")

            # Use minimal basis due to Fe complexity
            basis = 'sto-3g'
            logger.info(f"\nUsing {basis} basis...")

            representation = LCAORepresentation(feh)
            hamiltonian = CovalentHamiltonian(feh, representation, basis_name=basis)
            n_qubits = 2 * hamiltonian.n_orbitals

            logger.info(f"  Orbitals: {hamiltonian.n_orbitals}")
            logger.info(f"  Qubits: {n_qubits}")

            # Use hardware-efficient ansatz for large systems
            ansatz = HardwareEfficientAnsatz(
                n_qubits=n_qubits,
                n_electrons=feh.n_electrons,
                n_layers=2
            )

            mapper = JordanWignerMapper()

            logger.info(f"  Parameters: {ansatz.n_parameters}")

            # Build circuit
            qubit_ham = mapper.map_hamiltonian(hamiltonian)
            circuit = ansatz.build_circuit(np.zeros(ansatz.n_parameters))

            logger.info(f"  Circuit depth: {circuit.depth()}")

            # Check qubit limit
            if n_qubits > 30:
                logger.warning(f"  Too many qubits ({n_qubits}), would need active space reduction")
                return False

            # Run on cloud
            if self.backend_type == 'bluequbit':
                result = self.backend.run_circuit(circuit, shots=2048)
                energy = self._compute_energy_from_result(result, qubit_ham)
            else:
                observable = mapper.map_hamiltonian_to_observable(hamiltonian)
                result = self.backend.run_batch([circuit], [observable], shots=2048)
                energy = result['values'][0]

            hf_energy = hamiltonian.compute_hf_energy()

            self.results.append({
                'material': 'FeH (Catalyst)',
                'basis': basis,
                'mapper': 'jordan_wigner',
                'n_qubits': n_qubits,
                'hf_energy': hf_energy,
                'vqe_energy': energy,
                'correlation': energy - hf_energy if energy else None
            })

            logger.info(f"  HF Energy: {hf_energy:.6f} Ha")
            if energy:
                logger.info(f"  VQE Energy: {energy:.6f} Ha")
                logger.info(f"  Correlation: {energy - hf_energy:.6f} Ha")

            return True

        except Exception as e:
            logger.error(f"Error testing FeH: {e}", exc_info=True)
            return False

    def test_hydrogen_storage_material(self):
        """
        Test H2 molecule - fundamental for hydrogen storage materials.

        H2 dissociation on metal surfaces is critical for fuel cells.
        """
        logger.info(f"\n{'='*80}")
        logger.info("TEST 3: Hydrogen Molecule (Storage Material)")
        logger.info(f"{'='*80}")

        try:
            # Test H2 at different bond distances (dissociation curve)
            distances = [0.5, 0.74, 1.0, 1.5, 2.0]  # Angstroms

            for dist in distances:
                logger.info(f"\nH2 at distance: {dist} Angstrom")

                h1 = Atom('H', [0.0, 0.0, 0.0])
                h2 = Atom('H', [0.0, 0.0, dist])
                h2_mol = Molecule([h1, h2], charge=0, spin=0)

                # Test with different ansatze
                for ansatz_type in ['ucc', 'hardware_efficient']:
                    logger.info(f"  Ansatz: {ansatz_type}")

                    representation = LCAORepresentation(h2_mol)
                    hamiltonian = CovalentHamiltonian(h2_mol, representation, basis_name='sto-3g')
                    n_qubits = 2 * hamiltonian.n_orbitals

                    if ansatz_type == 'ucc':
                        ansatz = UCCAnsatz(
                            n_qubits=n_qubits,
                            n_electrons=h2_mol.n_electrons,
                            include_singles=True,
                            include_doubles=True
                        )
                    else:
                        ansatz = HardwareEfficientAnsatz(
                            n_qubits=n_qubits,
                            n_electrons=h2_mol.n_electrons,
                            n_layers=3
                        )

                    mapper = JordanWignerMapper()
                    qubit_ham = mapper.map_hamiltonian(hamiltonian)
                    circuit = ansatz.build_circuit(np.zeros(ansatz.n_parameters))

                    # Run on cloud
                    if self.backend_type == 'bluequbit':
                        result = self.backend.run_circuit(circuit, shots=1024)
                        energy = self._compute_energy_from_result(result, qubit_ham)
                    else:
                        observable = mapper.map_hamiltonian_to_observable(hamiltonian)
                        result = self.backend.run_batch([circuit], [observable], shots=1024)
                        energy = result['values'][0]

                    hf_energy = hamiltonian.compute_hf_energy()

                    self.results.append({
                        'material': f'H2 (d={dist}Å)',
                        'basis': 'sto-3g',
                        'ansatz': ansatz_type,
                        'n_qubits': n_qubits,
                        'hf_energy': hf_energy,
                        'vqe_energy': energy,
                        'correlation': energy - hf_energy if energy else None
                    })

                    logger.info(f"    HF: {hf_energy:.6f} Ha, VQE: {energy:.6f} Ha" if energy else f"    HF: {hf_energy:.6f} Ha")

            return True

        except Exception as e:
            logger.error(f"Error testing H2: {e}", exc_info=True)
            return False

    def _compute_energy_from_result(self, result, hamiltonian):
        """Compute energy from backend result."""
        try:
            if 'statevector' in result:
                statevector = result['statevector']
                pauli_sum = hamiltonian.to_pauli_sum()
                H_matrix = pauli_sum.to_matrix()
                energy = np.vdot(statevector, H_matrix @ statevector).real
                return energy
            elif 'counts' in result:
                # Estimate energy from measurement counts
                return None  # Would need to implement expectation value estimation
            return None
        except:
            return None

    def run_all_tests(self):
        """Run all materials science tests."""
        logger.info(f"\n{'#'*80}")
        logger.info(f"# Materials Science Cloud Testing on {self.backend_type.upper()}")
        logger.info(f"# Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"{'#'*80}\n")

        self.initialize_backend()

        # Run tests
        tests = [
            ('Lithium Hydride (Battery)', self.test_lithium_hydride_battery_material),
            ('Iron Hydride (Catalyst)', self.test_transition_metal_catalyst),
            ('Hydrogen (Storage)', self.test_hydrogen_storage_material),
        ]

        for name, test_func in tests:
            logger.info(f"\nRunning: {name}")
            success = test_func()
            logger.info(f"{'✓' if success else '✗'} {name} {'completed' if success else 'failed'}")

        # Print summary
        self.print_summary()

    def print_summary(self):
        """Print summary of all results."""
        logger.info(f"\n{'#'*80}")
        logger.info(f"# SUMMARY - Materials Science Results")
        logger.info(f"{'#'*80}\n")

        logger.info(f"Total simulations: {len(self.results)}\n")

        for r in self.results:
            logger.info(f"{r['material']}:")
            logger.info(f"  Basis: {r.get('basis', 'N/A')}, Qubits: {r['n_qubits']}")
            logger.info(f"  HF Energy: {r['hf_energy']:.6f} Ha")
            if r.get('vqe_energy'):
                logger.info(f"  VQE Energy: {r['vqe_energy']:.6f} Ha")
                logger.info(f"  Correlation: {r['correlation']:.6f} Ha")
            logger.info("")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description='Materials Science Cloud Testing')
    parser.add_argument(
        '--backend',
        choices=['ibm', 'bluequbit'],
        default='bluequbit',
        help='Cloud backend to use'
    )

    args = parser.parse_args()

    # Check environment
    if args.backend == 'ibm' and not os.getenv('IBM_API'):
        logger.error("IBM_API not set! Export IBM_API=your_token")
        sys.exit(1)

    if args.backend == 'bluequbit' and not os.getenv('BLUE_TOKEN'):
        logger.error("BLUE_TOKEN not set! Export BLUE_TOKEN=your_token")
        sys.exit(1)

    # Run tests
    tester = MaterialsScienceTester(backend_type=args.backend)
    tester.run_all_tests()


if __name__ == '__main__':
    main()
