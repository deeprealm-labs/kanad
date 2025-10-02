"""
Ionic bond with governance protocol enforcement.

Models electron transfer between atoms with localized orbitals.
"""

from typing import Dict, Any, Optional
import numpy as np

from kanad.bonds.base_bond import BaseBond
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule
from kanad.core.representations.second_quantization import SecondQuantizationRepresentation
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.governance.protocols.ionic_protocol import IonicGovernanceProtocol
from kanad.ansatze.governance_aware_ansatz import IonicGovernanceAnsatz
from kanad.analysis.energy_analysis import EnergyAnalyzer, BondingAnalyzer


class IonicBond(BaseBond):
    """
    Ionic bond with automatic governance.

    Models:
    - Electron transfer from donor to acceptor
    - Localized atomic orbitals
    - Coulombic interactions
    - Minimal entanglement (charge transfer only)

    Governance:
    - Enforces localization
    - Validates charge transfer
    - Constructs minimal-entanglement ansatz
    """

    def __init__(
        self,
        atom_1: Atom,
        atom_2: Atom,
        distance: Optional[float] = None
    ):
        """
        Initialize ionic bond.

        Args:
            atom_1: First atom (typically donor/cation)
            atom_2: Second atom (typically acceptor/anion)
            distance: Bond distance in Angstroms (optional)
        """
        super().__init__([atom_1, atom_2], 'ionic', distance)

        # Identify donor and acceptor based on electronegativity
        if atom_1.electronegativity < atom_2.electronegativity:
            self.donor = atom_1
            self.acceptor = atom_2
        else:
            self.donor = atom_2
            self.acceptor = atom_1

        # Create molecule
        self.molecule = Molecule([atom_1, atom_2])

        # Set up representation (second quantization for ionic)
        self.representation = SecondQuantizationRepresentation(self.molecule)

        # Set up Hamiltonian
        self.hamiltonian = IonicHamiltonian(
            self.molecule,
            self.representation
        )

        # Governance protocol
        self.governance = IonicGovernanceProtocol()

        # Mapper (Jordan-Wigner for localized states)
        self.mapper = JordanWignerMapper()

    def compute_energy(
        self,
        method: str = 'VQE',
        max_iterations: int = 100,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Compute ionic bond energy.

        Args:
            method: Computational method
                - 'VQE': Variational Quantum Eigensolver
                - 'HF': Hartree-Fock (classical)
                - 'exact': Exact diagonalization (small systems)
            max_iterations: Maximum VQE iterations
            **kwargs: Additional method parameters

        Returns:
            Dictionary with results:
                - energy: Ground state energy (Hartree)
                - method: Method used
                - converged: Whether calculation converged
                - bond_analysis: Bond property analysis
        """
        result = {}

        if method.upper() == 'HF':
            # Hartree-Fock energy
            density_matrix, hf_energy = self.hamiltonian.solve_scf(
                max_iterations=max_iterations
            )
            result['energy'] = hf_energy
            result['method'] = 'Hartree-Fock'
            result['converged'] = True
            result['density_matrix'] = density_matrix

        elif method.upper() == 'VQE':
            # VQE with governance-aware ansatz
            from kanad.solvers.vqe_solver import VQESolver

            # Get number of qubits from representation
            n_qubits = self.representation.n_qubits
            n_electrons = self.molecule.n_electrons

            # Create ionic governance ansatz
            ansatz = IonicGovernanceAnsatz(
                n_qubits=n_qubits,
                n_electrons=n_electrons,
                n_layers=2
            )

            # Create VQE solver
            vqe = VQESolver(
                hamiltonian=self.hamiltonian,
                ansatz=ansatz,
                mapper=self.mapper,
                max_iterations=max_iterations
            )

            # Solve
            vqe_result = vqe.solve()

            result['energy'] = vqe_result['energy']
            result['method'] = 'VQE (Ionic Governance)'
            result['converged'] = vqe_result['converged']
            result['iterations'] = vqe_result['iterations']
            result['parameters'] = vqe_result['parameters']
            result['energy_history'] = vqe_result['energy_history']

            # Compute HF for comparison
            _, hf_energy = self.hamiltonian.solve_scf(max_iterations=50)
            result['hf_energy'] = hf_energy
            result['correlation_energy'] = result['energy'] - hf_energy

        elif method.upper() == 'EXACT':
            # Exact diagonalization (for small systems)
            H_matrix = self.hamiltonian.to_matrix()
            eigenvalues, eigenvectors = np.linalg.eigh(H_matrix)

            result['energy'] = eigenvalues[0]
            result['method'] = 'Exact'
            result['converged'] = True
            result['eigenvalues'] = eigenvalues
            result['ground_state'] = eigenvectors[:, 0]

        else:
            raise ValueError(f"Unknown method: {method}")

        # Add bond analysis
        result['bond_analysis'] = self.analyze()
        result['bond_length'] = self.get_bond_length()

        return result

    def analyze(self) -> Dict[str, Any]:
        """
        Analyze ionic bond properties.

        Returns:
            Dictionary with:
                - bond_type: 'ionic'
                - donor: Donor atom symbol
                - acceptor: Acceptor atom symbol
                - electronegativity_difference: ΔEN
                - charge_transfer: Estimated charge transfer
                - ionic_character: Percentage ionic character
                - coulomb_energy: Electrostatic energy
        """
        # Basic properties
        delta_en = abs(
            self.donor.electronegativity - self.acceptor.electronegativity
        )

        # Estimate charge transfer from EN difference
        # Empirical: q ≈ 1 - exp(-0.25 * ΔEN^2)
        charge_transfer = 1.0 - np.exp(-0.25 * delta_en**2)

        # Ionic character (Pauling formula)
        ionic_character = 1.0 - np.exp(-0.25 * delta_en**2)

        # Estimate Coulombic energy (simple point charge model)
        # E_coulomb = -k * q^2 / r (in atomic units, k = 1)
        bond_length_bohr = self.get_bond_length() / 0.529177  # Angstrom to Bohr
        if bond_length_bohr > 0:
            coulomb_energy = -charge_transfer**2 / bond_length_bohr
        else:
            coulomb_energy = 0.0

        analysis = {
            'bond_type': 'ionic',
            'donor': self.donor.symbol,
            'acceptor': self.acceptor.symbol,
            'electronegativity_difference': delta_en,
            'charge_transfer': charge_transfer,
            'ionic_character': ionic_character,
            'coulomb_energy': coulomb_energy,
            'bond_length': self.get_bond_length(),
            'entanglement_type': 'minimal (charge transfer only)',
            'governance_protocol': 'IonicGovernanceProtocol'
        }

        # Add Hamiltonian analysis if available
        if hasattr(self.hamiltonian, 'get_charge_distribution'):
            charges = self.hamiltonian.get_charge_distribution()
            analysis['atomic_charges'] = {
                self.atoms[0].symbol: charges[0],
                self.atoms[1].symbol: charges[1]
            }

        return analysis

    def get_charge_distribution(self) -> Dict[str, float]:
        """
        Get charge distribution on each atom.

        Returns:
            Dictionary mapping atom symbols to charges
        """
        if hasattr(self.hamiltonian, 'get_charge_distribution'):
            charges = self.hamiltonian.get_charge_distribution()
            return {
                self.donor.symbol: charges[0],
                self.acceptor.symbol: charges[1]
            }
        else:
            # Estimate from electronegativity
            analysis = self.analyze()
            q = analysis['charge_transfer']
            return {
                self.donor.symbol: +q,  # Loses electron (positive)
                self.acceptor.symbol: -q  # Gains electron (negative)
            }

    def __repr__(self) -> str:
        """String representation."""
        return (f"IonicBond({self.donor.symbol}+ — {self.acceptor.symbol}-, "
                f"{self.get_bond_length():.2f} Å)")
