"""
Metallic bond with governance protocol enforcement.

Models delocalized electron systems and band structure.
"""

from typing import Dict, Any, Optional, List
import numpy as np

from kanad.bonds.base_bond import BaseBond
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian
from kanad.governance.protocols.metallic_protocol import MetallicGovernanceProtocol
from kanad.core.temperature import Temperature
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper


class MetallicBond(BaseBond):
    """
    Metallic bond with automatic governance and quantum framework.

    Models:
    - Delocalized electrons across multiple atoms
    - Quantum tight-binding Hamiltonian
    - Band structure and Fermi surface
    - GHZ-like collective entanglement
    - Temperature-dependent properties

    Governance:
    - Enforces delocalization
    - Validates metallic character
    - Constructs band-structure ansatz

    Quantum Framework:
    - MetallicHamiltonian (tight-binding + Hubbard U)
    - MetallicGovernanceProtocol (validation)
    - Temperature (thermal effects)
    - VQE support for ground state

    Note:
        Metallic bonding requires multiple atoms (typically > 2)
        and shows collective quantum behavior.
    """

    def __init__(
        self,
        atoms: List[Atom],
        lattice_type: str = '1d_chain',
        hopping_parameter: Optional[float] = None,
        hubbard_u: float = 0.0,
        temperature: Optional[float] = None,
        periodic: bool = True
    ):
        """
        Initialize metallic bond with quantum framework.

        Args:
            atoms: List of atoms (typically all same element)
            lattice_type: Lattice structure ('1d_chain', '2d_square', etc.)
            hopping_parameter: Electron hopping strength t (eV), default -1.0
            hubbard_u: Coulomb repulsion U (eV), 0 for non-interacting
            temperature: Temperature in Kelvin (None for T=0)
            periodic: Use periodic boundary conditions
        """
        super().__init__(atoms, 'metallic', distance=None)

        self.lattice_type = lattice_type
        self.n_atoms = len(atoms)
        self.hopping_parameter = hopping_parameter if hopping_parameter is not None else -1.0
        self.hubbard_u = hubbard_u
        self.periodic = periodic

        # Create molecule
        self.molecule = Molecule(atoms)

        # Temperature (None = T=0)
        self.temperature = Temperature(temperature) if temperature is not None else Temperature.zero()

        # Create quantum Hamiltonian
        self.hamiltonian = MetallicHamiltonian(
            molecule=self.molecule,
            lattice_type=lattice_type,
            hopping_parameter=self.hopping_parameter,
            onsite_energy=0.0,
            hubbard_u=hubbard_u,
            periodic=periodic,
            temperature=temperature
        )

        # Governance protocol
        self.governance = MetallicGovernanceProtocol()

        # Mapper for qubit operations
        self.mapper = JordanWignerMapper()

    def compute_energy(
        self,
        method: str = 'tight_binding',
        **kwargs
    ) -> Dict[str, Any]:
        """
        Compute metallic system energy.

        Args:
            method: Computational method
                - 'tight_binding': Classical tight-binding (fast)
                - 'quantum': Quantum Hamiltonian (exact diagonalization)
                - 'VQE': Variational Quantum Eigensolver
            **kwargs: Method parameters
                use_temperature: Include thermal effects (bool)
                n_layers: VQE layers (int)
                max_iterations: VQE max iterations (int)

        Returns:
            Dictionary with results including:
                - energy: Ground state or thermal energy
                - fermi_energy: Fermi level
                - band_energies: Eigenvalues
                - thermal_properties: If temperature > 0
        """
        result = {}

        if method.lower() == 'tight_binding':
            # Classical tight-binding (original implementation)
            eigenvalues = np.linalg.eigvalsh(self.hamiltonian.h_tight_binding)
            n_electrons = self.hamiltonian.n_electrons
            n_bands_occupied = int(np.ceil(n_electrons / 2.0))

            # Ensure n_bands_occupied doesn't exceed available eigenvalues
            n_available_bands = len(eigenvalues)
            n_bands_occupied = min(n_bands_occupied, n_available_bands)

            # Compute energy (T=0 or thermal)
            if self.temperature.T > 0 and kwargs.get('use_temperature', False):
                # Thermal energy with Fermi-Dirac distribution
                fermi_energy = self.hamiltonian.get_fermi_energy(eigenvalues)
                total_energy = self.temperature.thermal_energy(
                    eigenvalues, fermi_energy, degeneracy=2
                )
                result['thermal'] = True
                result['entropy'] = self.temperature.entropy(eigenvalues, fermi_energy)
                result['free_energy'] = self.temperature.free_energy(eigenvalues, fermi_energy)
            else:
                # T=0: fill bands
                # For tight-binding, we need to account for the fact that each eigenvalue
                # represents a band that can hold 2 electrons (spin up and down)

                # Build density matrix from occupied orbitals
                eigenvectors = np.linalg.eigh(self.hamiltonian.h_tight_binding)[1]
                density_matrix = np.zeros_like(self.hamiltonian.h_tight_binding)

                # Fill occupied orbitals (accounting for spin)
                for i in range(n_bands_occupied):
                    density_matrix += 2.0 * np.outer(eigenvectors[:, i], eigenvectors[:, i])

                # Compute energy including Hubbard U
                total_energy = self.hamiltonian.compute_energy(density_matrix)
                result['thermal'] = False

            result['energy'] = total_energy
            result['method'] = 'Tight-Binding'
            result['converged'] = True
            result['band_energies'] = eigenvalues
            result['fermi_energy'] = self.hamiltonian.get_fermi_energy(eigenvalues)
            result['n_electrons'] = n_electrons
            result['is_metallic'] = self.hamiltonian.is_metallic()

        elif method.lower() == 'quantum':
            # Quantum Hamiltonian (exact diagonalization)
            # Use quantum many-body Hamiltonian
            eigenvalues = np.linalg.eigvalsh(self.hamiltonian.h_tight_binding)
            fermi_energy = self.hamiltonian.get_fermi_energy(eigenvalues)

            # For quantum treatment, we'd need to diagonalize full many-body Hamiltonian
            # For now, use tight-binding result
            result['energy'] = self.temperature.thermal_energy(eigenvalues, fermi_energy, 2)
            result['method'] = 'Quantum (Tight-Binding)'
            result['converged'] = True
            result['band_energies'] = eigenvalues
            result['fermi_energy'] = fermi_energy
            result['is_metallic'] = self.hamiltonian.is_metallic()

        elif method.lower() == 'vqe':
            # VQE for metallic system
            result = self._compute_vqe(**kwargs)

        else:
            raise ValueError(f"Method '{method}' not implemented")

        # Validate with governance
        validation = self.governance.validate_physical_constraints(self.hamiltonian)
        result['governance_validation'] = validation

        # Add bond analysis
        result['bond_analysis'] = self.analyze(result)
        return result

    def _compute_vqe(self, **kwargs) -> Dict[str, Any]:
        """
        Run VQE for metallic system.

        Args:
            **kwargs: VQE parameters

        Returns:
            VQE results
        """
        # Get suggested parameters from governance
        params = self.governance.suggest_parameters(self.hamiltonian)
        params.update(kwargs)  # Override with user params

        n_layers = params.get('n_layers', 2)
        entanglement = params.get('entanglement_type', 'ghz')
        max_iter = params.get('max_iterations', 500)

        # Build ansatz
        n_qubits = 2 * self.hamiltonian.n_orbitals  # spin up + down
        ansatz = self.governance.construct_ansatz(
            n_qubits, n_layers, entanglement
        )

        # Convert Hamiltonian to qubit operator
        from kanad.core.hamiltonians.pauli_converter import PauliConverter
        qubit_op = PauliConverter.to_sparse_pauli_op(
            self.hamiltonian, self.mapper, use_qiskit_nature=False
        )

        # Run VQE (simplified - would use actual VQE solver)
        # For now, return exact result
        H_matrix = qubit_op.to_matrix()
        eigenvalues = np.linalg.eigvalsh(H_matrix)
        ground_state_energy = eigenvalues[0]

        return {
            'energy': ground_state_energy,
            'method': 'VQE (exact diagonalization)',
            'converged': True,
            'n_qubits': n_qubits,
            'n_layers': n_layers,
            'entanglement': entanglement,
            'iterations': 0,  # Exact diagonalization, no iterations
            'band_energies': eigenvalues
        }

    def analyze(self, energy_data: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Analyze metallic bond properties.

        Args:
            energy_data: Optional energy computation results

        Returns:
            Dictionary with bond properties
        """
        analysis = {
            'bond_type': 'metallic',
            'lattice_type': self.lattice_type,
            'n_atoms': self.n_atoms,
            'hopping_parameter': self.hopping_parameter,
            'entanglement_type': 'GHZ-like (collective)',
            'governance_protocol': 'MetallicGovernanceProtocol'
        }

        if energy_data and 'band_energies' in energy_data:
            band_energies = energy_data['band_energies']
            analysis['bandwidth'] = band_energies[-1] - band_energies[0]
            analysis['fermi_energy'] = energy_data.get('fermi_energy', None)

            # DOS at Fermi level (simplified: count states near Fermi energy)
            if 'fermi_energy' in energy_data:
                E_F = energy_data['fermi_energy']
                delta_E = 0.1  # eV window
                dos_at_fermi = np.sum(np.abs(band_energies - E_F) < delta_E)
                analysis['dos_at_fermi'] = dos_at_fermi

        return analysis

    def get_band_structure(self) -> Dict[str, np.ndarray]:
        """
        Compute band structure.

        Returns:
            Dictionary with k-points and energies
        """
        # Use the Hamiltonian's band structure method
        return self.hamiltonian.get_band_structure(n_k=50)

    def __repr__(self) -> str:
        """String representation."""
        atom_symbol = self.atoms[0].symbol if self.atoms else '?'
        return (f"MetallicBond({atom_symbol}_{self.n_atoms}, "
                f"{self.lattice_type})")
