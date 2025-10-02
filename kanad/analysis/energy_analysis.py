"""
Energy analysis and decomposition tools.

Provides utilities for analyzing VQE results and decomposing
molecular energies into physical components.
"""

from typing import Dict, List, Optional, Tuple
import numpy as np
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian


class EnergyAnalyzer:
    """
    Analyzes molecular energies and their components.

    Decomposes total energy into:
    - Nuclear repulsion
    - Kinetic energy
    - Electron-nuclear attraction
    - Electron-electron repulsion
    - Exchange energy
    - Correlation energy
    """

    def __init__(self, hamiltonian: MolecularHamiltonian):
        """
        Initialize energy analyzer.

        Args:
            hamiltonian: Molecular Hamiltonian
        """
        self.hamiltonian = hamiltonian

    def decompose_energy(self, density_matrix: np.ndarray) -> Dict[str, float]:
        """
        Decompose total energy into components.

        Args:
            density_matrix: One-particle density matrix

        Returns:
            Dictionary with energy components (Hartree)
        """
        decomposition = {}

        # Nuclear repulsion
        decomposition['nuclear_repulsion'] = self.hamiltonian.nuclear_repulsion

        # One-electron terms
        h_core = self.hamiltonian.h_core
        E_core = np.sum(density_matrix * h_core)
        decomposition['one_electron'] = E_core

        # Two-electron terms (if available)
        if hasattr(self.hamiltonian, 'eri'):
            eri = self.hamiltonian.eri
            n = len(h_core)

            # Coulomb energy
            J = 0.0
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        for l in range(n):
                            J += 0.5 * density_matrix[i, j] * density_matrix[k, l] * eri[i, j, k, l]

            decomposition['coulomb'] = J

            # Exchange energy
            K = 0.0
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        for l in range(n):
                            K += 0.25 * density_matrix[i, k] * density_matrix[j, l] * eri[i, j, k, l]

            decomposition['exchange'] = -K  # Note: exchange lowers energy

            decomposition['two_electron'] = J - K
        else:
            decomposition['two_electron'] = 0.0

        # Total energy
        decomposition['total'] = sum(decomposition.values())

        return decomposition

    def compute_binding_energy(
        self,
        molecular_energy: float,
        atomic_energies: List[float]
    ) -> float:
        """
        Compute binding energy (bond dissociation energy).

        BE = E_atoms - E_molecule

        Args:
            molecular_energy: Total molecular energy
            atomic_energies: List of atomic energies

        Returns:
            Binding energy (positive for stable molecule)
        """
        total_atomic_energy = sum(atomic_energies)
        binding_energy = total_atomic_energy - molecular_energy

        return binding_energy

    def compute_ionization_energy(
        self,
        neutral_energy: float,
        cation_energy: float
    ) -> float:
        """
        Compute ionization energy.

        IE = E_cation - E_neutral

        Args:
            neutral_energy: Neutral molecule energy
            cation_energy: Cation energy

        Returns:
            Ionization energy (positive)
        """
        return cation_energy - neutral_energy

    def compute_electron_affinity(
        self,
        neutral_energy: float,
        anion_energy: float
    ) -> float:
        """
        Compute electron affinity.

        EA = E_neutral - E_anion

        Args:
            neutral_energy: Neutral molecule energy
            anion_energy: Anion energy

        Returns:
            Electron affinity (positive for stable anion)
        """
        return neutral_energy - anion_energy

    def analyze_convergence(self, energy_history: np.ndarray) -> Dict:
        """
        Analyze VQE convergence.

        Args:
            energy_history: Energy at each iteration

        Returns:
            Convergence metrics
        """
        analysis = {
            'initial_energy': energy_history[0],
            'final_energy': energy_history[-1],
            'energy_change': energy_history[-1] - energy_history[0],
            'iterations': len(energy_history),
            'converged': self._check_convergence(energy_history),
        }

        # Convergence rate
        if len(energy_history) > 1:
            energy_diffs = np.diff(energy_history)
            analysis['mean_energy_change'] = np.mean(np.abs(energy_diffs))
            analysis['final_gradient'] = np.abs(energy_diffs[-1])

        return analysis

    def _check_convergence(
        self,
        energy_history: np.ndarray,
        threshold: float = 1e-6
    ) -> bool:
        """Check if VQE converged."""
        if len(energy_history) < 2:
            return False

        # Check last few iterations
        n_check = min(5, len(energy_history) - 1)
        recent_changes = np.abs(np.diff(energy_history[-n_check:]))

        return np.all(recent_changes < threshold)


class BondingAnalyzer:
    """
    Analyzes chemical bonding characteristics.

    Provides tools for:
    - Bond order analysis
    - Orbital population analysis
    - Hybridization analysis
    - Charge transfer analysis
    """

    def __init__(self, hamiltonian: MolecularHamiltonian):
        """
        Initialize bonding analyzer.

        Args:
            hamiltonian: Molecular Hamiltonian
        """
        self.hamiltonian = hamiltonian

    def analyze_bonding_type(self) -> Dict:
        """
        Determine bonding type (ionic, covalent, metallic).

        Returns:
            Bonding analysis
        """
        analysis = {'bonding_type': 'unknown'}

        if isinstance(self.hamiltonian, IonicHamiltonian):
            analysis['bonding_type'] = 'ionic'
            analysis['characteristics'] = [
                'Localized electrons',
                'Charge transfer',
                'Electrostatic interactions'
            ]

        elif isinstance(self.hamiltonian, CovalentHamiltonian):
            analysis['bonding_type'] = 'covalent'
            analysis['characteristics'] = [
                'Shared electrons',
                'Orbital hybridization',
                'Bonding/antibonding MOs'
            ]

            # Get HOMO-LUMO gap
            if hasattr(self.hamiltonian, 'get_homo_lumo_gap'):
                gap = self.hamiltonian.get_homo_lumo_gap()
                analysis['homo_lumo_gap'] = gap
                analysis['homo_lumo_gap_ev'] = gap * 27.211  # Convert to eV

        return analysis

    def compute_mulliken_charges(
        self,
        density_matrix: np.ndarray,
        overlap_matrix: np.ndarray
    ) -> np.ndarray:
        """
        Compute Mulliken atomic charges.

        Q_A = Z_A - Σ_μ∈A Σ_ν P_μν S_νμ

        Args:
            density_matrix: Density matrix
            overlap_matrix: Overlap matrix

        Returns:
            Atomic charges
        """
        # Mulliken population matrix
        PS = density_matrix @ overlap_matrix

        # Compute gross atomic populations
        # Simplified: assume equal orbitals per atom
        n_orbitals = len(density_matrix)

        if isinstance(self.hamiltonian, CovalentHamiltonian):
            n_atoms = len(self.hamiltonian.atoms)
            orbitals_per_atom = n_orbitals // n_atoms

            charges = np.zeros(n_atoms)
            for atom_idx in range(n_atoms):
                start = atom_idx * orbitals_per_atom
                end = (atom_idx + 1) * orbitals_per_atom

                # Sum diagonal blocks
                population = np.trace(PS[start:end, start:end])

                # Charge = Z - population
                Z = self.hamiltonian.atoms[atom_idx].atomic_number
                charges[atom_idx] = Z - population

            return charges
        else:
            # For ionic systems, use different approach
            return np.array([])

    def analyze_bond_orders(
        self,
        density_matrix: np.ndarray
    ) -> Dict:
        """
        Analyze bond orders between atoms.

        Args:
            density_matrix: Density matrix

        Returns:
            Bond order analysis
        """
        analysis = {}

        if isinstance(self.hamiltonian, CovalentHamiltonian):
            n_atoms = len(self.hamiltonian.atoms)
            bond_orders = np.zeros((n_atoms, n_atoms))

            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    bo = self.hamiltonian.compute_bond_order(density_matrix, i, j)
                    bond_orders[i, j] = bo
                    bond_orders[j, i] = bo

            analysis['bond_orders'] = bond_orders
            analysis['bond_classification'] = self._classify_bonds(bond_orders)

        return analysis

    def _classify_bonds(self, bond_orders: np.ndarray) -> Dict:
        """
        Classify bonds as single, double, triple, etc.

        Args:
            bond_orders: Bond order matrix

        Returns:
            Bond classification
        """
        classification = {}
        n_atoms = len(bond_orders)

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                bo = bond_orders[i, j]

                if bo > 0.5:  # Significant bond
                    if bo < 1.3:
                        bond_type = 'single'
                    elif bo < 2.3:
                        bond_type = 'double'
                    elif bo < 3.3:
                        bond_type = 'triple'
                    else:
                        bond_type = 'multiple'

                    classification[f'atom_{i}_atom_{j}'] = {
                        'type': bond_type,
                        'order': bo
                    }

        return classification

    def compute_overlap_populations(
        self,
        density_matrix: np.ndarray,
        overlap_matrix: np.ndarray
    ) -> np.ndarray:
        """
        Compute overlap population matrix.

        OP_μν = P_μν S_μν

        Args:
            density_matrix: Density matrix
            overlap_matrix: Overlap matrix

        Returns:
            Overlap population matrix
        """
        return density_matrix * overlap_matrix


class CorrelationAnalyzer:
    """
    Analyzes electron correlation effects.

    Compares VQE (correlated) vs Hartree-Fock (mean-field) energies.
    """

    def __init__(self, hamiltonian: MolecularHamiltonian):
        """
        Initialize correlation analyzer.

        Args:
            hamiltonian: Molecular Hamiltonian
        """
        self.hamiltonian = hamiltonian

    def compute_correlation_energy(
        self,
        vqe_energy: float,
        hf_energy: float
    ) -> float:
        """
        Compute correlation energy.

        E_corr = E_VQE - E_HF

        Args:
            vqe_energy: VQE ground state energy
            hf_energy: Hartree-Fock energy

        Returns:
            Correlation energy (negative for stable correlation)
        """
        return vqe_energy - hf_energy

    def compute_percent_correlation(
        self,
        vqe_energy: float,
        hf_energy: float,
        exact_energy: Optional[float] = None
    ) -> float:
        """
        Compute percentage of correlation energy recovered.

        % = (E_HF - E_VQE) / (E_HF - E_exact) × 100

        Args:
            vqe_energy: VQE energy
            hf_energy: Hartree-Fock energy
            exact_energy: Exact (FCI) energy (if known)

        Returns:
            Percentage of correlation recovered
        """
        if exact_energy is None:
            return 0.0

        E_corr_total = exact_energy - hf_energy
        E_corr_vqe = vqe_energy - hf_energy

        if abs(E_corr_total) < 1e-10:
            return 100.0

        return (E_corr_vqe / E_corr_total) * 100.0

    def analyze_electron_correlation(
        self,
        vqe_result: Dict,
        hf_energy: float
    ) -> Dict:
        """
        Comprehensive correlation analysis.

        Args:
            vqe_result: VQE result dictionary
            hf_energy: Hartree-Fock energy

        Returns:
            Correlation analysis
        """
        vqe_energy = vqe_result['energy']
        correlation_energy = self.compute_correlation_energy(vqe_energy, hf_energy)

        analysis = {
            'hf_energy': hf_energy,
            'vqe_energy': vqe_energy,
            'correlation_energy': correlation_energy,
            'correlation_energy_ev': correlation_energy * 27.211,  # eV
        }

        # Correlation strength classification
        if abs(correlation_energy) < 0.001:
            analysis['correlation_strength'] = 'negligible'
        elif abs(correlation_energy) < 0.01:
            analysis['correlation_strength'] = 'weak'
        elif abs(correlation_energy) < 0.1:
            analysis['correlation_strength'] = 'moderate'
        else:
            analysis['correlation_strength'] = 'strong'

        return analysis
