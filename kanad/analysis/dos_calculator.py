"""
Density of States (DOS) Calculator for Periodic Systems.

Computes electronic density of states from band structure with
Gaussian broadening, projected DOS, and integration capabilities.
"""

import numpy as np
from typing import Dict, Any, Optional, List, Tuple
import logging

logger = logging.getLogger(__name__)


class DOSCalculator:
    """
    Calculate density of states for periodic systems.

    Supports:
    - Total DOS with Gaussian/Lorentzian broadening
    - Projected DOS (PDOS) on atoms or orbitals
    - Integrated DOS (IDOS)
    - Fermi level calculation
    - Van Hove singularities identification
    """

    def __init__(self, periodic_hamiltonian):
        """
        Initialize DOS calculator.

        Args:
            periodic_hamiltonian: PeriodicHamiltonian object with solved SCF

        Examples:
            >>> from kanad.analysis import DOSCalculator
            >>> dos_calc = DOSCalculator(crystal.hamiltonian)
            >>> result = dos_calc.compute_dos(energy_range=(-15, 15))
        """
        self.hamiltonian = periodic_hamiltonian

        if not hasattr(periodic_hamiltonian, 'band_energies') or \
           periodic_hamiltonian.band_energies is None:
            raise ValueError("Must run solve_scf() on PeriodicHamiltonian first")

        self.band_energies = periodic_hamiltonian.band_energies  # (nk, n_bands) in Ha
        self.k_weights = periodic_hamiltonian.k_weights  # (nk,)
        self.n_k = len(self.k_weights)
        self.n_bands = self.band_energies.shape[1]

        # Constants
        self.Ha_to_eV = 27.2114
        self.eV_to_Ha = 1.0 / self.Ha_to_eV

        logger.debug(f"DOSCalculator: {self.n_k} k-points, {self.n_bands} bands")

    def compute_dos(self,
                    energy_range: Tuple[float, float] = (-10, 10),
                    n_points: int = 1000,
                    sigma: float = 0.1,
                    method: str = 'gaussian',
                    units: str = 'eV') -> Dict[str, Any]:
        """
        Compute total density of states.

        Formula:
            DOS(E) = Σ_{nk} w_k × δ(E - E_nk)

        With Gaussian broadening:
            DOS(E) ≈ Σ_{nk} w_k / (σ√2π) exp(-(E - E_nk)² / 2σ²)

        Args:
            energy_range: (E_min, E_max) in eV (or Hartree if units='Ha')
            n_points: Number of energy grid points
            sigma: Broadening parameter in eV
            method: 'gaussian' or 'lorentzian'
            units: 'eV' or 'Ha' for energy units

        Returns:
            result: Dictionary with:
                - energies: (n_points,) energy grid
                - dos: (n_points,) density of states (states/eV or states/Ha)
                - idos: (n_points,) integrated DOS
                - fermi_energy: Fermi level
                - n_electrons: Number of electrons (from integration)
        """
        # Convert energy range to Hartree if needed
        if units == 'eV':
            E_min_Ha = energy_range[0] * self.eV_to_Ha
            E_max_Ha = energy_range[1] * self.eV_to_Ha
            sigma_Ha = sigma * self.eV_to_Ha
        else:
            E_min_Ha, E_max_Ha = energy_range
            sigma_Ha = sigma

        # Create energy grid
        energies_Ha = np.linspace(E_min_Ha, E_max_Ha, n_points)
        dos = np.zeros(n_points)

        # Sum over all k-points and bands
        for k_idx in range(self.n_k):
            for band_idx in range(self.n_bands):
                E_nk = self.band_energies[k_idx, band_idx]  # Hartree
                weight = self.k_weights[k_idx]

                # Apply broadening
                if method == 'gaussian':
                    dos += weight * self._gaussian_broadening(energies_Ha, E_nk, sigma_Ha)
                elif method == 'lorentzian':
                    dos += weight * self._lorentzian_broadening(energies_Ha, E_nk, sigma_Ha)
                else:
                    raise ValueError(f"Unknown method: {method}")

        # Normalize DOS
        # Total integral should equal number of states (nk × n_bands)
        # But we want states/eV, so we need to account for k-point sampling
        # DOS is already normalized by k_weights (which sum to 1)

        # Integrated DOS
        idos = self._integrate_dos(energies_Ha, dos)

        # Fermi energy
        fermi_energy_Ha = self._compute_fermi_energy_from_dos(energies_Ha, dos)

        # Convert back to desired units
        if units == 'eV':
            energies_out = energies_Ha * self.Ha_to_eV
            dos_out = dos / self.Ha_to_eV  # states/eV (not states/Ha)
            fermi_energy_out = fermi_energy_Ha * self.Ha_to_eV
        else:
            energies_out = energies_Ha
            dos_out = dos
            fermi_energy_out = fermi_energy_Ha

        return {
            'energies': energies_out,
            'dos': dos_out,
            'idos': idos,
            'fermi_energy': fermi_energy_out,
            'n_electrons_from_integration': idos[np.argmin(np.abs(energies_Ha - fermi_energy_Ha))],
            'n_electrons_actual': self.hamiltonian.n_electrons,
            'units': units
        }

    def _gaussian_broadening(self, energies: np.ndarray, E_center: float, sigma: float) -> np.ndarray:
        """
        Gaussian broadening function.

        Args:
            energies: Energy grid
            E_center: Center energy
            sigma: Width parameter

        Returns:
            broadened: Gaussian centered at E_center
        """
        prefactor = 1.0 / (sigma * np.sqrt(2 * np.pi))
        exponent = -0.5 * ((energies - E_center) / sigma)**2
        return prefactor * np.exp(exponent)

    def _lorentzian_broadening(self, energies: np.ndarray, E_center: float, gamma: float) -> np.ndarray:
        """
        Lorentzian (Cauchy) broadening function.

        Args:
            energies: Energy grid
            E_center: Center energy
            gamma: Half-width at half-maximum

        Returns:
            broadened: Lorentzian centered at E_center
        """
        prefactor = gamma / np.pi
        denominator = (energies - E_center)**2 + gamma**2
        return prefactor / denominator

    def _integrate_dos(self, energies: np.ndarray, dos: np.ndarray) -> np.ndarray:
        """
        Integrate DOS to get number of states up to each energy.

        Args:
            energies: Energy grid
            dos: DOS values

        Returns:
            idos: Cumulative integral (number of states)
        """
        from scipy.integrate import cumulative_trapezoid

        idos = cumulative_trapezoid(dos, energies, initial=0.0)
        return idos

    def _compute_fermi_energy_from_dos(self, energies: np.ndarray, dos: np.ndarray) -> float:
        """
        Compute Fermi energy by finding where IDOS = n_electrons.

        Args:
            energies: Energy grid
            dos: DOS values

        Returns:
            E_F: Fermi energy
        """
        idos = self._integrate_dos(energies, dos)
        n_electrons = self.hamiltonian.n_electrons

        # Find energy where IDOS crosses n_electrons
        idx = np.argmin(np.abs(idos - n_electrons))
        E_F = energies[idx]

        return E_F

    def compute_pdos(self,
                    atom_indices: Optional[List[int]] = None,
                    energy_range: Tuple[float, float] = (-10, 10),
                    n_points: int = 1000,
                    sigma: float = 0.1,
                    units: str = 'eV') -> Dict[str, Any]:
        """
        Compute projected density of states (PDOS).

        Projects DOS onto specified atoms to identify their contributions
        to electronic states.

        Args:
            atom_indices: List of atom indices to project onto (None = all)
            energy_range: (E_min, E_max) in eV
            n_points: Energy grid points
            sigma: Broadening in eV
            units: 'eV' or 'Ha'

        Returns:
            result: Dictionary with:
                - energies: Energy grid
                - pdos: Dictionary mapping atom_idx → PDOS array
                - total_dos: Total DOS for comparison
        """
        # This requires projection matrices from PySCF
        # For now, return placeholder

        logger.warning("PDOS not fully implemented yet - returning total DOS")

        # Compute total DOS
        total_result = self.compute_dos(energy_range, n_points, sigma, 'gaussian', units)

        if atom_indices is None:
            atom_indices = list(range(len(self.hamiltonian.atoms)))

        # Placeholder: split DOS equally among atoms
        n_atoms = len(atom_indices)
        pdos_dict = {}
        for idx in atom_indices:
            pdos_dict[idx] = total_result['dos'] / n_atoms

        return {
            'energies': total_result['energies'],
            'pdos': pdos_dict,
            'total_dos': total_result['dos'],
            'units': units
        }

    def find_band_gap(self) -> Dict[str, float]:
        """
        Find band gap from DOS.

        Identifies energy window with zero DOS near Fermi level.

        Returns:
            result: Dictionary with:
                - gap: Band gap in eV
                - vbm: Valence band maximum
                - cbm: Conduction band minimum
                - gap_type: 'direct' or 'indirect'
        """
        # Use Hamiltonian's built-in method
        return self.hamiltonian.get_band_gap()

    def find_van_hove_singularities(self,
                                    energy_range: Tuple[float, float] = (-10, 10),
                                    threshold: float = 2.0) -> List[Dict[str, float]]:
        """
        Identify Van Hove singularities in DOS.

        Van Hove singularities are peaks in DOS where ∇_k E(k) = 0.

        Args:
            energy_range: Energy window to search (eV)
            threshold: Minimum DOS value to qualify as singularity

        Returns:
            singularities: List of dictionaries with:
                - energy: Energy of singularity
                - dos_value: DOS at singularity
                - type: 'M0' (minimum), 'M1' (saddle), 'M2' (saddle), 'M3' (maximum)
        """
        # Compute DOS
        result = self.compute_dos(energy_range, n_points=2000, sigma=0.05)
        energies = result['energies']
        dos = result['dos']

        # Find local maxima
        from scipy.signal import find_peaks

        peaks, properties = find_peaks(dos, height=threshold, prominence=0.5)

        singularities = []
        for peak_idx in peaks:
            singularities.append({
                'energy': energies[peak_idx],
                'dos_value': dos[peak_idx],
                'type': 'unknown'  # Would need ∇²E analysis to classify
            })

        logger.info(f"Found {len(singularities)} Van Hove singularities")

        return singularities

    def plot_dos(self,
                 result: Optional[Dict[str, Any]] = None,
                 energy_range: Tuple[float, float] = (-10, 10),
                 show_fermi: bool = True,
                 show_gap: bool = True,
                 save_path: Optional[str] = None):
        """
        Plot density of states.

        Args:
            result: Pre-computed DOS result (None = compute now)
            energy_range: Energy range for plot
            show_fermi: Show Fermi level line
            show_gap: Shade band gap region
            save_path: Path to save figure (None = show)
        """
        import matplotlib.pyplot as plt

        if result is None:
            result = self.compute_dos(energy_range)

        energies = result['energies']
        dos = result['dos']
        E_F = result['fermi_energy']

        fig, ax = plt.subplots(figsize=(8, 6))

        ax.plot(energies, dos, 'b-', linewidth=2, label='DOS')
        ax.fill_between(energies, 0, dos, alpha=0.3)

        if show_fermi:
            ax.axvline(E_F, color='r', linestyle='--', linewidth=2,
                      label=f'E_F = {E_F:.2f} eV')

        if show_gap:
            gap_info = self.find_band_gap()
            if gap_info['gap'] > 0.1:  # Only show if significant gap
                vbm = gap_info['vbm']
                cbm = gap_info['cbm']
                ax.axvspan(vbm, cbm, alpha=0.2, color='gray',
                          label=f'Gap = {gap_info["gap"]:.2f} eV')

        ax.set_xlabel('Energy (eV)', fontsize=14)
        ax.set_ylabel('DOS (states/eV)', fontsize=14)
        ax.set_title('Density of States', fontsize=16)
        ax.legend()
        ax.grid(alpha=0.3)
        ax.set_xlim(energy_range)
        ax.set_ylim(bottom=0)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"DOS plot saved to {save_path}")
        else:
            plt.show()

    def plot_band_structure_with_dos(self,
                                     band_result: Dict[str, Any],
                                     dos_result: Optional[Dict[str, Any]] = None,
                                     save_path: Optional[str] = None):
        """
        Plot band structure and DOS side-by-side.

        Args:
            band_result: Result from compute_band_structure()
            dos_result: Result from compute_dos() (None = compute)
            save_path: Path to save figure
        """
        import matplotlib.pyplot as plt

        if dos_result is None:
            # Determine energy range from bands
            E_min = np.min(band_result['band_energies']) - 2
            E_max = np.max(band_result['band_energies']) + 2
            dos_result = self.compute_dos((E_min, E_max))

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                                       gridspec_kw={'width_ratios': [2, 1]})

        # Band structure
        k_dist = band_result['k_distances']
        bands = band_result['band_energies']

        for band_idx in range(bands.shape[1]):
            ax1.plot(k_dist, bands[:, band_idx], 'b-', linewidth=1.5)

        # High-symmetry points
        if 'labels' in band_result and 'label_positions' in band_result:
            for pos in band_result['label_positions']:
                ax1.axvline(pos, color='gray', linestyle='--', alpha=0.5)
            ax1.set_xticks(band_result['label_positions'])
            ax1.set_xticklabels(band_result['labels'])

        ax1.axhline(dos_result['fermi_energy'], color='r', linestyle='--',
                   linewidth=2, label='E_F')
        ax1.set_xlabel('k-path', fontsize=14)
        ax1.set_ylabel('Energy (eV)', fontsize=14)
        ax1.set_title('Band Structure', fontsize=16)
        ax1.legend()
        ax1.grid(alpha=0.3)

        # DOS
        ax2.plot(dos_result['dos'], dos_result['energies'], 'b-', linewidth=2)
        ax2.fill_betweenx(dos_result['energies'], 0, dos_result['dos'], alpha=0.3)
        ax2.axhline(dos_result['fermi_energy'], color='r', linestyle='--', linewidth=2)

        # Match y-axis
        ax2.set_ylim(ax1.get_ylim())
        ax2.set_xlabel('DOS (states/eV)', fontsize=14)
        ax2.set_ylabel('Energy (eV)', fontsize=14)
        ax2.set_title('Density of States', fontsize=16)
        ax2.grid(alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Band+DOS plot saved to {save_path}")
        else:
            plt.show()

    def __repr__(self) -> str:
        """String representation."""
        return (f"DOSCalculator(n_k={self.n_k}, n_bands={self.n_bands}, "
                f"n_electrons={self.hamiltonian.n_electrons})")
