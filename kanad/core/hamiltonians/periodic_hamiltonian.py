"""
Periodic Hamiltonian for crystalline systems.

Uses PySCF PBC module for k-point sampling and band structure calculations.
"""

import numpy as np
from typing import List, Optional, Dict, Any, Tuple, Union
import logging

from kanad.core.atom import Atom
from kanad.core.lattice import Lattice

logger = logging.getLogger(__name__)


class PeriodicHamiltonian:
    """
    Hamiltonian for periodic systems with k-point sampling.

    Uses PySCF's Periodic Boundary Condition (PBC) module for:
    - K-point sampling (Monkhorst-Pack grids)
    - Band structure calculations
    - Density of states
    - Crystal orbital analysis

    Supports 1D, 2D, and 3D periodic systems.
    """

    def __init__(self,
                 atoms: List[Atom],
                 lattice: Lattice,
                 charge: int = 0,
                 spin: int = 0,
                 basis: str = 'gth-dzvp',
                 pseudo: str = 'gth-pade',
                 k_points: Optional[Union[Tuple[int, int, int], np.ndarray]] = None):
        """
        Initialize periodic Hamiltonian.

        Args:
            atoms: List of Atom objects in the unit cell
            lattice: Lattice object with periodic boundary conditions
            charge: Total charge of the unit cell
            spin: Spin multiplicity (2S)
            basis: Basis set (recommend 'gth-dzvp' for PBC)
            pseudo: Pseudopotential type ('gth-pade', 'gth-pbe', etc.)
            k_points: Either:
                - (nk_x, nk_y, nk_z) for Monkhorst-Pack grid
                - (N, 3) array of explicit k-points
                - None for Γ-point only

        Examples:
            >>> # Diamond crystal with 4x4x4 k-point grid
            >>> lattice = Lattice(...)
            >>> atoms = [Atom('C', [0, 0, 0]), Atom('C', [0.25, 0.25, 0.25])]
            >>> ham = PeriodicHamiltonian(atoms, lattice, k_points=(4, 4, 4))
        """
        self.atoms = atoms
        self.lattice = lattice
        self.charge = charge
        self.spin = spin
        self.basis = basis
        self.pseudo = pseudo

        # Import PySCF PBC modules
        try:
            from pyscf.pbc import gto, scf, tools
            self.pbc_gto = gto
            self.pbc_scf = scf
            self.pbc_tools = tools
        except ImportError:
            raise ImportError("PySCF PBC module not available. Install with: pip install pyscf[geomopt]")

        # Build PySCF Cell object
        self._build_pyscf_cell()

        # Generate k-points
        self.k_points, self.k_weights = self._generate_k_points(k_points)

        logger.info(f"PeriodicHamiltonian: {self.cell.nelectron} electrons, "
                    f"{self.k_points.shape[0]} k-points, PBC={self.lattice.pbc}")

        # These will be set after SCF
        self.mf = None
        self.hf_energy = None
        self._scf_converged = False
        self.band_energies = None
        self.mo_coefficients = None

    def _build_pyscf_cell(self):
        """Build PySCF Cell object for periodic system."""
        self.cell = self.pbc_gto.Cell()

        # Set lattice vectors (PySCF uses 'a' not 'lattice_vectors')
        self.cell.a = self.lattice.lattice_vectors.copy()

        # Build atom string
        atom_str = []
        for atom in self.atoms:
            symbol = atom.symbol
            pos = atom.position  # Already in Angstrom
            atom_str.append(f"{symbol} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}")

        self.cell.atom = '; '.join(atom_str)
        self.cell.unit = 'A'  # Angstrom

        # Basis and pseudopotential
        self.cell.basis = self.basis
        self.cell.pseudo = self.pseudo

        # Charge and spin
        self.cell.charge = self.charge
        self.cell.spin = self.spin

        # Build the cell
        self.cell.build()

        # Store basic properties
        self.n_electrons = self.cell.nelectron
        self.n_orbitals = self.cell.nao_nr()  # Number of atomic orbitals

        logger.debug(f"Built cell: {len(self.atoms)} atoms, "
                     f"{self.n_electrons} electrons, "
                     f"{self.n_orbitals} orbitals")

    def _generate_k_points(self,
                          k_input: Optional[Union[Tuple[int, int, int], np.ndarray]]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate k-points and weights.

        Args:
            k_input: Either Monkhorst-Pack grid size or explicit k-points

        Returns:
            k_points: (nk, 3) array in reciprocal space (units of 2π/a)
            k_weights: (nk,) array of weights (sum to 1.0)
        """
        if k_input is None:
            # Γ-point only
            k_points = np.array([[0.0, 0.0, 0.0]])
            k_weights = np.array([1.0])
            logger.debug("Using Γ-point only")

        elif isinstance(k_input, (tuple, list)) and len(k_input) == 3:
            # Monkhorst-Pack grid
            nk = k_input
            k_points = self.cell.make_kpts(nk, with_gamma_point=True)
            nk_total = k_points.shape[0]
            k_weights = np.ones(nk_total) / nk_total  # Uniform weights

            logger.debug(f"Generated {nk_total} k-points from {nk} grid")

        elif isinstance(k_input, np.ndarray):
            # Explicit k-points
            k_points = np.array(k_input, dtype=float)
            assert k_points.ndim == 2 and k_points.shape[1] == 3
            nk_total = k_points.shape[0]
            k_weights = np.ones(nk_total) / nk_total

            logger.debug(f"Using {nk_total} explicit k-points")

        else:
            raise ValueError(f"Invalid k_input: {k_input}")

        return k_points, k_weights

    def solve_scf(self,
                  max_iterations: int = 50,
                  conv_tol: float = 1e-7,
                  verbose: int = 0) -> Dict[str, Any]:
        """
        Solve self-consistent field equations with k-point sampling.

        Args:
            max_iterations: Maximum SCF iterations
            conv_tol: Energy convergence tolerance
            verbose: Verbosity level (0-9)

        Returns:
            result: Dictionary with:
                - energy: Total energy per unit cell (Ha)
                - converged: Whether SCF converged
                - n_iterations: Number of iterations
                - band_energies: (nk, n_bands) array
                - fermi_energy: Fermi level (Ha)
        """
        # Choose RHF or ROHF based on spin
        if self.spin == 0:
            self.mf = self.pbc_scf.KRHF(self.cell, self.k_points)
        else:
            self.mf = self.pbc_scf.KROHF(self.cell, self.k_points)

        # Set convergence parameters
        self.mf.max_cycle = max_iterations
        self.mf.conv_tol = conv_tol
        self.mf.verbose = verbose

        # Run SCF
        self.hf_energy = self.mf.kernel()
        self._scf_converged = self.mf.converged

        # Extract band energies (MO energies at each k-point)
        # mf.mo_energy is a list of arrays (one per k-point)
        self.band_energies = np.array([mo_e for mo_e in self.mf.mo_energy])  # (nk, n_orbitals)
        self.mo_coefficients = self.mf.mo_coeff  # List of (n_orbitals, n_orbitals) arrays

        # Compute Fermi energy
        self.fermi_energy = self._compute_fermi_energy()

        logger.info(f"SCF converged: {self._scf_converged}, "
                    f"E = {self.hf_energy:.6f} Ha, "
                    f"E_F = {self.fermi_energy:.4f} Ha")

        return {
            'energy': self.hf_energy,
            'converged': self._scf_converged,
            'band_energies': self.band_energies,
            'fermi_energy': self.fermi_energy,
            'n_iterations': getattr(self.mf, 'iter', 0)
        }

    def _compute_fermi_energy(self) -> float:
        """
        Compute Fermi energy from band occupation.

        Returns:
            E_F: Fermi energy in Hartree
        """
        if self.band_energies is None:
            raise ValueError("Must run solve_scf() first")

        # Flatten all band energies across k-points
        all_energies = self.band_energies.flatten()

        # Sort energies
        all_energies_sorted = np.sort(all_energies)

        # Number of occupied states (accounting for spin and k-point weights)
        # Each k-point contributes equally, each orbital holds 2 electrons (if RHF)
        n_occupied_per_k = self.n_electrons // 2  # For RHF

        # Total occupied states across all k-points
        nk = len(self.k_points)
        n_occupied_total = n_occupied_per_k * nk

        # Fermi energy is highest occupied level
        # (In metals, this is approximate; should use Fermi-Dirac smearing)
        if n_occupied_total <= len(all_energies_sorted):
            E_F = all_energies_sorted[n_occupied_total - 1]
        else:
            E_F = all_energies_sorted[-1]

        return E_F

    def compute_band_structure(self,
                               k_path: Union[str, np.ndarray],
                               n_bands: Optional[int] = None,
                               n_points: int = 50) -> Dict[str, Any]:
        """
        Compute band structure along high-symmetry path.

        Args:
            k_path: Either:
                - String like 'GXLG' (auto-generate path)
                - (N, 3) array of explicit k-points
            n_bands: Number of bands to compute (default: all)
            n_points: Points per segment (if k_path is string)

        Returns:
            result: Dictionary with:
                - k_points: (N, 3) k-points along path
                - k_distances: (N,) cumulative distance along path
                - band_energies: (N, n_bands) energies in eV
                - labels: High-symmetry point labels
                - label_positions: Positions of labels in k_distances
        """
        if self.mf is None:
            raise ValueError("Must run solve_scf() first")

        # Generate k-path
        if isinstance(k_path, str):
            k_points_band, k_labels, label_positions = self._generate_band_path(k_path, n_points)
        else:
            k_points_band = np.array(k_path)
            k_labels = []
            label_positions = []

        # Compute band energies at each k-point
        band_energies_list = []

        for k in k_points_band:
            # Get Fock matrix at this k-point
            # For band structure, we use the converged density
            # This is simplified - proper implementation uses get_bands()
            pass

        # Use PySCF's built-in band structure calculation
        from pyscf.pbc.tools import get_bandpath

        if isinstance(k_path, str):
            # Use PySCF's helper
            bands = self.mf.get_bands(k_points_band)
        else:
            bands = self.mf.get_bands(k_points_band)

        # bands is (nk, n_orbitals) array in Hartree
        # Select only n_bands if specified
        if n_bands is None:
            n_bands = self.n_orbitals

        band_energies_eV = bands[0][:, :n_bands] * 27.2114  # Ha -> eV

        # Compute cumulative k-distance
        k_distances = self._compute_k_distances(k_points_band)

        return {
            'k_points': k_points_band,
            'k_distances': k_distances,
            'band_energies': band_energies_eV,
            'labels': k_labels,
            'label_positions': label_positions
        }

    def _generate_band_path(self,
                           path_string: str,
                           n_points: int = 50) -> Tuple[np.ndarray, List[str], List[float]]:
        """
        Generate k-point path from string like 'GXLG'.

        Uses PySCF's lattice tools to get standard paths.

        Args:
            path_string: String of high-symmetry points
            n_points: Points per segment

        Returns:
            k_points: (N, 3) array
            labels: List of point labels
            label_positions: Cumulative distances of labels
        """
        try:
            from pyscf.pbc.tools import lattice
            k_points, k_path_segments, special_points = lattice.get_band(
                self.cell, path_string, n_points
            )
            return k_points, list(path_string), k_path_segments
        except:
            # Fallback: manual generation for common paths
            logger.warning(f"Could not auto-generate path '{path_string}', using Γ-X")

            # Simple fallback: Γ → X (works for cubic)
            gamma = np.array([0.0, 0.0, 0.0])
            X = np.array([0.5, 0.0, 0.0])

            k_points = np.linspace(gamma, X, n_points)
            labels = ['Γ', 'X']
            label_positions = [0.0, 1.0]

            return k_points, labels, label_positions

    def _compute_k_distances(self, k_points: np.ndarray) -> np.ndarray:
        """
        Compute cumulative distance along k-path.

        Args:
            k_points: (N, 3) array in reciprocal space

        Returns:
            distances: (N,) cumulative distances
        """
        distances = np.zeros(len(k_points))

        for i in range(1, len(k_points)):
            # Distance in reciprocal space
            dk = k_points[i] - k_points[i-1]
            # Convert to Cartesian reciprocal space
            dk_cart = dk @ self.lattice.reciprocal_vectors
            distances[i] = distances[i-1] + np.linalg.norm(dk_cart)

        return distances

    def compute_density_matrix(self) -> np.ndarray:
        """
        Compute density matrix from k-point averaged MOs.

        Returns:
            density_matrix: Real-space density matrix
        """
        if self.mf is None:
            raise ValueError("Must run solve_scf() first")

        # PySCF handles this internally
        return self.mf.make_rdm1()

    def get_band_gap(self) -> Dict[str, float]:
        """
        Compute band gap.

        Returns:
            result: Dictionary with:
                - gap: Band gap in eV
                - vbm: Valence band maximum (eV)
                - cbm: Conduction band minimum (eV)
                - type: 'direct' or 'indirect'
        """
        if self.band_energies is None:
            raise ValueError("Must run solve_scf() first")

        Ha_to_eV = 27.2114

        # Find highest occupied and lowest unoccupied
        n_occ = self.n_electrons // 2  # For RHF

        # VBM: maximum of all occupied bands
        vbm = np.max(self.band_energies[:, :n_occ]) * Ha_to_eV

        # CBM: minimum of all unoccupied bands
        if n_occ < self.n_orbitals:
            cbm = np.min(self.band_energies[:, n_occ:]) * Ha_to_eV
            gap = cbm - vbm
        else:
            # No unoccupied bands (shouldn't happen)
            cbm = np.nan
            gap = 0.0

        # Determine if direct or indirect
        # Direct: VBM and CBM at same k-point
        is_direct = False
        for k_idx in range(len(self.k_points)):
            k_vbm = np.max(self.band_energies[k_idx, :n_occ]) * Ha_to_eV
            k_cbm = np.min(self.band_energies[k_idx, n_occ:]) * Ha_to_eV if n_occ < self.n_orbitals else np.nan

            if np.isclose(k_vbm, vbm, atol=1e-6) and np.isclose(k_cbm, cbm, atol=1e-6):
                is_direct = True
                break

        return {
            'gap': gap,
            'vbm': vbm,
            'cbm': cbm,
            'type': 'direct' if is_direct else 'indirect'
        }

    def to_matrix(self) -> np.ndarray:
        """
        Get Hamiltonian matrix (Fock matrix).

        For periodic systems, this is k-point dependent.
        Returns Fock matrix at Γ-point.

        Returns:
            fock_matrix: (n_orbitals, n_orbitals) array
        """
        if self.mf is None:
            raise ValueError("Must run solve_scf() first")

        # Return Fock matrix at first k-point (usually Γ)
        fock = self.mf.get_fock()[0]  # List of Fock matrices
        return fock

    def __repr__(self) -> str:
        """String representation."""
        return (f"PeriodicHamiltonian(atoms={len(self.atoms)}, "
                f"n_electrons={self.n_electrons}, "
                f"n_k={len(self.k_points)}, "
                f"basis='{self.basis}')")
