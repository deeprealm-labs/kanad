"""
Bond Length Scanner for Potential Energy Surface (PES) Analysis

Systematically scans bond distances to explore potential energy surfaces.
Works for any bond in any molecule, preserving geometry of other atoms.
"""

import numpy as np
from typing import Dict, List, Any, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class BondLengthScanner:
    """
    Scan potential energy surface along a bond length coordinate.

    Supports:
    - Any bond in any molecule (diatomic or polyatomic)
    - Rigid scan (freeze all other atoms)
    - Relaxed scan (optimize other coordinates, future feature)
    - Smooth interpolation via cubic splines
    - Multiple electronic structure methods (HF, MP2)

    Example:
        >>> from kanad.io import from_smiles
        >>> from kanad.analysis import BondLengthScanner
        >>>
        >>> # Create water molecule
        >>> water = from_smiles("O")
        >>>
        >>> # Scan O-H bond (atoms 0 and 1)
        >>> scanner = BondLengthScanner(water, atom_i=0, atom_j=1)
        >>> result = scanner.scan(r_min=0.7, r_max=1.5, n_points=15)
        >>>
        >>> print(f"Equilibrium distance: {result['optimized_distance']:.4f} Å")
        >>> print(f"Minimum energy: {result['optimized_energy']:.6f} Ha")
    """

    def __init__(
        self,
        molecule: 'Molecule',
        atom_i: int,
        atom_j: int
    ):
        """
        Initialize bond length scanner.

        Args:
            molecule: Molecule object to scan
            atom_i: Index of first atom in bond (0-based)
            atom_j: Index of second atom in bond (0-based)

        Raises:
            ValueError: If atom indices are invalid
        """
        self.molecule = molecule
        self.atom_i = atom_i
        self.atom_j = atom_j

        # Validate indices
        n_atoms = len(molecule.atoms)
        if atom_i < 0 or atom_i >= n_atoms:
            raise ValueError(f"Invalid atom index {atom_i} (molecule has {n_atoms} atoms)")
        if atom_j < 0 or atom_j >= n_atoms:
            raise ValueError(f"Invalid atom index {atom_j} (molecule has {n_atoms} atoms)")
        if atom_i == atom_j:
            raise ValueError(f"Atom indices must be different (got {atom_i} and {atom_j})")

        # Store original geometry
        self.original_positions = [atom.position.copy() for atom in molecule.atoms]
        self.original_distance = self._compute_bond_length()

        logger.info(
            f"BondLengthScanner initialized for "
            f"{molecule.atoms[atom_i].symbol}-{molecule.atoms[atom_j].symbol} bond "
            f"(atoms {atom_i}-{atom_j})"
        )
        logger.info(f"Current bond length: {self.original_distance:.4f} Å")

    def _compute_bond_length(self) -> float:
        """Compute current bond length in Angstroms."""
        pos_i = self.molecule.atoms[self.atom_i].position
        pos_j = self.molecule.atoms[self.atom_j].position
        return np.linalg.norm(pos_j - pos_i)

    def _set_bond_length(self, target_distance: float) -> None:
        """
        Set bond length to target distance, preserving geometry of other atoms.

        Algorithm:
        1. Compute current bond vector: v = R_j - R_i
        2. Compute unit vector: u = v / |v|
        3. Set new position: R_j_new = R_i + target_distance * u
        4. Keep R_i and all other atoms fixed

        Args:
            target_distance: Desired bond length (Å)
        """
        atoms = self.molecule.atoms

        # Get current positions
        pos_i = atoms[self.atom_i].position.copy()
        pos_j = atoms[self.atom_j].position.copy()

        # Compute bond vector and unit vector
        bond_vector = pos_j - pos_i
        current_distance = np.linalg.norm(bond_vector)

        if current_distance < 1e-10:
            raise ValueError(f"Atoms {self.atom_i} and {self.atom_j} are too close")

        unit_vector = bond_vector / current_distance

        # Set new position for atom j (atom i remains fixed)
        atoms[self.atom_j].position = pos_i + target_distance * unit_vector

    def _rebuild_hamiltonian(self) -> None:
        """Rebuild molecular Hamiltonian after geometry change."""
        # Force rebuild of hamiltonian by setting it to None
        self.molecule._hamiltonian = None

    def scan(
        self,
        r_min: float,
        r_max: float,
        n_points: int = 20,
        method: str = 'HF',
        relax: bool = False,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Scan bond length from r_min to r_max.

        Args:
            r_min: Minimum bond length (Å)
            r_max: Maximum bond length (Å)
            n_points: Number of scan points (default: 20)
            method: Energy method ('HF', 'MP2') (default: 'HF')
            relax: If True, optimize other coordinates at each point (default: False)
                   Note: Relaxed scans not yet implemented
            verbose: Print progress (default: True)

        Returns:
            Dictionary with:
                distances: np.ndarray - Bond lengths scanned (Å)
                energies: np.ndarray - Energies at each point (Ha)
                geometries: List[List[np.ndarray]] - Geometries at each point
                optimized_distance: float - Minimum energy bond length (Å)
                optimized_energy: float - Minimum energy (Ha)
                original_distance: float - Starting bond length (Å)
                spline: callable - Cubic spline interpolation E(r)
                atom_i: int - First atom index
                atom_j: int - Second atom index
                atom_symbols: Tuple[str, str] - Atom symbols

        Raises:
            ValueError: If r_min >= r_max or n_points < 2
            NotImplementedError: If relax=True (not yet implemented)
        """
        # Validate inputs
        if r_min >= r_max:
            raise ValueError(f"r_min ({r_min}) must be less than r_max ({r_max})")
        if n_points < 2:
            raise ValueError(f"n_points ({n_points}) must be at least 2")
        if relax:
            raise NotImplementedError("Relaxed scans not yet implemented (coming soon)")

        # Create scan points
        distances = np.linspace(r_min, r_max, n_points)
        energies = []
        geometries = []

        # Get atom symbols for display
        symbol_i = self.molecule.atoms[self.atom_i].symbol
        symbol_j = self.molecule.atoms[self.atom_j].symbol

        if verbose:
            print(f"\nScanning {symbol_i}-{symbol_j} bond (atoms {self.atom_i}-{self.atom_j})")
            print(f"Range: {r_min:.3f} - {r_max:.3f} Å ({n_points} points)")
            print(f"Method: {method}")
            print(f"Original distance: {self.original_distance:.4f} Å")
            print("-" * 70)

        # Scan loop
        for i, r in enumerate(distances):
            # Set bond to target distance
            self._set_bond_length(r)

            # Rebuild Hamiltonian with new geometry
            self._rebuild_hamiltonian()

            # Compute energy
            try:
                if method.upper() == 'HF':
                    energy_ha = self.molecule.hamiltonian.hf_energy
                elif method.upper() == 'MP2':
                    from kanad.core.correlation import MP2Solver
                    mp2_solver = MP2Solver(self.molecule.hamiltonian)
                    result_mp2 = mp2_solver.compute_energy()
                    energy_ha = result_mp2['e_mp2']
                else:
                    raise ValueError(f"Unknown method: {method}")

                energies.append(energy_ha)

                # Save geometry
                geom = [atom.position.copy() for atom in self.molecule.atoms]
                geometries.append(geom)

                if verbose:
                    print(f"  Point {i+1:2d}/{n_points}: r={r:.4f} Å, E={energy_ha:.8f} Ha")

            except Exception as e:
                logger.error(f"Energy calculation failed at r={r:.4f} Å: {e}")
                energies.append(np.nan)
                geometries.append(None)

                if verbose:
                    print(f"  Point {i+1:2d}/{n_points}: r={r:.4f} Å - FAILED: {e}")

        energies = np.array(energies)

        # Restore original geometry
        for i, pos in enumerate(self.original_positions):
            self.molecule.atoms[i].position = pos.copy()
        self._rebuild_hamiltonian()

        # Find minimum (excluding NaN values)
        valid_mask = ~np.isnan(energies)
        if not np.any(valid_mask):
            logger.error("All energy calculations failed!")
            raise RuntimeError("All energy calculations failed")

        valid_distances = distances[valid_mask]
        valid_energies = energies[valid_mask]

        min_idx = np.argmin(valid_energies)
        opt_distance = valid_distances[min_idx]
        opt_energy = valid_energies[min_idx]

        # Create cubic spline interpolation
        spline = None
        if len(valid_distances) >= 4:  # Need at least 4 points for cubic spline
            from scipy.interpolate import CubicSpline
            spline = CubicSpline(valid_distances, valid_energies)

            # Refine minimum using spline
            r_fine = np.linspace(valid_distances[0], valid_distances[-1], 1000)
            e_fine = spline(r_fine)
            min_idx_fine = np.argmin(e_fine)
            opt_distance_spline = r_fine[min_idx_fine]
            opt_energy_spline = e_fine[min_idx_fine]

            if verbose:
                print("-" * 70)
                print(f"Grid minimum:   r={opt_distance:.4f} Å, E={opt_energy:.8f} Ha")
                print(f"Spline minimum: r={opt_distance_spline:.4f} Å, E={opt_energy_spline:.8f} Ha")

            # Use spline-refined values
            opt_distance = opt_distance_spline
            opt_energy = opt_energy_spline
        else:
            if verbose:
                print("-" * 70)
                print(f"Minimum: r={opt_distance:.4f} Å, E={opt_energy:.8f} Ha")
                print("(Not enough points for spline interpolation)")

        if verbose:
            print(f"\n✓ Scan complete!")

        return {
            'distances': distances,
            'energies': energies,
            'geometries': geometries,
            'optimized_distance': opt_distance,
            'optimized_energy': opt_energy,
            'original_distance': self.original_distance,
            'spline': spline,
            'atom_i': self.atom_i,
            'atom_j': self.atom_j,
            'atom_symbols': (symbol_i, symbol_j),
            'method': method,
            'n_points': n_points,
            'success': True
        }

    def plot(
        self,
        result: Dict[str, Any],
        show_spline: bool = True,
        save_path: Optional[str] = None
    ) -> None:
        """
        Plot potential energy curve.

        Args:
            result: Result dictionary from scan()
            show_spline: Show spline interpolation (default: True)
            save_path: Path to save figure (optional)
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed - cannot plot")
            return

        # Extract data
        distances = result['distances']
        energies = result['energies']
        spline = result.get('spline')
        opt_dist = result['optimized_distance']
        opt_energy = result['optimized_energy']
        symbol_i, symbol_j = result['atom_symbols']

        # Filter valid points
        valid_mask = ~np.isnan(energies)
        valid_distances = distances[valid_mask]
        valid_energies = energies[valid_mask]

        # Convert to kcal/mol relative to minimum
        e_kcal = (valid_energies - opt_energy) * 627.509  # Ha to kcal/mol

        # Create figure
        fig, ax = plt.subplots(figsize=(8, 6))

        # Plot computed points
        ax.plot(valid_distances, e_kcal, 'o', markersize=8, label='Computed', zorder=3)

        # Plot spline if available
        if show_spline and spline is not None:
            r_fine = np.linspace(valid_distances[0], valid_distances[-1], 500)
            e_fine = spline(r_fine)
            e_fine_kcal = (e_fine - opt_energy) * 627.509
            ax.plot(r_fine, e_fine_kcal, '-', linewidth=2, alpha=0.7,
                   label='Cubic spline', zorder=2)

        # Mark minimum
        ax.axvline(opt_dist, color='red', linestyle='--', alpha=0.5,
                  label=f'Minimum: {opt_dist:.4f} Å', zorder=1)

        # Labels and formatting
        ax.set_xlabel(f'{symbol_i}-{symbol_j} Bond Length (Å)', fontsize=12)
        ax.set_ylabel('Relative Energy (kcal/mol)', fontsize=12)
        ax.set_title(f'Potential Energy Surface: {symbol_i}-{symbol_j} Bond', fontsize=14)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Set y-axis to start at 0
        ax.set_ylim(bottom=0)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")

        plt.show()
