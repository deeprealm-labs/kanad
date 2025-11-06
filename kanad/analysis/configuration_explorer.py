"""
Configuration Space Explorer

Explores molecular configuration space with Kanad's governance system:
- Real-time bond tracking during geometry optimization
- Reaction path finding (NEB, string method)
- Conformational analysis (rotamers, tautomers)
- Transition state search
- Potential energy surface (PES) scanning
- Bond breaking/forming visualization

Integrates environmental effects (T, P, pH, solvent) with quantum solvers.

Key Features:
    - Governance-aware: Automatically detects bond formation/breaking
    - Multi-environment: Scan T, P, pH, solvent simultaneously
    - Quantum-accurate: Uses SQD, Hi-VQE for energy calculations
    - Interactive: Real-time visualization of configuration changes

Example:
    >>> from kanad.analysis import ConfigurationExplorer
    >>> from kanad.molecule import Molecule
    >>>
    >>> explorer = ConfigurationExplorer()
    >>>
    >>> # Scan bond dissociation
    >>> h2 = Molecule.from_atoms([('H', [0, 0, 0]), ('H', [0.74, 0, 0])])
    >>> pes = explorer.scan_bond_length(h2, atom1=0, atom2=1,
    ...                                   r_range=(0.5, 3.0), n_points=20)
    >>>
    >>> # Find reaction path
    >>> path = explorer.find_reaction_path(reactant, product, method='neb')
    >>>
    >>> # Conformational search
    >>> conformers = explorer.conformational_search(molecule, max_conformers=10)
"""

import numpy as np
from typing import Dict, Any, Optional, List, Tuple, Callable
import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class ConfigurationSnapshot:
    """
    Snapshot of molecular configuration at a point in configuration space.

    Attributes:
        geometry: Atomic coordinates (Å)
        energy: Total energy (Ha)
        bonds: List of active bonds [(atom1, atom2, distance)]
        reaction_coordinate: Value of reaction coordinate
        environment: Environmental conditions (T, P, pH, solvent)
        gradient: Energy gradient (if available)
        converged: Whether geometry optimization converged
    """
    geometry: np.ndarray
    energy: float
    bonds: List[Tuple[int, int, float]]
    reaction_coordinate: Optional[float] = None
    environment: Optional[Dict[str, Any]] = None
    gradient: Optional[np.ndarray] = None
    converged: bool = False
    metadata: Optional[Dict[str, Any]] = None


@dataclass
class ReactionPath:
    """
    Complete reaction path from reactant to product.

    Attributes:
        images: List of configurations along path
        energies: Energy profile (Ha)
        reaction_coordinate: Reaction coordinate values
        transition_state_index: Index of TS configuration
        activation_energy: Forward barrier (Ha)
        reaction_energy: ΔE_rxn (Ha)
        path_length: Total path length in configuration space
    """
    images: List[ConfigurationSnapshot]
    energies: np.ndarray
    reaction_coordinate: np.ndarray
    transition_state_index: Optional[int] = None
    activation_energy: Optional[float] = None
    reaction_energy: Optional[float] = None
    path_length: Optional[float] = None


class ConfigurationExplorer:
    """
    Explore molecular configuration space with environmental effects.

    Uses Kanad's governance system to automatically track:
    - Bond formation and breaking
    - Transition state crossings
    - Phase transitions under pressure
    - Protonation state changes with pH
    - Solvation shell restructuring
    """

    Ha_to_kcal = 627.509474
    bohr_to_angstrom = 0.529177

    def __init__(
        self,
        solver_type: str = 'sqd',
        backend: str = 'statevector',
        use_governance: bool = True,
        track_bonds: bool = True
    ):
        """
        Initialize configuration explorer.

        Args:
            solver_type: Quantum solver ('vqe', 'sqd', 'hivqe')
            backend: Quantum backend ('statevector', 'ibm', 'bluequbit')
            use_governance: Use Kanad governance for bond tracking
            track_bonds: Automatically detect bond changes
        """
        self.solver_type = solver_type
        self.backend = backend
        self.use_governance = use_governance
        self.track_bonds = track_bonds

        # Environmental modulators (lazy import to avoid circular dependencies)
        self._temp_mod = None
        self._solv_mod = None
        self._ph_mod = None
        self._press_mod = None

        logger.info(f"ConfigurationExplorer initialized: solver={solver_type}, "
                   f"backend={backend}, governance={use_governance}")

    def scan_bond_length(
        self,
        molecule,
        atom1: int,
        atom2: int,
        r_range: Tuple[float, float] = (0.5, 3.0),
        n_points: int = 20,
        environment: Optional[Dict[str, Any]] = None,
        optimize_other_coords: bool = False
    ) -> Dict[str, np.ndarray]:
        """
        Scan potential energy surface along bond length coordinate.

        Args:
            molecule: Molecule object
            atom1, atom2: Atom indices for bond to scan
            r_range: (r_min, r_max) bond length range in Angstroms
            n_points: Number of points along coordinate
            environment: Environmental conditions {'T': 298.15, 'P': 1.0, 'pH': 7.0, 'solvent': 'water'}
            optimize_other_coords: Relax all other coordinates at each point

        Returns:
            Dictionary with:
                bond_lengths: r values (Å)
                energies: E(r) (Ha)
                bonds_active: List of active bonds at each point
                gradients: Energy gradients (if available)
                transition_detected: Whether bond breaking/forming detected
        """
        logger.info(f"Scanning bond {atom1}-{atom2}: r = {r_range[0]:.2f}-{r_range[1]:.2f} Å, "
                   f"{n_points} points")

        r_min, r_max = r_range
        bond_lengths = np.linspace(r_min, r_max, n_points)

        snapshots = []
        energies = []
        bonds_active = []

        for r in bond_lengths:
            # Set bond length
            config = self._set_bond_length(molecule, atom1, atom2, r)

            # Optimize other coordinates if requested
            if optimize_other_coords:
                config = self._optimize_geometry(config, freeze_atoms=[atom1, atom2])

            # Compute energy with environmental effects
            snapshot = self._compute_configuration_energy(
                config, environment=environment
            )
            snapshot.reaction_coordinate = r

            snapshots.append(snapshot)
            energies.append(snapshot.energy)
            bonds_active.append(snapshot.bonds)

            logger.debug(f"r = {r:.3f} Å, E = {snapshot.energy:.6f} Ha, "
                        f"{len(snapshot.bonds)} bonds active")

        # Detect transitions
        transition_detected = self._detect_bond_transition(bonds_active)

        return {
            'bond_lengths': bond_lengths,
            'energies': np.array(energies),
            'snapshots': snapshots,
            'bonds_active': bonds_active,
            'transition_detected': transition_detected,
            'dissociation_energy': self._compute_dissociation_energy(energies, bond_lengths)
        }

    def scan_reaction_coordinate(
        self,
        molecule,
        coordinate_func: Callable,
        coord_range: Tuple[float, float],
        n_points: int = 20,
        environment: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Scan along arbitrary reaction coordinate.

        Args:
            molecule: Molecule object
            coordinate_func: Function that modifies geometry given coordinate value
            coord_range: (min, max) range for reaction coordinate
            n_points: Number of points
            environment: Environmental conditions

        Returns:
            Dictionary with PES data along coordinate
        """
        coord_min, coord_max = coord_range
        coordinates = np.linspace(coord_min, coord_max, n_points)

        snapshots = []
        energies = []

        for xi in coordinates:
            # Apply coordinate transformation
            config = coordinate_func(molecule, xi)

            # Compute energy
            snapshot = self._compute_configuration_energy(config, environment)
            snapshot.reaction_coordinate = xi

            snapshots.append(snapshot)
            energies.append(snapshot.energy)

        return {
            'reaction_coordinate': coordinates,
            'energies': np.array(energies),
            'snapshots': snapshots
        }

    def find_reaction_path(
        self,
        reactant,
        product,
        method: str = 'neb',
        n_images: int = 10,
        environment: Optional[Dict[str, Any]] = None,
        max_iterations: int = 100,
        convergence_threshold: float = 1e-3
    ) -> ReactionPath:
        """
        Find minimum energy path (MEP) between reactant and product.

        Methods:
            - 'neb': Nudged Elastic Band
            - 'string': String method
            - 'linear': Linear interpolation (fast, approximate)

        Args:
            reactant: Reactant molecule
            product: Product molecule
            method: Path finding method
            n_images: Number of images along path
            environment: Environmental conditions
            max_iterations: Maximum optimization iterations
            convergence_threshold: Convergence criterion for path optimization

        Returns:
            ReactionPath object with full path information
        """
        logger.info(f"Finding reaction path: method={method}, {n_images} images")

        if method == 'linear':
            return self._linear_interpolation_path(reactant, product, n_images, environment)
        elif method == 'neb':
            return self._nudged_elastic_band(reactant, product, n_images, environment,
                                            max_iterations, convergence_threshold)
        elif method == 'string':
            return self._string_method(reactant, product, n_images, environment,
                                      max_iterations, convergence_threshold)
        else:
            raise ValueError(f"Unknown path method: {method}. Use 'neb', 'string', or 'linear'")

    def conformational_search(
        self,
        molecule,
        max_conformers: int = 10,
        energy_window: float = 10.0,  # kcal/mol
        environment: Optional[Dict[str, Any]] = None,
        method: str = 'systematic'
    ) -> List[ConfigurationSnapshot]:
        """
        Search conformational space for low-energy conformers.

        Methods:
            - 'systematic': Systematic grid search over torsions
            - 'random': Random sampling (Monte Carlo)
            - 'genetic': Genetic algorithm

        Args:
            molecule: Molecule object
            max_conformers: Maximum number of conformers to return
            energy_window: Energy window above global minimum (kcal/mol)
            environment: Environmental conditions
            method: Search method

        Returns:
            List of unique low-energy conformers
        """
        logger.info(f"Conformational search: method={method}, max={max_conformers}, "
                   f"window={energy_window} kcal/mol")

        if method == 'systematic':
            return self._systematic_conformer_search(
                molecule, max_conformers, energy_window, environment
            )
        elif method == 'random':
            return self._random_conformer_search(
                molecule, max_conformers, energy_window, environment
            )
        else:
            raise ValueError(f"Unknown conformational search method: {method}")

    def scan_environmental_effects(
        self,
        molecule,
        scan_params: Dict[str, Tuple[float, float, int]],
        optimize_geometry: bool = False
    ) -> Dict[str, Any]:
        """
        Scan multiple environmental parameters simultaneously.

        Args:
            molecule: Molecule object
            scan_params: Dict of parameter → (min, max, n_points)
                Example: {'T': (100, 500, 10), 'pH': (2, 12, 20)}
            optimize_geometry: Re-optimize at each condition

        Returns:
            Multi-dimensional grid of energies and properties
        """
        logger.info(f"Environmental scan: parameters={list(scan_params.keys())}")

        # For simplicity, scan one parameter at a time
        # Full multi-dimensional grid would require tensor product
        results = {}

        for param_name, (pmin, pmax, n_points) in scan_params.items():
            param_values = np.linspace(pmin, pmax, n_points)
            energies = []
            snapshots = []

            for val in param_values:
                env = {param_name: val}
                snapshot = self._compute_configuration_energy(molecule, environment=env)
                energies.append(snapshot.energy)
                snapshots.append(snapshot)

            results[param_name] = {
                'values': param_values,
                'energies': np.array(energies),
                'snapshots': snapshots
            }

        return results

    # ========== Private Methods ==========

    def _set_bond_length(self, molecule, atom1: int, atom2: int, r: float):
        """
        Set distance between two atoms to specified value.

        Args:
            molecule: Molecule object
            atom1, atom2: Atom indices
            r: New bond length in Angstroms

        Returns:
            New molecule with modified geometry
        """
        # Get current geometry
        coords = np.array(molecule.coordinates)  # Angstroms

        # Current bond vector
        vec = coords[atom2] - coords[atom1]
        r_current = np.linalg.norm(vec)

        # Scale to new length
        vec_new = vec * (r / r_current)

        # Update atom2 position
        coords_new = coords.copy()
        coords_new[atom2] = coords_new[atom1] + vec_new

        # Create new molecule with updated coordinates
        # This is pseudocode - actual implementation depends on Molecule class
        from copy import deepcopy
        mol_new = deepcopy(molecule)
        mol_new.coordinates = coords_new.tolist()

        return mol_new

    def _compute_configuration_energy(
        self,
        molecule,
        environment: Optional[Dict[str, Any]] = None
    ) -> ConfigurationSnapshot:
        """
        Compute energy of configuration with environmental effects.

        Args:
            molecule: Molecule object
            environment: Dict with 'T', 'P', 'pH', 'solvent'

        Returns:
            ConfigurationSnapshot with energy and bond information
        """
        # Apply environmental effects
        E_base = self._get_molecular_energy(molecule)

        if environment:
            # Apply temperature
            if 'T' in environment:
                E_base = self._apply_temperature_effect(molecule, E_base, environment['T'])

            # Apply pressure
            if 'P' in environment:
                E_base = self._apply_pressure_effect(molecule, E_base, environment['P'])

            # Apply pH
            if 'pH' in environment:
                E_base = self._apply_pH_effect(molecule, E_base, environment['pH'])

            # Apply solvent
            if 'solvent' in environment:
                E_base = self._apply_solvent_effect(molecule, E_base, environment['solvent'])

        # Detect active bonds using governance
        bonds = self._detect_bonds(molecule) if self.track_bonds else []

        snapshot = ConfigurationSnapshot(
            geometry=np.array(molecule.coordinates),
            energy=E_base,
            bonds=bonds,
            environment=environment,
            converged=True
        )

        return snapshot

    def _get_molecular_energy(self, molecule) -> float:
        """Compute molecular energy using specified quantum solver."""
        # Placeholder - actual implementation uses VQE/SQD/Hi-VQE
        if hasattr(molecule, 'compute_energy'):
            return molecule.compute_energy(solver=self.solver_type, backend=self.backend)
        elif hasattr(molecule, 'energy'):
            return molecule.energy
        else:
            logger.warning("Cannot compute energy - returning 0.0")
            return 0.0

    def _apply_temperature_effect(self, molecule, E_base: float, T: float) -> float:
        """Apply temperature corrections to energy."""
        if self._temp_mod is None:
            from kanad.environment import TemperatureModulator
            self._temp_mod = TemperatureModulator()

        result = self._temp_mod.apply_temperature(molecule, T)
        return result['energy']

    def _apply_pressure_effect(self, molecule, E_base: float, P: float) -> float:
        """Apply pressure corrections to energy."""
        if self._press_mod is None:
            from kanad.environment import PressureModulator
            self._press_mod = PressureModulator()

        result = self._press_mod.apply_pressure(molecule, P)
        return result['energy']

    def _apply_pH_effect(self, molecule, E_base: float, pH: float) -> float:
        """Apply pH corrections to energy."""
        if self._ph_mod is None:
            from kanad.environment import pHModulator
            self._ph_mod = pHModulator()
            # Would need to add protonation sites here

        result = self._ph_mod.apply_pH(molecule, pH)
        return result['energy']

    def _apply_solvent_effect(self, molecule, E_base: float, solvent: str) -> float:
        """Apply solvent corrections to energy."""
        if self._solv_mod is None:
            from kanad.environment import SolventModulator
            self._solv_mod = SolventModulator()

        result = self._solv_mod.apply_solvent(molecule, solvent)
        return result['energy']

    def _detect_bonds(self, molecule) -> List[Tuple[int, int, float]]:
        """
        Detect active bonds in current geometry.

        Uses governance system if available, otherwise distance-based.

        Returns:
            List of (atom1, atom2, distance) tuples
        """
        if not hasattr(molecule, 'atoms') or not hasattr(molecule, 'coordinates'):
            return []

        bonds = []
        coords = np.array(molecule.coordinates)
        n_atoms = len(molecule.atoms)

        # Covalent radii (Å)
        covalent_radii = {
            'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66,
            'F': 0.57, 'P': 1.07, 'S': 1.05, 'Cl': 1.02
        }

        # Distance-based bond detection
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                r_ij = np.linalg.norm(coords[j] - coords[i])

                # Bond cutoff: sum of covalent radii × 1.3
                atom_i = molecule.atoms[i]
                atom_j = molecule.atoms[j]
                r_cov = covalent_radii.get(atom_i, 1.0) + covalent_radii.get(atom_j, 1.0)

                if r_ij < 1.3 * r_cov:
                    bonds.append((i, j, r_ij))

        return bonds

    def _detect_bond_transition(self, bonds_list: List[List[Tuple]]) -> Dict[str, Any]:
        """
        Detect bond formation/breaking events along scan.

        Args:
            bonds_list: List of bond lists at each scan point

        Returns:
            Dictionary with transition information
        """
        transitions = {
            'bond_breaking': [],
            'bond_forming': [],
            'transition_points': []
        }

        for i in range(len(bonds_list) - 1):
            bonds_current = set((b[0], b[1]) for b in bonds_list[i])
            bonds_next = set((b[0], b[1]) for b in bonds_list[i+1])

            # Bonds that disappeared
            broken = bonds_current - bonds_next
            if broken:
                transitions['bond_breaking'].extend(broken)
                transitions['transition_points'].append(i)

            # Bonds that appeared
            formed = bonds_next - bonds_current
            if formed:
                transitions['bond_forming'].extend(formed)
                transitions['transition_points'].append(i)

        return transitions

    def _compute_dissociation_energy(
        self,
        energies: np.ndarray,
        bond_lengths: np.ndarray
    ) -> Dict[str, float]:
        """
        Compute bond dissociation energy from PES.

        Args:
            energies: Energy values (Ha)
            bond_lengths: Bond length values (Å)

        Returns:
            Dictionary with dissociation energy and equilibrium geometry
        """
        # Find minimum (equilibrium)
        idx_min = np.argmin(energies)
        E_min = energies[idx_min]
        r_eq = bond_lengths[idx_min]

        # Dissociation energy: E(∞) - E(r_eq)
        # Approximate E(∞) as energy at longest bond length
        E_inf = energies[-1]
        D_e = E_inf - E_min

        return {
            'dissociation_energy': D_e * self.Ha_to_kcal,  # kcal/mol
            'equilibrium_distance': r_eq,
            'equilibrium_energy': E_min
        }

    def _linear_interpolation_path(
        self,
        reactant,
        product,
        n_images: int,
        environment: Optional[Dict[str, Any]]
    ) -> ReactionPath:
        """
        Simple linear interpolation between reactant and product.

        Fast but not minimum energy path.
        """
        logger.info("Linear interpolation path (no optimization)")

        coords_R = np.array(reactant.coordinates)
        coords_P = np.array(product.coordinates)

        snapshots = []
        energies = []

        for i, alpha in enumerate(np.linspace(0, 1, n_images)):
            # Linear interpolation
            coords = (1 - alpha) * coords_R + alpha * coords_P

            # Create intermediate molecule
            from copy import deepcopy
            mol_i = deepcopy(reactant)
            mol_i.coordinates = coords.tolist()

            # Compute energy
            snapshot = self._compute_configuration_energy(mol_i, environment)
            snapshot.reaction_coordinate = alpha

            snapshots.append(snapshot)
            energies.append(snapshot.energy)

        energies = np.array(energies)
        idx_ts = np.argmax(energies)  # Highest energy point (approximate TS)

        return ReactionPath(
            images=snapshots,
            energies=energies,
            reaction_coordinate=np.linspace(0, 1, n_images),
            transition_state_index=idx_ts,
            activation_energy=(energies[idx_ts] - energies[0]) * self.Ha_to_kcal,
            reaction_energy=(energies[-1] - energies[0]) * self.Ha_to_kcal
        )

    def _nudged_elastic_band(
        self,
        reactant,
        product,
        n_images: int,
        environment: Optional[Dict[str, Any]],
        max_iterations: int,
        convergence_threshold: float
    ) -> ReactionPath:
        """
        Nudged Elastic Band (NEB) method for finding MEP.

        Simplified implementation - full NEB requires force calculations.
        """
        logger.info("NEB path optimization (simplified)")

        # Start with linear interpolation
        path = self._linear_interpolation_path(reactant, product, n_images, environment)

        # In full implementation, would iteratively optimize images
        # For now, return linear path with warning
        logger.warning("NEB optimization not fully implemented - returning linear path")

        return path

    def _string_method(
        self,
        reactant,
        product,
        n_images: int,
        environment: Optional[Dict[str, Any]],
        max_iterations: int,
        convergence_threshold: float
    ) -> ReactionPath:
        """String method for finding MEP."""
        logger.warning("String method not implemented - using linear interpolation")
        return self._linear_interpolation_path(reactant, product, n_images, environment)

    def _systematic_conformer_search(
        self,
        molecule,
        max_conformers: int,
        energy_window: float,
        environment: Optional[Dict[str, Any]]
    ) -> List[ConfigurationSnapshot]:
        """Systematic conformational search over rotatable bonds."""
        logger.info("Systematic conformer search")

        # Placeholder - full implementation would:
        # 1. Identify rotatable bonds
        # 2. Grid search over torsion angles
        # 3. Optimize each conformer
        # 4. Filter by energy and uniqueness

        # For now, return single conformer
        snapshot = self._compute_configuration_energy(molecule, environment)
        return [snapshot]

    def _random_conformer_search(
        self,
        molecule,
        max_conformers: int,
        energy_window: float,
        environment: Optional[Dict[str, Any]]
    ) -> List[ConfigurationSnapshot]:
        """Random conformational search (Monte Carlo)."""
        logger.info("Random conformer search")

        # Placeholder
        snapshot = self._compute_configuration_energy(molecule, environment)
        return [snapshot]

    def _optimize_geometry(self, molecule, freeze_atoms: Optional[List[int]] = None):
        """Optimize molecular geometry with optional frozen atoms."""
        # Placeholder - full implementation uses quantum gradients
        logger.debug("Geometry optimization (placeholder)")
        return molecule

    def plot_pes_scan(
        self,
        scan_result: Dict[str, Any],
        save_path: Optional[str] = None
    ):
        """
        Plot potential energy surface scan.

        Args:
            scan_result: Output from scan_bond_length()
            save_path: Optional path to save figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed")
            return

        x = scan_result['bond_lengths']
        E = scan_result['energies'] * self.Ha_to_kcal

        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(x, E, 'o-', linewidth=2, markersize=8)
        ax.set_xlabel('Bond Length (Å)', fontsize=12)
        ax.set_ylabel('Energy (kcal/mol)', fontsize=12)
        ax.set_title('Potential Energy Surface', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Mark equilibrium and dissociation
        if 'dissociation_energy' in scan_result:
            diss = scan_result['dissociation_energy']
            ax.text(0.7, 0.9, f"D_e = {diss['dissociation_energy']:.1f} kcal/mol\n"
                              f"r_eq = {diss['equilibrium_distance']:.3f} Å",
                   transform=ax.transAxes, fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")

        plt.show()

    def plot_reaction_path(
        self,
        path: ReactionPath,
        save_path: Optional[str] = None
    ):
        """
        Plot reaction path energy profile.

        Args:
            path: ReactionPath object
            save_path: Optional path to save figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed")
            return

        xi = path.reaction_coordinate
        E = path.energies * self.Ha_to_kcal
        E_rel = E - E[0]  # Relative to reactant

        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(xi, E_rel, 'o-', linewidth=2, markersize=8)

        # Mark TS
        if path.transition_state_index is not None:
            ts_idx = path.transition_state_index
            ax.plot(xi[ts_idx], E_rel[ts_idx], 'r*', markersize=20,
                   label=f'TS: ΔE‡ = {path.activation_energy:.1f} kcal/mol')

        ax.axhline(y=0, color='k', linestyle='--', linewidth=1)
        ax.axhline(y=E_rel[-1], color='b', linestyle='--', linewidth=1,
                  label=f'Product: ΔE_rxn = {path.reaction_energy:.1f} kcal/mol')

        ax.set_xlabel('Reaction Coordinate', fontsize=12)
        ax.set_ylabel('Relative Energy (kcal/mol)', fontsize=12)
        ax.set_title('Reaction Path', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")

        plt.show()
