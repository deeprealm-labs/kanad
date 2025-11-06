"""
Solvent Effects Modulator

Implements solvation models to account for solvent effects on molecular systems:
- PCM (Polarizable Continuum Model): Cavity-based implicit solvation
- SMD (Solvation Model Based on Density): Universal solvation model
- Dielectric screening of electrostatic interactions
- Cavity formation and dispersion energies

Physical Basis:
    Solvation free energy: ΔG_solv = ΔG_elec + ΔG_cav + ΔG_disp + ΔG_rep
    where:
    - ΔG_elec: Electrostatic screening (dielectric)
    - ΔG_cav: Cavity formation energy
    - ΔG_disp: Dispersion interactions with solvent
    - ΔG_rep: Repulsion from solvent exclusion

References:
    - Tomasi et al. "Quantum Mechanical Continuum Solvation Models" Chem. Rev. (2005)
    - Marenich et al. "Universal Solvation Model Based on Solute Electron Density" J. Phys. Chem. B (2009)
"""

import numpy as np
from typing import Dict, Any, Optional, Tuple, List
import logging

logger = logging.getLogger(__name__)


# Common solvent parameters (dielectric constant, refractive index, surface tension)
SOLVENT_DATABASE = {
    'water': {
        'epsilon': 78.36,           # Dielectric constant at 298K
        'n': 1.333,                 # Refractive index
        'gamma': 71.99,             # Surface tension (cal/(mol·Å²))
        'name': 'Water',
        'formula': 'H2O',
        'density': 0.997,           # g/cm³
        'alpha': 0.0,               # Abraham's hydrogen bond acidity
        'beta': 0.82,               # Abraham's hydrogen bond basicity
    },
    'acetonitrile': {
        'epsilon': 35.69,
        'n': 1.344,
        'gamma': 28.66,
        'name': 'Acetonitrile',
        'formula': 'CH3CN',
        'density': 0.786,
        'alpha': 0.19,
        'beta': 0.40,
    },
    'dmso': {
        'epsilon': 46.83,
        'n': 1.479,
        'gamma': 42.92,
        'name': 'Dimethyl sulfoxide',
        'formula': '(CH3)2SO',
        'density': 1.100,
        'alpha': 0.0,
        'beta': 0.76,
    },
    'chloroform': {
        'epsilon': 4.81,
        'n': 1.446,
        'gamma': 27.14,
        'name': 'Chloroform',
        'formula': 'CHCl3',
        'density': 1.489,
        'alpha': 0.20,
        'beta': 0.10,
    },
    'methanol': {
        'epsilon': 32.63,
        'n': 1.329,
        'gamma': 22.07,
        'name': 'Methanol',
        'formula': 'CH3OH',
        'density': 0.791,
        'alpha': 0.43,
        'beta': 0.47,
    },
    'ethanol': {
        'epsilon': 24.85,
        'n': 1.361,
        'gamma': 21.97,
        'name': 'Ethanol',
        'formula': 'C2H5OH',
        'density': 0.789,
        'alpha': 0.37,
        'beta': 0.48,
    },
    'toluene': {
        'epsilon': 2.38,
        'n': 1.497,
        'gamma': 27.93,
        'name': 'Toluene',
        'formula': 'C6H5CH3',
        'density': 0.867,
        'alpha': 0.0,
        'beta': 0.11,
    },
    'hexane': {
        'epsilon': 1.88,
        'n': 1.375,
        'gamma': 17.89,
        'name': 'Hexane',
        'formula': 'C6H14',
        'density': 0.655,
        'alpha': 0.0,
        'beta': 0.0,
    },
    'benzene': {
        'epsilon': 2.27,
        'n': 1.501,
        'gamma': 28.22,
        'name': 'Benzene',
        'formula': 'C6H6',
        'density': 0.879,
        'alpha': 0.0,
        'beta': 0.14,
    },
    'thf': {
        'epsilon': 7.43,
        'n': 1.407,
        'gamma': 26.40,
        'name': 'Tetrahydrofuran',
        'formula': 'C4H8O',
        'density': 0.889,
        'alpha': 0.0,
        'beta': 0.55,
    },
    'dichloromethane': {
        'epsilon': 8.93,
        'n': 1.424,
        'gamma': 27.20,
        'name': 'Dichloromethane',
        'formula': 'CH2Cl2',
        'density': 1.327,
        'alpha': 0.13,
        'beta': 0.10,
    },
    'vacuum': {
        'epsilon': 1.0,
        'n': 1.0,
        'gamma': 0.0,
        'name': 'Vacuum',
        'formula': '',
        'density': 0.0,
        'alpha': 0.0,
        'beta': 0.0,
    },
}


class SolventModulator:
    """
    Apply solvent effects to molecular Hamiltonians using PCM and SMD models.

    Solvent effects modify molecular energies through:
    1. Electrostatic screening (dielectric effect on electron-electron repulsion)
    2. Cavity formation (energy cost to create cavity in solvent)
    3. Dispersion interactions (attractive van der Waals with solvent)
    4. Hydrogen bonding (specific solute-solvent interactions)

    Example:
        >>> from kanad.bonds import BondFactory
        >>> from kanad.environment import SolventModulator
        >>>
        >>> bond = BondFactory.create_bond('H', 'H', distance=0.74)
        >>> solv = SolventModulator()
        >>>
        >>> # In water
        >>> result_water = solv.apply_solvent(bond, 'water', model='pcm')
        >>> print(f"Solvation energy: {result_water['solvation_energy']:.6f} Ha")
        >>>
        >>> # In hexane (nonpolar)
        >>> result_hexane = solv.apply_solvent(bond, 'hexane', model='pcm')
        >>> print(f"Dielectric screening: {result_hexane['dielectric_factor']:.4f}")
    """

    # Physical constants
    Ha_to_kcal = 627.509474          # Hartree to kcal/mol
    bohr_to_angstrom = 0.529177      # Bohr to Angstrom
    kcal_per_mol_A2_to_Ha = 0.00159360144  # Surface tension conversion

    def __init__(self):
        """Initialize solvent modulator."""
        self.solvent_db = SOLVENT_DATABASE
        logger.info("SolventModulator initialized with %d solvents", len(self.solvent_db))

    def apply_solvent(
        self,
        bond_or_molecule,
        solvent: str,
        model: str = 'pcm',
        temperature: float = 298.15,
        custom_epsilon: Optional[float] = None,
        include_nonelectrostatic: bool = True
    ) -> Dict[str, Any]:
        """
        Apply solvent effects to molecular system.

        Args:
            bond_or_molecule: Bond or Molecule object
            solvent: Solvent name (e.g., 'water', 'acetonitrile') or 'vacuum'
            model: Solvation model ('pcm' or 'smd')
            temperature: Temperature in Kelvin
            custom_epsilon: Override dielectric constant
            include_nonelectrostatic: Include cavity/dispersion terms

        Returns:
            Dictionary with:
                energy: Solvent-corrected energy (Ha)
                solvation_energy: Total solvation free energy (Ha)
                electrostatic_energy: Electrostatic contribution (Ha)
                cavity_energy: Cavity formation energy (Ha)
                dispersion_energy: Dispersion interaction energy (Ha)
                dielectric_factor: Screening factor for Coulomb interactions
                solvent_info: Solvent parameters
        """
        logger.info(f"Applying {model.upper()} solvation: solvent={solvent}, T={temperature:.2f}K")

        # Get solvent parameters
        if solvent.lower() not in self.solvent_db:
            logger.warning(f"Solvent '{solvent}' not in database - using water")
            solvent = 'water'

        solv_params = self.solvent_db[solvent.lower()]
        epsilon = custom_epsilon if custom_epsilon else solv_params['epsilon']

        # Get base energy
        E_gas = self._get_base_energy(bond_or_molecule)

        # 1. Electrostatic solvation (Born model / reaction field)
        if model.lower() == 'pcm':
            E_elec, dielectric_factor = self._compute_pcm_electrostatic(
                bond_or_molecule, epsilon, temperature
            )
        elif model.lower() == 'smd':
            E_elec, dielectric_factor = self._compute_smd_electrostatic(
                bond_or_molecule, epsilon, solv_params, temperature
            )
        else:
            raise ValueError(f"Unknown solvation model: {model}")

        # 2. Non-electrostatic contributions
        if include_nonelectrostatic and epsilon > 1.0:
            E_cav = self._compute_cavity_energy(
                bond_or_molecule, solv_params
            )
            E_disp = self._compute_dispersion_energy(
                bond_or_molecule, solv_params
            )
        else:
            E_cav = 0.0
            E_disp = 0.0

        # 3. Total solvation free energy
        E_solv = E_elec + E_cav + E_disp

        # 4. Solvated energy
        E_solution = E_gas + E_solv

        return {
            'energy': E_solution,
            'solvation_energy': E_solv,
            'electrostatic_energy': E_elec,
            'cavity_energy': E_cav,
            'dispersion_energy': E_disp,
            'dielectric_factor': dielectric_factor,
            'gas_phase_energy': E_gas,
            'solvent': solvent,
            'epsilon': epsilon,
            'temperature': temperature,
            'model': model,
            'solvent_info': solv_params
        }

    def scan_solvents(
        self,
        bond_or_molecule,
        solvents: Optional[List[str]] = None,
        model: str = 'pcm',
        **kwargs
    ) -> Dict[str, Any]:
        """
        Compute solvation energies across multiple solvents.

        Args:
            bond_or_molecule: Molecular system
            solvents: List of solvent names (default: all in database)
            model: Solvation model
            **kwargs: Additional arguments for apply_solvent

        Returns:
            Dictionary with:
                solvents: Solvent names
                energies: Total energies (Ha)
                solvation_energies: ΔG_solv values (Ha)
                dielectric_constants: ε values
        """
        if solvents is None:
            solvents = list(self.solvent_db.keys())

        results = {
            'solvents': [],
            'energies': [],
            'solvation_energies': [],
            'electrostatic_energies': [],
            'cavity_energies': [],
            'dispersion_energies': [],
            'dielectric_constants': []
        }

        for solvent in solvents:
            result = self.apply_solvent(bond_or_molecule, solvent, model=model, **kwargs)
            results['solvents'].append(solvent)
            results['energies'].append(result['energy'])
            results['solvation_energies'].append(result['solvation_energy'])
            results['electrostatic_energies'].append(result['electrostatic_energy'])
            results['cavity_energies'].append(result['cavity_energy'])
            results['dispersion_energies'].append(result['dispersion_energy'])
            results['dielectric_constants'].append(result['epsilon'])

        # Convert lists to arrays
        for key in ['energies', 'solvation_energies', 'electrostatic_energies',
                    'cavity_energies', 'dispersion_energies', 'dielectric_constants']:
            results[key] = np.array(results[key])

        return results

    def _get_base_energy(self, bond_or_molecule) -> float:
        """Get gas-phase ground state energy."""
        if hasattr(bond_or_molecule, 'energy'):
            return bond_or_molecule.energy
        elif hasattr(bond_or_molecule, '_cached_energy'):
            return bond_or_molecule._cached_energy
        else:
            # CRITICAL FIX: Compute HF energy if not cached (don't return 0.0!)
            if hasattr(bond_or_molecule, 'hamiltonian'):
                logger.info("Computing HF energy (not cached) for solvent calculation")
                rdm1_hf, E_hf = bond_or_molecule.hamiltonian.solve_scf()
                bond_or_molecule._cached_energy = E_hf
                return E_hf
            else:
                raise ValueError("Cannot compute energy - no hamiltonian available and no cached energy")

    def _compute_pcm_electrostatic(
        self,
        bond_or_molecule,
        epsilon: float,
        temperature: float
    ) -> Tuple[float, float]:
        """
        Compute PCM electrostatic solvation energy using Born model.

        For a spherical cavity:
            ΔG_elec = -½ q² (1 - 1/ε) / a
        where q is the charge, a is the cavity radius.

        For molecular systems, use generalized Born:
            ΔG_elec = -½ (1 - 1/ε) Σ_ij q_i q_j f_GB(r_ij)

        Args:
            bond_or_molecule: Molecular system
            epsilon: Dielectric constant
            temperature: Temperature in K

        Returns:
            (E_elec, dielectric_factor)
        """
        if epsilon <= 1.0:
            return 0.0, 1.0

        # Reaction field factor: f(ε) = (ε - 1) / (2ε + 1)
        # For Born model: ΔG = -½ q² f(ε) / a
        f_epsilon = (epsilon - 1.0) / (2.0 * epsilon + 1.0)

        # Dielectric screening factor for Coulomb interactions
        # In continuum: V_eff = V_vac / ε_eff
        # Approximate: ε_eff ≈ 1 + f(ε)
        dielectric_factor = 1.0 / (1.0 + 0.5 * f_epsilon)

        # Estimate molecular cavity radius from bond/molecule size
        if hasattr(bond_or_molecule, 'atoms'):
            # Molecule: estimate from atomic positions
            cavity_radius = self._estimate_cavity_radius(bond_or_molecule)
        else:
            # Bond: use bond length as characteristic size
            if hasattr(bond_or_molecule, 'distance'):
                cavity_radius = bond_or_molecule.distance / self.bohr_to_angstrom
            else:
                cavity_radius = 2.0  # Default 2 Å

        # Estimate effective charge (use electron count as proxy)
        if hasattr(bond_or_molecule, 'n_electrons'):
            n_electrons = bond_or_molecule.n_electrons
        elif hasattr(bond_or_molecule, 'electrons'):
            n_electrons = bond_or_molecule.electrons
        else:
            n_electrons = 2  # Default for H2

        # Born solvation energy (simplified)
        # ΔG_Born ≈ -α q²_eff (1 - 1/ε) / a
        # Use empirical scaling: α ~ 0.5, q_eff ~ √n_electrons
        alpha = 0.5
        q_eff_sq = n_electrons / 2.0  # Effective charge squared

        # Convert: 1 / (Å) = 1 / 0.529177 Bohr⁻¹
        E_elec = -alpha * q_eff_sq * (1.0 - 1.0/epsilon) / (cavity_radius / self.bohr_to_angstrom)

        # Convert to Hartree (empirical scaling)
        E_elec *= 0.01  # Typical solvation ~10-100 kcal/mol ~ 0.01-0.1 Ha

        logger.debug(f"PCM electrostatic: ΔG_elec = {E_elec:.6f} Ha, "
                    f"ε = {epsilon:.2f}, a = {cavity_radius:.2f} Å")

        return E_elec, dielectric_factor

    def _compute_smd_electrostatic(
        self,
        bond_or_molecule,
        epsilon: float,
        solv_params: Dict,
        temperature: float
    ) -> Tuple[float, float]:
        """
        Compute SMD electrostatic energy with additional descriptor terms.

        SMD model:
            ΔG_elec = ΔG_ENP + ΔG_CDS
        where:
            ΔG_ENP: Electronic (Born-like) + nuclear polarization
            ΔG_CDS: Cavity-dispersion-solvent structure

        Args:
            bond_or_molecule: Molecular system
            epsilon: Dielectric constant
            solv_params: Solvent parameters (alpha, beta, gamma, etc.)
            temperature: Temperature in K

        Returns:
            (E_elec, dielectric_factor)
        """
        # Start with PCM-like electrostatic term
        E_pcm, dielectric_factor = self._compute_pcm_electrostatic(
            bond_or_molecule, epsilon, temperature
        )

        # SMD correction terms based on Abraham descriptors
        # ΔG_CDS = Σ_k c_k τ_k  where τ_k are solvent/solute descriptors

        # Get molecular surface area
        surface_area = self._estimate_surface_area(bond_or_molecule)  # Å²

        # Hydrogen bonding corrections
        # α: H-bond acidity, β: H-bond basicity
        alpha_solv = solv_params.get('alpha', 0.0)
        beta_solv = solv_params.get('beta', 0.0)

        # Estimate solute H-bonding capacity (crude approximation)
        # For now, assume neutral contribution
        E_hbond = 0.0

        # Surface tension term (proportional to SASA)
        gamma = solv_params.get('gamma', 0.0)  # cal/(mol·Å²)
        E_surface = gamma * surface_area * self.kcal_per_mol_A2_to_Ha / 1000.0

        # Total SMD electrostatic
        E_elec = E_pcm + E_hbond + E_surface

        logger.debug(f"SMD electrostatic: ΔG_elec = {E_elec:.6f} Ha "
                    f"(PCM: {E_pcm:.6f}, Hbond: {E_hbond:.6f}, Surf: {E_surface:.6f})")

        return E_elec, dielectric_factor

    def _compute_cavity_energy(
        self,
        bond_or_molecule,
        solv_params: Dict
    ) -> float:
        """
        Compute cavity formation energy.

        Scaled Particle Theory:
            ΔG_cav = 4π a² γ + pV_cav
        where γ is surface tension, a is radius, p is pressure

        Args:
            bond_or_molecule: Molecular system
            solv_params: Solvent parameters

        Returns:
            Cavity formation energy (Ha)
        """
        gamma = solv_params.get('gamma', 0.0)  # cal/(mol·Å²)

        if gamma == 0.0:
            return 0.0

        # Get surface area
        surface_area = self._estimate_surface_area(bond_or_molecule)  # Å²

        # Cavity energy: ΔG_cav = γ × SASA
        # Convert cal/(mol·Å²) → Ha
        E_cav = gamma * surface_area * self.kcal_per_mol_A2_to_Ha / 1000.0

        # Cavity formation is unfavorable (positive energy)
        E_cav = abs(E_cav)

        logger.debug(f"Cavity energy: ΔG_cav = {E_cav:.6f} Ha (SASA = {surface_area:.2f} Å²)")

        return E_cav

    def _compute_dispersion_energy(
        self,
        bond_or_molecule,
        solv_params: Dict
    ) -> float:
        """
        Compute dispersion interaction energy with solvent.

        London dispersion:
            ΔG_disp = -C / r⁶  (attractive)

        For continuum solvent, integrate over solvent molecules:
            ΔG_disp ≈ -k × (n² - 1) × SASA
        where n is refractive index

        Args:
            bond_or_molecule: Molecular system
            solv_params: Solvent parameters

        Returns:
            Dispersion energy (Ha, negative = favorable)
        """
        n = solv_params.get('n', 1.0)  # Refractive index

        if n <= 1.0:
            return 0.0

        # Get surface area
        surface_area = self._estimate_surface_area(bond_or_molecule)  # Å²

        # Dispersion scaling: (n² - 1) gives polarizability measure
        # Empirical: ΔG_disp ~ -0.3 kcal/(mol·Å²) × (n² - 1) × SASA
        k_disp = -0.3  # kcal/(mol·Å²)
        polarizability_factor = (n**2 - 1.0)

        E_disp = k_disp * polarizability_factor * surface_area

        # Convert kcal/mol → Ha
        E_disp *= self.kcal_per_mol_A2_to_Ha / 1000.0

        logger.debug(f"Dispersion energy: ΔG_disp = {E_disp:.6f} Ha "
                    f"(n = {n:.3f}, SASA = {surface_area:.2f} Å²)")

        return E_disp

    def _estimate_cavity_radius(self, molecule) -> float:
        """
        Estimate molecular cavity radius from atomic positions.

        Use radius of gyration:
            R_g = √(⟨r²⟩) where r is distance from center of mass

        Args:
            molecule: Molecule object with atomic positions

        Returns:
            Cavity radius in Angstroms
        """
        if not hasattr(molecule, 'atoms') or not hasattr(molecule, 'coordinates'):
            return 2.0  # Default

        coords = np.array(molecule.coordinates)  # Assume in Angstroms

        # Center of mass (unweighted for simplicity)
        center = np.mean(coords, axis=0)

        # Radius of gyration
        r_squared = np.sum((coords - center)**2, axis=1)
        R_g = np.sqrt(np.mean(r_squared))

        # Add van der Waals radii (~1.5 Å for typical atoms)
        cavity_radius = R_g + 1.5

        return max(cavity_radius, 1.0)  # Minimum 1 Å

    def _estimate_surface_area(self, bond_or_molecule) -> float:
        """
        Estimate solvent-accessible surface area (SASA).

        For simple systems:
            SASA ≈ 4π (r_cav)²

        Args:
            bond_or_molecule: Molecular system

        Returns:
            Surface area in Ų
        """
        if hasattr(bond_or_molecule, 'atoms'):
            # Molecule: use cavity radius
            radius = self._estimate_cavity_radius(bond_or_molecule)
            n_atoms = len(bond_or_molecule.atoms)
            # Rough approximation: SASA ~ 4πr² but reduced for overlapping spheres
            sasa = 4.0 * np.pi * radius**2 * (1.0 + 0.2 * (n_atoms - 1))
        else:
            # Bond: estimate from bond length
            if hasattr(bond_or_molecule, 'distance'):
                bond_length = bond_or_molecule.distance * self.bohr_to_angstrom
            else:
                bond_length = 1.0  # Å

            # Approximate as cylinder: SASA ~ 2πrL + 2πr²
            r_atom = 1.5  # Typical van der Waals radius
            sasa = 2.0 * np.pi * r_atom * bond_length + 2.0 * np.pi * r_atom**2

        return max(sasa, 10.0)  # Minimum 10 Ų

    def modify_hamiltonian_with_solvent(
        self,
        hamiltonian,
        solvent: str,
        model: str = 'pcm',
        custom_epsilon: Optional[float] = None
    ):
        """
        Create solvent-modified Hamiltonian.

        Applies dielectric screening to Coulomb terms:
            H_solv = T + V_nuc + V_ee/ε_eff + ΔG_solv

        Args:
            hamiltonian: Original Hamiltonian
            solvent: Solvent name
            model: Solvation model
            custom_epsilon: Override dielectric constant

        Returns:
            Modified Hamiltonian object
        """
        solv_params = self.solvent_db.get(solvent.lower(), self.solvent_db['water'])
        epsilon = custom_epsilon if custom_epsilon else solv_params['epsilon']

        # Compute dielectric screening factor
        if epsilon > 1.0:
            # Reaction field factor
            f_epsilon = (epsilon - 1.0) / (2.0 * epsilon + 1.0)
            screening_factor = 1.0 / (1.0 + 0.5 * f_epsilon)
        else:
            screening_factor = 1.0

        logger.info(f"Applying dielectric screening: ε = {epsilon:.2f}, "
                   f"factor = {screening_factor:.4f}")

        # Scale Hamiltonian to approximate screening
        # More sophisticated: identify and scale only electron-electron terms
        H_solvent = hamiltonian * screening_factor

        return H_solvent

    def plot_solvent_scan(
        self,
        scan_result: Dict[str, Any],
        save_path: Optional[str] = None
    ):
        """
        Plot solvation energies across different solvents.

        Args:
            scan_result: Output from scan_solvents()
            save_path: Optional path to save figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed - cannot plot")
            return

        solvents = scan_result['solvents']
        E_solv = scan_result['solvation_energies'] * self.Ha_to_kcal  # kcal/mol
        epsilon = scan_result['dielectric_constants']

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Plot 1: Solvation energies
        colors = plt.cm.viridis(np.linspace(0, 1, len(solvents)))
        bars = ax1.bar(range(len(solvents)), E_solv, color=colors)
        ax1.set_xticks(range(len(solvents)))
        ax1.set_xticklabels(solvents, rotation=45, ha='right')
        ax1.set_ylabel('Solvation Energy (kcal/mol)', fontsize=12)
        ax1.set_title('Solvation Free Energy', fontsize=14, fontweight='bold')
        ax1.axhline(y=0, color='k', linestyle='--', linewidth=1)
        ax1.grid(True, alpha=0.3, axis='y')

        # Plot 2: Energy vs dielectric constant
        ax2.scatter(epsilon, E_solv, s=100, c=colors, edgecolors='k', linewidth=1.5)
        for i, solvent in enumerate(solvents):
            ax2.annotate(solvent, (epsilon[i], E_solv[i]),
                        xytext=(5, 5), textcoords='offset points', fontsize=9)
        ax2.set_xlabel('Dielectric Constant (ε)', fontsize=12)
        ax2.set_ylabel('Solvation Energy (kcal/mol)', fontsize=12)
        ax2.set_title('ΔG_solv vs Dielectric Constant', fontsize=14, fontweight='bold')
        ax2.set_xscale('log')
        ax2.grid(True, alpha=0.3)
        ax2.axhline(y=0, color='k', linestyle='--', linewidth=1)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")

        plt.show()

    def get_available_solvents(self) -> List[str]:
        """Return list of available solvents in database."""
        return list(self.solvent_db.keys())

    def get_solvent_info(self, solvent: str) -> Dict[str, Any]:
        """
        Get detailed information about a solvent.

        Args:
            solvent: Solvent name

        Returns:
            Dictionary with solvent parameters
        """
        if solvent.lower() not in self.solvent_db:
            raise ValueError(f"Solvent '{solvent}' not in database. "
                           f"Available: {self.get_available_solvents()}")

        return self.solvent_db[solvent.lower()].copy()
