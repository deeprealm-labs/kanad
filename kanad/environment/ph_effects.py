"""
pH Effects Modulator

Handles pH-dependent protonation states and their effects on molecular energies:
- Protonation/deprotonation equilibria
- Henderson-Hasselbalch equation for ionization states
- pH-dependent charge distribution
- Electrostatic energy changes with protonation
- Tautomeric equilibria

Physical Basis:
    At pH, the fraction of protonated species follows Henderson-Hasselbalch:
        α_HA = 1 / (1 + 10^(pH - pKa))
        α_A⁻ = 1 / (1 + 10^(pKa - pH))

    Free energy of protonation:
        ΔG(pH) = ΔG°(pKa) + RT ln(10) × (pH - pKa)

Critical for:
    - Drug binding (charged vs neutral forms have vastly different affinities)
    - Enzyme catalysis (active site residues must be in correct protonation state)
    - Membrane permeability (neutral form crosses membranes, charged does not)
    - Protein stability (protonation affects folding)

References:
    - Jensen et al. "Very fast empirical prediction and rationalization of protein pKa values" Proteins (2005)
    - Ullmann et al. "Computational Simulations of Realistic Systems" J. Phys. Chem. B (2003)
"""

import numpy as np
from typing import Dict, Any, Optional, List, Tuple
import logging

logger = logging.getLogger(__name__)


# Common pKa values for functional groups (in water, 298K)
PKA_DATABASE = {
    # Acids
    'carboxylic_acid': {'pKa': 4.8, 'acidic_group': 'COOH', 'basic_group': 'COO-'},
    'phenol': {'pKa': 10.0, 'acidic_group': 'OH', 'basic_group': 'O-'},
    'thiol': {'pKa': 8.3, 'acidic_group': 'SH', 'basic_group': 'S-'},
    'phosphate': {'pKa': 7.2, 'acidic_group': 'PO4H2-', 'basic_group': 'PO4H2-'},

    # Bases
    'amine_primary': {'pKa': 10.6, 'acidic_group': 'NH3+', 'basic_group': 'NH2'},
    'amine_secondary': {'pKa': 10.7, 'acidic_group': 'NH2+', 'basic_group': 'NH'},
    'amine_tertiary': {'pKa': 9.8, 'acidic_group': 'NH+', 'basic_group': 'N'},
    'imidazole': {'pKa': 6.0, 'acidic_group': 'Im-H+', 'basic_group': 'Im'},  # Histidine
    'guanidinium': {'pKa': 12.5, 'acidic_group': 'Gu+', 'basic_group': 'Gu'},  # Arginine

    # Amino acids (protein residues)
    'asp': {'pKa': 3.9, 'acidic_group': 'COOH', 'basic_group': 'COO-', 'name': 'Aspartate'},
    'glu': {'pKa': 4.3, 'acidic_group': 'COOH', 'basic_group': 'COO-', 'name': 'Glutamate'},
    'his': {'pKa': 6.0, 'acidic_group': 'Im-H+', 'basic_group': 'Im', 'name': 'Histidine'},
    'cys': {'pKa': 8.3, 'acidic_group': 'SH', 'basic_group': 'S-', 'name': 'Cysteine'},
    'tyr': {'pKa': 10.1, 'acidic_group': 'OH', 'basic_group': 'O-', 'name': 'Tyrosine'},
    'lys': {'pKa': 10.5, 'acidic_group': 'NH3+', 'basic_group': 'NH2', 'name': 'Lysine'},
    'arg': {'pKa': 12.5, 'acidic_group': 'Gu+', 'basic_group': 'Gu', 'name': 'Arginine'},

    # N-terminus and C-terminus
    'n_terminus': {'pKa': 9.6, 'acidic_group': 'NH3+', 'basic_group': 'NH2'},
    'c_terminus': {'pKa': 2.3, 'acidic_group': 'COOH', 'basic_group': 'COO-'},
}


class ProtonationSite:
    """
    Represents a single protonatable site in a molecule.

    Attributes:
        atom_index: Index of the protonatable atom
        pKa: Acid dissociation constant
        group_type: Type of functional group
        is_acid: True if deprotonation (HA → A⁻), False if protonation (B + H⁺ → BH⁺)
        charge_change: Change in formal charge upon protonation (+1 or -1)
    """

    def __init__(
        self,
        atom_index: int,
        pKa: float,
        group_type: str,
        is_acid: bool = True
    ):
        self.atom_index = atom_index
        self.pKa = pKa
        self.group_type = group_type
        self.is_acid = is_acid
        self.charge_change = -1 if is_acid else +1  # Charge change upon deprotonation

    def get_protonated_fraction(self, pH: float) -> float:
        """
        Calculate fraction in protonated state using Henderson-Hasselbalch.

        For acid HA ⇌ A⁻ + H⁺:
            α_HA = 1 / (1 + 10^(pH - pKa))
        For base B + H⁺ ⇌ BH⁺:
            α_BH+ = 1 / (1 + 10^(pH - pKa))

        Args:
            pH: Solution pH

        Returns:
            Fraction in protonated form (0 to 1)
        """
        x = pH - self.pKa
        # Avoid overflow for extreme pH
        if x > 10:
            return 0.0
        elif x < -10:
            return 1.0
        else:
            return 1.0 / (1.0 + 10**x)

    def get_charge_state(self, pH: float) -> float:
        """
        Get average charge at given pH.

        Args:
            pH: Solution pH

        Returns:
            Average charge (continuous value between ionized and neutral)
        """
        f_prot = self.get_protonated_fraction(pH)

        if self.is_acid:
            # HA (neutral) ⇌ A⁻ (charge -1)
            charge_protonated = 0.0
            charge_deprotonated = -1.0
        else:
            # B (neutral) ⇌ BH⁺ (charge +1)
            charge_protonated = +1.0
            charge_deprotonated = 0.0

        return f_prot * charge_protonated + (1.0 - f_prot) * charge_deprotonated


class pHModulator:
    """
    Apply pH effects to molecular Hamiltonians through protonation state changes.

    pH affects:
    1. Charge distribution: Protonation changes formal charges
    2. Electrostatic energy: Charged species have different energies
    3. Hydrogen bonding: Protonation states affect H-bond donors/acceptors
    4. Solvation energy: Charged species have higher solvation
    5. Molecular geometry: Protonation can trigger conformational changes

    Example:
        >>> from kanad.environment import pHModulator
        >>>
        >>> # Create pH modulator
        >>> ph_mod = pHModulator()
        >>>
        >>> # Add protonatable sites (e.g., carboxylic acid)
        >>> ph_mod.add_site(atom_index=0, group_type='carboxylic_acid')
        >>>
        >>> # Compute pH-dependent properties
        >>> result_7 = ph_mod.apply_pH(molecule, pH=7.4)  # Physiological
        >>> result_2 = ph_mod.apply_pH(molecule, pH=2.0)  # Acidic
        >>>
        >>> print(f"Charge at pH 7.4: {result_7['net_charge']:.2f}")
        >>> print(f"Charge at pH 2.0: {result_2['net_charge']:.2f}")
    """

    # Physical constants
    RT_298K_Ha = 3.1668115634556076e-6 * 298.15  # kT at 298K in Hartree
    Ha_to_kcal = 627.509474

    def __init__(self):
        """Initialize pH modulator."""
        self.sites: List[ProtonationSite] = []
        self.pka_db = PKA_DATABASE
        logger.info("pHModulator initialized")

    def add_site(
        self,
        atom_index: int,
        group_type: str,
        custom_pKa: Optional[float] = None
    ):
        """
        Add a protonatable site to the molecule.

        Args:
            atom_index: Index of protonatable atom
            group_type: Type of functional group (e.g., 'carboxylic_acid')
            custom_pKa: Override database pKa value
        """
        if group_type not in self.pka_db:
            logger.warning(f"Group type '{group_type}' not in database. Available: {list(self.pka_db.keys())}")
            pKa = 7.0  # Default neutral
            is_acid = True
        else:
            group_data = self.pka_db[group_type]
            pKa = custom_pKa if custom_pKa else group_data['pKa']
            # Determine if acid or base from group type
            is_acid = 'acid' in group_type or group_type in ['asp', 'glu', 'cys', 'tyr', 'c_terminus']

        site = ProtonationSite(atom_index, pKa, group_type, is_acid)
        self.sites.append(site)
        logger.info(f"Added protonation site: atom {atom_index}, pKa={pKa:.1f}, type={group_type}")

    def add_sites_from_molecule(self, molecule):
        """
        Automatically detect and add protonatable sites from molecule structure.

        Args:
            molecule: Molecule object with atom types and connectivity

        Note:
            Automatic site detection requires SMARTS pattern matching.
            For now, users must manually add sites using add_site().

        Example:
            >>> ph_model = pHEffectsModel()
            >>> ph_model.add_site(site_id=0, pka=4.76, site_type='carboxyl')
            >>> ph_model.add_site(site_id=1, pka=9.25, site_type='amine')
        """
        logger.info("Automatic site detection not yet implemented")
        logger.info("Please add protonatable sites manually using add_site()")
        logger.info("  Example: model.add_site(site_id=0, pka=4.76, site_type='carboxyl')")

        # Future implementation would:
        # 1. Parse molecular structure
        # 2. Identify functional groups using SMARTS patterns
        # 3. Assign pKa values from database
        # 4. Add sites automatically
        pass

    def determine_protonation_state(
        self,
        molecule,
        pH: float,
        return_detailed: bool = False
    ) -> Dict[int, bool]:
        """
        Determine protonation state based on pKa values using Henderson-Hasselbalch.

        For each protonatable site:
            pH = pKa + log([A⁻]/[HA])

        Rearranging:
            fraction_protonated = 1 / (1 + 10^(pH - pKa))

        If fraction > 0.5, site is protonated; otherwise deprotonated.

        Args:
            molecule: Molecular system (can be bond or molecule object)
            pH: Solution pH (0-14)
            return_detailed: If True, return detailed state info; if False, return bool only

        Returns:
            Dictionary mapping site index → protonation state
                If return_detailed=False: {site_idx: True/False}
                If return_detailed=True:  {site_idx: {'protonated': bool,
                                                       'fraction': float,
                                                       'pKa': float,
                                                       'group_type': str}}

        Example:
            >>> ph_mod = pHModulator()
            >>> ph_mod.add_site(atom_index=0, group_type='carboxylic_acid')  # pKa=4.8
            >>> state = ph_mod.determine_protonation_state(molecule, pH=7.0)
            >>> print(state)  # {0: False} - deprotonated at pH 7
            >>>
            >>> state = ph_mod.determine_protonation_state(molecule, pH=2.0)
            >>> print(state)  # {0: True} - protonated at pH 2
        """
        if not self.sites:
            logger.warning("No protonation sites defined - returning empty state")
            return {}

        protonation_state = {}

        for site in self.sites:
            # Compute protonated fraction using Henderson-Hasselbalch
            f_prot = site.get_protonated_fraction(pH)

            # Determine discrete state (>50% = protonated)
            is_protonated = f_prot > 0.5

            if return_detailed:
                protonation_state[site.atom_index] = {
                    'protonated': is_protonated,
                    'fraction': f_prot,
                    'pKa': site.pKa,
                    'group_type': site.group_type,
                    'is_acid': site.is_acid
                }
            else:
                protonation_state[site.atom_index] = is_protonated

        logger.debug(f"Protonation state at pH {pH:.1f}: {protonation_state}")

        return protonation_state

    def apply_pH(
        self,
        bond_or_molecule,
        pH: float,
        temperature: float = 298.15,
        solvent: str = 'water',
        include_solvation: bool = True
    ) -> Dict[str, Any]:
        """
        Apply pH effects to molecular system.

        Args:
            bond_or_molecule: Bond or Molecule object
            pH: Solution pH (0-14, typically 1-14)
            temperature: Temperature in Kelvin
            solvent: Solvent (affects solvation of charged species)
            include_solvation: Include Born solvation for charged states

        Returns:
            Dictionary with:
                energy: pH-corrected energy (Ha)
                protonation_state: Dict of site → protonation fraction
                net_charge: Average net charge
                microstate_populations: Populations of ionization microstates
                free_energy: ΔG(pH) including entropy
                major_species: Most populated protonation microstate
        """
        logger.info(f"Applying pH effects: pH = {pH:.2f}, T = {temperature:.2f}K")

        if not self.sites:
            logger.warning("No protonation sites defined - pH will have no effect")
            E_base = self._get_base_energy(bond_or_molecule)
            return {
                'energy': E_base,
                'protonation_state': {},
                'net_charge': 0.0,
                'microstate_populations': {},
                'free_energy': E_base,
                'major_species': 'neutral',
                'pH': pH,
                'temperature': temperature
            }

        # Get base energy (reference state)
        E_base = self._get_base_energy(bond_or_molecule)

        # 1. Compute protonation fractions for each site
        protonation_state = {}
        net_charge = 0.0
        for site in self.sites:
            f_prot = site.get_protonated_fraction(pH)
            charge = site.get_charge_state(pH)
            protonation_state[site.atom_index] = {
                'protonated_fraction': f_prot,
                'charge': charge,
                'pKa': site.pKa,
                'group_type': site.group_type
            }
            net_charge += charge

        # 2. Compute protonation free energy
        # ΔG(pH) = ΔG°(pKa) + RT ln(10) × (pH - pKa)
        # Sum over all sites
        RT_ln10 = self.RT_298K_Ha * (temperature / 298.15) * np.log(10)

        Delta_G_protonation = 0.0
        for site in self.sites:
            Delta_pH = pH - site.pKa
            # Free energy of ionization
            Delta_G_site = RT_ln10 * Delta_pH
            # Weight by population
            f_ionized = 1.0 - site.get_protonated_fraction(pH)
            Delta_G_protonation += f_ionized * Delta_G_site

        # 3. Electrostatic penalty for charged states
        # Charged species have higher energy (less stable in vacuum)
        # Born model: ΔG_charge ∝ q² / r
        if include_solvation and abs(net_charge) > 0.1:
            Delta_G_charge = self._compute_charging_energy(
                bond_or_molecule, net_charge, solvent
            )
        else:
            Delta_G_charge = 0.0

        # 4. Total free energy at this pH
        G_pH = E_base + Delta_G_protonation + Delta_G_charge

        # 5. Generate microstate populations (if multiple sites)
        if len(self.sites) > 1 and len(self.sites) <= 5:  # Only for small number of sites
            microstate_pops = self._compute_microstate_populations(pH)
            major_species = max(microstate_pops, key=microstate_pops.get)
        else:
            microstate_pops = {'averaged': 1.0}
            major_species = 'averaged'

        return {
            'energy': G_pH,
            'protonation_state': protonation_state,
            'net_charge': net_charge,
            'microstate_populations': microstate_pops,
            'free_energy': G_pH,
            'protonation_free_energy': Delta_G_protonation,
            'charging_energy': Delta_G_charge,
            'major_species': major_species,
            'pH': pH,
            'temperature': temperature,
            'solvent': solvent
        }

    def scan_pH(
        self,
        bond_or_molecule,
        pH_range: Tuple[float, float] = (1.0, 14.0),
        n_points: int = 50,
        **kwargs
    ) -> Dict[str, np.ndarray]:
        """
        Scan pH and compute properties at each point.

        Args:
            bond_or_molecule: Molecular system
            pH_range: (pH_min, pH_max)
            n_points: Number of pH points
            **kwargs: Additional arguments for apply_pH

        Returns:
            Dictionary with arrays:
                pH_values: pH values
                energies: Total energy vs pH (Ha)
                free_energies: Free energy vs pH (Ha)
                net_charges: Net charge vs pH
                protonation_fractions: Fraction protonated for each site
        """
        pH_min, pH_max = pH_range
        pH_values = np.linspace(pH_min, pH_max, n_points)

        energies = []
        free_energies = []
        net_charges = []
        # Track individual site protonation
        site_protonation = {i: [] for i in range(len(self.sites))}

        for pH in pH_values:
            result = self.apply_pH(bond_or_molecule, pH, **kwargs)
            energies.append(result['energy'])
            free_energies.append(result['free_energy'])
            net_charges.append(result['net_charge'])

            # Store individual site fractions
            for idx, site in enumerate(self.sites):
                f_prot = site.get_protonated_fraction(pH)
                site_protonation[idx].append(f_prot)

        return {
            'pH_values': pH_values,
            'energies': np.array(energies),
            'free_energies': np.array(free_energies),
            'net_charges': np.array(net_charges),
            'site_protonation': {i: np.array(vals) for i, vals in site_protonation.items()}
        }

    def _get_base_energy(self, bond_or_molecule) -> float:
        """Get base ground state energy."""
        if hasattr(bond_or_molecule, 'energy'):
            return bond_or_molecule.energy
        elif hasattr(bond_or_molecule, '_cached_energy'):
            return bond_or_molecule._cached_energy
        else:
            # CRITICAL FIX: Compute HF energy if not cached (don't return 0.0!)
            if hasattr(bond_or_molecule, 'hamiltonian'):
                logger.info("Computing HF energy (not cached) for pH calculation")
                rdm1_hf, E_hf = bond_or_molecule.hamiltonian.solve_scf()
                bond_or_molecule._cached_energy = E_hf
                return E_hf
            else:
                raise ValueError("Cannot compute energy - no hamiltonian available and no cached energy")

    def _compute_charging_energy(
        self,
        bond_or_molecule,
        net_charge: float,
        solvent: str
    ) -> float:
        """
        Compute electrostatic energy of charging using Born model.

        ΔG_charge = q² / (2 * (1 - 1/ε) * r)

        Args:
            bond_or_molecule: Molecular system
            net_charge: Net charge in elementary units
            solvent: Solvent name

        Returns:
            Charging energy in Ha
        """
        if abs(net_charge) < 0.01:
            return 0.0

        # Solvent dielectric constants
        epsilon_db = {
            'water': 78.36,
            'acetonitrile': 35.69,
            'dmso': 46.83,
            'methanol': 32.63,
            'ethanol': 24.85,
            'vacuum': 1.0
        }
        epsilon = epsilon_db.get(solvent.lower(), 78.36)

        if epsilon <= 1.0:
            # Vacuum: very high penalty for charging
            return 10.0 * net_charge**2  # Large penalty

        # Estimate molecular radius (Å)
        if hasattr(bond_or_molecule, 'atoms'):
            # Use approximate radius from molecule size
            radius = 3.0  # Å (typical small molecule)
        else:
            radius = 2.0  # Å

        # Born charging energy
        # ΔG_Born = q² / (8πε₀r) × (1 - 1/ε)
        # In atomic units and simplified:
        bohr_to_angstrom = 0.529177
        radius_bohr = radius / bohr_to_angstrom

        # Empirical scaling for Born energy
        # Typical: ~50-100 kcal/mol for single charge in water
        # 100 kcal/mol ~ 0.16 Ha
        born_prefactor = 0.15  # Ha per unit charge² per Å

        Delta_G_charge = born_prefactor * net_charge**2 * (1.0 - 1.0/epsilon) / radius

        # Charging is unfavorable (positive energy)
        Delta_G_charge = abs(Delta_G_charge)

        logger.debug(f"Charging energy: ΔG_charge = {Delta_G_charge:.6f} Ha "
                    f"(q = {net_charge:.2f}, ε = {epsilon:.1f}, r = {radius:.1f} Å)")

        return Delta_G_charge

    def _compute_microstate_populations(self, pH: float) -> Dict[str, float]:
        """
        Compute populations of all protonation microstates.

        For N sites, there are 2^N microstates (each site protonated or not).
        Use Boltzmann weighting based on free energies.

        Args:
            pH: Solution pH

        Returns:
            Dictionary mapping microstate label → population
        """
        n_sites = len(self.sites)
        if n_sites > 5:
            logger.warning(f"Too many sites ({n_sites}) for microstate enumeration")
            return {'averaged': 1.0}

        # Generate all 2^N microstates
        n_states = 2 ** n_sites
        microstates = []
        weights = []

        RT = self.RT_298K_Ha

        for state_idx in range(n_states):
            # Binary representation: 1 = protonated, 0 = deprotonated
            protonation_pattern = [(state_idx >> i) & 1 for i in range(n_sites)]

            # Compute free energy of this microstate
            G_state = 0.0
            label = ""
            for i, (site, is_prot) in enumerate(zip(self.sites, protonation_pattern)):
                if is_prot:
                    # Protonated form
                    label += "H"
                    G_site = 0.0  # Reference
                else:
                    # Deprotonated form
                    label += "D"
                    Delta_pH = pH - site.pKa
                    G_site = RT * np.log(10) * Delta_pH

                G_state += G_site

            microstates.append(label)
            weights.append(np.exp(-G_state / RT))

        # Normalize to get populations
        Z = sum(weights)
        populations = {label: w/Z for label, w in zip(microstates, weights)}

        return populations

    def get_titration_curve(
        self,
        bond_or_molecule,
        pH_range: Tuple[float, float] = (0.0, 14.0),
        n_points: int = 100
    ) -> Dict[str, Any]:
        """
        Generate full titration curve showing protonation vs pH.

        Args:
            bond_or_molecule: Molecular system
            pH_range: pH range to scan
            n_points: Number of pH points

        Returns:
            Dictionary with titration data for plotting
        """
        scan_result = self.scan_pH(bond_or_molecule, pH_range, n_points)

        # Compute total protonated sites at each pH
        n_sites = len(self.sites)
        total_protonated = np.zeros(n_points)

        for idx in range(n_sites):
            total_protonated += scan_result['site_protonation'][idx]

        return {
            'pH': scan_result['pH_values'],
            'total_protonated': total_protonated,
            'fraction_protonated': total_protonated / n_sites if n_sites > 0 else total_protonated,
            'net_charge': scan_result['net_charges'],
            'free_energy': scan_result['free_energies'],
            'site_protonation': scan_result['site_protonation']
        }

    def plot_titration_curve(
        self,
        titration_data: Dict[str, Any],
        save_path: Optional[str] = None
    ):
        """
        Plot titration curve (protonation and charge vs pH).

        Args:
            titration_data: Output from get_titration_curve()
            save_path: Optional path to save figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed - cannot plot")
            return

        pH = titration_data['pH']

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # Plot 1: Fraction protonated
        ax1.plot(pH, titration_data['fraction_protonated'], 'b-', linewidth=2)

        # Plot individual sites if available
        if 'site_protonation' in titration_data:
            for idx, site_data in titration_data['site_protonation'].items():
                site = self.sites[idx]
                ax1.plot(pH, site_data, '--', alpha=0.5,
                        label=f"Site {idx}: pKa={site.pKa:.1f}")
                # Mark pKa
                ax1.axvline(site.pKa, color='gray', linestyle=':', alpha=0.5)

        ax1.set_ylabel('Fraction Protonated', fontsize=12)
        ax1.set_ylim(-0.05, 1.05)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=9)
        ax1.set_title('Titration Curve', fontsize=14, fontweight='bold')

        # Plot 2: Net charge
        ax2.plot(pH, titration_data['net_charge'], 'r-', linewidth=2)
        ax2.axhline(y=0, color='k', linestyle='--', linewidth=1)
        ax2.set_xlabel('pH', fontsize=12)
        ax2.set_ylabel('Net Charge', fontsize=12)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")

        plt.show()

    def plot_pH_scan(
        self,
        scan_result: Dict[str, np.ndarray],
        save_path: Optional[str] = None
    ):
        """
        Plot pH-dependent energies and charge.

        Args:
            scan_result: Output from scan_pH()
            save_path: Optional path to save figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed - cannot plot")
            return

        pH = scan_result['pH_values']
        energies = scan_result['free_energies'] * self.Ha_to_kcal  # kcal/mol
        charges = scan_result['net_charges']

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Plot 1: Free energy vs pH
        ax1.plot(pH, energies, 'b-', linewidth=2, marker='o', markersize=4)
        ax1.set_xlabel('pH', fontsize=12)
        ax1.set_ylabel('Free Energy (kcal/mol)', fontsize=12)
        ax1.set_title('pH-Dependent Free Energy', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)

        # Mark pKa values
        for site in self.sites:
            ax1.axvline(site.pKa, color='gray', linestyle='--', alpha=0.5)
            ax1.text(site.pKa, ax1.get_ylim()[1], f'pKa={site.pKa:.1f}',
                    rotation=90, va='top', fontsize=9)

        # Plot 2: Net charge vs pH
        ax2.plot(pH, charges, 'r-', linewidth=2, marker='s', markersize=4)
        ax2.axhline(y=0, color='k', linestyle='--', linewidth=1)
        ax2.set_xlabel('pH', fontsize=12)
        ax2.set_ylabel('Net Charge', fontsize=12)
        ax2.set_title('pH-Dependent Charge State', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)

        # Mark pKa values
        for site in self.sites:
            ax2.axvline(site.pKa, color='gray', linestyle='--', alpha=0.5)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")

        plt.show()

    def get_available_groups(self) -> List[str]:
        """Return list of available functional groups in pKa database."""
        return list(self.pka_db.keys())

    def get_group_info(self, group_type: str) -> Dict[str, Any]:
        """
        Get detailed information about a functional group.

        Args:
            group_type: Group name

        Returns:
            Dictionary with pKa and group information
        """
        if group_type not in self.pka_db:
            raise ValueError(f"Group '{group_type}' not in database. "
                           f"Available: {self.get_available_groups()}")

        return self.pka_db[group_type].copy()
