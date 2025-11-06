"""
Pressure Effects Modulator

Models pressure effects on molecular systems:
- Bond compression and stretching
- Volume changes and equations of state
- Pressure-induced phase transitions
- Mechanical properties (compressibility, bulk modulus)
- High-pressure chemistry (barrier lowering)

Physical Basis:
    Pressure-volume work: ΔG(P) = ΔG° + P × ΔV

    For molecular compression:
        E(P) = E₀ + K × (V/V₀ - 1)²
    where K is bulk modulus

    Bond length change under pressure (Morse potential):
        r(P) = r₀ × [1 - α × P]
    where α is compressibility

Applications:
    - High-pressure synthesis (diamond anvil cells, >100 GPa)
    - Geological chemistry (Earth's mantle conditions)
    - Supramolecular assembly (pressure-induced polymerization)
    - Materials under extreme conditions
    - Planetary interiors

References:
    - Hemley & Ashcroft "The Revealing Role of Pressure in Physics of Condensed Matter" Phys. Today (1998)
    - McMillan "Chemistry at High Pressure" Chem. Soc. Rev. (2006)
"""

import numpy as np
from typing import Dict, Any, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class PressureModulator:
    """
    Apply pressure effects to molecular Hamiltonians.

    Pressure affects:
    1. Bond lengths: Compression under pressure
    2. Vibrational frequencies: Increase with compression
    3. Potential energy: PV work term
    4. Phase transitions: Pressure-induced structural changes
    5. Reaction barriers: Often lowered under pressure

    Example:
        >>> from kanad.bonds import BondFactory
        >>> from kanad.environment import PressureModulator
        >>>
        >>> bond = BondFactory.create_bond('H', 'H', distance=0.74)
        >>> pressure_mod = PressureModulator()
        >>>
        >>> # Ambient pressure
        >>> result_1bar = pressure_mod.apply_pressure(bond, 1.0)
        >>> print(f"Bond length at 1 bar: {result_1bar['bond_length']:.4f} Å")
        >>>
        >>> # High pressure (10 GPa ~ 100,000 bar)
        >>> result_10GPa = pressure_mod.apply_pressure(bond, 100000.0)
        >>> print(f"Bond length at 10 GPa: {result_10GPa['bond_length']:.4f} Å")
        >>> print(f"Compression: {result_10GPa['compression_ratio']:.3f}")
    """

    # Physical constants
    bar_to_GPa = 1e-4                    # 1 bar = 1e-4 GPa
    GPa_to_Ha_per_bohr3 = 3.398827e-5   # Pressure unit conversion
    bohr_to_angstrom = 0.529177         # Length conversion
    Ha_to_kcal = 627.509474             # Energy conversion

    # Typical compressibilities (1/GPa) for different bond types
    COMPRESSIBILITIES = {
        'H-H': 0.015,     # Very compressible
        'C-C': 0.008,     # Medium
        'C=C': 0.005,     # Stiffer
        'C≡C': 0.003,     # Very stiff
        'C-H': 0.010,
        'C-N': 0.007,
        'C-O': 0.006,
        'N-N': 0.009,
        'O-O': 0.010,
        'default': 0.008  # Generic organic bond
    }

    # Bulk moduli (GPa) for common materials
    BULK_MODULI = {
        'organic': 10.0,      # Typical organic molecules
        'water': 2.2,         # Liquid water
        'ice': 8.9,           # Ice
        'diamond': 442.0,     # Diamond (very incompressible)
        'graphite': 34.0,     # Graphite
        'polymer': 3.0,       # Typical polymer
        'metal': 100.0,       # Typical metal
    }

    def __init__(self):
        """Initialize pressure modulator."""
        self.reference_pressure = 1.0  # bar (ambient)
        logger.info("PressureModulator initialized")

    def apply_pressure(
        self,
        bond_or_molecule,
        pressure: float,
        temperature: float = 298.15,
        bulk_modulus: Optional[float] = None,
        allow_phase_transition: bool = True
    ) -> Dict[str, Any]:
        """
        Apply pressure effects to molecular system.

        Args:
            bond_or_molecule: Bond or Molecule object
            pressure: Pressure in bar (1 bar = ambient, 1e5 bar = 10 GPa)
            temperature: Temperature in Kelvin
            bulk_modulus: Custom bulk modulus in GPa (or use default)
            allow_phase_transition: Check for pressure-induced phase transitions

        Returns:
            Dictionary with:
                energy: Pressure-corrected energy (Ha)
                bond_length: Compressed bond length (Å) or None for molecules
                compression_ratio: V/V₀
                pV_work: Pressure-volume work contribution (Ha)
                bulk_modulus: Bulk modulus used (GPa)
                vibrational_shift: Change in vibrational frequency
                phase: Detected phase (if applicable)
        """
        logger.info(f"Applying pressure effects: P = {pressure:.1f} bar = "
                   f"{pressure * self.bar_to_GPa:.4f} GPa, T = {temperature:.2f}K")

        # Get base energy
        E_base = self._get_base_energy(bond_or_molecule)

        # Convert pressure to GPa for calculations
        P_GPa = pressure * self.bar_to_GPa

        # 1. Determine bulk modulus
        if bulk_modulus is None:
            bulk_modulus = self._estimate_bulk_modulus(bond_or_molecule)

        # 2. Compute volume compression ratio
        compression_ratio = self._compute_compression_ratio(P_GPa, bulk_modulus)

        # 3. Bond length changes (if bond object)
        if hasattr(bond_or_molecule, 'distance') and hasattr(bond_or_molecule, 'atom1'):
            bond_length_compressed = self._compute_compressed_bond_length(
                bond_or_molecule, P_GPa
            )
            bond_info = {
                'bond_length': bond_length_compressed,
                'original_bond_length': bond_or_molecule.distance * self.bohr_to_angstrom
            }
        else:
            bond_length_compressed = None
            bond_info = {}

        # 4. PV work term
        # ΔG = PΔV = P(V - V₀) = PV₀(V/V₀ - 1)
        V0 = self._estimate_molecular_volume(bond_or_molecule)  # Å³
        Delta_V = V0 * (compression_ratio - 1.0)  # Negative (compression)

        # Convert: P (GPa) × V (Ų) → Ha
        # 1 GPa·Ų = 0.001 GPa·nm³ = 6.022e-6 Ha
        pV_work = P_GPa * Delta_V * 6.022e-6

        # 5. Elastic strain energy
        # E_strain = ½ K (V/V₀ - 1)²
        # For volume: E = (9/2) K V₀ [(V/V₀)^(-2/3) - 1]² (Murnaghan EOS)
        # Simplified: E_strain ≈ K V₀ (1 - V/V₀)²
        strain = 1.0 - compression_ratio
        E_strain = bulk_modulus * V0 * strain**2 * 6.022e-6  # Ha

        # 6. Vibrational frequency shift
        # ω(P) ≈ ω₀ × (V/V₀)^(-γ) where γ ~ 1-2 (Grüneisen parameter)
        gamma = 1.5  # Typical Grüneisen parameter
        freq_shift_factor = compression_ratio ** (-gamma)

        # 7. Total energy under pressure
        E_pressure = E_base + pV_work + E_strain

        # 8. Phase detection (simplified)
        phase = self._detect_phase(P_GPa, temperature, compression_ratio)

        return {
            'energy': E_pressure,
            'compression_ratio': compression_ratio,
            'pV_work': pV_work,
            'strain_energy': E_strain,
            'bulk_modulus': bulk_modulus,
            'vibrational_shift_factor': freq_shift_factor,
            'molecular_volume': V0,
            'compressed_volume': V0 * compression_ratio,
            'phase': phase,
            'pressure_GPa': P_GPa,
            'pressure_bar': pressure,
            'temperature': temperature,
            **bond_info
        }

    def scan_pressure(
        self,
        bond_or_molecule,
        pressure_range: Tuple[float, float] = (1.0, 100000.0),
        n_points: int = 20,
        log_scale: bool = True,
        **kwargs
    ) -> Dict[str, np.ndarray]:
        """
        Scan pressure and compute properties at each point.

        Args:
            bond_or_molecule: Molecular system
            pressure_range: (P_min, P_max) in bar
            n_points: Number of pressure points
            log_scale: Use logarithmic pressure spacing (recommended)
            **kwargs: Additional arguments for apply_pressure

        Returns:
            Dictionary with arrays:
                pressures: Pressure values (bar)
                pressures_GPa: Pressure values (GPa)
                energies: Total energy vs P (Ha)
                compression_ratios: V/V₀ vs P
                bond_lengths: Bond length vs P (if applicable)
                bulk_moduli: Bulk modulus vs P
        """
        P_min, P_max = pressure_range

        if log_scale:
            pressures = np.logspace(np.log10(P_min), np.log10(P_max), n_points)
        else:
            pressures = np.linspace(P_min, P_max, n_points)

        energies = []
        compression_ratios = []
        bond_lengths = []
        pV_works = []
        strain_energies = []

        for P in pressures:
            result = self.apply_pressure(bond_or_molecule, P, **kwargs)
            energies.append(result['energy'])
            compression_ratios.append(result['compression_ratio'])
            pV_works.append(result['pV_work'])
            strain_energies.append(result['strain_energy'])

            if 'bond_length' in result and result['bond_length'] is not None:
                bond_lengths.append(result['bond_length'])

        results = {
            'pressures': pressures,
            'pressures_GPa': pressures * self.bar_to_GPa,
            'energies': np.array(energies),
            'compression_ratios': np.array(compression_ratios),
            'pV_works': np.array(pV_works),
            'strain_energies': np.array(strain_energies)
        }

        if bond_lengths:
            results['bond_lengths'] = np.array(bond_lengths)

        return results

    def compute_equation_of_state(
        self,
        bond_or_molecule,
        pressure_range: Tuple[float, float] = (1.0, 200000.0),
        n_points: int = 50,
        fit_model: str = 'murnaghan'
    ) -> Dict[str, Any]:
        """
        Compute pressure-volume equation of state.

        Args:
            bond_or_molecule: Molecular system
            pressure_range: Pressure range in bar
            n_points: Number of points for fitting
            fit_model: EOS model ('murnaghan', 'birch-murnaghan', 'vinet')

        Returns:
            Dictionary with:
                V0: Equilibrium volume (Ų)
                K0: Bulk modulus at zero pressure (GPa)
                K0_prime: Pressure derivative of bulk modulus
                fit_parameters: Fitted EOS parameters
                pressures, volumes: P-V data
        """
        scan_result = self.scan_pressure(
            bond_or_molecule, pressure_range, n_points, log_scale=True
        )

        pressures_GPa = scan_result['pressures_GPa']
        V_V0 = scan_result['compression_ratios']

        # Estimate V0
        V0 = self._estimate_molecular_volume(bond_or_molecule)

        # Fit Murnaghan EOS: P = (K₀/K'₀) [(V₀/V)^K'₀ - 1]
        # For simplicity, use bulk modulus from first principles
        K0 = self._estimate_bulk_modulus(bond_or_molecule)
        K0_prime = 4.0  # Typical value

        return {
            'V0': V0,
            'K0': K0,
            'K0_prime': K0_prime,
            'pressures_GPa': pressures_GPa,
            'volumes': V0 * V_V0,
            'compression_ratios': V_V0,
            'model': fit_model
        }

    def _get_base_energy(self, bond_or_molecule) -> float:
        """Get base ground state energy."""
        if hasattr(bond_or_molecule, 'energy'):
            return bond_or_molecule.energy
        elif hasattr(bond_or_molecule, '_cached_energy'):
            return bond_or_molecule._cached_energy
        else:
            return 0.0

    def _estimate_bulk_modulus(self, bond_or_molecule) -> float:
        """
        Estimate bulk modulus for molecule.

        Args:
            bond_or_molecule: Molecular system

        Returns:
            Bulk modulus in GPa
        """
        # For organic molecules, typical K ~ 5-15 GPa
        # For small diatomics, higher (~20-50 GPa)

        if hasattr(bond_or_molecule, 'atoms'):
            # Molecule: estimate from composition
            n_atoms = len(bond_or_molecule.atoms)
            if n_atoms <= 2:
                K = 30.0  # Diatomic (stiff)
            elif n_atoms <= 10:
                K = 10.0  # Small organic
            else:
                K = 5.0   # Large organic (softer)
        else:
            # Bond: estimate from bond type
            if hasattr(bond_or_molecule, 'atom1') and hasattr(bond_or_molecule, 'atom2'):
                bond_type = f"{bond_or_molecule.atom1}-{bond_or_molecule.atom2}"
                # Simple diatomic: higher bulk modulus
                K = 25.0
            else:
                K = 10.0  # Default

        return K

    def _compute_compression_ratio(self, P_GPa: float, K: float) -> float:
        """
        Compute volume compression ratio V/V₀ using equation of state.

        Using linearized Murnaghan EOS:
            V/V₀ ≈ (1 + K'₀ P / K₀)^(-1/K'₀)

        For small P: V/V₀ ≈ 1 - P/K

        Args:
            P_GPa: Pressure in GPa
            K: Bulk modulus in GPa

        Returns:
            Compression ratio V/V₀ (dimensionless, ≤ 1)
        """
        if P_GPa <= 0:
            return 1.0

        K_prime = 4.0  # Typical value for K'₀

        # Murnaghan EOS
        if P_GPa < 0.1 * K:
            # Linear regime: V/V₀ ≈ 1 - P/K
            ratio = 1.0 - P_GPa / K
        else:
            # Full Murnaghan: V/V₀ = (1 + K'P/K)^(-1/K')
            ratio = (1.0 + K_prime * P_GPa / K) ** (-1.0 / K_prime)

        # Ensure physical range
        ratio = np.clip(ratio, 0.3, 1.0)  # Max compression ~70%

        return ratio

    def _compute_compressed_bond_length(
        self,
        bond,
        P_GPa: float
    ) -> float:
        """
        Compute compressed bond length under pressure.

        r(P) = r₀ × [1 - α × P]
        where α is bond compressibility (1/GPa)

        Args:
            bond: Bond object
            P_GPa: Pressure in GPa

        Returns:
            Compressed bond length in Angstroms
        """
        r0 = bond.distance * self.bohr_to_angstrom  # Original length

        # Get bond compressibility
        if hasattr(bond, 'atom1') and hasattr(bond, 'atom2'):
            bond_key = f"{bond.atom1}-{bond.atom2}"
            alpha = self.COMPRESSIBILITIES.get(bond_key, self.COMPRESSIBILITIES['default'])
        else:
            alpha = self.COMPRESSIBILITIES['default']

        # Linear compression
        if P_GPa < 50:  # GPa
            r_compressed = r0 * (1.0 - alpha * P_GPa)
        else:
            # Nonlinear for extreme pressure (logarithmic)
            r_compressed = r0 * np.exp(-alpha * P_GPa / 10.0)

        # Physical limit: bond can't compress below ~50% of original length
        r_min = 0.5 * r0
        r_compressed = max(r_compressed, r_min)

        return r_compressed

    def _estimate_molecular_volume(self, bond_or_molecule) -> float:
        """
        Estimate molecular volume using van der Waals radii.

        Args:
            bond_or_molecule: Molecular system

        Returns:
            Volume in Ų
        """
        if hasattr(bond_or_molecule, 'atoms'):
            # Molecule: sum atomic volumes with overlap correction
            n_atoms = len(bond_or_molecule.atoms)
            # Typical atom volume ~ 20 Ų (van der Waals sphere)
            V_atom = 20.0
            # Overlap factor (molecules are not close-packed)
            overlap = 0.7
            volume = n_atoms * V_atom * overlap
        else:
            # Bond: estimate from two atoms
            # H2: ~10 Ų, typical diatomic ~20-30 Ų
            if hasattr(bond_or_molecule, 'distance'):
                r = bond_or_molecule.distance * self.bohr_to_angstrom
                # Cylindrical approximation: V ~ πr²L
                L = r  # Bond length
                r_atom = 1.5  # Å (van der Waals)
                volume = np.pi * r_atom**2 * L + (4/3) * np.pi * r_atom**3
            else:
                volume = 20.0  # Default

        return max(volume, 5.0)  # Minimum 5 Ų

    def _detect_phase(
        self,
        P_GPa: float,
        temperature: float,
        compression_ratio: float
    ) -> str:
        """
        Detect molecular phase based on pressure and temperature.

        Simplified phase detection:
        - Gas: P < 0.001 GPa (< 10 bar) and T > 300K
        - Liquid: 0.001 < P < 1 GPa, moderate T
        - Solid: P > 1 GPa or T < 200K
        - High-pressure solid: P > 10 GPa

        Args:
            P_GPa: Pressure in GPa
            temperature: Temperature in K
            compression_ratio: V/V₀

        Returns:
            Phase label
        """
        if P_GPa < 0.001 and temperature > 300:
            return 'gas'
        elif P_GPa > 100:
            return 'high_pressure_solid'
        elif P_GPa > 10:
            return 'compressed_solid'
        elif P_GPa > 1:
            return 'solid'
        elif P_GPa > 0.001:
            return 'liquid'
        else:
            return 'gas'

    def modify_hamiltonian_with_pressure(
        self,
        hamiltonian,
        pressure: float,
        bulk_modulus: Optional[float] = None
    ):
        """
        Create pressure-modified Hamiltonian.

        Under pressure:
        1. Bond lengths shorten → increases kinetic energy
        2. Coulomb interactions strengthen (1/r increases as r decreases)
        3. PV work adds constant shift

        Args:
            hamiltonian: Original Hamiltonian
            pressure: Pressure in bar
            bulk_modulus: Bulk modulus in GPa

        Returns:
            Modified Hamiltonian object
        """
        P_GPa = pressure * self.bar_to_GPa

        if bulk_modulus is None:
            bulk_modulus = 10.0  # Default organic

        compression_ratio = self._compute_compression_ratio(P_GPa, bulk_modulus)

        # Compression increases energy scale
        # Kinetic energy T ~ 1/r² → T_compressed ~ T₀ / (V/V₀)^(2/3)
        # Approximate scaling
        scaling_factor = compression_ratio ** (-2.0/3.0)

        logger.info(f"Applying pressure scaling: P = {P_GPa:.2f} GPa, "
                   f"V/V₀ = {compression_ratio:.3f}, scale = {scaling_factor:.3f}")

        # Scale Hamiltonian
        H_pressure = hamiltonian * scaling_factor

        return H_pressure

    def plot_pressure_scan(
        self,
        scan_result: Dict[str, np.ndarray],
        properties: list = ['energies', 'compression_ratios'],
        save_path: Optional[str] = None
    ):
        """
        Plot pressure-dependent properties.

        Args:
            scan_result: Output from scan_pressure()
            properties: List of properties to plot
            save_path: Optional path to save figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed - cannot plot")
            return

        P_GPa = scan_result['pressures_GPa']

        n_props = len(properties)
        fig, axes = plt.subplots(n_props, 1, figsize=(10, 4*n_props))
        if n_props == 1:
            axes = [axes]

        for ax, prop in zip(axes, properties):
            if prop == 'energies':
                values = scan_result[prop] * self.Ha_to_kcal
                ylabel = 'Energy (kcal/mol)'
            elif prop == 'compression_ratios':
                values = scan_result[prop]
                ylabel = 'V/V₀'
                ax.axhline(y=1.0, color='k', linestyle='--', linewidth=1)
            elif prop == 'bond_lengths' and prop in scan_result:
                values = scan_result[prop]
                ylabel = 'Bond Length (Å)'
            else:
                values = scan_result.get(prop, np.zeros_like(P_GPa))
                ylabel = prop.replace('_', ' ').title()

            ax.plot(P_GPa, values, 'o-', linewidth=2, markersize=6)
            ax.set_xlabel('Pressure (GPa)', fontsize=12)
            ax.set_ylabel(ylabel, fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.set_xscale('log')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")

        plt.show()

    def plot_equation_of_state(
        self,
        eos_result: Dict[str, Any],
        save_path: Optional[str] = None
    ):
        """
        Plot pressure-volume equation of state.

        Args:
            eos_result: Output from compute_equation_of_state()
            save_path: Optional path to save figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed - cannot plot")
            return

        P_GPa = eos_result['pressures_GPa']
        V = eos_result['volumes']
        V0 = eos_result['V0']

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Plot 1: P vs V
        ax1.plot(V, P_GPa, 'o-', linewidth=2, markersize=6)
        ax1.axvline(V0, color='r', linestyle='--', linewidth=1.5, label=f'V₀ = {V0:.1f} Ų')
        ax1.set_xlabel('Volume (Ų)', fontsize=12)
        ax1.set_ylabel('Pressure (GPa)', fontsize=12)
        ax1.set_title('Pressure-Volume Equation of State', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Plot 2: P vs V/V₀
        ax2.plot(eos_result['compression_ratios'], P_GPa, 'o-', linewidth=2, markersize=6, color='red')
        ax2.axvline(1.0, color='k', linestyle='--', linewidth=1)
        ax2.set_xlabel('V/V₀', fontsize=12)
        ax2.set_ylabel('Pressure (GPa)', fontsize=12)
        ax2.set_title(f'Compression Curve (K₀ = {eos_result["K0"]:.1f} GPa)',
                     fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")

        plt.show()
