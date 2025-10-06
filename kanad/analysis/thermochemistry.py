"""
Thermochemistry Calculator

Computes thermodynamic properties (H, S, G) at finite temperature
using statistical mechanics and the rigid rotor-harmonic oscillator
(RRHO) approximation.

Theory:
    H(T) = E_elec + ZPE + E_thermal + RT
    S(T) = S_trans + S_rot + S_vib
    G(T) = H(T) - T·S(T)
"""

import numpy as np
from typing import Dict, List, Optional, Any, Tuple
import logging

logger = logging.getLogger(__name__)


class ThermochemistryCalculator:
    """
    Calculate thermodynamic properties at finite temperature.

    Uses rigid rotor-harmonic oscillator (RRHO) approximation with:
    - Electronic energy from QM calculation
    - Translational partition function (ideal gas)
    - Rotational partition function (rigid rotor)
    - Vibrational partition function (harmonic oscillator)

    Standard state: 298.15 K, 1 atm (101325 Pa)

    Example:
        >>> from kanad.io import from_smiles
        >>> from kanad.analysis import ThermochemistryCalculator
        >>>
        >>> water = from_smiles("O")
        >>> thermo = ThermochemistryCalculator(
        ...     water,
        ...     frequencies=[3657.0, 1595.0, 3756.0]  # cm⁻¹
        ... )
        >>> result = thermo.compute_thermochemistry(T=298.15)
        >>> print(f"Entropy: {result['s']:.2f} cal/(mol·K)")
        >>> print(f"Gibbs free energy: {result['g']:.6f} Ha")
    """

    # Physical constants (CODATA 2018)
    k_B = 1.380649e-23      # Boltzmann constant (J/K)
    h = 6.62607015e-34      # Planck constant (J·s)
    c = 2.99792458e10       # Speed of light (cm/s)
    N_A = 6.02214076e23     # Avogadro's number (mol⁻¹)
    R = 8.314462618         # Gas constant (J/(mol·K))
    amu = 1.66053906660e-27 # Atomic mass unit (kg)

    # Conversion factors
    Ha_to_J = 4.3597447222071e-18       # Hartree to Joules
    Ha_to_kcal = 627.509474             # Hartree to kcal/mol
    cal_to_J = 4.184                    # cal to J
    R_cal = R / cal_to_J                # Gas constant (cal/(mol·K))

    # Known vibrational frequencies (cm⁻¹) for common molecules
    KNOWN_FREQUENCIES = {
        'H2': [4401.2],
        'H2O': [3657.0, 1595.0, 3756.0],
        'NH3': [3337.0, 950.0, 3414.0, 3414.0, 1627.0, 1627.0],
        'CH4': [2917.0, 1534.0, 1534.0, 1534.0, 3019.0, 3019.0, 3019.0, 1306.0, 1306.0],
        'CO': [2143.0],
        'N2': [2330.0],
        'O2': [1580.0],
    }

    def __init__(
        self,
        molecule: 'Molecule',
        frequencies: Optional[List[float]] = None
    ):
        """
        Initialize thermochemistry calculator.

        Args:
            molecule: Molecule object with optimized geometry
            frequencies: Vibrational frequencies (cm⁻¹), optional
                        If None, will look up in database or warn

        Raises:
            ValueError: If molecule has no atoms
        """
        self.molecule = molecule

        if len(molecule.atoms) == 0:
            raise ValueError("Molecule has no atoms")

        # Set frequencies
        if frequencies is not None:
            self.frequencies = np.array(frequencies)
        else:
            # Try to look up
            formula = molecule.formula
            if formula in self.KNOWN_FREQUENCIES:
                self.frequencies = np.array(self.KNOWN_FREQUENCIES[formula])
                logger.info(f"Using known frequencies for {formula}")
            else:
                self.frequencies = None
                logger.warning(
                    f"No frequency data for {formula}. "
                    f"Vibrational contributions will be zero. "
                    f"Provide frequencies explicitly or compute Hessian."
                )

        # Compute molecular properties
        self.mass = self._compute_mass()
        self.com = self._compute_center_of_mass()
        self.I_principal = self._compute_principal_moments()
        self.is_linear = self._check_linearity()
        self.symmetry_number = self._estimate_symmetry_number()
        self.θ_rot = self._compute_rotational_temperatures()

        logger.info(f"ThermochemistryCalculator initialized for {molecule.formula}")
        logger.info(f"  Mass: {self.mass:.4f} amu")
        logger.info(f"  Linear: {self.is_linear}")
        logger.info(f"  Symmetry number: {self.symmetry_number}")
        if self.frequencies is not None:
            logger.info(f"  Vibrational modes: {len(self.frequencies)}")

    def _compute_mass(self) -> float:
        """Compute total molecular mass in amu."""
        return sum(atom.atomic_mass for atom in self.molecule.atoms)

    def _compute_center_of_mass(self) -> np.ndarray:
        """Compute center of mass."""
        total_mass = 0.0
        com = np.zeros(3)

        for atom in self.molecule.atoms:
            total_mass += atom.atomic_mass
            com += atom.atomic_mass * atom.position

        return com / total_mass

    def _compute_inertia_tensor(self) -> np.ndarray:
        """
        Compute inertia tensor about center of mass.

        Returns:
            3x3 inertia tensor (amu·Å²)
        """
        I = np.zeros((3, 3))

        for atom in self.molecule.atoms:
            r = atom.position - self.com  # Position relative to COM
            m = atom.atomic_mass

            # I_xx, I_yy, I_zz
            I[0, 0] += m * (r[1]**2 + r[2]**2)
            I[1, 1] += m * (r[0]**2 + r[2]**2)
            I[2, 2] += m * (r[0]**2 + r[1]**2)

            # I_xy, I_xz, I_yz (symmetric)
            I[0, 1] -= m * r[0] * r[1]
            I[0, 2] -= m * r[0] * r[2]
            I[1, 2] -= m * r[1] * r[2]

        # Symmetrize
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]

        return I

    def _compute_principal_moments(self) -> np.ndarray:
        """
        Compute principal moments of inertia.

        Returns:
            Principal moments in increasing order (amu·Å²)
        """
        n_atoms = len(self.molecule.atoms)

        if n_atoms == 1:
            # Atom - no rotation
            return np.array([0.0, 0.0, 0.0])

        I_tensor = self._compute_inertia_tensor()
        I_principal = np.linalg.eigvalsh(I_tensor)  # Sorted
        return np.maximum(I_principal, 1e-10)  # Avoid division by zero

    def _check_linearity(self) -> bool:
        """
        Check if molecule is linear.

        Linear if smallest moment of inertia is negligible.
        """
        n_atoms = len(self.molecule.atoms)

        if n_atoms <= 2:
            return True

        # Check if smallest moment << others
        I_min = self.I_principal[0]
        I_max = self.I_principal[2]

        return I_min < 1e-3 * I_max

    def _estimate_symmetry_number(self) -> int:
        """
        Estimate rotational symmetry number.

        σ = number of indistinguishable orientations
        """
        formula = self.molecule.formula
        n_atoms = len(self.molecule.atoms)

        if n_atoms == 1:
            return 1  # Atom

        # Known molecules
        SYMMETRY = {
            'H2': 2, 'N2': 2, 'O2': 2, 'F2': 2, 'Cl2': 2,
            'H2O': 2, 'NH3': 3, 'CH4': 12, 'C2H6': 6,
            'CO': 1, 'HCl': 1, 'HF': 1,
        }

        if formula in SYMMETRY:
            return SYMMETRY[formula]

        # Default: assume no symmetry
        logger.warning(f"Unknown symmetry for {formula}, using σ=1")
        return 1

    def _compute_rotational_temperatures(self) -> np.ndarray:
        """
        Compute rotational temperatures θ_rot = h²/(8π²Ik_B).

        Returns:
            Rotational temperatures (K) for each principal axis
        """
        if len(self.molecule.atoms) == 1:
            return np.array([0.0, 0.0, 0.0])

        # Convert amu·Å² to kg·m²
        I_SI = self.I_principal * self.amu * 1e-20  # kg·m²

        # θ_rot = h² / (8π² I k_B)
        θ_rot = self.h**2 / (8 * np.pi**2 * I_SI * self.k_B)

        return θ_rot

    def _translational_thermochemistry(
        self,
        T: float,
        P: float
    ) -> Tuple[float, float]:
        """
        Compute translational contributions.

        Args:
            T: Temperature (K)
            P: Pressure (Pa)

        Returns:
            E_trans: Translational energy (Ha)
            S_trans: Translational entropy (cal/(mol·K))
        """
        # Mass in kg
        m_kg = self.mass * self.amu

        # Volume (ideal gas)
        V = self.k_B * T / P

        # Thermal de Broglie wavelength
        Λ = self.h / np.sqrt(2 * np.pi * m_kg * self.k_B * T)

        # Partition function
        q_trans = V / Λ**3

        # Translational energy: (3/2) RT
        E_trans = (3.0 / 2.0) * self.R * T / self.N_A  # J/molecule
        E_trans /= self.Ha_to_J  # Convert to Ha

        # Translational entropy: R[ln(q/N) + 5/2]
        S_trans = self.R_cal * (np.log(q_trans) + 5.0 / 2.0)

        return E_trans, S_trans

    def _rotational_thermochemistry(self, T: float) -> Tuple[float, float]:
        """
        Compute rotational contributions.

        Args:
            T: Temperature (K)

        Returns:
            E_rot: Rotational energy (Ha)
            S_rot: Rotational entropy (cal/(mol·K))
        """
        n_atoms = len(self.molecule.atoms)

        if n_atoms == 1:
            # Atom - no rotation
            return 0.0, 0.0

        if self.is_linear:
            # Linear molecule: q_rot = T / (σ θ_rot)
            θ_rot = self.θ_rot[1]  # Use middle moment (perpendicular to axis)
            q_rot = T / (self.symmetry_number * θ_rot)

            E_rot = self.R * T / self.N_A  # RT/molecule
            E_rot /= self.Ha_to_J

            S_rot = self.R_cal * (np.log(q_rot) + 1.0)

        else:
            # Nonlinear molecule: q_rot = sqrt(π)/σ × (T³/(θ_A θ_B θ_C))^(1/2)
            θ_A, θ_B, θ_C = self.θ_rot
            q_rot = (np.sqrt(np.pi) / self.symmetry_number) * \
                    (T**3 / (θ_A * θ_B * θ_C))**0.5

            E_rot = (3.0 / 2.0) * self.R * T / self.N_A
            E_rot /= self.Ha_to_J

            S_rot = self.R_cal * (np.log(q_rot) + 3.0 / 2.0)

        return E_rot, S_rot

    def _vibrational_thermochemistry(self, T: float) -> Tuple[float, float, float]:
        """
        Compute vibrational contributions.

        Args:
            T: Temperature (K)

        Returns:
            ZPE: Zero-point energy (Ha)
            E_vib: Vibrational thermal energy (Ha)
            S_vib: Vibrational entropy (cal/(mol·K))
        """
        if self.frequencies is None or len(self.frequencies) == 0:
            return 0.0, 0.0, 0.0

        ZPE = 0.0
        E_vib = 0.0
        S_vib = 0.0

        for ν_cm in self.frequencies:
            # Convert cm⁻¹ to J
            ν_J = ν_cm * self.c * self.h

            # Vibrational temperature: θ_vib = hν/k_B
            θ_vib = ν_J / self.k_B

            # Zero-point energy: hν/2 per mode
            ZPE += 0.5 * ν_J / self.Ha_to_J

            # Thermal vibrational energy
            x = θ_vib / T
            if x < 100:  # Avoid overflow
                E_vib_mode = self.R * θ_vib / (np.exp(x) - 1)
                E_vib += E_vib_mode / self.N_A / self.Ha_to_J

                # Vibrational entropy
                S_vib_mode = self.R_cal * (x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)))
                S_vib += S_vib_mode

        return ZPE, E_vib, S_vib

    def compute_thermochemistry(
        self,
        temperature: float = 298.15,
        pressure: float = 101325.0,
        method: str = 'HF'
    ) -> Dict[str, Any]:
        """
        Compute thermodynamic properties at given T and P.

        Args:
            temperature: Temperature (K), default 298.15
            pressure: Pressure (Pa), default 101325 (1 atm)
            method: Electronic structure method ('HF', 'MP2')

        Returns:
            Dictionary with:
                temperature: T (K)
                pressure: P (Pa)
                e_elec: Electronic energy (Ha)
                zpe: Zero-point energy (Ha)
                e_trans: Translational energy (Ha)
                e_rot: Rotational energy (Ha)
                e_vib: Vibrational thermal energy (Ha)
                e_thermal: Total thermal correction (Ha)
                h: Enthalpy H = E + ZPE + E_thermal + RT (Ha)
                s_trans: Translational entropy (cal/(mol·K))
                s_rot: Rotational entropy (cal/(mol·K))
                s_vib: Vibrational entropy (cal/(mol·K))
                s: Total entropy (cal/(mol·K))
                g: Gibbs free energy G = H - TS (Ha)
                cv: Heat capacity at const V (cal/(mol·K))
                cp: Heat capacity at const P (cal/(mol·K))
        """
        T = temperature
        P = pressure

        # Electronic energy
        if method.upper() == 'HF':
            E_elec = self.molecule.hamiltonian.hf_energy
        elif method.upper() == 'MP2':
            from kanad.core.correlation import MP2Solver
            mp2_solver = MP2Solver(self.molecule.hamiltonian)
            result_mp2 = mp2_solver.compute_energy()
            E_elec = result_mp2['e_mp2']
        else:
            raise ValueError(f"Unknown method: {method}")

        # Thermal contributions
        E_trans, S_trans = self._translational_thermochemistry(T, P)
        E_rot, S_rot = self._rotational_thermochemistry(T)
        ZPE, E_vib, S_vib = self._vibrational_thermochemistry(T)

        # Total thermal energy
        E_thermal = E_trans + E_rot + E_vib

        # Enthalpy: H = E_elec + ZPE + E_thermal + RT
        RT = self.R * T / self.N_A / self.Ha_to_J  # Ha
        H = E_elec + ZPE + E_thermal + RT

        # Total entropy
        S_total = S_trans + S_rot + S_vib  # cal/(mol·K)

        # Gibbs free energy: G = H - TS
        TS = T * S_total * self.cal_to_J / self.N_A / self.Ha_to_J  # Ha
        G = H - TS

        # Heat capacities (ideal gas)
        n_atoms = len(self.molecule.atoms)
        if n_atoms == 1:
            # Atom: Cv = (3/2)R
            Cv = 1.5 * self.R_cal
            Cp = 2.5 * self.R_cal
        elif self.is_linear:
            # Linear: Cv = (5/2)R + Cv_vib
            Cv = 2.5 * self.R_cal + self._vibrational_heat_capacity(T)
            Cp = Cv + self.R_cal
        else:
            # Nonlinear: Cv = 3R + Cv_vib
            Cv = 3.0 * self.R_cal + self._vibrational_heat_capacity(T)
            Cp = Cv + self.R_cal

        return {
            'temperature': T,
            'pressure': P,
            'method': method,
            'e_elec': E_elec,
            'zpe': ZPE,
            'e_trans': E_trans,
            'e_rot': E_rot,
            'e_vib': E_vib,
            'e_thermal': E_thermal,
            'h': H,
            's_trans': S_trans,
            's_rot': S_rot,
            's_vib': S_vib,
            's': S_total,
            'g': G,
            'cv': Cv,
            'cp': Cp,
        }

    def _vibrational_heat_capacity(self, T: float) -> float:
        """
        Compute vibrational contribution to heat capacity.

        Cv_vib = R Σ_i (θ_i/T)² exp(θ_i/T) / (exp(θ_i/T) - 1)²

        Args:
            T: Temperature (K)

        Returns:
            Cv_vib (cal/(mol·K))
        """
        if self.frequencies is None or len(self.frequencies) == 0:
            return 0.0

        Cv_vib = 0.0

        for ν_cm in self.frequencies:
            # Vibrational temperature
            θ_vib = self.h * self.c * ν_cm / self.k_B

            x = θ_vib / T
            if x < 100:  # Avoid overflow
                exp_x = np.exp(x)
                Cv_vib += self.R_cal * x**2 * exp_x / (exp_x - 1)**2

        return Cv_vib
