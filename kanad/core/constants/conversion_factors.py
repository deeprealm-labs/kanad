"""
Unit conversion factors for quantum chemistry calculations.

All conversions are relative to atomic units (Hartree atomic units).
"""

from typing import Final


class ConversionFactors:
    """
    Unit conversion factors for quantum chemistry.

    All factors convert TO atomic units unless specified otherwise.
    """

    # Energy conversions
    EV_TO_HARTREE: Final[float] = 0.036749322175655  # eV → Eₕ
    HARTREE_TO_EV: Final[float] = 27.211386245988  # Eₕ → eV
    KCAL_MOL_TO_HARTREE: Final[float] = 0.0015936011  # kcal/mol → Eₕ
    HARTREE_TO_KCAL_MOL: Final[float] = 627.5094740631  # Eₕ → kcal/mol
    KJ_MOL_TO_HARTREE: Final[float] = 0.00038087988  # kJ/mol → Eₕ
    HARTREE_TO_KJ_MOL: Final[float] = 2625.4996394799  # Eₕ → kJ/mol
    WAVENUMBER_TO_HARTREE: Final[float] = 4.556335e-6  # cm⁻¹ → Eₕ
    HARTREE_TO_WAVENUMBER: Final[float] = 219474.63  # Eₕ → cm⁻¹

    # Length conversions
    ANGSTROM_TO_BOHR: Final[float] = 1.8897261245650618  # Å → a₀
    BOHR_TO_ANGSTROM: Final[float] = 0.529177210903  # a₀ → Å
    NANOMETER_TO_BOHR: Final[float] = 18.897261245650618  # nm → a₀
    BOHR_TO_NANOMETER: Final[float] = 0.0529177210903  # a₀ → nm
    PICOMETER_TO_BOHR: Final[float] = 0.018897261245650618  # pm → a₀
    BOHR_TO_PICOMETER: Final[float] = 52.9177210903  # a₀ → pm

    # Time conversions
    FEMTOSECOND_TO_AU: Final[float] = 41.341373336561364  # fs → ℏ/Eₕ
    AU_TO_FEMTOSECOND: Final[float] = 0.024188843265857  # ℏ/Eₕ → fs

    # Mass conversions
    AMU_TO_ME: Final[float] = 1822.888486209  # amu → mₑ (electron mass)
    ME_TO_AMU: Final[float] = 0.0005485799090  # mₑ → amu

    # Electric field conversions
    V_PER_ANGSTROM_TO_AU: Final[float] = 0.0194469064916  # V/Å → Eₕ/(e·a₀)

    # Dipole moment conversions
    DEBYE_TO_AU: Final[float] = 0.393430307  # D → e·a₀
    AU_TO_DEBYE: Final[float] = 2.541746473  # e·a₀ → D

    # Temperature conversions
    KELVIN_TO_HARTREE: Final[float] = 3.166811563e-6  # K → Eₕ (via kB)
    HARTREE_TO_KELVIN: Final[float] = 315775.02480407  # Eₕ → K

    @staticmethod
    def energy_to_hartree(value: float, unit: str) -> float:
        """
        Convert energy to Hartree atomic units.

        Args:
            value: Energy value
            unit: Unit name ('ev', 'kcal/mol', 'kj/mol', 'cm-1', 'K')

        Returns:
            Energy in Hartree

        Raises:
            ValueError: If unit not recognized
        """
        unit = unit.lower().replace('/', '_').replace('-', '_')

        conversions = {
            'ev': ConversionFactors.EV_TO_HARTREE,
            'kcal_mol': ConversionFactors.KCAL_MOL_TO_HARTREE,
            'kj_mol': ConversionFactors.KJ_MOL_TO_HARTREE,
            'cm_1': ConversionFactors.WAVENUMBER_TO_HARTREE,
            'wavenumber': ConversionFactors.WAVENUMBER_TO_HARTREE,
            'k': ConversionFactors.KELVIN_TO_HARTREE,
            'kelvin': ConversionFactors.KELVIN_TO_HARTREE,
            'hartree': 1.0,
            'au': 1.0,
        }

        if unit not in conversions:
            raise ValueError(
                f"Unknown energy unit: {unit}. "
                f"Supported: {', '.join(conversions.keys())}"
            )

        return value * conversions[unit]

    @staticmethod
    def hartree_to_energy(value: float, unit: str) -> float:
        """
        Convert energy from Hartree to specified unit.

        Args:
            value: Energy in Hartree
            unit: Target unit ('ev', 'kcal/mol', 'kj/mol', 'cm-1', 'K')

        Returns:
            Energy in specified unit
        """
        unit = unit.lower().replace('/', '_').replace('-', '_')

        conversions = {
            'ev': ConversionFactors.HARTREE_TO_EV,
            'kcal_mol': ConversionFactors.HARTREE_TO_KCAL_MOL,
            'kj_mol': ConversionFactors.HARTREE_TO_KJ_MOL,
            'cm_1': ConversionFactors.HARTREE_TO_WAVENUMBER,
            'wavenumber': ConversionFactors.HARTREE_TO_WAVENUMBER,
            'k': ConversionFactors.HARTREE_TO_KELVIN,
            'kelvin': ConversionFactors.HARTREE_TO_KELVIN,
            'hartree': 1.0,
            'au': 1.0,
        }

        if unit not in conversions:
            raise ValueError(
                f"Unknown energy unit: {unit}. "
                f"Supported: {', '.join(conversions.keys())}"
            )

        return value * conversions[unit]

    @staticmethod
    def length_to_bohr(value: float, unit: str) -> float:
        """Convert length to Bohr radii."""
        unit = unit.lower()

        conversions = {
            'angstrom': ConversionFactors.ANGSTROM_TO_BOHR,
            'a': ConversionFactors.ANGSTROM_TO_BOHR,
            'nm': ConversionFactors.NANOMETER_TO_BOHR,
            'nanometer': ConversionFactors.NANOMETER_TO_BOHR,
            'pm': ConversionFactors.PICOMETER_TO_BOHR,
            'picometer': ConversionFactors.PICOMETER_TO_BOHR,
            'bohr': 1.0,
            'au': 1.0,
        }

        if unit not in conversions:
            raise ValueError(f"Unknown length unit: {unit}")

        return value * conversions[unit]

    @staticmethod
    def bohr_to_length(value: float, unit: str) -> float:
        """Convert Bohr radii to specified length unit."""
        unit = unit.lower()

        conversions = {
            'angstrom': ConversionFactors.BOHR_TO_ANGSTROM,
            'a': ConversionFactors.BOHR_TO_ANGSTROM,
            'nm': ConversionFactors.BOHR_TO_NANOMETER,
            'nanometer': ConversionFactors.BOHR_TO_NANOMETER,
            'pm': ConversionFactors.BOHR_TO_PICOMETER,
            'picometer': ConversionFactors.BOHR_TO_PICOMETER,
            'bohr': 1.0,
            'au': 1.0,
        }

        if unit not in conversions:
            raise ValueError(f"Unknown length unit: {unit}")

        return value * conversions[unit]
