"""
Fundamental physical constants in SI units.

All values are from CODATA 2018 recommended values.
"""

from typing import Final
from dataclasses import dataclass


@dataclass(frozen=True)
class PhysicalConstants:
    """
    Fundamental physical constants in SI units.

    All values are exact or from CODATA 2018 recommendations.
    """

    # Universal constants
    SPEED_OF_LIGHT: Final[float] = 299792458.0  # m/s (exact)
    PLANCK_CONSTANT: Final[float] = 6.62607015e-34  # J⋅s (exact)
    HBAR: Final[float] = 1.054571817e-34  # J⋅s (ℏ = h/2π)

    # Atomic constants
    ELEMENTARY_CHARGE: Final[float] = 1.602176634e-19  # C (exact)
    ELECTRON_MASS: Final[float] = 9.1093837015e-31  # kg
    PROTON_MASS: Final[float] = 1.67262192369e-27  # kg
    NEUTRON_MASS: Final[float] = 1.67492749804e-27  # kg
    BOHR_RADIUS: Final[float] = 5.29177210903e-11  # m (a₀)

    # Derived atomic units (Hartree atomic units)
    HARTREE_ENERGY: Final[float] = 4.3597447222071e-18  # J (Eₕ)
    HARTREE_TO_EV: Final[float] = 27.211386245988  # eV
    BOHR_TO_ANGSTROM: Final[float] = 0.529177210903  # Å

    # Energy conversions
    EV_TO_JOULE: Final[float] = 1.602176634e-19  # J
    KCAL_MOL_TO_HARTREE: Final[float] = 0.0015936011  # Eₕ
    KJ_MOL_TO_HARTREE: Final[float] = 0.00038087988  # Eₕ

    # Other useful constants
    AVOGADRO_CONSTANT: Final[float] = 6.02214076e23  # mol⁻¹ (exact)
    BOLTZMANN_CONSTANT: Final[float] = 1.380649e-23  # J/K (exact)
    GAS_CONSTANT: Final[float] = 8.314462618  # J/(mol⋅K)

    # Vacuum permittivity
    VACUUM_PERMITTIVITY: Final[float] = 8.8541878128e-12  # F/m (ε₀)

    # Fine structure constant
    FINE_STRUCTURE_CONSTANT: Final[float] = 7.2973525693e-3  # α ≈ 1/137

    def __post_init__(self):
        """Validate constants on initialization."""
        # Verify derived relationships
        assert abs(self.HBAR - self.PLANCK_CONSTANT / (2 * 3.141592653589793)) < 1e-40
        assert abs(self.HARTREE_TO_EV - self.HARTREE_ENERGY / self.EV_TO_JOULE) < 1e-6


# Singleton instance for global access
CONSTANTS = PhysicalConstants()


# Convenience accessors for common values
def hartree_to_ev(energy_hartree: float) -> float:
    """Convert energy from Hartree to eV."""
    return energy_hartree * CONSTANTS.HARTREE_TO_EV


def ev_to_hartree(energy_ev: float) -> float:
    """Convert energy from eV to Hartree."""
    return energy_ev / CONSTANTS.HARTREE_TO_EV


def bohr_to_angstrom(distance_bohr: float) -> float:
    """Convert distance from Bohr radii to Angstroms."""
    return distance_bohr * CONSTANTS.BOHR_TO_ANGSTROM


def angstrom_to_bohr(distance_angstrom: float) -> float:
    """Convert distance from Angstroms to Bohr radii."""
    return distance_angstrom / CONSTANTS.BOHR_TO_ANGSTROM


def kcal_mol_to_hartree(energy_kcal_mol: float) -> float:
    """Convert energy from kcal/mol to Hartree."""
    return energy_kcal_mol * CONSTANTS.KCAL_MOL_TO_HARTREE


def hartree_to_kcal_mol(energy_hartree: float) -> float:
    """Convert energy from Hartree to kcal/mol."""
    return energy_hartree / CONSTANTS.KCAL_MOL_TO_HARTREE
