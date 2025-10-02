"""
Unit tests for physical constants and atomic data.

Tests ensure scientific accuracy of fundamental constants and
atomic properties critical for quantum chemistry calculations.
"""

import pytest
import numpy as np
from kanad.core.constants.physical_constants import (
    CONSTANTS,
    PhysicalConstants,
    hartree_to_ev,
    ev_to_hartree,
    bohr_to_angstrom,
    angstrom_to_bohr,
)
from kanad.core.constants.atomic_data import PeriodicTable, AtomicProperties
from kanad.core.constants.conversion_factors import ConversionFactors


class TestPhysicalConstants:
    """Test fundamental physical constants."""

    def test_constants_immutable(self):
        """Constants should be immutable (frozen dataclass)."""
        with pytest.raises(Exception):  # FrozenInstanceError
            CONSTANTS.SPEED_OF_LIGHT = 1.0

    def test_planck_constant_value(self):
        """Verify Planck constant matches CODATA 2018."""
        assert CONSTANTS.PLANCK_CONSTANT == 6.62607015e-34

    def test_hbar_derived(self):
        """ℏ should equal h/2π."""
        expected_hbar = CONSTANTS.PLANCK_CONSTANT / (2 * np.pi)
        assert abs(CONSTANTS.HBAR - expected_hbar) < 1e-40

    def test_elementary_charge(self):
        """Verify elementary charge is exact (CODATA 2018)."""
        assert CONSTANTS.ELEMENTARY_CHARGE == 1.602176634e-19

    def test_hartree_to_ev_conversion(self):
        """Test Hartree to eV conversion factor."""
        # 1 Hartree should equal 27.211... eV
        assert abs(CONSTANTS.HARTREE_TO_EV - 27.211386245988) < 1e-10

    def test_bohr_to_angstrom_conversion(self):
        """Test Bohr to Angstrom conversion."""
        # 1 Bohr radius = 0.529177... Angstrom
        assert abs(CONSTANTS.BOHR_TO_ANGSTROM - 0.529177210903) < 1e-10

    def test_energy_conversion_functions(self):
        """Test energy conversion helper functions."""
        # Test round-trip conversion
        energy_hartree = 1.0
        energy_ev = hartree_to_ev(energy_hartree)
        back_to_hartree = ev_to_hartree(energy_ev)

        assert abs(energy_hartree - back_to_hartree) < 1e-12
        assert abs(energy_ev - 27.211386245988) < 1e-10

    def test_length_conversion_functions(self):
        """Test length conversion helper functions."""
        # Test round-trip conversion
        distance_bohr = 1.0
        distance_angstrom = bohr_to_angstrom(distance_bohr)
        back_to_bohr = angstrom_to_bohr(distance_angstrom)

        assert abs(distance_bohr - back_to_bohr) < 1e-12
        assert abs(distance_angstrom - 0.529177210903) < 1e-10


class TestAtomicData:
    """Test periodic table and atomic properties."""

    def test_hydrogen_properties(self):
        """Verify hydrogen atomic properties."""
        h = PeriodicTable.get_element('H')

        assert h.symbol == 'H'
        assert h.atomic_number == 1
        assert h.valence_electrons == 1
        assert abs(h.electronegativity - 2.20) < 0.01
        assert abs(h.ionization_energy - 13.598) < 0.01
        assert h.is_metal is False

    def test_carbon_properties(self):
        """Verify carbon atomic properties for organic chemistry."""
        c = PeriodicTable.get_element('C')

        assert c.symbol == 'C'
        assert c.atomic_number == 6
        assert c.valence_electrons == 4
        assert abs(c.electronegativity - 2.55) < 0.01
        assert c.is_metal is False
        assert c.group == 14
        assert c.period == 2

    def test_sodium_chlorine_ionic(self):
        """Test Na-Cl ionic bond criteria."""
        na = PeriodicTable.get_element('Na')
        cl = PeriodicTable.get_element('Cl')

        # Na is metal, Cl is not
        assert na.is_metal is True
        assert cl.is_metal is False

        # Large electronegativity difference (> 1.7 → ionic)
        delta_en = PeriodicTable.electronegativity_difference('Na', 'Cl')
        assert delta_en > 1.7  # Ionic bond threshold
        assert abs(delta_en - (3.16 - 0.93)) < 0.01

    def test_covalent_radius_h2(self):
        """Test H-H bond length estimation."""
        estimated_length = PeriodicTable.estimate_bond_length('H', 'H')

        # H2 bond length ≈ 0.74 Å
        # 2 × covalent radius should be close
        assert abs(estimated_length - 0.62) < 0.05  # 2 × 0.31 Å

    def test_element_not_found(self):
        """Should raise KeyError for unknown elements."""
        with pytest.raises(KeyError):
            PeriodicTable.get_element('Xx')

    def test_get_by_atomic_number(self):
        """Test lookup by atomic number."""
        oxygen = PeriodicTable.get_by_atomic_number(8)

        assert oxygen is not None
        assert oxygen.symbol == 'O'
        assert oxygen.atomic_number == 8

    def test_list_elements(self):
        """Test element list retrieval."""
        elements = PeriodicTable.list_elements()

        assert 'H' in elements
        assert 'C' in elements
        assert 'O' in elements
        assert len(elements) > 0


class TestConversionFactors:
    """Test unit conversion utilities."""

    def test_energy_conversions_hartree_ev(self):
        """Test Hartree ↔ eV conversions."""
        energy_hartree = 1.0

        # Convert to eV
        energy_ev = ConversionFactors.hartree_to_energy(energy_hartree, 'eV')
        assert abs(energy_ev - 27.211386245988) < 1e-10

        # Convert back to Hartree
        back = ConversionFactors.energy_to_hartree(energy_ev, 'eV')
        assert abs(back - energy_hartree) < 1e-12

    def test_energy_conversions_kcal_mol(self):
        """Test Hartree ↔ kcal/mol conversions."""
        energy_hartree = 1.0

        energy_kcal = ConversionFactors.hartree_to_energy(energy_hartree, 'kcal/mol')
        assert abs(energy_kcal - 627.5094740631) < 1e-6

        back = ConversionFactors.energy_to_hartree(energy_kcal, 'kcal/mol')
        assert abs(back - energy_hartree) < 1e-6  # Relaxed tolerance for round-trip

    def test_length_conversions_angstrom(self):
        """Test Bohr ↔ Angstrom conversions."""
        length_bohr = 1.0

        length_angstrom = ConversionFactors.bohr_to_length(length_bohr, 'angstrom')
        assert abs(length_angstrom - 0.529177210903) < 1e-10

        back = ConversionFactors.length_to_bohr(length_angstrom, 'angstrom')
        assert abs(back - length_bohr) < 1e-9  # Relaxed tolerance for round-trip

    def test_invalid_unit_raises_error(self):
        """Test that invalid units raise ValueError."""
        with pytest.raises(ValueError):
            ConversionFactors.energy_to_hartree(1.0, 'invalid_unit')

        with pytest.raises(ValueError):
            ConversionFactors.length_to_bohr(1.0, 'invalid_unit')

    def test_temperature_to_hartree(self):
        """Test temperature to energy conversion."""
        # Room temperature ~300 K
        temp_k = 300.0
        energy_hartree = ConversionFactors.energy_to_hartree(temp_k, 'K')

        # Should be very small energy (~0.001 Hartree)
        assert 0.0005 < energy_hartree < 0.002

    def test_wavenumber_conversions(self):
        """Test cm⁻¹ ↔ Hartree conversions."""
        wavenumber = 1000.0  # cm⁻¹

        energy_hartree = ConversionFactors.energy_to_hartree(wavenumber, 'cm-1')
        assert energy_hartree > 0

        back = ConversionFactors.hartree_to_energy(energy_hartree, 'cm-1')
        assert abs(back - wavenumber) < 0.001  # Relaxed tolerance for wavenumbers


class TestScientificAccuracy:
    """
    Integration tests for scientific accuracy.

    These tests verify that constants and conversions maintain
    the precision required for quantum chemistry calculations.
    """

    def test_hartree_energy_precision(self):
        """
        Verify Hartree energy is precise enough for chemistry.

        Chemical accuracy requires ~1.6 mHa (milli-Hartree) precision.
        """
        # Hartree energy should be known to at least 10 significant figures
        hartree_joules = CONSTANTS.HARTREE_ENERGY
        assert hartree_joules == 4.3597447222071e-18

        # Verify derived value matches
        derived = CONSTANTS.HARTREE_TO_EV * CONSTANTS.EV_TO_JOULE
        relative_error = abs(hartree_joules - derived) / hartree_joules
        assert relative_error < 1e-10

    def test_bohr_radius_precision(self):
        """Verify Bohr radius precision for accurate geometries."""
        # Bohr radius should be precise to picometer scale
        bohr_meters = CONSTANTS.BOHR_RADIUS
        assert bohr_meters == 5.29177210903e-11

        # Convert to Angstroms and verify
        bohr_angstrom = bohr_meters * 1e10
        assert abs(bohr_angstrom - CONSTANTS.BOHR_TO_ANGSTROM) < 1e-12

    def test_fine_structure_constant(self):
        """
        Verify fine structure constant.

        α ≈ 1/137.036... is crucial for relativistic corrections.
        """
        alpha = CONSTANTS.FINE_STRUCTURE_CONSTANT
        assert abs(alpha - 7.2973525693e-3) < 1e-12

        # α ≈ 1/137.036
        inverse = 1.0 / alpha
        assert abs(inverse - 137.036) < 0.001

    def test_h2_bond_length_calculation(self):
        """
        Test realistic bond length calculation for H₂.

        H₂ experimental bond length: 0.74 Å
        """
        # Estimate from covalent radii
        h_radius = PeriodicTable.get_covalent_radius('H')
        estimated_angstrom = 2 * h_radius

        # Should be reasonably close (covalent radii are approximate)
        experimental = 0.74
        assert abs(estimated_angstrom - experimental) < 0.2  # Within 0.2 Å

    def test_nacl_ionic_character(self):
        """
        Test ionic character determination for NaCl.

        NaCl is highly ionic (ΔEN = 2.23).
        """
        delta_en = PeriodicTable.electronegativity_difference('Na', 'Cl')

        # ΔEN > 1.7 indicates ionic bonding
        assert delta_en > 1.7
        assert abs(delta_en - 2.23) < 0.05

        # Ionic character ≈ 1 - exp(-0.25 × ΔEN²) ≈ 0.71 for NaCl
        ionic_character = 1.0 - np.exp(-0.25 * delta_en**2)
        assert 0.70 < ionic_character < 0.85


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
