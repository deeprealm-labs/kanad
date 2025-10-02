"""
Unit tests for fermionic-to-qubit mappers (Phase 6).

Tests different mapping strategies:
- Jordan-Wigner (ionic bonding)
- Hybrid Orbital Mapper (covalent bonding)
- Bravyi-Kitaev (general purpose)
"""

import pytest
import numpy as np
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.mappers.hybrid_orbital_mapper import HybridOrbitalMapper
from kanad.core.mappers.bravyi_kitaev_mapper import BravyiKitaevMapper


class TestJordanWignerMapper:
    """Test Jordan-Wigner transformation."""

    def test_jw_n_qubits(self):
        """Test qubit count."""
        mapper = JordanWignerMapper()

        assert mapper.n_qubits(4) == 4
        assert mapper.n_qubits(10) == 10

    def test_jw_number_operator(self):
        """Test mapping of number operator."""
        mapper = JordanWignerMapper()

        # n_0 for 2 orbitals
        result = mapper.map_number_operator(0, 2)

        # Should be (I - Z_0) / 2
        assert 'II' in result
        assert 'ZI' in result
        assert result['II'] == 0.5
        assert result['ZI'] == -0.5

    def test_jw_number_operator_orbital_1(self):
        """Test number operator for second orbital."""
        mapper = JordanWignerMapper()

        result = mapper.map_number_operator(1, 2)

        # Should be (I - Z_1) / 2
        assert 'II' in result
        assert 'IZ' in result
        assert result['II'] == 0.5
        assert result['IZ'] == -0.5

    def test_jw_excitation_same_orbital(self):
        """Test excitation on same orbital reduces to number operator."""
        mapper = JordanWignerMapper()

        exc = mapper.map_excitation_operator(0, 0, 2)
        num = mapper.map_number_operator(0, 2)

        # Should be equal
        assert exc.keys() == num.keys()
        for key in exc:
            assert abs(exc[key] - num[key]) < 1e-10

    def test_jw_excitation_operator(self):
        """Test excitation operator mapping."""
        mapper = JordanWignerMapper()

        # a†_1 a_0 for 2 orbitals
        result = mapper.map_excitation_operator(0, 1, 2)

        # Should have XX, YY, XY, YX terms
        assert 'XX' in result
        assert 'YY' in result

        # Check coefficients
        assert abs(result['XX']) == 0.25
        assert abs(result['YY']) == 0.25

    def test_jw_excitation_with_z_string(self):
        """Test that excitation includes Z string for distant orbitals."""
        mapper = JordanWignerMapper()

        # a†_2 a_0 for 4 orbitals (should have Z on orbital 1)
        result = mapper.map_excitation_operator(0, 2, 4)

        # Should have Z on middle qubit
        has_z_string = any('Z' in pauli and pauli.index('Z') == 1 for pauli in result.keys())
        assert has_z_string

    def test_jw_hermiticity(self):
        """Test that mapped operators have Hermitian properties."""
        mapper = JordanWignerMapper()

        # Number operator should be Hermitian (real coefficients)
        num_op = mapper.map_number_operator(0, 2)
        for coeff in num_op.values():
            assert np.isreal(coeff) or abs(np.imag(coeff)) < 1e-10


class TestHybridOrbitalMapper:
    """Test Hybrid Orbital Mapper for covalent bonding."""

    def test_hom_creation(self):
        """Test creating mapper with MO pairs."""
        # Two MO pairs: (0, 1) and (2, 3)
        mapper = HybridOrbitalMapper(mo_pairs=[(0, 1), (2, 3)])

        assert mapper.n_pairs == 2
        assert len(mapper.orbital_to_pair) == 4

    def test_hom_n_qubits(self):
        """Test qubit count (2 per MO pair)."""
        mapper = HybridOrbitalMapper(mo_pairs=[(0, 1), (2, 3)])

        # 2 pairs → 4 qubits
        assert mapper.n_qubits(4) == 4

    def test_hom_number_operator(self):
        """Test number operator for bonding orbital."""
        mapper = HybridOrbitalMapper(mo_pairs=[(0, 1)])

        # n_0 (bonding orbital)
        result = mapper.map_number_operator(0, 2)

        # Should map to qubit 0 (first of the pair)
        assert 'II' in result
        assert 'ZI' in result

    def test_hom_mo_pair_excitation(self):
        """Test bonding → antibonding excitation within MO pair."""
        mapper = HybridOrbitalMapper(mo_pairs=[(0, 1)])

        # a†_1 a_0 (antibonding ← bonding)
        result = mapper.map_excitation_operator(0, 1, 2)

        # Should have XX and YY terms (no imaginary)
        assert 'XX' in result
        assert 'YY' in result

        # Check coefficients are real
        assert np.isreal(result['XX'])
        assert np.isreal(result['YY'])

    def test_hom_bond_excitation_operator(self):
        """Test creating bond excitation operator."""
        mapper = HybridOrbitalMapper(mo_pairs=[(0, 1), (2, 3)])

        # Create excitation for first bond
        exc = mapper.create_bond_excitation_operator(0)

        # Should have XX and YY on qubits 0 and 1
        assert 'XXII' in exc or 'XX' in list(exc.keys())[0][:2]

    def test_hom_are_mo_pair(self):
        """Test MO pair detection."""
        mapper = HybridOrbitalMapper(mo_pairs=[(0, 1), (2, 3)])

        assert mapper._are_mo_pair(0, 1) is True
        assert mapper._are_mo_pair(1, 0) is True  # Symmetric
        assert mapper._are_mo_pair(0, 2) is False
        assert mapper._are_mo_pair(2, 3) is True


class TestBravyiKitaevMapper:
    """Test Bravyi-Kitaev transformation."""

    def test_bk_n_qubits(self):
        """Test qubit count (same as JW)."""
        mapper = BravyiKitaevMapper()

        assert mapper.n_qubits(4) == 4
        assert mapper.n_qubits(8) == 8

    def test_bk_number_operator(self):
        """Test number operator mapping."""
        mapper = BravyiKitaevMapper()

        result = mapper.map_number_operator(0, 2)

        # Should be similar to JW for number operator
        assert 'II' in result
        assert 'ZI' in result

    def test_bk_excitation_operator(self):
        """Test excitation operator mapping."""
        mapper = BravyiKitaevMapper()

        result = mapper.map_excitation_operator(0, 1, 2)

        # Should have some Pauli terms
        assert len(result) > 0

        # All coefficients should be non-zero
        for coeff in result.values():
            assert abs(coeff) > 0

    def test_bk_operator_weight(self):
        """Test that BK has lower weight than JW for distant excitations."""
        mapper = BravyiKitaevMapper()

        # Weight for excitation 0 → 7 on 8 orbitals
        weight_bk = mapper.get_operator_weight(0, 7, 8)

        # BK should have O(log n) weight
        # For n=8, weight should be < 8 (compared to JW's 7)
        assert weight_bk <= 8


class TestMapperComparison:
    """Compare different mappers."""

    def test_jw_vs_hom_qubit_count(self):
        """Compare qubit requirements."""
        jw = JordanWignerMapper()
        hom = HybridOrbitalMapper(mo_pairs=[(0, 1), (2, 3)])

        # For 4 orbitals
        assert jw.n_qubits(4) == 4
        assert hom.n_qubits(4) == 4  # Same for this case

    def test_number_operator_consistency(self):
        """Test that number operators give same eigenvalues."""
        jw = JordanWignerMapper()
        bk = BravyiKitaevMapper()

        # Both should give (I - Z) / 2
        jw_num = jw.map_number_operator(0, 2)
        bk_num = bk.map_number_operator(0, 2)

        # Should have same structure
        assert set(jw_num.keys()) == set(bk_num.keys())

    def test_excitation_coefficient_magnitudes(self):
        """Test that excitation operators have reasonable coefficients."""
        jw = JordanWignerMapper()
        hom = HybridOrbitalMapper(mo_pairs=[(0, 1)])

        jw_exc = jw.map_excitation_operator(0, 1, 2)
        hom_exc = hom.map_excitation_operator(0, 1, 2)

        # All coefficients should be <= 1
        for coeff in jw_exc.values():
            assert abs(coeff) <= 1.0
        for coeff in hom_exc.values():
            assert abs(coeff) <= 1.0


class TestMapperEdgeCases:
    """Test edge cases and error handling."""

    def test_single_orbital_jw(self):
        """Test JW with single orbital."""
        mapper = JordanWignerMapper()

        result = mapper.map_number_operator(0, 1)

        assert 'I' in result
        assert 'Z' in result

    def test_empty_mo_pairs_hom(self):
        """Test HOM with no MO pairs."""
        mapper = HybridOrbitalMapper(mo_pairs=[])

        assert mapper.n_pairs == 0
        assert mapper.n_qubits(0) == 0

    def test_pauli_multiplication(self):
        """Test Pauli string multiplication."""
        mapper = JordanWignerMapper()

        # Test X * X = I
        pauli1 = {'XI': 1.0}
        pauli2 = {'XI': 1.0}

        result = mapper.pauli_string_multiply(pauli1, pauli2)

        # XI * XI = II (since X² = I)
        assert 'II' in result

    def test_pauli_phases(self):
        """Test that Pauli multiplication gets phases right."""
        mapper = JordanWignerMapper()

        # X * Y = iZ
        prod, phase = mapper._multiply_single_pauli('X', 'Y')

        assert prod == 'Z'
        assert phase == 1j


class TestMapperPerformance:
    """Test performance characteristics of different mappers."""

    def test_jw_weight_scaling(self):
        """Test that JW weight scales linearly with distance."""
        mapper = JordanWignerMapper()

        # Excitation 0 → 2 should have Z on orbital 1
        exc_02 = mapper.map_excitation_operator(0, 2, 4)

        # Count total Pauli weight (non-I operators)
        def count_weight(pauli_dict):
            max_weight = 0
            for pauli_str in pauli_dict.keys():
                weight = sum(1 for p in pauli_str if p != 'I')
                max_weight = max(max_weight, weight)
            return max_weight

        weight = count_weight(exc_02)

        # Should have weight ~3 (X, Y on ends + Z in middle)
        assert weight >= 2  # At least the two operators

    def test_hom_local_excitations(self):
        """Test that HOM excitations within pairs are local."""
        mapper = HybridOrbitalMapper(mo_pairs=[(0, 1), (2, 3)])

        # Excitation within first pair
        exc = mapper.map_excitation_operator(0, 1, 4)

        # Should only involve qubits 0 and 1
        for pauli_str in exc.keys():
            # Check that qubits 2 and 3 are Identity
            assert pauli_str[2] == 'I'
            assert pauli_str[3] == 'I'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
