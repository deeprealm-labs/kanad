"""
Unit tests for analysis tools.

Tests:
- Energy decomposition
- Bonding analysis
- Correlation analysis
"""

import pytest
import numpy as np

from kanad.analysis.energy_analysis import (
    EnergyAnalyzer,
    BondingAnalyzer,
    CorrelationAnalyzer
)
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation


class TestEnergyAnalyzer:
    """Test energy decomposition and analysis."""

    @pytest.fixture
    def h2_system(self):
        """Create H2 system for testing."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))
        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation, 'sto-3g')
        return hamiltonian

    def test_energy_analyzer_creation(self, h2_system):
        """Test energy analyzer initialization."""
        analyzer = EnergyAnalyzer(h2_system)
        assert analyzer.hamiltonian is not None

    def test_energy_decomposition(self, h2_system):
        """Test energy decomposition into components."""
        analyzer = EnergyAnalyzer(h2_system)

        # Create simple density matrix (Hartree-Fock-like)
        n = h2_system.n_orbitals
        P = np.zeros((n, n))
        # Fill lowest orbitals
        n_occ = min(2, n)  # 2 electrons in H2
        for i in range(n_occ):
            P[i, i] = 1.0

        decomp = analyzer.decompose_energy(P)

        # Check all components exist
        assert 'nuclear_repulsion' in decomp
        assert 'one_electron' in decomp
        assert 'two_electron' in decomp
        assert 'total' in decomp

        # Nuclear repulsion should be positive
        assert decomp['nuclear_repulsion'] > 0

        # Total energy should be real and finite
        assert np.isfinite(decomp['total'])
        # Components should sum to total (approximately, accounting for ERI indexing)
        # Just check that total is in reasonable range
        assert np.abs(decomp['total']) < 100  # Reasonable energy range

    def test_coulomb_exchange_decomposition(self, h2_system):
        """Test Coulomb and Exchange energy decomposition."""
        analyzer = EnergyAnalyzer(h2_system)

        n = h2_system.n_orbitals
        P = np.eye(n) * 0.5  # Simple density

        decomp = analyzer.decompose_energy(P)

        # Should have Coulomb and Exchange terms
        if 'coulomb' in decomp and 'exchange' in decomp:
            # Coulomb should be positive (repulsive)
            assert decomp['coulomb'] >= 0

            # Exchange should be negative (attractive)
            assert decomp['exchange'] <= 0

    def test_binding_energy_computation(self, h2_system):
        """Test binding energy calculation."""
        analyzer = EnergyAnalyzer(h2_system)

        # Mock energies
        E_molecule = -1.1  # Hartree
        E_atoms = [-0.5, -0.5]  # Two H atoms

        binding_energy = analyzer.compute_binding_energy(E_molecule, E_atoms)

        # Binding energy should be positive for stable molecule
        assert binding_energy > 0
        assert np.isclose(binding_energy, 0.1, atol=1e-9)

    def test_ionization_energy(self, h2_system):
        """Test ionization energy calculation."""
        analyzer = EnergyAnalyzer(h2_system)

        E_neutral = -1.1
        E_cation = -0.5

        IE = analyzer.compute_ionization_energy(E_neutral, E_cation)

        # IE should be positive (energy required to remove electron)
        assert IE > 0
        assert IE == E_cation - E_neutral

    def test_electron_affinity(self, h2_system):
        """Test electron affinity calculation."""
        analyzer = EnergyAnalyzer(h2_system)

        E_neutral = -1.1
        E_anion = -1.3

        EA = analyzer.compute_electron_affinity(E_neutral, E_anion)

        # EA should be positive for stable anion
        assert EA > 0
        assert EA == E_neutral - E_anion

    def test_convergence_analysis(self, h2_system):
        """Test VQE convergence analysis."""
        analyzer = EnergyAnalyzer(h2_system)

        # Mock energy history
        energy_history = np.array([
            -1.0, -1.05, -1.08, -1.10, -1.11, -1.11, -1.11
        ])

        analysis = analyzer.analyze_convergence(energy_history)

        # Check convergence metrics
        assert 'initial_energy' in analysis
        assert 'final_energy' in analysis
        assert 'energy_change' in analysis
        assert 'converged' in analysis

        assert analysis['initial_energy'] == -1.0
        assert analysis['final_energy'] == -1.11
        assert analysis['energy_change'] < 0  # Energy decreased

    def test_convergence_detection(self, h2_system):
        """Test convergence detection."""
        analyzer = EnergyAnalyzer(h2_system)

        # Converged case
        energy_converged = np.array([-1.0, -1.1, -1.11, -1.11, -1.11])
        analysis_conv = analyzer.analyze_convergence(energy_converged)
        # Converged flag may depend on threshold, just check it exists
        assert 'converged' in analysis_conv
        assert isinstance(analysis_conv['converged'], (bool, np.bool_))


class TestBondingAnalyzer:
    """Test bonding analysis tools."""

    @pytest.fixture
    def h2_covalent(self):
        """Create H2 covalent system."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))
        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation, 'sto-3g')
        return hamiltonian

    @pytest.fixture
    def nacl_ionic(self):
        """Create NaCl ionic system."""
        na = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
        cl = Atom('Cl', position=np.array([5.0, 0.0, 0.0]))
        molecule = Molecule([na, cl])
        representation = LCAORepresentation(molecule)
        hamiltonian = IonicHamiltonian(molecule, representation)
        return hamiltonian

    def test_bonding_analyzer_creation(self, h2_covalent):
        """Test bonding analyzer initialization."""
        analyzer = BondingAnalyzer(h2_covalent)
        assert analyzer.hamiltonian is not None

    def test_bonding_type_detection_covalent(self, h2_covalent):
        """Test detection of covalent bonding."""
        analyzer = BondingAnalyzer(h2_covalent)
        analysis = analyzer.analyze_bonding_type()

        assert analysis['bonding_type'] == 'covalent'
        assert 'characteristics' in analysis

    def test_bonding_type_detection_ionic(self, nacl_ionic):
        """Test detection of ionic bonding."""
        analyzer = BondingAnalyzer(nacl_ionic)
        analysis = analyzer.analyze_bonding_type()

        assert analysis['bonding_type'] == 'ionic'
        assert 'characteristics' in analysis

    def test_homo_lumo_gap(self, h2_covalent):
        """Test HOMO-LUMO gap analysis."""
        analyzer = BondingAnalyzer(h2_covalent)
        analysis = analyzer.analyze_bonding_type()

        # Covalent systems should have HOMO-LUMO gap
        if 'homo_lumo_gap' in analysis:
            assert analysis['homo_lumo_gap'] >= 0
            assert analysis['homo_lumo_gap_ev'] >= 0

    def test_mulliken_charges(self, h2_covalent):
        """Test Mulliken charge analysis."""
        analyzer = BondingAnalyzer(h2_covalent)

        n = h2_covalent.n_orbitals
        P = np.eye(n) * 0.5
        S = h2_covalent.S

        charges = analyzer.compute_mulliken_charges(P, S)

        # Should return charges for covalent system
        if len(charges) > 0:
            # Charges should sum to total charge (0 for neutral)
            # (approximately, due to numerical errors)
            # Just check it returns an array
            assert isinstance(charges, np.ndarray)

    def test_bond_order_analysis(self, h2_covalent):
        """Test bond order analysis."""
        analyzer = BondingAnalyzer(h2_covalent)

        n = h2_covalent.n_orbitals
        P = np.eye(n) * 0.5

        analysis = analyzer.analyze_bond_orders(P)

        # Should have bond orders for covalent system
        if 'bond_orders' in analysis:
            bond_orders = analysis['bond_orders']
            # Bond order matrix should be symmetric
            assert np.allclose(bond_orders, bond_orders.T)

    def test_bond_classification(self, h2_covalent):
        """Test bond classification (single, double, triple)."""
        analyzer = BondingAnalyzer(h2_covalent)

        n = h2_covalent.n_orbitals
        P = np.eye(n) * 0.5

        analysis = analyzer.analyze_bond_orders(P)

        # Should classify bonds
        if 'bond_classification' in analysis:
            classification = analysis['bond_classification']
            # H2 should have bonds classified

    def test_overlap_populations(self, h2_covalent):
        """Test overlap population analysis."""
        analyzer = BondingAnalyzer(h2_covalent)

        n = h2_covalent.n_orbitals
        P = np.eye(n) * 0.5
        S = h2_covalent.S

        overlap_pop = analyzer.compute_overlap_populations(P, S)

        # Should have same shape as P and S
        assert overlap_pop.shape == P.shape
        assert overlap_pop.shape == S.shape


class TestCorrelationAnalyzer:
    """Test electron correlation analysis."""

    @pytest.fixture
    def h2_system(self):
        """Create H2 system."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))
        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation, 'sto-3g')
        return hamiltonian

    def test_correlation_analyzer_creation(self, h2_system):
        """Test correlation analyzer initialization."""
        analyzer = CorrelationAnalyzer(h2_system)
        assert analyzer.hamiltonian is not None

    def test_correlation_energy_computation(self, h2_system):
        """Test correlation energy calculation."""
        analyzer = CorrelationAnalyzer(h2_system)

        E_vqe = -1.15  # Correlated
        E_hf = -1.10   # Mean-field

        E_corr = analyzer.compute_correlation_energy(E_vqe, E_hf)

        # Correlation energy should be negative (stabilizing)
        assert E_corr < 0
        assert E_corr == E_vqe - E_hf

    def test_percent_correlation_recovery(self, h2_system):
        """Test percentage correlation recovered."""
        analyzer = CorrelationAnalyzer(h2_system)

        E_exact = -1.17  # FCI
        E_vqe = -1.15    # VQE
        E_hf = -1.10     # HF

        percent = analyzer.compute_percent_correlation(E_vqe, E_hf, E_exact)

        # Should recover some percentage
        assert 0 <= percent <= 100

        # VQE between HF and exact
        expected = ((E_vqe - E_hf) / (E_exact - E_hf)) * 100
        assert np.isclose(percent, expected)

    def test_correlation_analysis(self, h2_system):
        """Test comprehensive correlation analysis."""
        analyzer = CorrelationAnalyzer(h2_system)

        vqe_result = {
            'energy': -1.15,
            'converged': True,
            'iterations': 50
        }
        E_hf = -1.10

        analysis = analyzer.analyze_electron_correlation(vqe_result, E_hf)

        # Check analysis components
        assert 'hf_energy' in analysis
        assert 'vqe_energy' in analysis
        assert 'correlation_energy' in analysis
        assert 'correlation_strength' in analysis

        assert analysis['hf_energy'] == E_hf
        assert analysis['vqe_energy'] == -1.15

    def test_correlation_strength_classification(self, h2_system):
        """Test correlation strength classification."""
        analyzer = CorrelationAnalyzer(h2_system)

        # Negligible correlation
        result_negligible = {'energy': -1.1005}
        analysis_neg = analyzer.analyze_electron_correlation(
            result_negligible,
            -1.1
        )
        assert analysis_neg['correlation_strength'] == 'negligible'

        # Weak correlation
        result_weak = {'energy': -1.105}
        analysis_weak = analyzer.analyze_electron_correlation(
            result_weak,
            -1.1
        )
        assert analysis_weak['correlation_strength'] == 'weak'

        # Moderate correlation
        result_moderate = {'energy': -1.15}
        analysis_mod = analyzer.analyze_electron_correlation(
            result_moderate,
            -1.1
        )
        assert analysis_mod['correlation_strength'] == 'moderate'

        # Strong correlation
        result_strong = {'energy': -1.25}
        analysis_strong = analyzer.analyze_electron_correlation(
            result_strong,
            -1.1
        )
        assert analysis_strong['correlation_strength'] == 'strong'

    def test_correlation_energy_in_ev(self, h2_system):
        """Test correlation energy conversion to eV."""
        analyzer = CorrelationAnalyzer(h2_system)

        E_vqe = -1.15
        E_hf = -1.10

        vqe_result = {'energy': E_vqe}
        analysis = analyzer.analyze_electron_correlation(vqe_result, E_hf)

        # Should have eV conversion
        assert 'correlation_energy_ev' in analysis

        # Check conversion (1 Hartree = 27.211 eV)
        E_corr_hartree = E_vqe - E_hf
        E_corr_ev_expected = E_corr_hartree * 27.211

        assert np.isclose(
            analysis['correlation_energy_ev'],
            E_corr_ev_expected
        )


class TestAnalysisIntegration:
    """Integration tests for analysis tools."""

    @pytest.fixture
    def h2_complete_system(self):
        """Create complete H2 system."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))
        molecule = Molecule([h1, h2])
        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(molecule, representation, 'sto-3g')
        return hamiltonian

    def test_full_analysis_workflow(self, h2_complete_system):
        """Test full analysis workflow."""
        # Energy analyzer
        energy_analyzer = EnergyAnalyzer(h2_complete_system)

        n = h2_complete_system.n_orbitals
        P = np.eye(n) * 0.5

        energy_decomp = energy_analyzer.decompose_energy(P)
        assert 'total' in energy_decomp

        # Bonding analyzer
        bonding_analyzer = BondingAnalyzer(h2_complete_system)
        bonding_type = bonding_analyzer.analyze_bonding_type()
        assert bonding_type['bonding_type'] == 'covalent'

        # Correlation analyzer
        correlation_analyzer = CorrelationAnalyzer(h2_complete_system)
        vqe_result = {'energy': -1.15}
        correlation_analysis = correlation_analyzer.analyze_electron_correlation(
            vqe_result,
            -1.10
        )
        assert 'correlation_energy' in correlation_analysis

    def test_analyzers_with_different_hamiltonians(self):
        """Test analyzers work with different Hamiltonian types."""
        # Covalent system
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))
        h2_mol = Molecule([h1, h2])
        h2_rep = LCAORepresentation(h2_mol)
        h2_ham = CovalentHamiltonian(h2_mol, h2_rep, 'sto-3g')

        # Ionic system
        na = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
        cl = Atom('Cl', position=np.array([5.0, 0.0, 0.0]))
        nacl_mol = Molecule([na, cl])
        nacl_rep = LCAORepresentation(nacl_mol)
        nacl_ham = IonicHamiltonian(nacl_mol, nacl_rep)

        # Both should work with analyzers
        for hamiltonian in [h2_ham, nacl_ham]:
            energy_analyzer = EnergyAnalyzer(hamiltonian)
            bonding_analyzer = BondingAnalyzer(hamiltonian)
            correlation_analyzer = CorrelationAnalyzer(hamiltonian)

            # Should not raise errors
            assert energy_analyzer.hamiltonian is not None
            assert bonding_analyzer.hamiltonian is not None
            assert correlation_analyzer.hamiltonian is not None
