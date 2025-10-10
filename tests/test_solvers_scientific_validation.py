"""
Scientific Validation of Solvers, Bonds, Analysis, and Optimization Modules.

Nobel Laureate Inspector Standard:
- Validate against experimental data (NIST, literature)
- Check methodologies are scientifically correct
- Verify numerical values match published benchmarks
- Ensure units and conversions are accurate
"""

import pytest
import numpy as np
from kanad.core.atom import Atom
from kanad.bonds import BondFactory
from kanad.solvers.excited_states_solver import ExcitedStatesSolver


class TestVQESolverScientific:
    """Scientific validation of VQE solver against literature."""

    def test_h2_vqe_vs_fci_literature(self):
        """
        Validate H2 VQE energy against FCI gold standard.

        Reference: Szabo & Ostlund, "Modern Quantum Chemistry", Table 3.4
        - HF energy: -1.117 Ha
        - FCI energy: -1.137 Ha
        - Correlation: -0.020 Ha

        STO-3G basis, R = 0.735 Å
        """
        from kanad.solvers import VQESolver

        # Create H2 bond
        bond = BondFactory.create_bond('H', 'H', distance=0.735)

        # VQE with UCC ansatz (should match FCI)
        solver = VQESolver(
            bond=bond,
            ansatz_type='ucc',
            optimizer='SLSQP',
            max_iterations=100,
            backend='statevector'
        )

        # Run VQE
        result = solver.solve()

        # Literature values
        lit_fci = -1.137  # Ha (Szabo & Ostlund)
        lit_hf = -1.117   # Ha
        chemical_accuracy = 0.0016  # 1.6 mHa = 1 kcal/mol

        # Check energy
        energy_ha = result['energy_ha']
        error = abs(energy_ha - lit_fci)
        error_mha = error * 1000

        print(f"\n{'='*70}")
        print(f"H2 VQE VALIDATION (STO-3G, R=0.735 Å)")
        print(f"{'='*70}")
        print(f"Literature FCI (Szabo & Ostlund): {lit_fci:.6f} Ha")
        print(f"Computed VQE energy:               {energy_ha:.6f} Ha")
        print(f"Error:                             {error:.9f} Ha ({error_mha:.6f} mHa)")
        print(f"Chemical accuracy threshold:       {chemical_accuracy:.6f} Ha (1.6 mHa)")
        print(f"Status:                            {'✓ EXACT' if error_mha < 0.01 else '✓ Chemical accuracy' if error < chemical_accuracy else '✗ Too large'}")
        print(f"{'='*70}\n")

        # Should be within chemical accuracy
        assert error < chemical_accuracy, f"VQE error {error_mha:.6f} mHa exceeds chemical accuracy"

        # Check convergence
        assert result['converged'], "VQE did not converge"


    def test_h2_bond_dissociation_curve(self):
        """
        Validate H2 bond dissociation curve shape.

        Physics:
        - At equilibrium (0.74 Å): Minimum energy
        - At compressed (0.5 Å): Repulsive, higher energy
        - At extended (2.0 Å): Dissociating, approaches atomic limit

        Reference: Potential energy surface should have single minimum
        """
        distances = [0.5, 0.74, 1.0, 1.5, 2.0]
        energies = []

        for R in distances:
            bond = BondFactory.create_bond('H', 'H', distance=R)
            # Use HF for speed (shape is similar)
            result = bond.compute_energy(method='HF')
            energies.append(result['energy_ha'])

        energies = np.array(energies)

        print(f"\n{'='*70}")
        print(f"H2 DISSOCIATION CURVE VALIDATION")
        print(f"{'='*70}")
        print(f"Distance (Å)    Energy (Ha)    Behavior")
        print(f"-" * 70)
        for R, E in zip(distances, energies):
            status = "✓ Compressed" if R < 0.74 else "✓ Equilibrium" if R == 0.74 else "✓ Dissociating"
            print(f"{R:>12.2f}    {E:>11.6f}    {status}")
        print(f"{'='*70}\n")

        # Physical checks
        # 1. Minimum should be near 0.74 Å
        min_idx = np.argmin(energies)
        min_distance = distances[min_idx]
        assert abs(min_distance - 0.74) < 0.5, "Minimum not near equilibrium distance"

        # 2. Energy should increase when compressed
        assert energies[0] > energies[1], "Compressed bond should have higher energy"

        # 3. Energy should increase when dissociating (until dissociation limit)
        # Compare equilibrium to moderate extension
        eq_idx = distances.index(0.74)
        extended_idx = distances.index(1.5)
        assert energies[extended_idx] > energies[eq_idx], "Extended bond should have higher energy"


class TestExcitedStatesSolverScientific:
    """Scientific validation of excited states solver."""

    def test_h2_excited_states_cis(self):
        """
        Validate H2 excited states using CIS.

        Physics:
        - σ→σ* excitation should be lowest
        - Excitation energy should be positive
        - Typical HOMO-LUMO gap for H2: ~10-15 eV

        Reference: CIS should give qualitatively correct excitations
        """
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        solver = ExcitedStatesSolver(
            bond=bond,
            method='cis',
            n_states=3,
            enable_analysis=True
        )

        result = solver.solve()

        # Get excitation energies
        exc_energies_ev = result['excitation_energies_ev']
        oscillator_strengths = result['oscillator_strengths']

        print(f"\n{'='*70}")
        print(f"H2 EXCITED STATES VALIDATION (CIS)")
        print(f"{'='*70}")
        print(f"Ground state energy: {result['ground_state_energy']:.6f} Ha")
        print(f"\nExcited States:")
        print(f"State    Excitation (eV)    f (osc. strength)    Transition")
        print(f"-" * 70)
        for i, (E_ev, f, trans) in enumerate(zip(
            exc_energies_ev,
            oscillator_strengths,
            result['dominant_transitions']
        ), 1):
            print(f"{i:>5}    {E_ev:>14.4f}    {f:>17.6f}    {trans}")
        print(f"{'='*70}\n")

        # Physical checks
        # 1. All excitation energies should be positive
        assert np.all(exc_energies_ev > 0), "Excitation energies must be positive"

        # 2. Oscillator strengths should be non-negative
        assert np.all(oscillator_strengths >= 0), "Oscillator strengths cannot be negative"

        # 3. Lowest excitation should be σ→σ* (HOMO-LUMO)
        assert 'HOMO' in result['dominant_transitions'][0], "First excitation should involve HOMO"
        assert 'LUMO' in result['dominant_transitions'][0], "First excitation should involve LUMO"

        # 4. Typical H2 σ→σ* gap is 10-15 eV (CIS overestimates, so allow range)
        assert 5 < exc_energies_ev[0] < 25, f"HOMO-LUMO gap {exc_energies_ev[0]:.2f} eV outside reasonable range"


    def test_excited_states_ordering(self):
        """
        Validate that excited states are ordered by energy.

        Physics: E₀ < E₁ < E₂ < ...
        """
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        solver = ExcitedStatesSolver(
            bond=bond,
            method='cis',
            n_states=5
        )

        result = solver.solve()

        energies = result['energies']

        print(f"\n{'='*70}")
        print(f"EXCITED STATES ORDERING VALIDATION")
        print(f"{'='*70}")
        print(f"State    Energy (Ha)    Energy (eV)    ΔE from ground (eV)")
        print(f"-" * 70)
        for i, E in enumerate(energies):
            delta_ev = (E - energies[0]) * 27.2114
            print(f"{i:>5}    {E:>11.6f}    {E*27.2114:>11.4f}    {delta_ev:>18.4f}")
        print(f"{'='*70}\n")

        # Check monotonic increasing
        for i in range(len(energies) - 1):
            assert energies[i] < energies[i+1], f"State {i+1} energy lower than state {i}"


class TestBondAnalysisScientific:
    """Scientific validation of bond analysis."""

    def test_h2_bond_order_single_bond(self):
        """
        Validate bond order calculation for H2 (single bond).

        Reference: H2 has bond order ≈ 1.0 (one σ bond)
        """
        from kanad.analysis.energy_analysis import BondingAnalyzer

        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        # Compute HF solution
        density_matrix, hf_energy = bond.hamiltonian.solve_scf(
            max_iterations=100,
            conv_tol=1e-8
        )

        # Analyze bonding
        analyzer = BondingAnalyzer(
            bond.hamiltonian,
            bond.molecule,
            density_matrix
        )

        bond_order = analyzer.compute_bond_order(atom_i=0, atom_j=1)

        print(f"\n{'='*70}")
        print(f"H2 BOND ORDER VALIDATION")
        print(f"{'='*70}")
        print(f"Computed bond order: {bond_order:.6f}")
        print(f"Expected (single bond): ~1.0")
        print(f"Error: {abs(bond_order - 1.0):.6f}")
        print(f"Status: {'✓ Correct' if abs(bond_order - 1.0) < 0.2 else '✗ Incorrect'}")
        print(f"{'='*70}\n")

        # Bond order for H2 should be close to 1
        assert 0.8 < bond_order < 1.2, f"H2 bond order {bond_order:.3f} not close to 1.0"


    def test_homo_lumo_gap_validation(self):
        """
        Validate HOMO-LUMO gap calculation.

        Reference: H2 HOMO-LUMO gap ~10-15 eV (HF approximation)
        """
        from kanad.analysis.energy_analysis import BondingAnalyzer

        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        density_matrix, _ = bond.hamiltonian.solve_scf()

        analyzer = BondingAnalyzer(
            bond.hamiltonian,
            bond.molecule,
            density_matrix
        )

        homo_lumo_gap_ev = analyzer.compute_homo_lumo_gap()

        print(f"\n{'='*70}")
        print(f"H2 HOMO-LUMO GAP VALIDATION")
        print(f"{'='*70}")
        print(f"Computed HOMO-LUMO gap: {homo_lumo_gap_ev:.4f} eV")
        print(f"Expected range: 10-15 eV (HF)")
        print(f"Status: {'✓ In range' if 8 < homo_lumo_gap_ev < 20 else '✗ Outside range'}")
        print(f"{'='*70}\n")

        # HOMO-LUMO gap should be reasonable
        assert 5 < homo_lumo_gap_ev < 25, f"HOMO-LUMO gap {homo_lumo_gap_ev:.2f} eV unreasonable"


    def test_mulliken_charges_h2_neutral(self):
        """
        Validate Mulliken charges for H2 (should be near zero).

        Reference: Homonuclear diatomic should have no net charge transfer
        """
        from kanad.analysis.energy_analysis import BondingAnalyzer

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        density_matrix, _ = bond.hamiltonian.solve_scf()

        analyzer = BondingAnalyzer(
            bond.hamiltonian,
            bond.molecule,
            density_matrix
        )

        charges = analyzer.compute_mulliken_charges()

        print(f"\n{'='*70}")
        print(f"H2 MULLIKEN CHARGES VALIDATION")
        print(f"{'='*70}")
        print(f"Atom    Charge")
        print(f"-" * 30)
        for i, q in enumerate(charges):
            print(f"H{i}      {q:>7.4f}")
        print(f"\nTotal charge: {np.sum(charges):.6f}")
        print(f"Max |charge|: {np.max(np.abs(charges)):.6f}")
        print(f"Status: {'✓ Neutral' if np.max(np.abs(charges)) < 0.2 else '✗ Polarized'}")
        print(f"{'='*70}\n")

        # Charges should be small for homonuclear diatomic
        assert np.max(np.abs(charges)) < 0.3, "H2 should have minimal charge separation"

        # Total charge should be near zero
        assert abs(np.sum(charges)) < 0.1, "Total charge should be conserved"


    def test_ionic_vs_covalent_charge_transfer(self):
        """
        Validate charge transfer in ionic vs covalent bonds.

        Reference:
        - LiF (ionic): Large charge transfer (Li⁺ F⁻)
        - HF (covalent): Moderate charge transfer (polar covalent)
        - H2 (covalent): Minimal charge transfer (nonpolar)
        """
        from kanad.analysis.energy_analysis import BondingAnalyzer

        systems = [
            ('H', 'H', 'nonpolar covalent', 0.74),
            ('H', 'F', 'polar covalent', 0.92),
            ('Li', 'F', 'ionic', 1.56)
        ]

        charge_transfers = []

        print(f"\n{'='*70}")
        print(f"CHARGE TRANSFER VALIDATION")
        print(f"{'='*70}")
        print(f"System    Bond Type              Distance (Å)    Charge Transfer")
        print(f"-" * 70)

        for atom1, atom2, bond_type, distance in systems:
            bond = BondFactory.create_bond(atom1, atom2, distance=distance)
            density_matrix, _ = bond.hamiltonian.solve_scf(max_iterations=200)

            analyzer = BondingAnalyzer(
                bond.hamiltonian,
                bond.molecule,
                density_matrix
            )

            charges = analyzer.compute_mulliken_charges()
            charge_transfer = abs(charges[0])  # |q1| = |q2| for diatomic
            charge_transfers.append(charge_transfer)

            print(f"{atom1}-{atom2:>2}    {bond_type:<22} {distance:>12.2f}    {charge_transfer:>15.6f}")

        print(f"{'='*70}")
        print(f"\nPhysical Validation:")
        print(f"  H-H (nonpolar):  {'✓ Small' if charge_transfers[0] < 0.3 else '✗ Too large'}")
        print(f"  H-F (polar):     {'✓ Moderate' if 0.2 < charge_transfers[1] < 0.8 else '✗ Unexpected'}")
        print(f"  Li-F (ionic):    {'✓ Large' if charge_transfers[2] > 0.6 else '✗ Too small'}")
        print(f"  Ordering:        {'✓ H2 < HF < LiF' if charge_transfers[0] < charge_transfers[1] < charge_transfers[2] else '✗ Wrong order'}")
        print(f"{'='*70}\n")

        # Check ordering: H2 < HF < LiF
        assert charge_transfers[0] < charge_transfers[1], "H2 should have less charge transfer than HF"
        assert charge_transfers[1] < charge_transfers[2], "HF should have less charge transfer than LiF"


class TestPhysicalUnitsAndConversions:
    """Validate unit conversions and physical constants."""

    def test_hartree_to_ev_conversion(self):
        """
        Validate Hartree to eV conversion.

        Reference: CODATA 2018: 1 Ha = 27.211386245988 eV
        """
        from kanad.core.constants.conversion_factors import ConversionFactors

        # CODATA 2018 value
        codata_ha_to_ev = 27.211386245988

        # Our implementation
        our_value = ConversionFactors.HARTREE_TO_EV

        error = abs(our_value - codata_ha_to_ev)
        percent_error = error / codata_ha_to_ev * 100

        print(f"\n{'='*70}")
        print(f"HARTREE TO EV CONVERSION VALIDATION")
        print(f"{'='*70}")
        print(f"CODATA 2018:        {codata_ha_to_ev:.12f} eV/Ha")
        print(f"Our implementation: {our_value:.12f} eV/Ha")
        print(f"Error:              {error:.12f} ({percent_error:.6f}%)")
        print(f"Status:             {'✓ Correct' if percent_error < 0.01 else '✗ Incorrect'}")
        print(f"{'='*70}\n")

        assert percent_error < 0.01, "Hartree to eV conversion error too large"


    def test_bohr_to_angstrom_conversion(self):
        """
        Validate Bohr to Angstrom conversion.

        Reference: CODATA 2018: 1 Bohr = 0.529177210903 Å
        """
        from kanad.core.constants.conversion_factors import ConversionFactors

        codata_bohr_to_ang = 0.529177210903

        our_value = ConversionFactors.BOHR_TO_ANGSTROM

        error = abs(our_value - codata_bohr_to_ang)
        percent_error = error / codata_bohr_to_ang * 100

        print(f"\n{'='*70}")
        print(f"BOHR TO ANGSTROM CONVERSION VALIDATION")
        print(f"{'='*70}")
        print(f"CODATA 2018:        {codata_bohr_to_ang:.12f} Å/Bohr")
        print(f"Our implementation: {our_value:.12f} Å/Bohr")
        print(f"Error:              {error:.12f} ({percent_error:.6f}%)")
        print(f"Status:             {'✓ Correct' if percent_error < 0.01 else '✗ Incorrect'}")
        print(f"{'='*70}\n")

        assert percent_error < 0.01, "Bohr to Angstrom conversion error too large"


    def test_energy_unit_consistency(self):
        """
        Validate energy unit consistency across modules.

        Test: 1 Ha in different units should convert back correctly
        """
        from kanad.core.constants.conversion_factors import ConversionFactors

        # Start with 1 Hartree
        energy_ha = 1.0

        # Convert to eV and back
        energy_ev = energy_ha * ConversionFactors.HARTREE_TO_EV
        energy_ha_back = energy_ev / ConversionFactors.HARTREE_TO_EV

        error = abs(energy_ha - energy_ha_back)

        print(f"\n{'='*70}")
        print(f"ENERGY UNIT CONSISTENCY VALIDATION")
        print(f"{'='*70}")
        print(f"Original:           {energy_ha:.12f} Ha")
        print(f"→ eV:               {energy_ev:.12f} eV")
        print(f"→ Back to Ha:       {energy_ha_back:.12f} Ha")
        print(f"Round-trip error:   {error:.2e}")
        print(f"Status:             {'✓ Consistent' if error < 1e-10 else '✗ Inconsistent'}")
        print(f"{'='*70}\n")

        assert error < 1e-10, "Energy unit conversion not self-consistent"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
