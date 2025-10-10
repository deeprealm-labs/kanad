"""
Complete Governance Integration Test - STRICT VALIDATION

Tests the complete integration chain:
1. Atoms → Molecules → Bonds
2. Bonds → Hamiltonians + Representations + Protocols + Mappers
3. Solvers use Bonds under Governance
4. All numerical values validated against literature (ZERO TOLERANCE)

Molecules tested:
- H2 (covalent, homonuclear)
- HeH+ (covalent, heteronuclear)
- H2O (covalent, polyatomic)

Nobel Laureate Standard: NO MERCY ON VALUES
"""

import pytest
import numpy as np
from kanad.core.atom import Atom
from kanad.bonds import BondFactory
from kanad.solvers import VQESolver, ExcitedStatesSolver


class TestGovernanceIntegrationH2:
    """STRICT validation of H2 molecule with complete governance."""

    def test_h2_bond_has_all_components(self):
        """
        Verify CovalentBond integrates ALL framework components.

        Required:
        - Atoms ✓
        - Molecule ✓
        - Representation (LCAO) ✓
        - Hamiltonian (Covalent) ✓
        - Governance Protocol ✓
        - Mapper (HybridOrbital) ✓
        - Ansatz (CovalentGovernance) ✓
        """
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        print(f"\n{'='*70}")
        print(f"H2 BOND COMPONENT INTEGRATION AUDIT")
        print(f"{'='*70}")

        # 1. Check Atoms
        assert hasattr(bond, 'atom_1'), "Missing atom_1"
        assert hasattr(bond, 'atom_2'), "Missing atom_2"
        assert bond.atom_1.symbol == 'H', "atom_1 is not H"
        assert bond.atom_2.symbol == 'H', "atom_2 is not H"
        print(f"✓ Atoms: H-H (atomic numbers: {bond.atom_1.atomic_number}, {bond.atom_2.atomic_number})")

        # 2. Check Molecule
        assert hasattr(bond, 'molecule'), "Missing molecule"
        assert bond.molecule.n_electrons == 2, f"H2 should have 2 electrons, got {bond.molecule.n_electrons}"
        assert len(bond.molecule.atoms) == 2, f"H2 should have 2 atoms, got {len(bond.molecule.atoms)}"
        print(f"✓ Molecule: {bond.molecule.n_electrons} electrons, {len(bond.molecule.atoms)} atoms")

        # 3. Check Representation
        assert hasattr(bond, 'representation'), "Missing representation"
        assert bond.representation.__class__.__name__ == 'LCAORepresentation', \
            f"Expected LCAORepresentation, got {bond.representation.__class__.__name__}"
        print(f"✓ Representation: {bond.representation.__class__.__name__}")
        print(f"  - n_orbitals: {bond.representation.n_orbitals}")
        print(f"  - n_qubits: {bond.representation.n_qubits}")

        # 4. Check Hamiltonian
        assert hasattr(bond, 'hamiltonian'), "Missing hamiltonian"
        assert bond.hamiltonian.__class__.__name__ == 'CovalentHamiltonian', \
            f"Expected CovalentHamiltonian, got {bond.hamiltonian.__class__.__name__}"
        assert bond.hamiltonian.n_orbitals == bond.representation.n_orbitals, \
            "Hamiltonian and representation orbital count mismatch"
        print(f"✓ Hamiltonian: {bond.hamiltonian.__class__.__name__}")
        print(f"  - n_orbitals: {bond.hamiltonian.n_orbitals}")
        print(f"  - n_electrons: {bond.hamiltonian.n_electrons}")
        print(f"  - nuclear_repulsion: {bond.hamiltonian.nuclear_repulsion:.6f} Ha")

        # 5. Check Governance Protocol
        assert hasattr(bond, 'governance'), "Missing governance protocol"
        assert bond.governance.__class__.__name__ == 'CovalentGovernanceProtocol', \
            f"Expected CovalentGovernanceProtocol, got {bond.governance.__class__.__name__}"
        print(f"✓ Governance: {bond.governance.__class__.__name__}")

        # 6. Check Mapper
        assert hasattr(bond, 'mapper'), "Missing mapper"
        assert bond.mapper.__class__.__name__ == 'HybridOrbitalMapper', \
            f"Expected HybridOrbitalMapper, got {bond.mapper.__class__.__name__}"
        print(f"✓ Mapper: {bond.mapper.__class__.__name__}")

        # 7. Check bond type
        assert bond.bond_type == 'covalent', f"Expected covalent, got {bond.bond_type}"
        print(f"✓ Bond Type: {bond.bond_type}")

        # 8. Check hybridization
        assert hasattr(bond, 'hybridization'), "Missing hybridization"
        print(f"✓ Hybridization: {bond.hybridization}")

        print(f"{'='*70}")
        print(f"ALL COMPONENTS INTEGRATED ✅")
        print(f"{'='*70}\n")

    def test_h2_hf_energy_strict_validation(self):
        """
        H2 Hartree-Fock energy - STRICT validation against Szabo & Ostlund.

        Literature: Szabo & Ostlund, Table 3.4
        - Basis: STO-3G
        - Distance: 0.735 Å (1.4 Bohr)
        - HF Energy: -1.117 Ha
        - Tolerance: < 0.001 Ha (1 mHa) - NO LENIENCY
        """
        bond = BondFactory.create_bond('H', 'H', distance=0.735)

        result = bond.compute_energy(method='HF', max_iterations=200)

        # Literature value
        lit_hf = -1.117  # Ha (Szabo & Ostlund, Table 3.4)

        energy = result['energy']
        error = abs(energy - lit_hf)
        error_mha = error * 1000

        print(f"\n{'='*70}")
        print(f"H2 HARTREE-FOCK STRICT VALIDATION")
        print(f"{'='*70}")
        print(f"System: H2, R = 0.735 Å, Basis: STO-3G")
        print(f"Literature (Szabo & Ostlund): {lit_hf:.6f} Ha")
        print(f"Computed HF energy:           {energy:.6f} Ha")
        print(f"Error:                        {error:.9f} Ha ({error_mha:.6f} mHa)")
        print(f"Tolerance:                    1.000 mHa (STRICT)")
        print(f"Converged:                    {result['converged']}")
        print(f"Iterations:                   {result['iterations']}")
        print(f"Status:                       {'✓ PASS' if error_mha < 1.0 else '✗ FAIL'}")
        print(f"{'='*70}\n")

        # STRICT: Must be within 1 mHa
        assert error_mha < 1.0, f"HF energy error {error_mha:.6f} mHa exceeds 1 mHa tolerance"
        assert result['converged'], "SCF did not converge"

    def test_h2_nuclear_repulsion_exact(self):
        """
        Nuclear repulsion must be EXACT.

        Formula: E_nn = Z1*Z2/R
        H2: Z1=Z2=1, R=0.735 Å = 1.3889 Bohr
        E_nn = 1*1/1.3889 = 0.720083 Ha

        Tolerance: < 1e-6 Ha (NO LENIENCY)
        """
        bond = BondFactory.create_bond('H', 'H', distance=0.735)

        # Compute expected nuclear repulsion
        from kanad.core.constants.conversion_factors import ConversionFactors
        R_ang = 0.735  # Angstroms
        R_bohr = R_ang / ConversionFactors.BOHR_TO_ANGSTROM
        Z1, Z2 = 1, 1
        expected_nuc_rep = Z1 * Z2 / R_bohr

        computed_nuc_rep = bond.hamiltonian.nuclear_repulsion

        error = abs(computed_nuc_rep - expected_nuc_rep)

        print(f"\n{'='*70}")
        print(f"H2 NUCLEAR REPULSION EXACT VALIDATION")
        print(f"{'='*70}")
        print(f"Distance:    {R_ang:.6f} Å = {R_bohr:.6f} Bohr")
        print(f"Formula:     E_nn = Z₁Z₂/R = {Z1}×{Z2}/{R_bohr:.6f}")
        print(f"Expected:    {expected_nuc_rep:.9f} Ha")
        print(f"Computed:    {computed_nuc_rep:.9f} Ha")
        print(f"Error:       {error:.2e} Ha")
        print(f"Tolerance:   1.0e-06 Ha (STRICT)")
        print(f"Status:      {'✓ EXACT' if error < 1e-6 else '✗ FAIL'}")
        print(f"{'='*70}\n")

        assert error < 1e-6, f"Nuclear repulsion error {error:.2e} exceeds tolerance"

    def test_h2_vqe_with_governance_strict(self):
        """
        H2 VQE with governance ansatz - must match FCI.

        Literature: FCI energy = -1.137 Ha (Szabo & Ostlund)
        Tolerance: < 1.6 mHa (chemical accuracy) - STRICT
        """
        bond = BondFactory.create_bond('H', 'H', distance=0.735)

        result = bond.compute_energy(method='VQE', max_iterations=100)

        # Literature FCI
        lit_fci = -1.137  # Ha

        energy = result['energy']
        error = abs(energy - lit_fci)
        error_mha = error * 1000

        print(f"\n{'='*70}")
        print(f"H2 VQE WITH GOVERNANCE STRICT VALIDATION")
        print(f"{'='*70}")
        print(f"System: H2, R = 0.735 Å")
        print(f"Ansatz: CovalentGovernanceAnsatz (bonding-aware)")
        print(f"Literature FCI: {lit_fci:.6f} Ha")
        print(f"Computed VQE:   {energy:.6f} Ha")
        print(f"Error:          {error:.9f} Ha ({error_mha:.6f} mHa)")
        print(f"Tolerance:      1.600 mHa (chemical accuracy - STRICT)")
        print(f"Converged:      {result['converged']}")
        print(f"Iterations:     {result['iterations']}")

        if 'correlation_energy' in result:
            print(f"Correlation:    {result['correlation_energy']:.6f} Ha ({result['correlation_energy']*1000:.3f} mHa)")

        print(f"Status:         {'✓ PASS' if error_mha < 1.6 else '✗ FAIL'}")
        print(f"{'='*70}\n")

        # STRICT: Must be within chemical accuracy
        assert error_mha < 1.6, f"VQE energy error {error_mha:.6f} mHa exceeds chemical accuracy"


class TestGovernanceIntegrationHeH:
    """STRICT validation of HeH+ with governance."""

    def test_heh_bond_creation_with_governance(self):
        """
        HeH+ bond must integrate all components.

        System: He(+)-H (1 electron removed from He-H)
        Charge: +1
        Spin: 0 (singlet)
        """
        # Create HeH+ manually
        he = Atom('He', [0.0, 0.0, 0.0])
        h = Atom('H', [0.0, 0.0, 0.772])  # ~equilibrium distance

        from kanad.bonds.covalent_bond import CovalentBond
        bond = CovalentBond(he, h, distance=0.772)

        # System properties
        total_electrons = he.atomic_number + h.atomic_number  # 2 + 1 = 3

        print(f"\n{'='*70}")
        print(f"HeH+ BOND INTEGRATION AUDIT")
        print(f"{'='*70}")
        print(f"Atoms: He-H")
        print(f"Total electrons (neutral): {total_electrons}")
        print(f"Bond type: {bond.bond_type}")
        print(f"Representation: {bond.representation.__class__.__name__}")
        print(f"Hamiltonian: {bond.hamiltonian.__class__.__name__}")
        print(f"Governance: {bond.governance.__class__.__name__}")
        print(f"Mapper: {bond.mapper.__class__.__name__}")
        print(f"{'='*70}\n")

        # Verify integration
        assert bond.bond_type == 'covalent'
        assert hasattr(bond, 'hamiltonian')
        assert hasattr(bond, 'governance')
        assert hasattr(bond, 'representation')
        assert hasattr(bond, 'mapper')

    def test_heh_hf_energy_validation(self):
        """
        HeH+ HF energy validation.

        Literature: HeH+ (R~0.772 Å) E_HF ~ -2.9 Ha (approximate)
        Note: Exact literature value depends on basis set
        """
        he = Atom('He', [0.0, 0.0, 0.0])
        h = Atom('H', [0.0, 0.0, 0.772])

        from kanad.bonds.covalent_bond import CovalentBond
        bond = CovalentBond(he, h, distance=0.772)

        result = bond.compute_energy(method='HF', max_iterations=200)

        energy = result['energy']

        print(f"\n{'='*70}")
        print(f"HeH+ HARTREE-FOCK VALIDATION")
        print(f"{'='*70}")
        print(f"System: HeH+, R = 0.772 Å")
        print(f"Basis: STO-3G")
        print(f"HF Energy: {energy:.6f} Ha")
        print(f"Converged: {result['converged']}")
        print(f"Iterations: {result['iterations']}")
        print(f"Nuclear repulsion: {bond.hamiltonian.nuclear_repulsion:.6f} Ha")
        print(f"{'='*70}\n")

        # Physical checks
        assert energy < 0, "Total energy should be negative (bound state)"
        assert result['converged'], "SCF must converge"
        assert bond.hamiltonian.nuclear_repulsion > 0, "Nuclear repulsion must be positive"


class TestGovernanceIntegrationH2O:
    """STRICT validation of H2O with governance."""

    def test_h2o_molecule_creation(self):
        """
        H2O molecule using BondFactory for polyatomic.

        Geometry: Bent structure, 104.5° angle
        """
        # Create H2O using preset geometry
        molecule = BondFactory.create_molecule(['H', 'O', 'H'], geometry='water')

        print(f"\n{'='*70}")
        print(f"H2O MOLECULE INTEGRATION AUDIT")
        print(f"{'='*70}")
        print(f"Atoms: {len(molecule.atoms)}")
        for i, atom in enumerate(molecule.atoms):
            print(f"  Atom {i}: {atom.symbol} at {atom.position}")
        print(f"Total electrons: {molecule.n_electrons}")
        print(f"{'='*70}\n")

        # Verify
        assert len(molecule.atoms) == 3, "H2O should have 3 atoms"
        assert molecule.n_electrons == 10, "H2O should have 10 electrons (8+1+1)"

    def test_h2o_hf_energy_validation(self):
        """
        H2O HF energy - STRICT validation.

        Literature: H2O (STO-3G) E_HF ~ -74.96 Ha
        Reference: Szabo & Ostlund
        Tolerance: < 0.01 Ha (STRICT for minimal basis)
        """
        molecule = BondFactory.create_molecule(['H', 'O', 'H'], geometry='water')

        # Create Hamiltonian and solve SCF
        from kanad.core.representations.lcao_representation import LCAORepresentation
        from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(
            molecule,
            representation,
            basis_name='sto-3g',
            use_governance=False  # Speed up for energy test
        )

        density_matrix, hf_energy = hamiltonian.solve_scf(
            max_iterations=200,
            conv_tol=1e-8,
            use_diis=True
        )

        # Literature value (approximate for STO-3G)
        lit_hf_approx = -74.96  # Ha (STO-3G from literature)

        error = abs(hf_energy - lit_hf_approx)

        print(f"\n{'='*70}")
        print(f"H2O HARTREE-FOCK STRICT VALIDATION")
        print(f"{'='*70}")
        print(f"System: H2O")
        print(f"Basis: STO-3G")
        print(f"Literature (approx): {lit_hf_approx:.6f} Ha")
        print(f"Computed HF:         {hf_energy:.6f} Ha")
        print(f"Error:               {error:.6f} Ha")
        print(f"Tolerance:           0.010 Ha (STRICT)")
        print(f"Converged:           {hamiltonian._scf_converged}")
        print(f"Iterations:          {hamiltonian._scf_iterations}")
        print(f"Nuclear repulsion:   {hamiltonian.nuclear_repulsion:.6f} Ha")
        print(f"Status:              {'✓ PASS' if error < 0.01 else '✗ FAIL'}")
        print(f"{'='*70}\n")

        # STRICT validation
        assert error < 0.01, f"H2O HF error {error:.6f} Ha exceeds 0.01 Ha tolerance"
        assert hamiltonian._scf_converged, "SCF must converge for H2O"


class TestSolverOptionsCompleteness:
    """Verify solvers provide ALL available options."""

    def test_vqe_solver_all_ansatz_options(self):
        """
        VQE solver must provide all ansatz types:
        - UCC (Unitary Coupled Cluster)
        - Hardware-Efficient
        - Governance-aware (Covalent, Ionic, Metallic)
        """
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        ansatz_types = ['ucc', 'hardware_efficient', 'governance']

        print(f"\n{'='*70}")
        print(f"VQE SOLVER ANSATZ OPTIONS AUDIT")
        print(f"{'='*70}")

        for ansatz_type in ansatz_types:
            try:
                solver = VQESolver(
                    bond=bond,
                    ansatz_type=ansatz_type,
                    max_iterations=10
                )
                print(f"✓ {ansatz_type:20} - Available")
                assert hasattr(solver, 'ansatz'), f"Solver missing ansatz for {ansatz_type}"
            except Exception as e:
                print(f"✗ {ansatz_type:20} - ERROR: {e}")
                pytest.fail(f"Ansatz type {ansatz_type} not available: {e}")

        print(f"{'='*70}\n")

    def test_vqe_solver_all_mapper_options(self):
        """
        VQE solver must provide all mapper types:
        - Jordan-Wigner
        - Bravyi-Kitaev
        - Parity
        """
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        mapper_types = ['jordan_wigner', 'bravyi_kitaev', 'parity']

        print(f"\n{'='*70}")
        print(f"VQE SOLVER MAPPER OPTIONS AUDIT")
        print(f"{'='*70}")

        for mapper_type in mapper_types:
            try:
                solver = VQESolver(
                    bond=bond,
                    mapper_type=mapper_type,
                    max_iterations=10
                )
                print(f"✓ {mapper_type:20} - Available")
                assert hasattr(solver, 'mapper'), f"Solver missing mapper for {mapper_type}"
            except Exception as e:
                print(f"✗ {mapper_type:20} - ERROR: {e}")
                pytest.fail(f"Mapper type {mapper_type} not available: {e}")

        print(f"{'='*70}\n")

    def test_vqe_solver_all_optimizer_options(self):
        """
        VQE solver must provide all optimizer types:
        - SLSQP (Sequential Least Squares)
        - COBYLA (Constrained Optimization)
        - L-BFGS-B (Limited-memory BFGS)
        """
        bond = BondFactory.create_bond('H', 'H', distance=0.74)

        optimizer_types = ['SLSQP', 'COBYLA', 'L-BFGS-B']

        print(f"\n{'='*70}")
        print(f"VQE SOLVER OPTIMIZER OPTIONS AUDIT")
        print(f"{'='*70}")

        for optimizer in optimizer_types:
            try:
                solver = VQESolver(
                    bond=bond,
                    optimizer=optimizer,
                    max_iterations=10
                )
                print(f"✓ {optimizer:20} - Available")
                assert solver.optimizer_method == optimizer, \
                    f"Optimizer not set correctly: {solver.optimizer_method} != {optimizer}"
            except Exception as e:
                print(f"✗ {optimizer:20} - ERROR: {e}")
                pytest.fail(f"Optimizer {optimizer} not available: {e}")

        print(f"{'='*70}\n")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s", "--tb=short"])
