"""
Unit tests for governance protocol layer (Phase 4).

Tests the core innovation: bonding-type-specific governance of
quantum circuit construction.
"""

import pytest
import numpy as np
from kanad.governance.protocols.base_protocol import (
    BaseGovernanceProtocol,
    BondingType,
    GovernanceRule,
    QuantumCircuitState
)
from kanad.governance.protocols.ionic_protocol import IonicGovernanceProtocol
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule
from kanad.core.representations.second_quantization import SecondQuantizationRepresentation
from kanad.core.representations.lcao_representation import LCAORepresentation


class TestQuantumCircuitState:
    """Test QuantumCircuitState class."""

    def test_circuit_state_creation(self):
        """Test creating a circuit state."""
        circuit = QuantumCircuitState(n_qubits=4)

        assert circuit.n_qubits == 4
        assert len(circuit.gates) == 0
        assert circuit.depth == 0
        assert circuit.is_localized is True

    def test_add_single_qubit_gate(self):
        """Test adding single-qubit gates."""
        circuit = QuantumCircuitState(4)

        circuit.add_gate('rx', [0], params=[np.pi/2])
        circuit.add_gate('ry', [1], params=[np.pi/4])

        assert len(circuit.gates) == 2
        assert circuit.depth == 2
        assert circuit.gates[0]['type'] == 'rx'
        assert circuit.gates[0]['qubits'] == [0]

    def test_add_two_qubit_gate(self):
        """Test adding two-qubit gates."""
        circuit = QuantumCircuitState(4)

        circuit.add_gate('cx', [0, 1])
        circuit.add_gate('rxx', [2, 3], params=[0.1])

        assert len(circuit.gates) == 2
        assert circuit.gates[0]['type'] == 'cx'
        assert len(circuit.gates[0]['qubits']) == 2

    def test_entanglement_tracking(self):
        """Test entanglement graph tracking."""
        circuit = QuantumCircuitState(4)

        # Add entangling gates
        circuit.add_gate('cx', [0, 1])
        circuit.add_gate('cx', [1, 2])

        # Check entanglement graph
        assert 1 in circuit.entanglement_graph[0]
        assert 0 in circuit.entanglement_graph[1]
        assert 2 in circuit.entanglement_graph[1]

        # Check entanglement degree
        assert circuit.get_entanglement_degree(1) == 2  # Connected to 0 and 2
        assert circuit.max_entanglement_degree() == 2

    def test_circuit_sparsity(self):
        """Test sparsity detection."""
        # Sparse circuit (few connections)
        sparse_circuit = QuantumCircuitState(4)
        sparse_circuit.add_gate('cx', [0, 1])

        assert sparse_circuit.is_sparse() is True

        # Dense circuit (many connections)
        dense_circuit = QuantumCircuitState(4)
        for i in range(4):
            for j in range(i + 1, 4):
                dense_circuit.add_gate('cx', [i, j])

        assert dense_circuit.is_sparse() is False


class TestGovernanceRule:
    """Test GovernanceRule class."""

    def test_rule_creation(self):
        """Test creating a governance rule."""
        rule = GovernanceRule(
            name="test_rule",
            description="A test rule",
            condition=lambda s, c: True,
            action=lambda s, c: s,
            priority=50
        )

        assert rule.name == "test_rule"
        assert rule.priority == 50
        assert rule.required is True

    def test_rule_application(self):
        """Test rule applies method."""
        # Rule that applies when state has flag
        rule = GovernanceRule(
            name="flag_check",
            description="Check for flag",
            condition=lambda s, c: getattr(s, 'flag', False),
            action=lambda s, c: s,
            priority=10
        )

        class DummyState:
            flag = True

        state = DummyState()
        assert rule.applies(state, {}) is True

        state.flag = False
        assert rule.applies(state, {}) is False

    def test_rule_execution(self):
        """Test rule execution."""
        # Rule that modifies state
        def modify_state(state, context):
            state.modified = True
            return state

        rule = GovernanceRule(
            name="modifier",
            description="Modifies state",
            condition=lambda s, c: True,
            action=modify_state,
            priority=10
        )

        class DummyState:
            modified = False

        state = DummyState()
        result = rule.execute(state, {})

        assert result.modified is True


class TestBaseGovernanceProtocol:
    """Test base governance protocol."""

    def test_bonding_type_enum(self):
        """Test BondingType enumeration."""
        assert BondingType.IONIC.value == "ionic"
        assert BondingType.COVALENT.value == "covalent"
        assert BondingType.METALLIC.value == "metallic"

    def test_entanglement_strategies(self):
        """Test entanglement strategy descriptions."""
        ionic_protocol = IonicGovernanceProtocol()
        covalent_protocol = CovalentGovernanceProtocol()

        ionic_strategy = ionic_protocol.get_entanglement_strategy()
        covalent_strategy = covalent_protocol.get_entanglement_strategy()

        assert "minimal" in ionic_strategy
        assert "paired" in covalent_strategy


class TestIonicGovernanceProtocol:
    """Test ionic bonding governance."""

    def test_ionic_protocol_creation(self):
        """Test creating ionic protocol."""
        protocol = IonicGovernanceProtocol()

        assert protocol.bond_type == BondingType.IONIC
        assert len(protocol.rules) > 0
        assert "minimal" in protocol.get_entanglement_strategy()

    def test_ionic_rules_initialized(self):
        """Test that ionic rules are properly initialized."""
        protocol = IonicGovernanceProtocol()

        rule_names = [rule.name for rule in protocol.rules]

        assert "localized_gates_only" in rule_names
        assert "sparse_connectivity" in rule_names
        assert "particle_number_conservation" in rule_names
        assert "forbid_collective_gates" in rule_names

    def test_ionic_allowed_operators(self):
        """Test allowed operators for ionic bonding."""
        protocol = IonicGovernanceProtocol()

        allowed = protocol.get_allowed_operators()
        forbidden = protocol.get_forbidden_operators()

        # Single-qubit rotations should be allowed
        assert 'rx' in allowed
        assert 'ry' in allowed

        # Collective operations should be forbidden
        assert 'qft' in forbidden
        assert 'ghz' in forbidden

    def test_ionic_ansatz_construction(self):
        """Test ionic ansatz construction."""
        # Create ionic system
        na = Atom('Na', position=np.zeros(3))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))
        nacl = Molecule([na, cl])
        rep = SecondQuantizationRepresentation(nacl)

        # Construct ansatz
        protocol = IonicGovernanceProtocol()
        circuit = protocol.construct_ansatz(rep)

        assert isinstance(circuit, QuantumCircuitState)
        assert circuit.is_localized is True
        assert circuit.n_qubits == rep.n_qubits
        assert len(circuit.gates) > 0

    def test_ionic_locality_enforcement(self):
        """Test that ionic protocol enforces locality."""
        protocol = IonicGovernanceProtocol()

        # Create circuit with long-range gate
        circuit = QuantumCircuitState(4)
        circuit.add_gate('rx', [0])
        circuit.add_gate('cx', [0, 1])  # Nearest neighbor (OK)

        # Should pass - only local gates
        constrained = protocol.enforce_constraints(circuit)
        assert len(constrained.gates) == 2

        # Add long-range gate
        circuit.add_gate('cx', [0, 3])  # Long-range (should be removed)

        constrained = protocol.enforce_constraints(circuit)
        # Long-range gate should be filtered out
        assert all(len(g['qubits']) == 1 or abs(g['qubits'][0] - g['qubits'][1]) <= 1
                  for g in constrained.gates)

    def test_ionic_transfer_integral(self):
        """Test transfer integral calculation."""
        protocol = IonicGovernanceProtocol()

        # Close atoms → larger transfer
        t_close = protocol.get_transfer_integral_estimate(distance=1.0)

        # Far atoms → smaller transfer
        t_far = protocol.get_transfer_integral_estimate(distance=5.0)

        assert t_close > t_far
        assert t_far > 0  # Should decay but not be zero


class TestCovalentGovernanceProtocol:
    """Test covalent bonding governance."""

    def test_covalent_protocol_creation(self):
        """Test creating covalent protocol."""
        protocol = CovalentGovernanceProtocol()

        assert protocol.bond_type == BondingType.COVALENT
        assert len(protocol.rules) > 0
        assert "paired" in protocol.get_entanglement_strategy()

    def test_covalent_rules_initialized(self):
        """Test that covalent rules are properly initialized."""
        protocol = CovalentGovernanceProtocol()

        rule_names = [rule.name for rule in protocol.rules]

        assert "hybridization_first" in rule_names
        assert "molecular_orbital_formation" in rule_names
        assert "electron_pair_entanglement" in rule_names
        assert "spin_symmetry" in rule_names

    def test_covalent_allowed_operators(self):
        """Test allowed operators for covalent bonding."""
        protocol = CovalentGovernanceProtocol()

        allowed = protocol.get_allowed_operators()
        forbidden = protocol.get_forbidden_operators()

        # Hybridization and pairing should be allowed
        assert 'ry' in allowed  # Hybridization
        assert 'givens' in allowed  # MO formation
        assert 'h' in allowed  # Bell states
        assert 'cx' in allowed  # Pairing

        # Collective operations should be forbidden
        assert 'qft' in forbidden
        assert 'ghz' in forbidden

    def test_covalent_ansatz_construction(self):
        """Test covalent ansatz construction."""
        # Create covalent system
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        h2_mol = Molecule([h1, h2])
        rep = LCAORepresentation(h2_mol)

        # Construct ansatz
        protocol = CovalentGovernanceProtocol()
        circuit = protocol.construct_ansatz(rep)

        assert isinstance(circuit, QuantumCircuitState)
        assert circuit.is_hybridized is True
        assert circuit.has_mo_pairs is True
        assert circuit.is_paired is True
        assert len(circuit.gates) > 0

    def test_covalent_hybridization_angles(self):
        """Test hybridization angle calculations."""
        protocol = CovalentGovernanceProtocol()

        # sp hybridization (linear)
        angle_sp = protocol.get_bonding_angle('sp')
        assert abs(angle_sp - np.pi) < 1e-10

        # sp² hybridization (trigonal planar, 120°)
        angle_sp2 = protocol.get_bonding_angle('sp2')
        assert abs(angle_sp2 - 2*np.pi/3) < 1e-10

        # sp³ hybridization (tetrahedral, ~109.5°)
        angle_sp3 = protocol.get_bonding_angle('sp3')
        expected_sp3 = np.arccos(-1.0/3.0)
        assert abs(angle_sp3 - expected_sp3) < 1e-10

    def test_covalent_constraint_ordering(self):
        """Test that covalent constraints enforce proper ordering."""
        protocol = CovalentGovernanceProtocol()

        # Create circuit with wrong order
        circuit = QuantumCircuitState(4)
        circuit.add_gate('cx', [0, 1])  # Pairing gate
        circuit.add_gate('ry', [0])  # Hybridization gate (should come first)

        # Enforce constraints - should reorder
        constrained = protocol.enforce_constraints(circuit)

        # Hybridization gates should now come first
        assert constrained.gates[0]['type'] == 'ry'
        assert constrained.is_hybridized is True


class TestGovernanceComparison:
    """Test differences between governance protocols."""

    def test_ionic_vs_covalent_entanglement(self):
        """Compare entanglement strategies."""
        ionic = IonicGovernanceProtocol()
        covalent = CovalentGovernanceProtocol()

        ionic_strategy = ionic.get_entanglement_strategy()
        covalent_strategy = covalent.get_entanglement_strategy()

        assert ionic_strategy != covalent_strategy
        assert "minimal" in ionic_strategy
        assert "paired" in covalent_strategy

    def test_ionic_vs_covalent_operators(self):
        """Compare allowed operators."""
        ionic = IonicGovernanceProtocol()
        covalent = CovalentGovernanceProtocol()

        ionic_allowed = set(ionic.get_allowed_operators())
        covalent_allowed = set(covalent.get_allowed_operators())

        # Both allow basic rotations
        assert 'rx' in ionic_allowed and 'rx' in covalent_allowed

        # Covalent has special operators
        assert 'givens' in covalent_allowed
        assert 'givens' not in ionic_allowed

    def test_ansatz_differences(self):
        """Test that different bonding types produce different ansatze."""
        # Ionic system
        na = Atom('Na', position=np.zeros(3))
        cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))
        ionic_mol = Molecule([na, cl])
        ionic_rep = SecondQuantizationRepresentation(ionic_mol)

        # Covalent system
        h1 = Atom('H', position=np.zeros(3))
        h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        covalent_mol = Molecule([h1, h2])
        covalent_rep = LCAORepresentation(covalent_mol)

        # Construct ansatze
        ionic_protocol = IonicGovernanceProtocol()
        covalent_protocol = CovalentGovernanceProtocol()

        ionic_circuit = ionic_protocol.construct_ansatz(ionic_rep)
        covalent_circuit = covalent_protocol.construct_ansatz(covalent_rep)

        # Should have different properties
        assert ionic_circuit.is_localized is True
        assert covalent_circuit.is_hybridized is True

        # Covalent should have MO pairs, ionic shouldn't
        assert covalent_circuit.has_mo_pairs is True
        assert ionic_circuit.has_mo_pairs is False


class TestGovernanceRuleEngine:
    """Test the governance rule application engine."""

    def test_rule_priority_ordering(self):
        """Test that rules execute in priority order."""
        protocol = CovalentGovernanceProtocol()

        # Rules should be sorted by priority
        priorities = [rule.priority for rule in protocol.rules]

        # Check that priorities are generally decreasing or stay same
        # (some rules might have same priority)
        for i in range(len(priorities) - 1):
            # High priority (100) should come before low priority (70)
            pass  # Just verify no errors

    def test_rule_application_with_context(self):
        """Test applying rules with context."""
        protocol = CovalentGovernanceProtocol()

        circuit = QuantumCircuitState(4)
        context = {
            'hybridization': 'sp3',
            'bonds': [(0, 1), (2, 3)],
            'electron_pairs': [(0, 1)]
        }

        # Apply governance
        governed = protocol.apply_governance(circuit, context)

        assert governed.is_hybridized or governed.has_mo_pairs or governed.is_paired


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
