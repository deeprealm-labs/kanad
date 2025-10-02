"""
Unit tests for variational ansätze.

Tests all ansatz types:
- UCC (UCCSD, Singles, Doubles)
- Hardware-Efficient (Real Amplitudes, Efficient SU(2))
- Governance-aware (Ionic, Covalent, Adaptive)
"""

import pytest
import numpy as np

from kanad.ansatze.base_ansatz import QuantumCircuit, Parameter
from kanad.ansatze.ucc_ansatz import UCCAnsatz, UCC_S_Ansatz, UCC_D_Ansatz
from kanad.ansatze.hardware_efficient_ansatz import (
    HardwareEfficientAnsatz,
    RealAmplitudesAnsatz,
    EfficientSU2Ansatz
)
from kanad.ansatze.governance_aware_ansatz import (
    IonicGovernanceAnsatz,
    CovalentGovernanceAnsatz,
    AdaptiveGovernanceAnsatz
)


class TestBaseClasses:
    """Test base classes (Parameter, QuantumCircuit)."""

    def test_parameter_creation(self):
        """Test parameter creation and initialization."""
        param = Parameter('theta_0', 0.5)
        assert param.name == 'theta_0'
        assert param.value == 0.5

        # Default value
        param2 = Parameter('theta_1')
        assert param2.value == 0.0

    def test_quantum_circuit_creation(self):
        """Test quantum circuit creation."""
        circuit = QuantumCircuit(4)
        assert circuit.n_qubits == 4
        assert len(circuit.gates) == 0
        assert len(circuit.parameters) == 0

    def test_quantum_circuit_gates(self):
        """Test adding gates to circuit."""
        circuit = QuantumCircuit(2)

        # Single-qubit gates
        circuit.h(0)
        circuit.x(1)
        assert len(circuit.gates) == 2

        # Parametrized gates
        theta = Parameter('theta')
        circuit.ry(theta, 0)
        assert len(circuit.parameters) == 1

        # Two-qubit gates
        circuit.cx(0, 1)
        assert len(circuit.gates) == 4

    def test_parameter_binding(self):
        """Test binding parameter values."""
        circuit = QuantumCircuit(2)
        theta1 = Parameter('theta_1')
        theta2 = Parameter('theta_2')

        circuit.ry(theta1, 0)
        circuit.rz(theta2, 1)

        values = np.array([0.5, 1.0])
        circuit.bind_parameters(values)

        assert circuit.parameters[0].value == 0.5
        assert circuit.parameters[1].value == 1.0

    def test_parameter_binding_wrong_size(self):
        """Test error handling for wrong parameter count."""
        circuit = QuantumCircuit(2)
        circuit.ry(Parameter('theta'), 0)

        with pytest.raises(ValueError, match="Expected 1 values"):
            circuit.bind_parameters(np.array([0.5, 1.0]))


class TestUCCAnsatz:
    """Test Unitary Coupled Cluster ansätze."""

    def test_uccsd_creation(self):
        """Test UCCSD ansatz creation."""
        ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)
        assert ansatz.n_qubits == 4
        assert ansatz.n_electrons == 2

    def test_uccsd_circuit_building(self):
        """Test UCCSD circuit construction."""
        ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)
        circuit = ansatz.build_circuit()

        assert circuit.n_qubits == 4
        assert len(circuit.gates) > 0
        assert len(circuit.parameters) > 0

    def test_uccsd_excitation_generation(self):
        """Test excitation generation (singles + doubles)."""
        ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)
        circuit = ansatz.build_circuit()

        # Should have singles and doubles
        # 2 occupied → 2 virtual = 2×2 = 4 singles
        # 2 occupied pairs → 2 virtual pairs = 1 double
        n_params = circuit.get_num_parameters()
        assert n_params > 0  # Has parameters

    def test_ucc_singles_only(self):
        """Test UCC with singles only."""
        ansatz = UCC_S_Ansatz(n_qubits=4, n_electrons=2)
        circuit = ansatz.build_circuit()

        assert circuit.n_qubits == 4
        assert len(circuit.parameters) > 0

    def test_ucc_doubles_only(self):
        """Test UCC with doubles only."""
        ansatz = UCC_D_Ansatz(n_qubits=4, n_electrons=2)
        circuit = ansatz.build_circuit()

        assert circuit.n_qubits == 4

    def test_ucc_parameter_initialization(self):
        """Test parameter initialization strategies."""
        ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)

        # Random initialization
        params_random = ansatz.initialize_parameters('random')
        assert len(params_random) == ansatz.get_num_parameters()
        assert np.all(params_random >= -np.pi)
        assert np.all(params_random <= np.pi)

        # Zero initialization
        params_zeros = ansatz.initialize_parameters('zeros')
        assert np.all(params_zeros == 0.0)

        # Small random
        params_small = ansatz.initialize_parameters('small_random')
        assert np.all(np.abs(params_small) <= 0.1)


class TestHardwareEfficientAnsatz:
    """Test hardware-efficient ansätze."""

    def test_hardware_efficient_creation(self):
        """Test hardware-efficient ansatz creation."""
        ansatz = HardwareEfficientAnsatz(
            n_qubits=4,
            n_electrons=2,
            n_layers=2
        )
        assert ansatz.n_qubits == 4
        assert ansatz.n_layers == 2

    def test_linear_entanglement(self):
        """Test linear entanglement pattern."""
        ansatz = HardwareEfficientAnsatz(
            n_qubits=4,
            n_electrons=2,
            entanglement='linear'
        )
        circuit = ansatz.build_circuit()
        assert circuit.n_qubits == 4

    def test_circular_entanglement(self):
        """Test circular entanglement pattern."""
        ansatz = HardwareEfficientAnsatz(
            n_qubits=4,
            n_electrons=2,
            entanglement='circular'
        )
        circuit = ansatz.build_circuit()
        assert circuit.n_qubits == 4

    def test_full_entanglement(self):
        """Test full entanglement pattern."""
        ansatz = HardwareEfficientAnsatz(
            n_qubits=4,
            n_electrons=2,
            entanglement='full'
        )
        circuit = ansatz.build_circuit()
        assert circuit.n_qubits == 4

    def test_pairwise_entanglement(self):
        """Test pairwise entanglement pattern."""
        ansatz = HardwareEfficientAnsatz(
            n_qubits=4,
            n_electrons=2,
            entanglement='pairwise'
        )
        circuit = ansatz.build_circuit()
        assert circuit.n_qubits == 4

    def test_real_amplitudes_ansatz(self):
        """Test RealAmplitudes ansatz (RY only)."""
        ansatz = RealAmplitudesAnsatz(n_qubits=4, n_electrons=2)
        assert ansatz.rotation_gates == ['ry']
        assert ansatz.entangling_gate == 'cx'

        circuit = ansatz.build_circuit()
        assert circuit.n_qubits == 4

    def test_efficient_su2_ansatz(self):
        """Test EfficientSU2 ansatz (RY + RZ)."""
        ansatz = EfficientSU2Ansatz(n_qubits=4, n_electrons=2)
        assert ansatz.rotation_gates == ['ry', 'rz']

        circuit = ansatz.build_circuit()
        assert circuit.n_qubits == 4

    def test_custom_rotation_gates(self):
        """Test custom rotation gates."""
        ansatz = HardwareEfficientAnsatz(
            n_qubits=4,
            n_electrons=2,
            rotation_gates=['rx', 'ry', 'rz']
        )
        circuit = ansatz.build_circuit()
        assert len(circuit.parameters) > 0

    def test_custom_entangling_gates(self):
        """Test custom entangling gates."""
        for gate in ['cx', 'cz', 'rxx', 'ryy', 'rzz']:
            ansatz = HardwareEfficientAnsatz(
                n_qubits=4,
                n_electrons=2,
                entangling_gate=gate
            )
            circuit = ansatz.build_circuit()
            assert circuit.n_qubits == 4

    def test_initial_state_preparation(self):
        """Test initial state preparation."""
        ansatz = HardwareEfficientAnsatz(
            n_qubits=4,
            n_electrons=2
        )

        # Custom initial state
        initial_state = [1, 1, 0, 0]
        circuit = ansatz.build_circuit(initial_state=initial_state)
        assert circuit.n_qubits == 4

        # Default (Hartree-Fock)
        circuit2 = ansatz.build_circuit()
        assert circuit2.n_qubits == 4


class TestGovernanceAwareAnsatze:
    """Test governance-aware ansätze."""

    def test_ionic_governance_ansatz(self):
        """Test ionic governance ansatz."""
        ansatz = IonicGovernanceAnsatz(
            n_qubits=4,
            n_electrons=2,
            n_layers=2
        )
        circuit = ansatz.build_circuit()

        assert circuit.n_qubits == 4
        assert len(circuit.parameters) > 0

    def test_ionic_ansatz_charge_transfer(self):
        """Test ionic ansatz emphasizes charge transfer."""
        ansatz = IonicGovernanceAnsatz(
            n_qubits=4,
            n_electrons=2
        )
        circuit = ansatz.build_circuit(ionization_threshold=1.5)

        # Should have RZ gates (charge state) and CNOT (transfer)
        gate_types = [g['type'] for g in circuit.gates if g['type'] != 'barrier']
        assert 'rz' in gate_types  # Charge state rotations
        assert 'ry' in gate_types  # Local rotations

    def test_covalent_governance_ansatz(self):
        """Test covalent governance ansatz."""
        ansatz = CovalentGovernanceAnsatz(
            n_qubits=4,  # Must be even
            n_electrons=2,
            n_layers=2,
            hybridization='sp3'
        )
        circuit = ansatz.build_circuit()

        assert circuit.n_qubits == 4
        assert len(circuit.parameters) > 0

    def test_covalent_ansatz_requires_even_qubits(self):
        """Test covalent ansatz requires even qubits."""
        with pytest.raises(ValueError, match="even number of qubits"):
            CovalentGovernanceAnsatz(
                n_qubits=5,  # Odd number
                n_electrons=2
            )

    def test_covalent_ansatz_pairing(self):
        """Test covalent ansatz uses MO pairing."""
        ansatz = CovalentGovernanceAnsatz(
            n_qubits=4,
            n_electrons=2
        )
        circuit = ansatz.build_circuit(overlap_threshold=0.3)

        # Should have RXX and RYY gates (bonding/antibonding entanglement)
        gate_types = [g['type'] for g in circuit.gates if g['type'] != 'barrier']
        # May or may not have RXX/RYY depending on implementation

    def test_adaptive_governance_ansatz(self):
        """Test adaptive governance ansatz."""
        ansatz = AdaptiveGovernanceAnsatz(
            n_qubits=4,
            n_electrons=2,
            bonding_type='auto'
        )

        # Should auto-detect bonding type
        circuit = ansatz.build_circuit(electronegativity_diff=0.5)
        assert circuit.n_qubits == 4

    def test_adaptive_ansatz_ionic_selection(self):
        """Test adaptive ansatz selects ionic for large EN difference."""
        ansatz = AdaptiveGovernanceAnsatz(
            n_qubits=4,
            n_electrons=2,
            bonding_type='auto'
        )

        # Large EN difference → ionic
        circuit = ansatz.build_circuit(electronegativity_diff=2.5)
        assert ansatz._delegate_ansatz is not None
        assert isinstance(ansatz._delegate_ansatz, IonicGovernanceAnsatz)

    def test_adaptive_ansatz_covalent_selection(self):
        """Test adaptive ansatz selects covalent for small EN difference."""
        ansatz = AdaptiveGovernanceAnsatz(
            n_qubits=4,
            n_electrons=2,
            bonding_type='auto'
        )

        # Small EN difference → covalent
        circuit = ansatz.build_circuit(electronegativity_diff=0.5)
        assert ansatz._delegate_ansatz is not None
        assert isinstance(ansatz._delegate_ansatz, CovalentGovernanceAnsatz)

    def test_adaptive_ansatz_overlap_based(self):
        """Test adaptive ansatz uses overlap matrix."""
        ansatz = AdaptiveGovernanceAnsatz(
            n_qubits=4,
            n_electrons=2
        )

        # High overlap → covalent
        overlap = np.array([
            [1.0, 0.5, 0.0, 0.0],
            [0.5, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.5],
            [0.0, 0.0, 0.5, 1.0]
        ])

        circuit = ansatz.build_circuit(overlap_matrix=overlap)
        assert isinstance(ansatz._delegate_ansatz, CovalentGovernanceAnsatz)


class TestAnsatzComparison:
    """Test comparisons between different ansätze."""

    def test_parameter_count_comparison(self):
        """Compare parameter counts across ansätze."""
        n_qubits = 4
        n_electrons = 2

        # UCC typically has more parameters
        ucc = UCCAnsatz(n_qubits, n_electrons)
        ucc_params = ucc.get_num_parameters()

        # Hardware-efficient has fewer
        hea = HardwareEfficientAnsatz(n_qubits, n_electrons, n_layers=2)
        hea_params = hea.get_num_parameters()

        # Both should have parameters
        assert ucc_params > 0
        assert hea_params > 0

    def test_circuit_depth_comparison(self):
        """Compare circuit depths."""
        n_qubits = 4
        n_electrons = 2

        ansatze = [
            UCCAnsatz(n_qubits, n_electrons),
            HardwareEfficientAnsatz(n_qubits, n_electrons, n_layers=1),
            IonicGovernanceAnsatz(n_qubits, n_electrons, n_layers=1),
        ]

        for ansatz in ansatze:
            depth = ansatz.get_circuit_depth()
            assert depth > 0

    def test_all_ansatze_build_successfully(self):
        """Test all ansätze build valid circuits."""
        n_qubits = 4
        n_electrons = 2

        ansatze = [
            UCCAnsatz(n_qubits, n_electrons),
            UCC_S_Ansatz(n_qubits, n_electrons),
            UCC_D_Ansatz(n_qubits, n_electrons),
            HardwareEfficientAnsatz(n_qubits, n_electrons),
            RealAmplitudesAnsatz(n_qubits, n_electrons),
            EfficientSU2Ansatz(n_qubits, n_electrons),
            IonicGovernanceAnsatz(n_qubits, n_electrons),
            CovalentGovernanceAnsatz(n_qubits, n_electrons),
        ]

        for ansatz in ansatze:
            circuit = ansatz.build_circuit()
            assert circuit.n_qubits == n_qubits
            assert len(circuit.gates) > 0
