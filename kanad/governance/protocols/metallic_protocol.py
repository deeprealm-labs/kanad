"""
Metallic Governance Protocol.

Enforces physics rules for metallic bonding:
- Delocalization validation
- Fermi surface validation
- GHZ-like collective entanglement
- Periodic boundary conditions
"""

from typing import Dict, Any, List, Optional
import numpy as np

from kanad.governance.protocols.base_protocol import BaseGovernanceProtocol


class MetallicGovernanceProtocol(BaseGovernanceProtocol):
    """
    Governance protocol for metallic bonds.

    Physical Rules:
        1. Electrons must be delocalized across lattice
        2. Fermi surface must exist (DOS(E_F) > 0)
        3. Wavefunction should show GHZ-like entanglement
        4. Band structure must show continuous bands
        5. No localized charge accumulation

    Ansatz Strategy:
        - GHZ-like states for delocalized electrons
        - Collective entanglement (all qubits entangled)
        - Parameterized rotations preserving translational symmetry
    """

    def __init__(self):
        """Initialize metallic governance protocol."""
        super().__init__(bond_type='metallic')
        self._initialize_rules()

    def _initialize_rules(self):
        """Initialize metallic bonding rules."""
        self.rules = {
            'delocalization': {
                'min_score': 0.3,
                'description': 'Electrons must be delocalized across lattice'
            },
            'fermi_surface': {
                'dos_at_fermi': '> 0',
                'description': 'Metallic character requires DOS(E_F) > 0'
            },
            'periodic_bc': {
                'min_atoms': 3,
                'description': 'Periodic BC require >= 3 atoms'
            },
            'mott_transition': {
                'max_u_over_t': 4.0,
                'description': 'Large U/t may cause insulator transition'
            }
        }

    def validate_operator(self, operator) -> bool:
        """
        Validate if operator is appropriate for metallic system.

        For metallic bonds, operators should preserve translational symmetry
        and not localize charge.

        Args:
            operator: Quantum operator

        Returns:
            True if valid
        """
        # For now, accept all operators
        # Could add checks for:
        # - Translational symmetry
        # - Momentum conservation
        # - No charge localization
        return True

    def enforce_constraints(self, circuit):
        """
        Apply constraints to circuit for metallic system.

        Ensures circuit respects metallic bonding physics.

        Args:
            circuit: Quantum circuit

        Returns:
            Constrained circuit
        """
        # For metallic systems, ensure GHZ-like entanglement
        # Already handled in construct_ansatz
        return circuit

    def validate_physical_constraints(
        self,
        hamiltonian,
        wavefunction: Optional[np.ndarray] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Validate that system satisfies metallic bonding physics.

        Args:
            hamiltonian: MetallicHamiltonian instance
            wavefunction: Optional quantum state
            **kwargs: Additional parameters

        Returns:
            Validation results dictionary
        """
        results = {
            'valid': True,
            'violations': [],
            'warnings': []
        }

        # 1. Check metallic character (DOS at Fermi level)
        if hasattr(hamiltonian, 'is_metallic'):
            is_metallic = hamiltonian.is_metallic()
            if not is_metallic:
                results['violations'].append("System is not metallic (DOS(E_F) = 0)")
                results['valid'] = False

        # 2. Validate delocalization
        if wavefunction is not None:
            delocalization_score = self._compute_delocalization(wavefunction, hamiltonian.n_sites)

            # For metallic systems, electrons should be spread across sites
            # Delocalization score should be > 0.3 (empirical threshold)
            if delocalization_score < 0.3:
                results['warnings'].append(
                    f"Low delocalization score: {delocalization_score:.3f} (expect > 0.3)"
                )

            results['delocalization_score'] = delocalization_score

        # 3. Check periodic boundary conditions (if applicable)
        if hasattr(hamiltonian, 'periodic') and hamiltonian.periodic:
            if hamiltonian.n_sites < 3:
                results['warnings'].append(
                    "Periodic BC with < 3 sites may cause artifacts"
                )

        # 4. Validate number of electrons
        if hamiltonian.n_electrons < hamiltonian.n_sites:
            results['warnings'].append(
                f"Partially filled bands: {hamiltonian.n_electrons} electrons, "
                f"{hamiltonian.n_sites} sites (good for metallic behavior)"
            )

        # 5. Check for localized states (strong Hubbard U)
        if hasattr(hamiltonian, 'hubbard_u'):
            if hamiltonian.hubbard_u > abs(4 * hamiltonian.hopping_parameter):
                results['warnings'].append(
                    f"Large U/t ratio ({hamiltonian.hubbard_u / abs(hamiltonian.hopping_parameter):.2f}) "
                    "may cause Mott insulator transition"
                )

        return results

    def _compute_delocalization(self, wavefunction: np.ndarray, n_sites: int) -> float:
        """
        Compute delocalization score from wavefunction.

        Measures how spread out the wavefunction is across sites.
        Higher score = more delocalized (metallic).

        Args:
            wavefunction: Quantum state vector
            n_sites: Number of lattice sites

        Returns:
            Delocalization score (0 = fully localized, 1 = fully delocalized)
        """
        # Compute participation ratio: PR = 1 / Σ_i |ψ_i|^4
        # PR = 1 for fully localized, PR = N for fully delocalized

        probabilities = np.abs(wavefunction) ** 2
        participation_ratio = 1.0 / np.sum(probabilities ** 2)

        # Normalize to [0, 1]: score = (PR - 1) / (N - 1)
        n_basis = len(wavefunction)
        if n_basis > 1:
            score = (participation_ratio - 1.0) / (n_basis - 1.0)
        else:
            score = 0.0

        return min(score, 1.0)  # Cap at 1.0

    def construct_ansatz(
        self,
        n_qubits: int,
        n_layers: int = 1,
        entanglement_type: str = 'ghz'
    ) -> 'QuantumCircuit':
        """
        Construct ansatz for metallic system.

        Creates GHZ-like states for collective delocalization.

        Args:
            n_qubits: Number of qubits (2 * n_sites for spin up/down)
            n_layers: Number of variational layers
            entanglement_type: 'ghz' or 'chain' or 'all-to-all'

        Returns:
            Parameterized quantum circuit
        """
        try:
            from qiskit import QuantumCircuit
            from qiskit.circuit import Parameter
        except ImportError:
            raise ImportError("Qiskit required for circuit construction")

        qc = QuantumCircuit(n_qubits)

        param_index = 0

        for layer in range(n_layers):
            # 1. Single-qubit rotations (local chemistry)
            for q in range(n_qubits):
                theta = Parameter(f'θ_{layer}_{q}')
                phi = Parameter(f'φ_{layer}_{q}')
                qc.ry(theta, q)
                qc.rz(phi, q)
                param_index += 2

            # 2. Entanglement layer
            if entanglement_type == 'ghz':
                # GHZ-like entanglement: H on first qubit, then CNOTs
                if layer == 0:
                    qc.h(0)
                for q in range(n_qubits - 1):
                    qc.cx(q, q + 1)

            elif entanglement_type == 'chain':
                # Linear chain entanglement
                for q in range(n_qubits - 1):
                    qc.cx(q, q + 1)

            elif entanglement_type == 'all-to-all':
                # Full entanglement (expensive but thorough)
                for q1 in range(n_qubits):
                    for q2 in range(q1 + 1, n_qubits):
                        qc.cx(q1, q2)

            # 3. Periodic boundary (optional)
            if n_qubits > 2 and entanglement_type in ['chain', 'ghz']:
                qc.cx(n_qubits - 1, 0)

        return qc

    def validate_convergence(
        self,
        energies: List[float],
        tolerance: float = 1e-6,
        min_iterations: int = 10
    ) -> Dict[str, Any]:
        """
        Validate VQE convergence for metallic system.

        Args:
            energies: Energy history
            tolerance: Convergence tolerance
            min_iterations: Minimum iterations required

        Returns:
            Convergence analysis
        """
        if len(energies) < min_iterations:
            return {
                'converged': False,
                'reason': f'Need at least {min_iterations} iterations'
            }

        # Check last 5 iterations
        recent = energies[-5:]
        if len(recent) < 5:
            return {'converged': False, 'reason': 'Not enough iterations'}

        energy_std = np.std(recent)
        energy_range = max(recent) - min(recent)

        converged = energy_std < tolerance and energy_range < tolerance

        return {
            'converged': converged,
            'final_energy': energies[-1],
            'energy_std': energy_std,
            'energy_range': energy_range,
            'iterations': len(energies)
        }

    def suggest_parameters(self, hamiltonian) -> Dict[str, Any]:
        """
        Suggest optimal parameters for metallic system VQE.

        Args:
            hamiltonian: MetallicHamiltonian instance

        Returns:
            Suggested parameters
        """
        n_sites = hamiltonian.n_sites
        n_qubits = 2 * n_sites  # spin up + spin down

        # More qubits → more layers needed
        suggested_layers = max(2, int(np.log2(n_sites)) + 1)

        # Choose entanglement based on system size
        if n_sites <= 4:
            entanglement = 'all-to-all'
        elif n_sites <= 8:
            entanglement = 'ghz'
        else:
            entanglement = 'chain'

        return {
            'n_layers': suggested_layers,
            'entanglement_type': entanglement,
            'optimizer': 'COBYLA',  # Good for noisy landscapes
            'max_iterations': 500,
            'conv_tolerance': 1e-6,
            'n_qubits': n_qubits
        }

    def validate_result(
        self,
        computed_energy: float,
        hamiltonian,
        reference_energy: Optional[float] = None
    ) -> Dict[str, Any]:
        """
        Validate computed energy result.

        Args:
            computed_energy: VQE or computed energy
            hamiltonian: MetallicHamiltonian
            reference_energy: Optional reference (e.g., from tight-binding)

        Returns:
            Validation results
        """
        results = {
            'valid': True,
            'energy': computed_energy,
            'warnings': []
        }

        # Check if energy is reasonable
        # For metallic systems with hopping t, typical energies ~ t * n_sites
        expected_scale = abs(hamiltonian.hopping_parameter) * hamiltonian.n_sites

        if abs(computed_energy) > 10 * expected_scale:
            results['warnings'].append(
                f"Energy seems too large: {computed_energy:.3f} eV "
                f"(expected scale: {expected_scale:.3f} eV)"
            )

        # Compare with reference if available
        if reference_energy is not None:
            error = abs(computed_energy - reference_energy)
            relative_error = error / abs(reference_energy) if reference_energy != 0 else np.inf

            results['reference_energy'] = reference_energy
            results['absolute_error'] = error
            results['relative_error'] = relative_error

            if relative_error > 0.1:  # 10% error
                results['warnings'].append(
                    f"Large deviation from reference: {relative_error*100:.1f}%"
                )

        return results

    def __repr__(self) -> str:
        """String representation."""
        return "MetallicGovernanceProtocol(delocalized, GHZ-like)"
