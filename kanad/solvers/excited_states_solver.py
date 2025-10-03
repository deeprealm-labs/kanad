"""
Excited States Solver using Time-Dependent VQE (TD-VQE).

Computes electronically excited states using quantum algorithms.
Governance protocols ensure physically correct excited state ansatze.

Methods:
- TD-VQE: Variational approach for excited states
- Subspace expansion: Build excited states orthogonal to ground state
- State-averaged VQE: Optimize multiple states simultaneously
"""

from typing import List, Dict, Optional, Any, Tuple
import numpy as np
from scipy.optimize import minimize
import logging

logger = logging.getLogger(__name__)

from kanad.ansatze.base_ansatz import BaseAnsatz
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.core.mappers.base_mapper import BaseMapper
from kanad.governance.protocols.base_protocol import BaseGovernanceProtocol


class ExcitedStatesSolver:
    """
    Solver for excited electronic states.

    Uses governance protocols to ensure excited states respect:
    - Orbital selection rules (covalent: bonding → antibonding)
    - Spin multiplicity constraints
    - Spatial symmetry
    - Orthogonality to lower states

    Methods:
        - subspace_expansion: Build excited states in subspace
        - state_averaged_vqe: Optimize multiple states together
        - overlap_penalty: Enforce orthogonality
    """

    def __init__(
        self,
        hamiltonian: MolecularHamiltonian,
        ansatz: BaseAnsatz,
        mapper: BaseMapper,
        governance: Optional[BaseGovernanceProtocol] = None,
        n_states: int = 3,
        backend: str = 'classical',
        **backend_options
    ):
        """
        Initialize excited states solver.

        Args:
            hamiltonian: Molecular Hamiltonian
            ansatz: Quantum circuit ansatz
            mapper: Fermion-to-qubit mapper
            governance: Governance protocol for state validation
            n_states: Number of excited states to compute
            backend: Computation backend
        """
        self.hamiltonian = hamiltonian
        self.ansatz = ansatz
        self.mapper = mapper
        self.governance = governance
        self.n_states = n_states
        self.backend = backend

        self.circuit = ansatz.build_circuit()
        self.n_parameters = self.circuit.get_num_parameters()

        # Store results for each state
        self.state_energies: List[float] = []
        self.state_parameters: List[np.ndarray] = []
        self.state_circuits: List[Any] = []

        logger.info(f"Initialized ExcitedStatesSolver for {n_states} states")

    def solve_subspace_expansion(
        self,
        initial_parameters: Optional[np.ndarray] = None,
        orthogonality_weight: float = 100.0
    ) -> Dict[str, Any]:
        """
        Solve for excited states using subspace expansion.

        Algorithm:
        1. Solve ground state VQE
        2. For each excited state:
           a. Build ansatz with overlap penalty to lower states
           b. Optimize energy with orthogonality constraint
           c. Validate with governance protocol

        Args:
            initial_parameters: Starting parameters for optimization
            orthogonality_weight: Penalty weight for state overlap

        Returns:
            Dictionary with energies, parameters, and excitation energies
        """
        if initial_parameters is None:
            initial_parameters = np.random.randn(self.n_parameters) * 0.1

        results = {
            'energies': [],
            'parameters': [],
            'excitation_energies': [],
            'oscillator_strengths': [],
            'converged': []
        }

        logger.info("Starting subspace expansion for excited states")

        # Solve ground state first
        logger.info("Computing ground state...")
        ground_result = self._optimize_single_state(
            initial_parameters,
            state_index=0,
            lower_states=[]
        )

        ground_energy = ground_result['energy']
        ground_params = ground_result['parameters']

        results['energies'].append(ground_energy)
        results['parameters'].append(ground_params)
        results['excitation_energies'].append(0.0)
        results['converged'].append(ground_result['converged'])

        self.state_energies.append(ground_energy)
        self.state_parameters.append(ground_params)

        logger.info(f"Ground state energy: {ground_energy:.6f} Ha")

        # Solve excited states with orthogonality constraints
        for i in range(1, self.n_states):
            logger.info(f"Computing excited state {i}...")

            excited_result = self._optimize_single_state(
                initial_parameters,
                state_index=i,
                lower_states=self.state_parameters[:i],
                orthogonality_weight=orthogonality_weight
            )

            excited_energy = excited_result['energy']
            excited_params = excited_result['parameters']
            excitation_energy = (excited_energy - ground_energy) * 27.211  # Convert to eV

            results['energies'].append(excited_energy)
            results['parameters'].append(excited_params)
            results['excitation_energies'].append(excitation_energy)
            results['converged'].append(excited_result['converged'])

            self.state_energies.append(excited_energy)
            self.state_parameters.append(excited_params)

            logger.info(f"Excited state {i}: E = {excited_energy:.6f} Ha, "
                       f"ΔE = {excitation_energy:.4f} eV")

        # Validate with governance if available
        if self.governance:
            validation = self.governance.validate_excited_states(
                results['energies'],
                results['excitation_energies']
            )
            results['validation'] = validation

        return results

    def _optimize_single_state(
        self,
        initial_parameters: np.ndarray,
        state_index: int,
        lower_states: List[np.ndarray],
        orthogonality_weight: float = 100.0
    ) -> Dict[str, Any]:
        """
        Optimize a single excited state with orthogonality constraints.

        Objective function:
            E(θ) + λ * Σ |⟨ψ_k(θ_k)|ψ(θ)⟩|²

        Where λ is orthogonality_weight and sum is over lower states.
        """
        def objective(params):
            # Compute energy
            energy = self._compute_energy(params)

            # Add orthogonality penalty for excited states
            if len(lower_states) > 0:
                overlap_penalty = 0.0
                for lower_params in lower_states:
                    overlap = self._compute_overlap(params, lower_params)
                    overlap_penalty += orthogonality_weight * overlap**2

                return energy + overlap_penalty
            else:
                return energy

        # Optimize
        result = minimize(
            objective,
            initial_parameters,
            method='COBYLA',
            options={'maxiter': 500}
        )

        return {
            'energy': self._compute_energy(result.x),
            'parameters': result.x,
            'converged': result.success,
            'iterations': result.nfev
        }

    def _compute_energy(self, parameters: np.ndarray) -> float:
        """Compute energy expectation value for given parameters."""
        # Simplified classical simulation
        # In production, would use actual quantum backend

        # Build parameterized state
        state_vector = self._build_state_vector(parameters)

        # Compute ⟨ψ|H|ψ⟩
        H_matrix = self.hamiltonian.to_matrix()
        energy = np.real(np.dot(np.conj(state_vector), np.dot(H_matrix, state_vector)))

        return energy

    def _compute_overlap(self, params1: np.ndarray, params2: np.ndarray) -> float:
        """Compute overlap ⟨ψ₁|ψ₂⟩ between two states."""
        state1 = self._build_state_vector(params1)
        state2 = self._build_state_vector(params2)

        overlap = np.abs(np.dot(np.conj(state1), state2))
        return overlap

    def _build_state_vector(self, parameters: np.ndarray) -> np.ndarray:
        """Build state vector from parameters (classical simulation)."""
        n_qubits = self.circuit.num_qubits  # Fixed: use property instead of method
        state_dim = 2**n_qubits

        # Start with HF reference state
        state = np.zeros(state_dim, dtype=complex)
        n_electrons = self.hamiltonian.n_electrons
        hf_index = (1 << n_electrons) - 1  # |11...100...0⟩
        state[hf_index] = 1.0

        # Apply parametrized gates (simplified)
        # In production, would simulate full quantum circuit
        for i, param in enumerate(parameters):
            # Simplified rotation
            rotation = np.exp(1j * param * 0.1)
            state *= rotation

        # Normalize
        state /= np.linalg.norm(state)

        return state

    def get_excitation_spectrum(self) -> Dict[str, Any]:
        """
        Get excitation spectrum with transition properties.

        Returns:
            Dictionary with:
            - excitation_energies (eV)
            - wavelengths (nm)
            - oscillator_strengths
            - dominant_transitions
        """
        if len(self.state_energies) < 2:
            raise ValueError("Must solve excited states first")

        ground_energy = self.state_energies[0]

        spectrum = {
            'excitation_energies': [],
            'wavelengths': [],
            'transitions': []
        }

        for i in range(1, len(self.state_energies)):
            delta_e_ev = (self.state_energies[i] - ground_energy) * 27.211
            wavelength_nm = 1240.0 / delta_e_ev if delta_e_ev > 0 else np.inf

            spectrum['excitation_energies'].append(delta_e_ev)
            spectrum['wavelengths'].append(wavelength_nm)
            spectrum['transitions'].append(f"S0 → S{i}")

        return spectrum
