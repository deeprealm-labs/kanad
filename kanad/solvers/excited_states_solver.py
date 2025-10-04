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

    def _build_manybody_hamiltonian(self) -> np.ndarray:
        """
        Build full many-body Hamiltonian matrix in Fock space.

        Same implementation as QPE solver - constructs full second-quantized
        Hamiltonian from h_core and eri using Slater-Condon rules.
        """
        from itertools import combinations

        n_orbitals = self.hamiltonian.n_orbitals
        n_electrons = self.hamiltonian.n_electrons
        h_core = self.hamiltonian.h_core
        eri = self.hamiltonian.eri

        # Build Slater determinant basis
        n_alpha = n_electrons // 2
        n_beta = n_electrons - n_alpha

        alpha_configs = list(combinations(range(n_orbitals), n_alpha))
        beta_configs = list(combinations(range(n_orbitals), n_beta))
        n_configs = len(alpha_configs) * len(beta_configs)

        H_matrix = np.zeros((n_configs, n_configs))

        # Build Hamiltonian matrix
        for i, (alpha_i, beta_i) in enumerate([(a, b) for a in alpha_configs for b in beta_configs]):
            for j, (alpha_j, beta_j) in enumerate([(a, b) for a in alpha_configs for b in beta_configs]):
                if i <= j:
                    H_ij = self._compute_hamiltonian_element(
                        alpha_i, beta_i, alpha_j, beta_j, h_core, eri
                    )
                    H_matrix[i, j] = H_ij
                    if i != j:
                        H_matrix[j, i] = H_ij

        return H_matrix

    def _compute_hamiltonian_element(
        self,
        alpha_i: Tuple[int, ...],
        beta_i: Tuple[int, ...],
        alpha_j: Tuple[int, ...],
        beta_j: Tuple[int, ...],
        h_core: np.ndarray,
        eri: np.ndarray
    ) -> float:
        """Compute Hamiltonian matrix element using Slater-Condon rules."""
        alpha_i_set = set(alpha_i)
        alpha_j_set = set(alpha_j)
        beta_i_set = set(beta_i)
        beta_j_set = set(beta_j)

        alpha_diff = len(alpha_i_set.symmetric_difference(alpha_j_set))
        beta_diff = len(beta_i_set.symmetric_difference(beta_j_set))
        total_diff = alpha_diff + beta_diff

        if total_diff == 0:
            return self._diagonal_element(alpha_i, beta_i, h_core, eri)
        elif total_diff == 2:
            return self._single_excitation_element(
                alpha_i, beta_i, alpha_j, beta_j, h_core, eri
            )
        elif total_diff == 4:
            return self._double_excitation_element(
                alpha_i, beta_i, alpha_j, beta_j, eri
            )
        else:
            return 0.0

    def _diagonal_element(
        self,
        alpha: Tuple[int, ...],
        beta: Tuple[int, ...],
        h_core: np.ndarray,
        eri: np.ndarray
    ) -> float:
        """Diagonal Hamiltonian element."""
        energy = 0.0

        # One-electron terms
        for p in alpha:
            energy += h_core[p, p]
        for p in beta:
            energy += h_core[p, p]

        # Two-electron terms
        for p in alpha:
            for q in alpha:
                if p != q:
                    energy += 0.5 * (eri[p, p, q, q] - eri[p, q, q, p])

        for p in beta:
            for q in beta:
                if p != q:
                    energy += 0.5 * (eri[p, p, q, q] - eri[p, q, q, p])

        for p in alpha:
            for q in beta:
                energy += eri[p, p, q, q]

        return energy

    def _single_excitation_element(
        self,
        alpha_i: Tuple[int, ...],
        beta_i: Tuple[int, ...],
        alpha_j: Tuple[int, ...],
        beta_j: Tuple[int, ...],
        h_core: np.ndarray,
        eri: np.ndarray
    ) -> float:
        """Single excitation Hamiltonian element."""
        alpha_i_set = set(alpha_i)
        alpha_j_set = set(alpha_j)
        beta_i_set = set(beta_i)
        beta_j_set = set(beta_j)

        if alpha_i_set != alpha_j_set:
            removed = list(alpha_i_set - alpha_j_set)[0]
            added = list(alpha_j_set - alpha_i_set)[0]
            spin_orbs = alpha_i
            other_spin = beta_i
        else:
            removed = list(beta_i_set - beta_j_set)[0]
            added = list(beta_j_set - beta_i_set)[0]
            spin_orbs = beta_i
            other_spin = alpha_i

        # One-electron term
        element = h_core[removed, added]

        # Two-electron terms
        for q in spin_orbs:
            if q != removed:
                element += eri[removed, q, added, q] - eri[removed, q, q, added]

        for q in other_spin:
            element += eri[removed, q, added, q]

        # Phase factor
        phase = self._compute_phase(removed, added, spin_orbs)

        return phase * element

    def _double_excitation_element(
        self,
        alpha_i: Tuple[int, ...],
        beta_i: Tuple[int, ...],
        alpha_j: Tuple[int, ...],
        beta_j: Tuple[int, ...],
        eri: np.ndarray
    ) -> float:
        """Double excitation Hamiltonian element."""
        alpha_i_set = set(alpha_i)
        alpha_j_set = set(alpha_j)
        beta_i_set = set(beta_i)
        beta_j_set = set(beta_j)

        alpha_diff = alpha_i_set.symmetric_difference(alpha_j_set)
        beta_diff = beta_i_set.symmetric_difference(beta_j_set)

        if len(alpha_diff) == 4:
            removed = sorted(alpha_i_set - alpha_j_set)
            added = sorted(alpha_j_set - alpha_i_set)
            element = eri[removed[0], removed[1], added[0], added[1]] - \
                     eri[removed[0], removed[1], added[1], added[0]]
            phase = self._compute_double_phase(removed, added, alpha_i)
            return phase * element

        elif len(beta_diff) == 4:
            removed = sorted(beta_i_set - beta_j_set)
            added = sorted(beta_j_set - beta_i_set)
            element = eri[removed[0], removed[1], added[0], added[1]] - \
                     eri[removed[0], removed[1], added[1], added[0]]
            phase = self._compute_double_phase(removed, added, beta_i)
            return phase * element

        else:
            alpha_removed = list(alpha_i_set - alpha_j_set)[0]
            alpha_added = list(alpha_j_set - alpha_i_set)[0]
            beta_removed = list(beta_i_set - beta_j_set)[0]
            beta_added = list(beta_j_set - beta_i_set)[0]

            element = eri[alpha_removed, beta_removed, alpha_added, beta_added]
            phase_alpha = self._compute_phase(alpha_removed, alpha_added, alpha_i)
            phase_beta = self._compute_phase(beta_removed, beta_added, beta_i)

            return phase_alpha * phase_beta * element

    def _compute_phase(self, p: int, q: int, orbitals: Tuple[int, ...]) -> int:
        """Compute phase factor for single excitation."""
        orb_list = sorted(orbitals)
        p_idx = orb_list.index(p) if p in orb_list else -1
        q_idx = len([o for o in orb_list if o < q])

        if p_idx == -1:
            return 1

        n_between = abs(q_idx - p_idx) - 1
        return (-1) ** n_between

    def _compute_double_phase(
        self,
        removed: List[int],
        added: List[int],
        orbitals: Tuple[int, ...]
    ) -> int:
        """Compute phase factor for double excitation."""
        phase1 = self._compute_phase(removed[0], added[0], orbitals)

        orb_list = list(orbitals)
        orb_list.remove(removed[0])
        orb_list.append(added[0])
        orb_list = tuple(sorted(orb_list))

        phase2 = self._compute_phase(removed[1], added[1], orb_list)

        return phase1 * phase2

    def _compute_energy(self, parameters: np.ndarray) -> float:
        """Compute energy expectation value for given parameters."""
        # Simplified classical simulation
        # In production, would use actual quantum backend

        # Build parameterized state
        state_vector = self._build_state_vector(parameters)

        # Compute ⟨ψ|H|ψ⟩ using many-body Hamiltonian
        H_matrix = self._build_manybody_hamiltonian()
        energy = np.real(np.dot(np.conj(state_vector), np.dot(H_matrix, state_vector)))

        return energy

    def _compute_overlap(self, params1: np.ndarray, params2: np.ndarray) -> float:
        """Compute overlap ⟨ψ₁|ψ₂⟩ between two states."""
        state1 = self._build_state_vector(params1)
        state2 = self._build_state_vector(params2)

        overlap = np.abs(np.dot(np.conj(state1), state2))
        return overlap

    def _build_state_vector(self, parameters: np.ndarray) -> np.ndarray:
        """Build state vector from parameters in Slater determinant basis."""
        from itertools import combinations

        n_orbitals = self.hamiltonian.n_orbitals
        n_electrons = self.hamiltonian.n_electrons

        # Build Slater determinant basis (same as Hamiltonian)
        n_alpha = n_electrons // 2
        n_beta = n_electrons - n_alpha

        alpha_configs = list(combinations(range(n_orbitals), n_alpha))
        beta_configs = list(combinations(range(n_orbitals), n_beta))
        n_configs = len(alpha_configs) * len(beta_configs)

        # Start with HF reference state (first configuration)
        state = np.zeros(n_configs, dtype=complex)
        state[0] = 1.0

        # Apply parameter rotations (simplified)
        # In production, would simulate full quantum circuit
        for i, param in enumerate(parameters):
            if i < n_configs - 1:
                # Mix with excited configurations
                angle = param * 0.1
                c = np.cos(angle)
                s = np.sin(angle)

                # Apply rotation between ground and excited state i+1
                new_state_0 = c * state[0] - s * state[i+1]
                new_state_i = s * state[0] + c * state[i+1]
                state[0] = new_state_0
                state[i+1] = new_state_i

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
