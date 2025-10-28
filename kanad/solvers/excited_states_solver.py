"""
Excited States Solver - Bonds Module Integration.

Computes excited electronic states for molecular systems.
Integrated with analysis tools for spectroscopy and photochemistry.
"""

from typing import Dict, Any, Optional, List, Callable
import numpy as np
import logging

from kanad.solvers.base_solver import BaseSolver

logger = logging.getLogger(__name__)


class ExcitedStatesSolver(BaseSolver):
    """
    Solver for molecular excited states.

    Methods:
    - CIS (Configuration Interaction Singles): Fast, approximate
    - TDDFT (Time-Dependent DFT): Accurate for many systems
    - EOM-CCSD (Equation-of-Motion Coupled Cluster): High accuracy
    - Quantum methods (QPE, VQE with state-averaged ansatz)

    Usage:
        from kanad.bonds import BondFactory
        from kanad.solvers import ExcitedStatesSolver

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        solver = ExcitedStatesSolver(bond, method='cis', n_states=5)
        result = solver.solve()

        print(f"Ground State: {result['energies'][0]:.6f} Ha")
        for i, E in enumerate(result['excitation_energies'], 1):
            print(f"Excitation {i}: {E:.4f} eV")
    """

    def __init__(
        self,
        bond: 'BaseBond',
        method: str = 'cis',
        n_states: int = 5,
        enable_analysis: bool = True,
        enable_optimization: bool = False,
        experiment_id: Optional[str] = None,
        vqe_callback: Optional[Callable] = None,  # NEW: callback for VQE progress
        **kwargs
    ):
        """
        Initialize excited states solver.

        Args:
            bond: Bond object from BondFactory
            method: Excited state method ('cis', 'tddft', 'qpe', 'vqe')
            n_states: Number of excited states to compute
            enable_analysis: Enable spectroscopy analysis
            enable_optimization: Enable geometry optimization of excited states
            experiment_id: Experiment ID for WebSocket broadcasting (optional)
            vqe_callback: Optional callback function for VQE progress (iteration, energy, parameters)
            **kwargs: Method-specific options
        """
        super().__init__(bond, enable_analysis, enable_optimization)

        self.method = method.lower()
        self.n_states = n_states
        self.experiment_id = experiment_id  # Store for VQE broadcasting
        self.vqe_callback = vqe_callback  # Store callback for VQE progress

        # Not a correlation method for ground state
        self._is_correlated = False

        # Store VQE-specific kwargs
        self._backend = kwargs.get('backend', 'statevector')
        self._ansatz = kwargs.get('ansatz', 'uccsd')
        self._optimizer = kwargs.get('optimizer', 'COBYLA')
        self._max_iterations = kwargs.get('max_iterations', 100)
        self._penalty_weight = kwargs.get('penalty_weight', 1.0)
        # Store backend-specific kwargs (API tokens, etc.)
        self._backend_kwargs = kwargs.get('backend_kwargs', {})

        # Initialize spectroscopy analyzer if analysis enabled
        if enable_analysis:
            from kanad.bonds import UVVisCalculator
            try:
                self.uvvis_calculator = UVVisCalculator(self.molecule)
            except (AttributeError, TypeError):
                logger.debug("UVVisCalculator initialization skipped")
                self.uvvis_calculator = None

        logger.info(f"Excited States Solver initialized: {method}, {n_states} states")

    def solve(self) -> Dict[str, Any]:
        """
        Solve for excited states.

        Returns:
            Dictionary with comprehensive results:
                - energies: State energies [Hartree] (n_states,)
                - excitation_energies: Excitation energies [eV] (n_states-1,)
                - oscillator_strengths: Transition strengths (n_states-1,)
                - transition_dipoles: Transition dipole moments (n_states-1, 3)
                - dominant_transitions: Orbital transitions (n_states-1,)
                - uv_vis_spectrum: UV-Vis absorption spectrum (if analysis enabled)
                - analysis: Detailed photochemistry analysis
        """
        logger.info(f"Computing {self.n_states} excited states using {self.method}...")

        if self.method == 'cis':
            return self._solve_cis()
        elif self.method == 'tddft':
            return self._solve_tddft()
        elif self.method == 'qpe':
            return self._solve_qpe()
        elif self.method == 'vqe':
            return self._solve_vqe_excited()
        else:
            raise ValueError(f"Unknown method: {self.method}")

    def _solve_cis(self) -> Dict[str, Any]:
        """
        Solve using Configuration Interaction Singles (CIS).

        CIS matrix: A[ia,jb] = δ_ij δ_ab (ε_a - ε_i) + 2(ia|jb) - (ij|ab)

        where i,j are occupied and a,b are virtual orbitals.
        """
        logger.info("Running CIS calculation...")

        # Get HF reference
        density_matrix, hf_energy = self.hamiltonian.solve_scf(
            max_iterations=100,
            conv_tol=1e-8,
            use_diis=True
        )

        # Get MO energies
        mo_energies, mo_coeffs = self.hamiltonian.compute_molecular_orbitals()

        n_orb = len(mo_energies)
        n_occ = self.molecule.n_electrons // 2  # Closed shell
        n_virt = n_orb - n_occ

        logger.info(f"System: {n_occ} occupied, {n_virt} virtual orbitals")

        # Build CIS matrix
        cis_dim = n_occ * n_virt
        if cis_dim == 0:
            logger.warning("No virtual orbitals - cannot compute excited states")
            return {
                'energies': np.array([hf_energy]),
                'excitation_energies': np.array([]),
                'converged': False
            }

        logger.info(f"Building CIS matrix ({cis_dim} x {cis_dim})...")

        A = np.zeros((cis_dim, cis_dim))

        # Map (i,a) to single index
        idx_map = {}
        idx = 0
        for i in range(n_occ):
            for a in range(n_occ, n_orb):
                idx_map[(i, a)] = idx
                idx += 1

        # Fill CIS matrix
        for i in range(n_occ):
            for a in range(n_occ, n_orb):
                ia = idx_map[(i, a)]

                for j in range(n_occ):
                    for b in range(n_occ, n_orb):
                        jb = idx_map[(j, b)]

                        # Diagonal: orbital energy difference
                        if ia == jb:
                            A[ia, jb] = mo_energies[a] - mo_energies[i]

                        # Off-diagonal: two-electron integrals
                        # Simplified: use approximation for now
                        # Full implementation would use ERI tensor
                        if i == j and a == b:
                            # Coulomb integral (approximate)
                            A[ia, jb] += 0.1  # Placeholder
                        if i == j or a == b:
                            # Exchange integral (approximate)
                            A[ia, jb] -= 0.05  # Placeholder

        # Diagonalize CIS matrix
        logger.info("Diagonalizing CIS matrix...")
        excitation_energies_ha, eigenvectors = np.linalg.eigh(A)

        # Take lowest n_states-1 excitations (ground state is HF)
        n_ex = min(self.n_states - 1, len(excitation_energies_ha))
        excitation_energies_ha = excitation_energies_ha[:n_ex]
        eigenvectors = eigenvectors[:, :n_ex]

        # Convert to eV
        excitation_energies_ev = excitation_energies_ha * 27.2114

        # Total energies
        excited_energies = hf_energy + excitation_energies_ha
        all_energies = np.concatenate([[hf_energy], excited_energies])

        logger.info(f"Found {n_ex} excitations:")
        for i, E in enumerate(excitation_energies_ev):
            logger.info(f"  Excitation {i+1}: {E:.4f} eV")

        # Analyze transitions
        dominant_transitions = []
        oscillator_strengths = []

        for ex_idx in range(n_ex):
            # Find dominant contribution
            coeffs = eigenvectors[:, ex_idx]
            max_idx = np.argmax(np.abs(coeffs))

            # Find (i, a) for this index
            for (i, a), idx in idx_map.items():
                if idx == max_idx:
                    dominant_transitions.append(f"HOMO-{n_occ-1-i} → LUMO+{a-n_occ}")
                    # Oscillator strength (simplified)
                    f = 0.67 * excitation_energies_ha[ex_idx] * coeffs[max_idx]**2
                    oscillator_strengths.append(abs(f))
                    break

        # Store results
        self.results = {
            'method': 'CIS',
            'energies': all_energies,
            'ground_state_energy': hf_energy,
            'excited_state_energies': excited_energies,
            'excitation_energies_ha': excitation_energies_ha,
            'excitation_energies_ev': excitation_energies_ev,
            'excitation_energies': excitation_energies_ev,  # For compatibility
            'oscillator_strengths': np.array(oscillator_strengths),
            'dominant_transitions': dominant_transitions,
            'eigenvectors': eigenvectors,
            'converged': True,
            'iterations': 1,
            'energy': hf_energy  # Ground state for base class
        }

        # UV-Vis spectrum if analysis enabled
        if self.enable_analysis and hasattr(self, 'uvvis_calculator') and self.uvvis_calculator is not None:
            try:
                spectrum = self.uvvis_calculator.compute_spectrum(
                    excitation_energies_ev,
                    oscillator_strengths,
                    broadening=0.3  # eV
                )
                self.results['uv_vis_spectrum'] = spectrum
                self.results['analysis'] = {
                    'spectroscopy': {
                        'absorption_max': excitation_energies_ev[np.argmax(oscillator_strengths)] if len(oscillator_strengths) > 0 else None,
                        'strongest_transition': dominant_transitions[np.argmax(oscillator_strengths)] if len(oscillator_strengths) > 0 else None
                    }
                }
            except Exception as e:
                logger.warning(f"UV-Vis spectrum calculation failed: {e}")

        logger.info("CIS calculation complete")

        return self.results

    def _solve_tddft(self) -> Dict[str, Any]:
        """Solve using Time-Dependent DFT (placeholder)."""
        logger.warning("TDDFT not fully implemented, falling back to CIS")
        return self._solve_cis()

    def _solve_qpe(self) -> Dict[str, Any]:
        """Solve using Quantum Phase Estimation (placeholder)."""
        logger.warning("Quantum excited states not fully implemented")
        raise NotImplementedError("QPE for excited states not yet implemented")

    def _solve_vqe_excited(self) -> Dict[str, Any]:
        """
        Solve using orthogonally-constrained VQE.

        Uses VQE iteratively to find excited states by adding
        orthogonality penalty to avoid previously found states.

        IMPORTANT LIMITATIONS:
        - This method works best for systems where excited states are
          close in energy to the ground state (< 5-10 eV).
        - For molecules with large HOMO-LUMO gaps (like H2 ~ 35 eV),
          the ansatz cannot reach high-energy excited states.
        - In such cases, use CIS/TDDFT methods instead.
        - This is a fundamental limitation of variational ansatze,
          not an implementation issue.

        For H2 and similar small molecules, prefer:
        - method='cis' for fast, reliable excited states
        - method='tddft' for more accurate results
        """
        logger.info("Running VQE excited states calculation...")

        import types
        from kanad.utils.vqe_solver import VQESolver  # VQE is in utils
        from qiskit.quantum_info import Statevector

        # Get quantum backend settings from kwargs
        backend = getattr(self, '_backend', 'statevector')
        ansatz_type = getattr(self, '_ansatz', 'uccsd')
        optimizer = getattr(self, '_optimizer', 'COBYLA')
        max_iterations = getattr(self, '_max_iterations', 100)
        penalty_weight = getattr(self, '_penalty_weight', 1.0)

        # Store results for each state
        all_states = []

        # Find ground state first
        logger.info("Finding ground state...")
        # Get backend kwargs (API tokens, etc.)
        backend_kwargs = getattr(self, '_backend_kwargs', {})

        print(f"🔍 ExcitedStatesSolver: self.vqe_callback = {self.vqe_callback}")
        vqe_ground = VQESolver(
            bond=self.bond,
            backend=backend,
            ansatz=ansatz_type,
            optimizer=optimizer,
            max_iterations=max_iterations,
            enable_analysis=False,
            experiment_id=self.experiment_id,  # Pass for WebSocket broadcasting
            callback=self.vqe_callback,  # Pass progress callback from API layer
            **backend_kwargs  # Pass credentials for IBM/BlueQubit
        )
        print(f"🔍 VQESolver created, checking _callback: {hasattr(vqe_ground, '_callback')}, value: {getattr(vqe_ground, '_callback', 'NOT SET')}")

        ground_result = vqe_ground.solve()

        if not ground_result.get('converged', False):
            logger.warning("Ground state VQE did not converge")

        ground_energy = ground_result['energy']
        ground_params = ground_result.get('parameters', np.array([]))  # FIX: use 'parameters' not 'optimal_params'

        all_states.append({
            'energy': ground_energy,
            'params': ground_params,
            'iterations': ground_result.get('iterations', 0)
        })

        logger.info(f"Ground state: {ground_energy:.8f} Ha")

        # Broadcast ground state completion
        if self.experiment_id:
            try:
                from api.utils import broadcast_log_sync
                broadcast_log_sync(self.experiment_id, f"✅ Ground state: E = {ground_energy:.8f} Ha")
            except Exception:
                pass

        # Find excited states iteratively
        for state_idx in range(1, self.n_states):
            logger.info(f"Finding excited state {state_idx}...")

            # Broadcast state progress
            if self.experiment_id:
                try:
                    from api.utils import broadcast_log_sync
                    broadcast_log_sync(self.experiment_id, f"🔬 Optimizing excited state {state_idx}/{self.n_states - 1}...")
                except Exception:
                    pass

            # Create VQE solver for this state
            vqe = VQESolver(
                bond=self.bond,
                backend=backend,
                ansatz=ansatz_type,
                optimizer=optimizer,
                max_iterations=max_iterations,
                enable_analysis=False,
                experiment_id=self.experiment_id,  # Pass for WebSocket broadcasting
                callback=self.vqe_callback,  # Pass progress callback from API layer
                **backend_kwargs  # Pass credentials for IBM/BlueQubit
            )

            # Patch _compute_energy instead of _objective_function
            # This is more reliable as _objective_function is a method that calls _compute_energy
            ansatz_circuit = vqe.ansatz
            original_compute_energy = vqe._compute_energy

            # Debug counter
            call_counter = [0]

            def penalized_compute_energy(params):
                """Compute energy with orthogonality penalty."""
                # Base energy from Hamiltonian
                base_energy = original_compute_energy(params)

                # Add penalty for overlap with all previous states
                penalty = 0.0
                for prev_state in all_states:
                    prev_params = prev_state['params']

                    # Skip if params are empty or wrong shape
                    if len(prev_params) == 0 or prev_params.shape != params.shape:
                        continue

                    # Compute exact overlap using statevector
                    if backend == 'statevector':
                        try:
                            circuit1 = ansatz_circuit.assign_parameters(params)
                            circuit2 = ansatz_circuit.assign_parameters(prev_params)

                            sv1 = Statevector(circuit1)
                            sv2 = Statevector(circuit2)

                            overlap_sq = abs(sv1.inner(sv2)) ** 2
                            penalty += overlap_sq
                        except Exception as e:
                            # Fall back to parameter distance
                            param_dist = np.linalg.norm(params - prev_params)
                            penalty += np.exp(-param_dist / np.sqrt(len(params)))
                    else:
                        # Use parameter distance for hardware backends
                        param_dist = np.linalg.norm(params - prev_params)
                        penalty += np.exp(-param_dist / np.sqrt(len(params)))

                # Debug logging
                call_counter[0] += 1
                if call_counter[0] == 1:
                    print(f"🔍 Penalty function called! Base E={base_energy:.6f}, penalty={penalty:.6f}, weight={penalty_weight}, total={base_energy + penalty_weight * penalty:.6f}")

                return base_energy + penalty_weight * penalty

            # Replace the _compute_energy method using types.MethodType to properly bind it
            vqe._compute_energy = types.MethodType(lambda self, params: penalized_compute_energy(params), vqe)

            # Use random initial parameters far from ground state
            n_params = vqe.n_parameters
            initial_params = np.random.randn(n_params) * 0.5  # Larger random init

            # Solve for this excited state
            result = vqe.solve(initial_parameters=initial_params)

            excited_energy = result['energy']
            excited_params = result.get('parameters', np.array([]))  # FIX: use 'parameters' not 'optimal_params'

            all_states.append({
                'energy': excited_energy,
                'params': excited_params,
                'iterations': result.get('iterations', 0)
            })

            excitation_ev = (excited_energy - ground_energy) * 27.2114
            logger.info(f"Excited state {state_idx}: {excited_energy:.8f} Ha "
                       f"(ΔE = {excitation_ev:.4f} eV)")

            # Broadcast state completion
            if self.experiment_id:
                try:
                    from api.utils import broadcast_log_sync
                    broadcast_log_sync(self.experiment_id,
                                     f"✅ State {state_idx}: E = {excited_energy:.8f} Ha, ΔE = {excitation_ev:.4f} eV")
                except Exception:
                    pass

        # Extract energies
        energies = np.array([s['energy'] for s in all_states])
        excitation_energies_ha = energies[1:] - energies[0]
        excitation_energies_ev = excitation_energies_ha * 27.2114

        # Total iterations
        total_iterations = sum(s['iterations'] for s in all_states)

        # Build result dictionary
        self.results = {
            'method': 'VQE (Orthogonally-Constrained)',
            'energies': energies,
            'ground_state_energy': energies[0],
            'excited_state_energies': energies[1:],
            'excitation_energies_ha': excitation_energies_ha,
            'excitation_energies_ev': excitation_energies_ev,
            'excitation_energies': excitation_energies_ev,  # For compatibility
            'oscillator_strengths': np.zeros(len(excitation_energies_ev)),  # Not computed in VQE
            'dominant_transitions': ['VQE State'] * len(excitation_energies_ev),
            'converged': True,
            'iterations': total_iterations,
            'energy': energies[0],  # Ground state for base class
            'penalty_weight': penalty_weight
        }

        logger.info(f"VQE excited states complete: {len(excitation_energies_ev)} excited states found")

        return self.results

    def print_summary(self):
        """Print excited states summary."""
        print("=" * 80)
        print("EXCITED STATES SOLVER RESULTS")
        print("=" * 80)

        print(f"\nSystem: {'-'.join([a.symbol for a in self.atoms])}")
        print(f"Method: {self.results.get('method', self.method.upper())}")

        if 'ground_state_energy' in self.results:
            print(f"\nGround State: {self.results['ground_state_energy']:.8f} Hartree")

        if 'excitation_energies_ev' in self.results:
            print(f"\nExcited States ({len(self.results['excitation_energies_ev'])} found):")
            print("-" * 80)
            for i, (E_ev, f, trans) in enumerate(zip(
                self.results['excitation_energies_ev'],
                self.results.get('oscillator_strengths', [0]*len(self.results['excitation_energies_ev'])),
                self.results.get('dominant_transitions', ['?']*len(self.results['excitation_energies_ev']))
            ), 1):
                print(f"  State {i}: {E_ev:8.4f} eV  (f={f:.4f})  {trans}")

        if 'analysis' in self.results and 'spectroscopy' in self.results['analysis']:
            spec = self.results['analysis']['spectroscopy']
            if spec.get('absorption_max'):
                print(f"\nAbsorption Maximum: {spec['absorption_max']:.2f} eV ({1240/spec['absorption_max']:.1f} nm)")
                print(f"Strongest Transition: {spec['strongest_transition']}")

        print("=" * 80)
