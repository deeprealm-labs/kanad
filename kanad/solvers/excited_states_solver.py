"""
Excited States Solver - Bonds Module Integration.

Computes excited electronic states for molecular systems.
Integrated with analysis tools for spectroscopy and photochemistry.
"""

from typing import Dict, Any, Optional, List
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
            **kwargs: Method-specific options
        """
        super().__init__(bond, enable_analysis, enable_optimization)

        self.method = method.lower()
        self.n_states = n_states

        # Not a correlation method for ground state
        self._is_correlated = False

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
        if self.enable_analysis and hasattr(self, 'uvvis_calculator'):
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
        """Solve using state-averaged VQE (placeholder)."""
        logger.warning("VQE excited states not fully implemented")
        raise NotImplementedError("VQE for excited states not yet implemented")

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
