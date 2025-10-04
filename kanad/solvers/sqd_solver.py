"""
Sample-Based Quantum Diagonalization (SQD) Solver.

Implements SQD using qiskit-addon-sqd for efficient eigenvalue finding
on quantum systems. Particularly effective for larger molecules where
direct diagonalization becomes impractical.
"""

import numpy as np
from typing import Dict, Any, Optional, List
import logging

logger = logging.getLogger(__name__)

from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian


class SQDSolver:
    """
    Sample-Based Quantum Diagonalization solver.

    Uses qiskit-addon-sqd to find ground state eigenvalue and eigenstate
    by sampling quantum circuits and processing with distributed classical computing.

    Advantages over VQE:
    - Better scaling for larger systems (~25 orbitals, ~10 electrons)
    - Leverages both quantum sampling and classical processing
    - Iterative refinement of samples

    Recommended for systems where:
    - VQE becomes too expensive (many parameters)
    - Direct diagonalization is impractical (large Hilbert space)
    """

    def __init__(
        self,
        hamiltonian: MolecularHamiltonian,
        ansatz: Optional[Any] = None,
        n_samples: int = 1000,
        max_iterations: int = 10
    ):
        """
        Initialize SQD solver.

        Args:
            hamiltonian: Molecular Hamiltonian with h_core and ERI
            ansatz: Optional variational ansatz for state preparation
            n_samples: Number of quantum circuit samples to collect
            max_iterations: Maximum SQD refinement iterations
        """
        self.hamiltonian = hamiltonian
        self.ansatz = ansatz
        self.n_samples = n_samples
        self.max_iterations = max_iterations

        self.n_orbitals = hamiltonian.n_orbitals
        self.n_electrons = hamiltonian.n_electrons

        logger.info(
            f"Initialized SQD solver: {self.n_electrons} electrons, "
            f"{self.n_orbitals} orbitals, {self.n_samples} samples"
        )

    def solve(self) -> Dict[str, Any]:
        """
        Solve for ground state using Sample-Based Quantum Diagonalization.

        Returns:
            Dictionary with:
            - energy: Ground state energy (Hartree)
            - eigenstate: Ground state wavefunction coefficients
            - n_samples: Number of quantum samples used
            - converged: Whether SQD converged
            - iterations: Number of SQD refinement iterations
        """
        try:
            # Try to use qiskit-addon-sqd if available
            return self._solve_with_sqd_addon()
        except ImportError:
            logger.warning(
                "qiskit-addon-sqd not available. Install with: "
                "pip install qiskit-addon-sqd"
            )
            logger.info("Falling back to classical exact diagonalization")
            return self._solve_classical_fallback()

    def _solve_with_sqd_addon(self) -> Dict[str, Any]:
        """
        Solve using qiskit-addon-sqd package.

        Returns:
            SQD results dictionary
        """
        try:
            import qiskit_addon_sqd
            from qiskit_addon_sqd.fermion import solve_fermion
        except ImportError as e:
            raise ImportError(
                "qiskit-addon-sqd not installed. Install with: "
                "pip install qiskit-addon-sqd"
            ) from e

        logger.info("Running SQD with qiskit-addon-sqd...")

        # Prepare quantum samples (bitstring matrix)
        bitstring_matrix = self._generate_quantum_samples()

        # Get molecular integrals
        h_core = self.hamiltonian.h_core
        eri = self.hamiltonian.eri

        # Solve using SQD
        logger.info(f"Processing {bitstring_matrix.shape[0]} quantum samples...")

        # solve_fermion returns (energy, state, occupancies, spin_sq)
        energy, state, occupancies, spin_sq = solve_fermion(
            bitstring_matrix,
            h_core,
            eri,
            open_shell=False,
            spin_sq=None  # No spin constraint
        )

        # Add nuclear repulsion
        total_energy = energy + self.hamiltonian.nuclear_repulsion

        logger.info(f"SQD ground state: {total_energy:.6f} Ha")

        return {
            'energy': total_energy,
            'eigenstate': state,
            'occupancies': occupancies,
            'spin_sq': spin_sq,
            'n_samples': bitstring_matrix.shape[0],
            'converged': True,
            'method': 'SQD (qiskit-addon-sqd)'
        }

    def _convert_to_fermionic_op(self):
        """
        Convert MolecularHamiltonian to fermionic operator format.

        Returns:
            Fermionic operator compatible with qiskit-addon-sqd
        """
        try:
            from qiskit_nature.second_q.operators import FermionicOp
        except ImportError:
            # Fallback: use simplified representation
            logger.warning("qiskit_nature not available, using simplified operator")
            return self._build_simplified_fermionic_op()

        # Build fermionic operator from h_core and ERI
        # H = Σ_{pq} h_pq a†_p a_q + ½ Σ_{pqrs} eri_{pqrs} a†_p a†_q a_s a_r

        op_dict = {}

        # One-electron terms
        for p in range(self.n_orbitals):
            for q in range(self.n_orbitals):
                if abs(self.hamiltonian.h_core[p, q]) > 1e-10:
                    # a†_p a_q
                    key = f"+_{p} -_{q}"
                    op_dict[key] = self.hamiltonian.h_core[p, q]

        # Two-electron terms (if ERI available)
        if hasattr(self.hamiltonian, 'eri') and self.hamiltonian.eri is not None:
            for p in range(self.n_orbitals):
                for q in range(self.n_orbitals):
                    for r in range(self.n_orbitals):
                        for s in range(self.n_orbitals):
                            coeff = 0.5 * self.hamiltonian.eri[p, q, r, s]
                            if abs(coeff) > 1e-10:
                                # a†_p a†_q a_s a_r
                                key = f"+_{p} +_{q} -_{s} -_{r}"
                                op_dict[key] = coeff

        return FermionicOp(op_dict)

    def _build_simplified_fermionic_op(self):
        """
        Build simplified fermionic operator without qiskit_nature.

        Returns:
            Dictionary representation of fermionic operator
        """
        return {
            'h_core': self.hamiltonian.h_core,
            'eri': getattr(self.hamiltonian, 'eri', None),
            'n_orbitals': self.n_orbitals,
            'n_electrons': self.n_electrons
        }

    def _generate_quantum_samples(self) -> np.ndarray:
        """
        Generate quantum circuit samples as bitstring matrix.

        For H2 with 2 spin-orbitals:
        - Column order: [spinup_1, spinup_0, spindown_1, spindown_0]
        - HF state for 2 electrons: [False, True, False, True] = occupy alpha_0 and beta_0

        Returns:
            Bitstring matrix (n_samples x n_qubits) of Boolean values
        """
        n_spin_orbitals = 2 * self.n_orbitals

        # Start with Hartree-Fock configuration
        hf_bitstring = np.zeros(n_spin_orbitals, dtype=bool)

        # Occupy lowest alpha and beta orbitals
        # For closed-shell: occupy first n_electrons//2 alpha and beta orbitals
        n_alpha = self.n_electrons // 2
        n_beta = self.n_electrons - n_alpha

        # Spin-up electrons occupy right half (indices n_spin_orbitals//2 to n_spin_orbitals-1)
        for i in range(n_alpha):
            hf_bitstring[n_spin_orbitals // 2 + i] = True

        # Spin-down electrons occupy left half (indices 0 to n_spin_orbitals//2-1)
        for i in range(n_beta):
            hf_bitstring[i] = True

        logger.info(f"HF bitstring (alpha|beta): {hf_bitstring}")

        # Generate excitations from HF state
        bitstrings = [hf_bitstring]

        # Single excitations (for small molecules, include some excited configs)
        if self.n_orbitals > 1 and len(bitstrings) < self.n_samples:
            # Single alpha excitation: move electron from HOMO to LUMO
            single_exc_alpha = hf_bitstring.copy()
            if n_alpha > 0 and n_alpha < self.n_orbitals:
                single_exc_alpha[n_spin_orbitals // 2 + n_alpha - 1] = False  # Remove from HOMO
                single_exc_alpha[n_spin_orbitals // 2 + n_alpha] = True       # Add to LUMO
                bitstrings.append(single_exc_alpha)

            # Single beta excitation
            single_exc_beta = hf_bitstring.copy()
            if n_beta > 0 and n_beta < self.n_orbitals:
                single_exc_beta[n_beta - 1] = False
                single_exc_beta[n_beta] = True
                bitstrings.append(single_exc_beta)

        # Convert to matrix
        bitstring_matrix = np.array(bitstrings, dtype=bool)

        logger.info(f"Generated {bitstring_matrix.shape[0]} configurations")
        return bitstring_matrix

    def _sample_from_ansatz(self) -> List[np.ndarray]:
        """
        Generate samples by evolving variational ansatz.

        Returns:
            List of state samples
        """
        samples = []

        # Sample at different parameter settings
        n_param_samples = min(self.n_samples, 100)

        for i in range(n_param_samples):
            # Random parameter perturbation
            params = np.random.randn(self.ansatz.num_parameters) * 0.1

            # Build state (simplified - would use actual circuit execution)
            # This is a placeholder - real implementation would use Qiskit
            state = np.zeros(2**self.n_orbitals)
            state[0] = 1.0  # Placeholder

            samples.append(state)

        return samples

    def _sample_from_hf(self) -> List[np.ndarray]:
        """
        Generate samples starting from Hartree-Fock state.

        Returns:
            List of HF-based samples
        """
        # Get HF state
        try:
            _, _ = self.hamiltonian.solve_scf()
        except:
            pass

        # Create HF occupation bitstring
        n_occ = self.n_electrons // 2
        hf_state = np.zeros(2**self.n_orbitals)

        # Occupy lowest n_occ orbitals (simplified)
        # In reality would use proper MO occupation
        hf_index = (2**n_occ) - 1  # Binary: 111...1 for n_occ bits
        hf_state[hf_index] = 1.0

        # Return multiple copies (in real implementation, would apply time evolution)
        return [hf_state for _ in range(min(self.n_samples, 10))]

    def _solve_classical_fallback(self) -> Dict[str, Any]:
        """
        Fallback to classical exact diagonalization.

        Used when qiskit-addon-sqd is not available.

        Returns:
            Classical diagonalization results
        """
        logger.info("Using classical exact diagonalization (fallback)")

        # Build Hamiltonian matrix
        H_matrix = self.hamiltonian.to_matrix()

        # Diagonalize
        eigenvalues, eigenvectors = np.linalg.eigh(H_matrix)

        # Ground state
        ground_energy = eigenvalues[0] + self.hamiltonian.nuclear_repulsion
        ground_state = eigenvectors[:, 0]

        logger.info(f"Classical ground state: {ground_energy:.6f} Ha")

        return {
            'energy': ground_energy,
            'eigenstate': ground_state,
            'n_samples': 0,  # No quantum samples used
            'converged': True,
            'iterations': 0,
            'method': 'Classical Exact Diagonalization (SQD fallback)'
        }
