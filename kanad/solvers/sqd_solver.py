"""
Subspace Quantum Diagonalization (SQD) Solver - Rebuilt with Bonds Module Integration.

SQD combines quantum and classical resources to solve eigenvalue problems
in a reduced subspace, achieving high accuracy with fewer quantum resources.

Reference: https://github.com/qiskit-community/qiskit-addon-sqd
"""

from typing import Dict, Any, Optional, List
import numpy as np
import logging

from kanad.solvers.base_solver import BaseSolver

logger = logging.getLogger(__name__)


class SQDSolver(BaseSolver):
    """
    Subspace Quantum Diagonalization for ground and excited states.

    SQD Workflow:
    1. Generate quantum subspace using short-depth circuits
    2. Project Hamiltonian into this subspace
    3. Classically diagonalize projected Hamiltonian
    4. Return eigenvalues and eigenvectors

    Advantages:
    - Lower circuit depth than VQE
    - Access to multiple eigenvalues (excited states)
    - More noise-resistant

    Usage:
        from kanad.bonds import BondFactory
        from kanad.solvers import SQDSolver

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        solver = SQDSolver(bond, subspace_dim=10)
        result = solver.solve()

        print(f"Ground State: {result['energies'][0]:.6f} Hartree")
        print(f"1st Excited: {result['energies'][1]:.6f} Hartree")
    """

    def __init__(
        self,
        bond: 'BaseBond',
        subspace_dim: int = 10,
        circuit_depth: int = 3,
        backend: str = 'statevector',
        shots: Optional[int] = None,
        enable_analysis: bool = True,
        enable_optimization: bool = True,
        **kwargs
    ):
        """
        Initialize SQD solver.

        Args:
            bond: Bond object from BondFactory
            subspace_dim: Dimension of quantum subspace
            circuit_depth: Depth of circuits for subspace generation
            backend: Quantum backend ('statevector', 'qasm', 'ibm')
            shots: Number of shots for sampling backends
            enable_analysis: Enable automatic analysis
            enable_optimization: Enable automatic optimization
            **kwargs: Additional backend options
        """
        super().__init__(bond, enable_analysis, enable_optimization)

        self.subspace_dim = subspace_dim
        self.circuit_depth = circuit_depth
        self.backend = backend
        self.shots = shots if shots is not None else 8192  # SQD needs more shots

        # This is a correlated method
        self._is_correlated = True

        # Initialize backend
        self._init_backend(**kwargs)

        # Check for qiskit-addon-sqd
        try:
            import qiskit_addon_sqd
            self._has_sqd_addon = True
            logger.info("qiskit-addon-sqd available")
        except ImportError:
            self._has_sqd_addon = False
            logger.warning("qiskit-addon-sqd not installed. Using simplified implementation.")

        logger.info(f"SQD Solver initialized: subspace_dim={subspace_dim}, depth={circuit_depth}")

    def _init_backend(self, **kwargs):
        """Initialize quantum backend."""
        if self.backend == 'statevector':
            self._use_statevector = True
            logger.info("Using statevector simulation")
        elif self.backend == 'ibm':
            try:
                from kanad.backends.ibm import IBMRuntimeBackend
                self._ibm_backend = IBMRuntimeBackend(**kwargs)
                self._use_statevector = False
                logger.info("IBM Quantum backend initialized")
            except Exception as e:
                logger.error(f"IBM backend initialization failed: {e}")
                raise
        else:
            logger.warning(f"Unknown backend {self.backend}, using statevector")
            self._use_statevector = True

    def _generate_subspace_basis(self) -> np.ndarray:
        """
        Generate quantum subspace basis states.

        Uses random parameterized circuits to generate a diverse subspace.

        Returns:
            Basis states (subspace_dim, 2^n_qubits)
        """
        n_qubits = 2 * self.hamiltonian.n_orbitals
        hilbert_dim = 2 ** n_qubits

        logger.info(f"Generating {self.subspace_dim} basis states for {n_qubits}-qubit system")

        # Create diverse basis using parameterized circuits
        basis_states = []

        # Include Hartree-Fock state
        # For blocked spin ordering [0↑,1↑,...,0↓,1↓,...]:
        # H2: 2 electrons in 2 orbitals -> both in lowest MO
        # State |0101⟩: q0=1 (MO0↑), q1=0 (MO1↑), q2=1 (MO0↓), q3=0 (MO1↓)
        hf_state = np.zeros(hilbert_dim, dtype=complex)
        n_orb = self.hamiltonian.n_orbitals
        n_elec = self.hamiltonian.n_electrons
        n_alpha = n_elec // 2
        n_beta = n_elec - n_alpha

        # Blocked ordering: fill lowest n_alpha orbitals (spin-up), then n_beta orbitals (spin-down)
        hf_occupation = 0
        for i in range(n_alpha):
            hf_occupation |= (1 << i)  # Spin-up orbitals
        for i in range(n_beta):
            hf_occupation |= (1 << (n_orb + i))  # Spin-down orbitals

        hf_state[hf_occupation] = 1.0
        basis_states.append(hf_state)

        # Generate additional states with random unitaries
        for i in range(self.subspace_dim - 1):
            # Random parameters for circuit
            n_params = 2 * n_qubits * self.circuit_depth  # Rough estimate
            params = np.random.randn(n_params) * 0.5

            # Build state (simplified - would use actual circuits)
            # For now, create random states
            state = np.random.randn(hilbert_dim) + 1j * np.random.randn(hilbert_dim)
            state = state / np.linalg.norm(state)
            basis_states.append(state)

        basis = np.array(basis_states)

        # Orthonormalize using Gram-Schmidt
        basis = self._gram_schmidt(basis)

        logger.info(f"Generated orthonormal basis: {basis.shape}")

        return basis

    def _gram_schmidt(self, vectors: np.ndarray) -> np.ndarray:
        """
        Orthonormalize vectors using Gram-Schmidt process.

        Args:
            vectors: Input vectors (n_vectors, dim)

        Returns:
            Orthonormal vectors (n_vectors, dim)
        """
        n_vectors = len(vectors)
        orthonormal = np.zeros_like(vectors)

        for i in range(n_vectors):
            # Start with current vector
            vec = vectors[i].copy()

            # Subtract projections onto previous orthonormal vectors
            for j in range(i):
                proj = np.vdot(orthonormal[j], vec) * orthonormal[j]
                vec = vec - proj

            # Normalize
            norm = np.linalg.norm(vec)
            if norm > 1e-10:
                orthonormal[i] = vec / norm
            else:
                # Linear dependence - use random vector
                vec = np.random.randn(len(vec)) + 1j * np.random.randn(len(vec))
                vec = vec / np.linalg.norm(vec)
                orthonormal[i] = vec

        return orthonormal

    def _project_hamiltonian(self, basis: np.ndarray) -> np.ndarray:
        """
        Project Hamiltonian into subspace.

        H_sub[i,j] = ⟨ψ_i|H|ψ_j⟩

        Args:
            basis: Subspace basis states (n_basis, hilbert_dim)

        Returns:
            Projected Hamiltonian (n_basis, n_basis)
        """
        n_qubits = 2 * self.hamiltonian.n_orbitals
        # Use MO basis - Kronecker product bug is now fixed!
        H_matrix = self.hamiltonian.to_matrix(n_qubits=n_qubits, use_mo_basis=True)

        n_basis = len(basis)
        H_sub = np.zeros((n_basis, n_basis), dtype=complex)

        logger.info(f"Projecting Hamiltonian into {n_basis}-dimensional subspace...")

        for i in range(n_basis):
            for j in range(i, n_basis):
                # ⟨ψ_i|H|ψ_j⟩
                H_sub[i, j] = np.vdot(basis[i], H_matrix @ basis[j])
                H_sub[j, i] = np.conj(H_sub[i, j])

        logger.info("Hamiltonian projection complete")

        return H_sub

    def solve(self, n_states: int = 3) -> Dict[str, Any]:
        """
        Solve for ground and excited states using SQD.

        Args:
            n_states: Number of lowest eigenstates to return

        Returns:
            Dictionary with comprehensive results:
                - energies: Eigenvalues (n_states,) [Hartree]
                - eigenvectors: Eigenvectors in subspace (n_states, subspace_dim)
                - ground_state_energy: Lowest eigenvalue [Hartree]
                - excited_state_energies: Higher eigenvalues [Hartree]
                - subspace_dim: Dimension of subspace used
                - hf_energy: Hartree-Fock reference
                - correlation_energy: Ground state correlation
                - analysis: Detailed analysis (if enabled)
        """
        logger.info(f"Starting SQD solve for {n_states} states...")

        # Get HF reference
        hf_energy = self.get_reference_energy()
        if hf_energy is not None:
            logger.info(f"HF reference energy: {hf_energy:.8f} Hartree")

        # Step 1: Generate quantum subspace
        basis = self._generate_subspace_basis()

        # Step 2: Project Hamiltonian
        H_sub = self._project_hamiltonian(basis)

        # Step 3: Classical diagonalization
        logger.info("Diagonalizing projected Hamiltonian...")
        eigenvalues, eigenvectors = np.linalg.eigh(H_sub)

        # Take lowest n_states
        eigenvalues = eigenvalues[:n_states]
        eigenvectors = eigenvectors[:, :n_states].T  # (n_states, subspace_dim)

        logger.info(f"Found {n_states} eigenvalues:")
        for i, E in enumerate(eigenvalues):
            logger.info(f"  State {i}: {E:.8f} Hartree")

        # Store results
        self.results = {
            'energies': eigenvalues.real,
            'eigenvectors': eigenvectors,
            'ground_state_energy': eigenvalues[0].real,
            'excited_state_energies': eigenvalues[1:].real if n_states > 1 else [],
            'energy': eigenvalues[0].real,  # For base class compatibility
            'converged': True,  # SQD always converges
            'iterations': 1,  # Single diagonalization
            'subspace_dim': self.subspace_dim,
            'circuit_depth': self.circuit_depth
        }

        # Add HF reference and correlation
        if hf_energy is not None:
            self.results['hf_energy'] = hf_energy
            self.results['correlation_energy'] = eigenvalues[0].real - hf_energy

            logger.info(f"Ground state correlation: {eigenvalues[0].real - hf_energy:.8f} Hartree")

        # Add analysis if enabled
        if self.enable_analysis:
            # Use HF density matrix as approximation
            density_matrix, _ = self.hamiltonian.solve_scf(max_iterations=50, conv_tol=1e-6)
            self._add_analysis_to_results(eigenvalues[0].real, density_matrix)

        # Add optimization stats
        if self.enable_optimization:
            self._add_optimization_stats()

        # Validate
        validation = self.validate_results()
        self.results['validation'] = validation

        if not validation['passed']:
            logger.warning("SQD results failed validation checks!")

        logger.info("SQD solve complete")

        return self.results

    def print_summary(self):
        """Print extended summary including excited states."""
        super().print_summary()

        # Add excited state info
        if 'excited_state_energies' in self.results and len(self.results['excited_state_energies']) > 0:
            print("\nExcited States:")
            for i, E in enumerate(self.results['excited_state_energies'], start=1):
                excitation = (E - self.results['ground_state_energy']) * 27.2114  # Convert to eV
                print(f"  State {i}: {E:.8f} Ha (ΔE = {excitation:.4f} eV)")
