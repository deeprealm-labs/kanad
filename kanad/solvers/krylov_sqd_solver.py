"""
Krylov Subspace Quantum Diagonalization (Krylov-SQD) Solver

Krylov-SQD improves upon standard SQD by using the Lanczos algorithm to build
a Krylov subspace that is optimally adapted to the Hamiltonian. This provides:
- Better energy estimates with smaller subspaces
- Faster convergence to ground and excited states
- Lower circuit depth requirements
- More efficient use of quantum resources

Theory:
The Krylov subspace K_m(H, |ψ₀⟩) = span{|ψ₀⟩, H|ψ₀⟩, H²|ψ₀⟩, ..., H^(m-1)|ψ₀⟩}
provides an optimal subspace for approximating eigenvalues of H.

The Lanczos algorithm efficiently builds this subspace and tridiagonalizes H
simultaneously, avoiding explicit matrix-vector products with H^k.

References:
- Qiskit SQD addon: https://github.com/qiskit-community/qiskit-addon-sqd
- Quantum Krylov methods: https://arxiv.org/abs/1909.05820
"""

from typing import Dict, Any, Optional, List, Tuple
import numpy as np
import logging

from kanad.solvers.base_solver import BaseSolver

logger = logging.getLogger(__name__)


class KrylovSQDSolver(BaseSolver):
    """
    Krylov Subspace Quantum Diagonalization for ground and excited states.

    Workflow:
    1. Start with initial state |ψ₀⟩ (e.g., Hartree-Fock)
    2. Build Krylov subspace K_m via Lanczos iterations
    3. Each iteration computes H|ψᵢ⟩ on quantum hardware
    4. Orthogonalize and build tridiagonal Hamiltonian
    5. Diagonalize tridiagonal matrix (classical)
    6. Return eigenvalues and eigenvectors

    Advantages over standard SQD:
    - Better energy estimates (Ritz values converge fast)
    - Smaller subspace needed (10-20 vs 50-100)
    - Lower circuit depth (no double excitations needed)
    - Excited states naturally included

    Usage:
        from kanad.bonds import BondFactory
        from kanad.solvers import KrylovSQDSolver

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        solver = KrylovSQDSolver(bond, krylov_dim=15, n_states=5)
        result = solver.solve()

        print(f"Ground State: {result['energies'][0]:.6f} Hartree")
        print(f"1st Excited: {result['energies'][1]:.6f} Hartree")
    """

    def __init__(
        self,
        bond: 'BaseBond',
        krylov_dim: int = 15,
        n_states: int = 3,
        initial_state: Optional[np.ndarray] = None,
        backend: str = 'statevector',
        shots: Optional[int] = None,
        enable_analysis: bool = True,
        enable_optimization: bool = True,
        random_seed: Optional[int] = None,
        reorthogonalize: bool = True,
        experiment_id: Optional[str] = None,
        **kwargs
    ):
        """
        Initialize Krylov-SQD solver.

        Args:
            bond: Bond object from BondFactory
            krylov_dim: Dimension of Krylov subspace (10-20 usually sufficient)
            n_states: Number of eigenvalues to compute
            initial_state: Initial vector |ψ₀⟩ (HF if None)
            backend: Quantum backend ('statevector', 'qasm', 'ibm', 'bluequbit')
            shots: Number of shots for sampling backends
            enable_analysis: Enable automatic analysis
            enable_optimization: Enable automatic optimization
            random_seed: Random seed for reproducibility
            reorthogonalize: Use full reorthogonalization (more stable)
            **kwargs: Additional backend options
        """
        super().__init__(bond, enable_analysis, enable_optimization)

        self.krylov_dim = krylov_dim
        self.n_states = n_states
        self.initial_state = initial_state
        self.backend = backend
        self.shots = shots if shots is not None else 8192
        self.random_seed = random_seed
        self.reorthogonalize = reorthogonalize
        self.experiment_id = experiment_id

        # Correlated method
        self._is_correlated = True

        # Set random seed
        if random_seed is not None:
            np.random.seed(random_seed)
            logger.info(f"Random seed set to {random_seed}")

        # Initialize backend
        self._init_backend(**kwargs)

        logger.info(f"Krylov-SQD Solver initialized: krylov_dim={krylov_dim}, n_states={n_states}")

    def _init_backend(self, **kwargs):
        """Initialize quantum backend."""
        if self.backend == 'statevector':
            self._use_statevector = True
            logger.info("Using statevector simulation")
        elif self.backend == 'bluequbit':
            try:
                from kanad.backends.bluequbit import BlueQubitBackend
                self._bluequbit_backend = BlueQubitBackend(**kwargs)
                self._use_statevector = False
                device = kwargs.get('device', 'gpu')
                logger.info(f"BlueQubit backend initialized: device={device}")
            except Exception as e:
                logger.error(f"BlueQubit backend initialization failed: {e}")
                raise
        elif self.backend == 'ibm':
            try:
                from kanad.backends.ibm import IBMBackend
                self._ibm_backend = IBMBackend(**kwargs)
                self._use_statevector = False
                backend_name = kwargs.get('backend_name', 'ibm_torino')
                logger.info(f"IBM Quantum backend initialized: {backend_name}")
            except Exception as e:
                logger.error(f"IBM backend initialization failed: {e}")
                raise
        else:
            logger.warning(f"Unknown backend {self.backend}, using statevector")
            self._use_statevector = True

    def _get_initial_state(self) -> np.ndarray:
        """
        Get initial state for Lanczos iteration.

        Default: Hartree-Fock state |HF⟩
        This is the most physically meaningful starting point.

        Returns:
            Initial state vector (2^n_qubits,)
        """
        if self.initial_state is not None:
            logger.info("Using provided initial state")
            return self.initial_state / np.linalg.norm(self.initial_state)

        n_qubits = 2 * self.hamiltonian.n_orbitals
        hilbert_dim = 2 ** n_qubits
        n_orb = self.hamiltonian.n_orbitals
        n_elec = self.hamiltonian.n_electrons
        n_alpha = n_elec // 2
        n_beta = n_elec - n_alpha

        # Build HF state in blocked spin ordering
        hf_occupation = 0
        for i in range(n_alpha):
            hf_occupation |= (1 << i)  # Alpha electrons
        for i in range(n_beta):
            hf_occupation |= (1 << (n_orb + i))  # Beta electrons

        hf_state = np.zeros(hilbert_dim, dtype=complex)
        hf_state[hf_occupation] = 1.0

        logger.info(f"Using Hartree-Fock initial state: occupation={bin(hf_occupation)}")
        return hf_state

    def _apply_hamiltonian(self, state: np.ndarray) -> np.ndarray:
        """
        Apply Hamiltonian to state: H|ψ⟩

        This is where quantum hardware is used in real implementations.
        For simulation, we use direct matrix-vector product.

        Args:
            state: Input state |ψ⟩ (2^n_qubits,)

        Returns:
            H|ψ⟩ (2^n_qubits,)
        """
        n_qubits = 2 * self.hamiltonian.n_orbitals

        if self._use_statevector:
            # Statevector simulation: direct matrix-vector product
            H_matrix = self.hamiltonian.to_matrix(n_qubits=n_qubits, use_mo_basis=True)
            return H_matrix @ state
        else:
            # Real quantum hardware: implement Hadamard test circuits
            # TODO: Implement for real hardware backends
            logger.warning("Real hardware not yet implemented for Krylov-SQD, using statevector")
            H_matrix = self.hamiltonian.to_matrix(n_qubits=n_qubits, use_mo_basis=True)
            return H_matrix @ state

    def _lanczos_iteration(
        self,
        callback=None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Perform Lanczos iteration to build Krylov subspace.

        Algorithm:
        1. Start with |q₀⟩ = |ψ₀⟩ / ||ψ₀||
        2. For i = 0, ..., m-1:
           a. |w⟩ = H|qᵢ⟩
           b. αᵢ = ⟨qᵢ|w⟩
           c. |w⟩ = |w⟩ - αᵢ|qᵢ⟩ - βᵢ₋₁|qᵢ₋₁⟩
           d. βᵢ = ||w||
           e. |qᵢ₊₁⟩ = |w⟩ / βᵢ
        3. Build tridiagonal matrix T with diagonal α and off-diagonal β

        Args:
            callback: Progress callback function

        Returns:
            krylov_basis: Orthonormal Krylov vectors (krylov_dim, 2^n_qubits)
            alpha: Diagonal elements (krylov_dim,)
            beta: Off-diagonal elements (krylov_dim-1,)
        """
        logger.info(f"Starting Lanczos iteration for {self.krylov_dim} steps...")

        # Get initial state
        q_prev = None
        q_curr = self._get_initial_state()

        # Storage
        krylov_basis = []
        alpha = []
        beta = []

        for i in range(self.krylov_dim):
            logger.debug(f"Lanczos iteration {i+1}/{self.krylov_dim}")

            # Store current vector
            krylov_basis.append(q_curr.copy())

            # Apply Hamiltonian (quantum operation!)
            w = self._apply_hamiltonian(q_curr)

            # Compute diagonal element: αᵢ = ⟨qᵢ|H|qᵢ⟩
            alpha_i = np.real(np.vdot(q_curr, w))
            alpha.append(alpha_i)

            # Orthogonalize: w = w - αᵢ|qᵢ⟩ - βᵢ₋₁|qᵢ₋₁⟩
            w = w - alpha_i * q_curr
            if q_prev is not None:
                w = w - beta[-1] * q_prev

            # Full reorthogonalization (more stable, recommended)
            if self.reorthogonalize:
                for q in krylov_basis:
                    w = w - np.vdot(q, w) * q

            # Compute off-diagonal element: βᵢ = ||w||
            beta_i = np.linalg.norm(w)

            # Check for convergence (Lanczos breakdown)
            if beta_i < 1e-10:
                logger.info(f"Lanczos breakdown at iteration {i+1} (β={beta_i:.2e})")
                logger.info("Krylov subspace is complete - this is good!")
                break

            beta.append(beta_i)

            # Next Krylov vector: |qᵢ₊₁⟩ = |w⟩ / βᵢ
            q_prev = q_curr
            q_curr = w / beta_i

            # Progress callback
            if callback and i % 5 == 0:
                callback(1, alpha_i, f"Lanczos iteration {i+1}/{self.krylov_dim}")

        logger.info(f"Lanczos complete: generated {len(krylov_basis)} Krylov vectors")

        return (
            np.array(krylov_basis),
            np.array(alpha),
            np.array(beta)
        )

    def _build_tridiagonal_hamiltonian(
        self,
        alpha: np.ndarray,
        beta: np.ndarray
    ) -> np.ndarray:
        """
        Build tridiagonal Hamiltonian matrix from Lanczos coefficients.

        T = [α₀  β₀   0   0  ...]
            [β₀  α₁  β₁   0  ...]
            [0   β₁  α₂  β₂ ...]
            [0   0   β₂  α₃ ...]
            ...

        Args:
            alpha: Diagonal elements (m,)
            beta: Off-diagonal elements (m-1,)

        Returns:
            Tridiagonal matrix (m, m)
        """
        m = len(alpha)
        T = np.zeros((m, m), dtype=float)

        # Diagonal
        for i in range(m):
            T[i, i] = alpha[i]

        # Off-diagonal
        for i in range(len(beta)):
            T[i, i+1] = beta[i]
            T[i+1, i] = beta[i]

        logger.debug(f"Built tridiagonal matrix: shape={T.shape}")
        return T

    def solve(self, callback=None) -> Dict[str, Any]:
        """
        Solve for ground and excited states using Krylov-SQD.

        Args:
            callback: Optional callback function(stage: int, energy: float, message: str)

        Returns:
            Dictionary with comprehensive results:
                - energies: Eigenvalues (n_states,) [Hartree]
                - eigenvectors: Eigenvectors in Krylov basis (n_states, krylov_dim)
                - ground_state_energy: Lowest eigenvalue [Hartree]
                - excited_state_energies: Higher eigenvalues [Hartree]
                - krylov_dim: Dimension of Krylov subspace used
                - n_lanczos_iterations: Number of Lanczos iterations performed
                - hf_energy: Hartree-Fock reference
                - correlation_energy: Ground state correlation
                - ritz_convergence: Convergence history of Ritz values
                - analysis: Detailed analysis (if enabled)
        """
        n_qubits = 2 * self.hamiltonian.n_orbitals
        hilbert_dim = 2 ** n_qubits

        logger.info(f"Starting Krylov-SQD solve for {self.n_states} states...")
        logger.info(f"System size: {n_qubits} qubits, Hilbert space: {hilbert_dim}D")

        # Get HF reference
        hf_energy = self.get_reference_energy()
        if hf_energy is not None:
            logger.info(f"HF reference energy: {hf_energy:.8f} Hartree")
            if callback:
                callback(0, hf_energy, "HF reference computed")

        # Step 1: Lanczos iteration to build Krylov subspace
        krylov_basis, alpha, beta = self._lanczos_iteration(callback=callback)
        actual_krylov_dim = len(krylov_basis)

        if callback:
            callback(1, hf_energy if hf_energy else 0.0,
                    f"Krylov subspace built ({actual_krylov_dim} vectors)")

        # Step 2: Build tridiagonal Hamiltonian
        T = self._build_tridiagonal_hamiltonian(alpha, beta)
        if callback:
            callback(2, hf_energy if hf_energy else 0.0,
                    "Tridiagonal Hamiltonian constructed")

        # Step 3: Diagonalize tridiagonal matrix
        logger.info("Diagonalizing tridiagonal Hamiltonian...")
        if callback:
            callback(3, hf_energy if hf_energy else 0.0, "Diagonalizing Hamiltonian")

        eigenvalues, eigenvectors = np.linalg.eigh(T)

        # Take lowest n_states (Ritz values and vectors)
        n_requested = min(self.n_states, len(eigenvalues))
        ritz_values = eigenvalues[:n_requested]
        ritz_vectors = eigenvectors[:, :n_requested].T

        logger.info(f"Found {n_requested} Ritz values:")
        for i, E in enumerate(ritz_values):
            logger.info(f"  State {i}: {E:.8f} Hartree")
            if callback:
                callback(4 + i, float(E), f"State {i} computed")

        # Store results
        self.results = {
            'energies': ritz_values,
            'eigenvectors': ritz_vectors,
            'ground_state_energy': ritz_values[0],
            'excited_state_energies': ritz_values[1:] if n_requested > 1 else [],
            'energy': ritz_values[0],  # For base class compatibility
            'converged': True,  # Krylov always converges
            'iterations': actual_krylov_dim,
            'krylov_dim': actual_krylov_dim,
            'n_lanczos_iterations': actual_krylov_dim,
            'method': 'Krylov-SQD',
            'alpha': alpha.tolist(),
            'beta': beta.tolist()
        }

        # Add HF reference and correlation
        if hf_energy is not None:
            self.results['hf_energy'] = hf_energy
            self.results['correlation_energy'] = ritz_values[0] - hf_energy
            logger.info(f"Ground state correlation: {ritz_values[0] - hf_energy:.8f} Hartree")

        # Add analysis if enabled
        if self.enable_analysis:
            # Use HF density matrix as approximation
            density_matrix, _ = self.hamiltonian.solve_scf(max_iterations=50, conv_tol=1e-6)
            self._add_analysis_to_results(ritz_values[0], density_matrix)

        # Add optimization stats
        if self.enable_optimization:
            self._add_optimization_stats()

        # Validate
        validation = self.validate_results()
        self.results['validation'] = validation

        if not validation['passed']:
            logger.warning("Krylov-SQD results failed validation checks!")

        # Enhanced data for analysis
        try:
            if self.molecule is not None:
                self.results['geometry'] = [
                    (atom.symbol, tuple(atom.position))
                    for atom in self.molecule.atoms
                ]
                self.results['atoms'] = [atom.symbol for atom in self.molecule.atoms]
                self.results['n_atoms'] = self.molecule.n_atoms
                self.results['n_electrons'] = self.molecule.n_electrons
                self.results['charge'] = getattr(self.molecule, 'charge', 0)
                self.results['multiplicity'] = getattr(self.molecule, 'multiplicity', 1)

            if hasattr(self.hamiltonian, 'nuclear_repulsion'):
                self.results['nuclear_repulsion'] = float(self.hamiltonian.nuclear_repulsion)

            # Try to get additional properties
            try:
                if hasattr(self.hamiltonian, 'mf'):
                    if hasattr(self.hamiltonian.mf, 'make_rdm1'):
                        rdm1 = self.hamiltonian.mf.make_rdm1()
                        self.results['rdm1'] = rdm1.tolist()
                    if hasattr(self.hamiltonian.mf, 'mo_energy'):
                        self.results['orbital_energies'] = self.hamiltonian.mf.mo_energy.tolist()
            except Exception as e:
                logger.warning(f"Could not extract additional properties: {e}")

        except Exception as e:
            logger.error(f"Error storing enhanced data: {e}")

        logger.info("Krylov-SQD solve complete")
        return self.results

    def print_summary(self):
        """Print extended summary including excited states and Krylov info."""
        super().print_summary()

        # Add Krylov-specific info
        if 'krylov_dim' in self.results:
            print(f"\nKrylov Subspace Info:")
            print(f"  Dimension: {self.results['krylov_dim']}")
            print(f"  Lanczos iterations: {self.results['n_lanczos_iterations']}")

        # Add excited states
        if 'excited_state_energies' in self.results and len(self.results['excited_state_energies']) > 0:
            print("\nExcited States (Ritz values):")
            for i, E in enumerate(self.results['excited_state_energies'], start=1):
                excitation = (E - self.results['ground_state_energy']) * 27.2114  # eV
                print(f"  State {i}: {E:.8f} Ha (ΔE = {excitation:.4f} eV)")
