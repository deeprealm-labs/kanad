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
        random_seed: Optional[int] = None,
        experiment_id: Optional[str] = None,  # For WebSocket broadcasting
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
            random_seed: Random seed for reproducible subspace generation
            **kwargs: Additional backend options
        """
        super().__init__(bond, enable_analysis, enable_optimization)

        self.subspace_dim = subspace_dim
        self.circuit_depth = circuit_depth
        self.backend = backend
        self.shots = shots if shots is not None else 8192  # SQD needs more shots
        self.random_seed = random_seed

        # Store experiment_id for WebSocket broadcasting
        self.experiment_id = experiment_id

        # This is a correlated method
        self._is_correlated = True

        # Set random seed for reproducibility
        if random_seed is not None:
            np.random.seed(random_seed)
            logger.info(f"Random seed set to {random_seed}")

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
        """Initialize quantum backend (uses base class implementation)."""
        super()._init_backend(self.backend, **kwargs)

    def _generate_subspace_basis(self) -> np.ndarray:
        """
        Generate quantum subspace basis states.

        Uses physically meaningful excited determinants (singles, doubles)
        from the HF reference to build a correlation-aware subspace.

        Returns:
            Basis states (subspace_dim, 2^n_qubits)
        """
        n_qubits = 2 * self.hamiltonian.n_orbitals
        hilbert_dim = 2 ** n_qubits
        n_orb = self.hamiltonian.n_orbitals
        n_elec = self.hamiltonian.n_electrons
        n_alpha = n_elec // 2
        n_beta = n_elec - n_alpha

        logger.info(f"Generating {self.subspace_dim} basis states for {n_qubits}-qubit system")

        # Create diverse basis using excited determinants
        basis_states = []

        # 1. Include Hartree-Fock state (most important!)
        # For blocked spin ordering [0↑,1↑,...,0↓,1↓,...]:
        hf_occupation = 0
        for i in range(n_alpha):
            hf_occupation |= (1 << i)  # Spin-up orbitals
        for i in range(n_beta):
            hf_occupation |= (1 << (n_orb + i))  # Spin-down orbitals

        hf_state = np.zeros(hilbert_dim, dtype=complex)
        hf_state[hf_occupation] = 1.0
        basis_states.append(hf_state)
        logger.debug(f"Added HF state: occupation={bin(hf_occupation)}")

        # 2. Add single excitations (capture orbital relaxation)
        single_excitations = []
        for i in range(n_alpha):  # Occupied alpha
            for a in range(n_alpha, n_orb):  # Virtual alpha
                # Alpha single: i→a
                occ = hf_occupation ^ (1 << i) ^ (1 << a)
                single_excitations.append(occ)

        for i in range(n_beta):  # Occupied beta
            for a in range(n_beta, n_orb):  # Virtual beta
                # Beta single: i→a (in beta space: n_orb+i → n_orb+a)
                occ = hf_occupation ^ (1 << (n_orb + i)) ^ (1 << (n_orb + a))
                single_excitations.append(occ)

        # Add single excitations to basis
        for occ in single_excitations[:min(len(single_excitations), self.subspace_dim - 1)]:
            state = np.zeros(hilbert_dim, dtype=complex)
            state[occ] = 1.0
            basis_states.append(state)

        logger.debug(f"Added {min(len(single_excitations), self.subspace_dim - 1)} single excitations")

        # 3. Add double excitations (capture correlation)
        if len(basis_states) < self.subspace_dim:
            double_excitations = []
            for i in range(n_alpha):
                for j in range(i + 1, n_alpha):
                    for a in range(n_alpha, n_orb):
                        for b in range(a + 1, n_orb):
                            # Alpha-alpha double: i,j→a,b
                            occ = hf_occupation ^ (1 << i) ^ (1 << j) ^ (1 << a) ^ (1 << b)
                            double_excitations.append(occ)

            # Alpha-beta doubles (most important for correlation!)
            for i in range(n_alpha):  # Occ alpha
                for j in range(n_beta):  # Occ beta
                    for a in range(n_alpha, n_orb):  # Virt alpha
                        for b in range(n_beta, n_orb):  # Virt beta
                            occ = hf_occupation ^ (1 << i) ^ (1 << (n_orb + j)) ^ (1 << a) ^ (1 << (n_orb + b))
                            double_excitations.append(occ)

            # Add double excitations to fill remaining subspace
            remaining = self.subspace_dim - len(basis_states)
            for occ in double_excitations[:remaining]:
                state = np.zeros(hilbert_dim, dtype=complex)
                state[occ] = 1.0
                basis_states.append(state)

            logger.debug(f"Added {min(len(double_excitations), remaining)} double excitations")

        # 4. If subspace_dim exceeds available determinants, cap it
        max_determinants = len(basis_states)
        actual_dim = min(self.subspace_dim, max_determinants)

        if actual_dim < self.subspace_dim:
            logger.info(f"Subspace auto-adjusted: requested {self.subspace_dim}, "
                       f"using {actual_dim} available determinants (HF + singles + doubles)")

        # 5. If still need more states (shouldn't happen often), add carefully constructed random states
        attempts = 0
        while len(basis_states) < actual_dim and attempts < 100:
            # Create random state in particle-conserving subspace
            state = np.zeros(hilbert_dim, dtype=complex)
            # Randomly weight existing determinants (this preserves particle number)
            weights = np.random.randn(len(basis_states)) + 1j * np.random.randn(len(basis_states))
            for i, bs in enumerate(basis_states):
                state += weights[i] * bs
            state = state / np.linalg.norm(state)

            # Check if linearly independent
            if len(basis_states) > 0:
                overlap = max(abs(np.vdot(bs, state)) for bs in basis_states)
                if overlap < 0.99:  # Not too similar to existing states
                    basis_states.append(state)
            else:
                basis_states.append(state)
            attempts += 1

        basis = np.array(basis_states[:actual_dim])

        # Orthonormalize using Gram-Schmidt
        basis = self._gram_schmidt(basis)

        logger.info(f"Generated orthonormal basis: {basis.shape} (HF + {len(single_excitations)} singles + {len(double_excitations) if 'double_excitations' in locals() else 0} doubles)")

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
        hilbert_dim = 2 ** n_qubits
        n_basis = len(basis)

        logger.info(f"Projecting Hamiltonian into {n_basis}-dimensional subspace...")

        # IMPORTANT: Use dense matrix construction for projection
        # SparsePauliOp.to_matrix() has qubit ordering issues that cause wrong eigenvalues
        # For small systems (< 8 qubits), dense matrix is fast and correct
        if hilbert_dim <= 256:  # 8 qubits or less
            logger.info(f"Using dense Hamiltonian matrix ({hilbert_dim}×{hilbert_dim}) for accurate projection")
            H_matrix = self.hamiltonian.to_matrix(n_qubits=n_qubits, use_mo_basis=True)
            H_sub = np.zeros((n_basis, n_basis), dtype=complex)

            for i in range(n_basis):
                for j in range(i, n_basis):
                    H_sub[i, j] = np.vdot(basis[i], H_matrix @ basis[j])
                    H_sub[j, i] = np.conj(H_sub[i, j])

                    # Progress logging
                    if n_basis > 5 and (i * n_basis + j) % 10 == 0:
                        progress = ((i * n_basis + j) / (n_basis * (n_basis + 1) // 2)) * 100
                        logger.debug(f"Projection progress: {progress:.1f}%")
        else:
            # For large systems, warn and use sparse (may have accuracy issues)
            logger.warning(f"Large system detected ({n_qubits} qubits, {hilbert_dim}D Hilbert space)")
            logger.warning(f"⚠️  SQD may not work correctly for systems > 8 qubits due to sparse Hamiltonian issues")
            logger.warning(f"⚠️  Consider using VQE instead")

            H_matrix = self.hamiltonian.to_matrix(n_qubits=n_qubits, use_mo_basis=True)
            H_sub = np.zeros((n_basis, n_basis), dtype=complex)

            for i in range(n_basis):
                for j in range(i, n_basis):
                    H_sub[i, j] = np.vdot(basis[i], H_matrix @ basis[j])
                    H_sub[j, i] = np.conj(H_sub[i, j])

        logger.info("Hamiltonian projection complete")
        return H_sub

    def solve(self, n_states: int = 3, callback=None) -> Dict[str, Any]:
        """
        Solve for ground and excited states using SQD.

        Args:
            n_states: Number of lowest eigenstates to return
            callback: Optional callback function(stage: int, energy: float, message: str)
                     Called at different stages: 0=init, 1=basis, 2=projection, 3=diag, 4+=states

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
        n_qubits = 2 * self.hamiltonian.n_orbitals
        hilbert_dim = 2 ** n_qubits

        logger.info(f"Starting SQD solve for {n_states} states...")
        logger.info(f"System size: {n_qubits} qubits, Hilbert space: {hilbert_dim}D")

        # Warn and auto-adjust for large systems
        if hilbert_dim > 256:
            logger.warning(f"⚠️  Large system detected! SQD may be very slow.")
            logger.warning(f"⚠️  Consider using VQE instead for systems > 8 qubits.")

            # Auto-reduce subspace dimension for large systems
            original_subspace_dim = self.subspace_dim
            if self.subspace_dim > 4:
                self.subspace_dim = min(4, self.subspace_dim)
                logger.warning(f"⚠️  Auto-reducing subspace_dim: {original_subspace_dim} → {self.subspace_dim}")

        # Get HF reference
        hf_energy = self.get_reference_energy()
        if hf_energy is not None:
            logger.info(f"HF reference energy: {hf_energy:.8f} Hartree")
            if callback:
                callback(0, hf_energy, "HF reference computed")

        # Step 1: Generate quantum subspace
        basis = self._generate_subspace_basis()
        if callback:
            callback(1, hf_energy if hf_energy else 0.0, f"Subspace basis generated ({len(basis)} states)")

        # Step 2: Project Hamiltonian
        H_sub = self._project_hamiltonian(basis)
        if callback:
            callback(2, hf_energy if hf_energy else 0.0, "Hamiltonian projection complete")

        # Step 3: Classical diagonalization
        logger.info("Diagonalizing projected Hamiltonian...")
        if callback:
            callback(3, hf_energy if hf_energy else 0.0, "Diagonalizing Hamiltonian")
        eigenvalues, eigenvectors = np.linalg.eigh(H_sub)

        # Take lowest n_states
        eigenvalues = eigenvalues[:n_states]
        eigenvectors = eigenvectors[:, :n_states].T  # (n_states, subspace_dim)

        logger.info(f"Found {n_states} eigenvalues:")
        for i, E in enumerate(eigenvalues):
            logger.info(f"  State {i}: {E:.8f} Hartree")
            if callback:
                callback(4 + i, float(E), f"State {i} computed")

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

        # ADD ENHANCED DATA FOR ANALYSIS SERVICE
        try:
            # Store molecule geometry for ADME and other analyses
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
                logger.info(f"✅ Stored molecule geometry for analysis")

            # Store nuclear repulsion energy
            if hasattr(self.hamiltonian, 'nuclear_repulsion'):
                self.results['nuclear_repulsion'] = float(self.hamiltonian.nuclear_repulsion)

            # Try to get density matrix from HF reference
            try:
                if hasattr(self.hamiltonian, 'mf'):
                    # PySCF mean-field object has density matrix
                    if hasattr(self.hamiltonian.mf, 'make_rdm1'):
                        rdm1 = self.hamiltonian.mf.make_rdm1()
                        self.results['rdm1'] = rdm1.tolist()
                        logger.info(f"✅ Stored RDM1 for bonding analysis")
            except Exception as e:
                logger.warning(f"Could not extract RDM1: {e}")

            # Try to get orbital energies
            try:
                logger.info(f"🔍 Checking hamiltonian for orbital energies: hasattr(mf)={hasattr(self.hamiltonian, 'mf')}")
                if hasattr(self.hamiltonian, 'mf'):
                    logger.info(f"🔍 Hamiltonian has mf attribute, checking mo_energy: hasattr(mo_energy)={hasattr(self.hamiltonian.mf, 'mo_energy')}")
                    if hasattr(self.hamiltonian.mf, 'mo_energy'):
                        orb_energies = self.hamiltonian.mf.mo_energy
                        logger.info(f"🔍 Found orbital energies: shape={orb_energies.shape}, dtype={orb_energies.dtype}")
                        self.results['orbital_energies'] = orb_energies.tolist()
                        logger.info(f"✅ Stored orbital energies for DOS analysis")
                    else:
                        logger.warning(f"⚠️  Hamiltonian.mf does not have mo_energy attribute")
                else:
                    logger.warning(f"⚠️  Hamiltonian does not have mf attribute (type: {type(self.hamiltonian).__name__})")
            except Exception as e:
                logger.warning(f"Could not extract orbital energies: {e}")

            # Try to get dipole moment
            try:
                if hasattr(self.hamiltonian, 'mf'):
                    from pyscf import scf
                    if hasattr(scf, 'hf') and hasattr(scf.hf, 'dip_moment'):
                        dipole = scf.hf.dip_moment(self.hamiltonian.mf.mol, self.hamiltonian.mf.make_rdm1())
                        self.results['dipole'] = dipole.tolist()
                        logger.info(f"✅ Stored dipole moment")
            except Exception as e:
                logger.warning(f"Could not calculate dipole: {e}")

        except Exception as e:
            logger.error(f"Error storing enhanced data: {e}")

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
