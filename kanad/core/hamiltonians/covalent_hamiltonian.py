"""
Covalent Hamiltonian for orbital hybridization systems.

Models covalent bonding via hybrid orbitals and molecular orbital formation.
"""

from typing import List, Dict, Tuple, Optional
import numpy as np
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.core.atom import Atom
from kanad.core.integrals.basis_sets import BasisSet
from kanad.core.integrals.overlap import OverlapIntegrals
from kanad.core.integrals.one_electron import OneElectronIntegrals
from kanad.core.integrals.two_electron import TwoElectronIntegrals


class CovalentHamiltonian(MolecularHamiltonian):
    """
    Hamiltonian for covalent bonding.

    Physical Model:
        H_covalent = Σ_μν h_μν c†_μ c_ν + ½ Σ_μνλσ (μν|λσ) c†_μ c†_ν c_σ c_λ

    where μ,ν run over atomic or hybrid orbitals.

    KEY PHYSICS:
        - Orbital overlap → bonding/antibonding splitting
        - Hybridization (sp, sp², sp³)
        - Shared electron pairs
        - Bond order from MO occupation
    """

    def __init__(
        self,
        molecule: 'Molecule',
        representation: 'LCAORepresentation',
        basis_name: str = 'sto-3g'
    ):
        """
        Initialize covalent Hamiltonian.

        Args:
            molecule: Molecule object
            representation: LCAO representation with hybridization
            basis_name: Basis set name
        """
        self.molecule = molecule
        self.representation = representation
        self.atoms = molecule.atoms
        self.basis_name = basis_name

        # Build basis set
        self.basis = BasisSet(basis_name)
        self.basis.build_basis(self.atoms)

        # Compute nuclear repulsion
        nuclear_rep = self._compute_nuclear_repulsion()

        super().__init__(
            n_orbitals=self.basis.n_basis_functions,
            n_electrons=molecule.n_electrons,
            nuclear_repulsion=nuclear_rep
        )

        # Build Hamiltonian
        self._build_hamiltonian()

    def _compute_nuclear_repulsion(self) -> float:
        """
        Compute nuclear-nuclear repulsion energy in atomic units.

        E_nn = Σ_{i<j} Z_i Z_j / |R_i - R_j|

        Returns:
            Nuclear repulsion energy in Hartree
        """
        from kanad.core.constants.conversion_factors import ConversionFactors

        energy = 0.0
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                Z_i = self.atoms[i].atomic_number
                Z_j = self.atoms[j].atomic_number
                r_ij_angstrom = self.atoms[i].distance_to(self.atoms[j])
                if r_ij_angstrom > 1e-10:
                    # Convert distance to Bohr for atomic units
                    r_ij_bohr = r_ij_angstrom * ConversionFactors.ANGSTROM_TO_BOHR
                    energy += Z_i * Z_j / r_ij_bohr

        return energy

    def _build_hamiltonian(self):
        """
        Build covalent Hamiltonian using full integral calculation.

        Uses:
        - Overlap integrals
        - Kinetic energy integrals
        - Nuclear attraction integrals
        - Electron repulsion integrals
        """
        # Compute one-electron integrals
        one_electron_ints = OneElectronIntegrals(self.atoms, self.basis.basis_functions)

        # Core Hamiltonian = T + V_ne
        T = one_electron_ints.compute_kinetic()
        V = one_electron_ints.compute_nuclear_attraction()
        self.h_core = T + V

        # Two-electron integrals
        two_electron_ints = TwoElectronIntegrals(self.basis.basis_functions)
        self.eri = two_electron_ints.compute_eri_tensor()

        # Overlap matrix (for analysis)
        self.S = OverlapIntegrals.build_overlap_matrix(self.basis.basis_functions)

    def to_matrix(self) -> np.ndarray:
        """
        Convert to matrix form (one-body part only).

        Returns:
            Core Hamiltonian matrix
        """
        return self.h_core.copy()

    def compute_energy(self, density_matrix: np.ndarray) -> float:
        """
        Compute total energy from density matrix.

        E = Σ_μν P_μν h_μν + ½ Σ_μνλσ P_μν P_λσ [(μν|λσ) - ½(μλ|νσ)] + E_nn

        Args:
            density_matrix: One-particle density matrix

        Returns:
            Total electronic energy (Hartree)
        """
        # One-electron contribution
        E_core = np.sum(density_matrix * self.h_core)

        # Two-electron contribution
        E_ee = 0.0
        for i in range(self.n_orbitals):
            for j in range(self.n_orbitals):
                for k in range(self.n_orbitals):
                    for l in range(self.n_orbitals):
                        # Coulomb
                        E_ee += 0.5 * density_matrix[i, j] * density_matrix[k, l] * self.eri[i, j, k, l]
                        # Exchange (closed-shell)
                        E_ee -= 0.25 * density_matrix[i, k] * density_matrix[j, l] * self.eri[i, j, k, l]

        # Total energy
        E_total = E_core + E_ee + self.nuclear_repulsion

        return E_total

    def compute_molecular_orbitals(self, use_hf=True) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute molecular orbitals.

        If use_hf=True (default), solves Hartree-Fock equations self-consistently.
        If use_hf=False, just diagonalizes core Hamiltonian (for testing).

        Returns:
            (mo_energies, mo_coefficients)
        """
        if use_hf:
            energies, coeffs, converged, iterations = self._solve_hartree_fock()
            return energies, coeffs
        else:
            # Just diagonalize core Hamiltonian (not physically meaningful for energies)
            from scipy.linalg import eigh
            energies, coefficients = eigh(self.h_core, self.S)
            return energies, coefficients

    def solve_scf(
        self,
        max_iterations=100,
        conv_tol=1e-8,
        use_diis=True,
        level_shift=0.0,
        damping_factor=0.0
    ) -> Tuple[np.ndarray, float]:
        """
        Solve Hartree-Fock self-consistently and return density matrix and energy.

        This is the public interface for SCF calculations.

        Args:
            max_iterations: Maximum SCF iterations
            conv_tol: Convergence tolerance for energy
            use_diis: Use DIIS convergence acceleration
            level_shift: Level shift for virtual orbitals (Ha) - helps difficult convergence
            damping_factor: Density damping (0-1) - helps oscillatory convergence

        Returns:
            (density_matrix, total_energy)
        """
        from kanad.core.scf_solver import SCFSolver

        # Create SCF solver
        solver = SCFSolver(
            h_core=self.h_core,
            S=self.S,
            eri=self.eri,
            n_electrons=self.n_electrons,
            nuclear_repulsion=self.nuclear_repulsion
        )

        # Solve SCF
        density_matrix, mo_energies, total_energy, converged, iterations = solver.solve(
            max_iterations=max_iterations,
            energy_tol=conv_tol,
            use_diis=use_diis,
            level_shift=level_shift,
            damping_factor=damping_factor
        )

        # Store convergence info for bond classes to access
        self._scf_converged = converged
        self._scf_iterations = iterations
        self._mo_energies = mo_energies

        return density_matrix, total_energy

    def _solve_hartree_fock(self, max_iter=100, conv_tol=1e-8, use_diis=True) -> Tuple[np.ndarray, np.ndarray, bool, int]:
        """
        Solve Hartree-Fock equations self-consistently.

        Restricted Hartree-Fock for closed-shell systems.

        Returns:
            (mo_energies, mo_coefficients, converged, iterations)
        """
        from kanad.core.scf_solver import SCFSolver

        # Create SCF solver
        solver = SCFSolver(
            h_core=self.h_core,
            S=self.S,
            eri=self.eri,
            n_electrons=self.n_electrons,
            nuclear_repulsion=self.nuclear_repulsion
        )

        # Solve SCF
        density_matrix, mo_energies, total_energy, converged, iterations = solver.solve(
            max_iterations=max_iter,
            energy_tol=conv_tol,
            use_diis=use_diis
        )

        # Compute MO coefficients from final diagonalization
        from scipy.linalg import eigh
        F = solver._build_fock_matrix(density_matrix)
        mo_energies, C = eigh(F, self.S)

        if not converged:
            print(f"Warning: HF did not converge in {max_iter} iterations")

        return mo_energies, C, converged, iterations

    def get_bonding_antibonding_split(self) -> Dict[str, float]:
        """
        Compute bonding/antibonding energy splitting.

        For a simple diatomic like H2:
            E_bonding = (h_aa + h_bb - 2S*h_ab) / (2(1 - S²))
            E_antibonding = (h_aa + h_bb + 2S*h_ab) / (2(1 + S²))

        Returns:
            Dictionary with bonding/antibonding info
        """
        energies, coeffs = self.compute_molecular_orbitals()

        # Bonding MOs have lower energies
        n_occ = self.n_electrons // 2

        bonding_energies = energies[:n_occ]
        antibonding_energies = energies[n_occ:]

        if len(antibonding_energies) > 0:
            splitting = antibonding_energies[0] - bonding_energies[-1]
        else:
            splitting = 0.0

        return {
            'bonding_energies': bonding_energies,
            'antibonding_energies': antibonding_energies,
            'homo_lumo_gap': splitting,
            'homo_energy': bonding_energies[-1],
            'lumo_energy': antibonding_energies[0] if len(antibonding_energies) > 0 else np.nan
        }

    def compute_bond_order(self, density_matrix: np.ndarray, atom_i: int, atom_j: int) -> float:
        """
        Compute bond order between two atoms using Mulliken population analysis.

        BO_ij = Σ_μ∈i Σ_ν∈j P_μν S_μν

        Args:
            density_matrix: Density matrix
            atom_i: Index of first atom
            atom_j: Index of second atom

        Returns:
            Bond order
        """
        # Get basis function indices for each atom
        # Simplified: assume contiguous basis functions per atom
        orbitals_per_atom = self.n_orbitals // len(self.atoms)
        start_i = atom_i * orbitals_per_atom
        end_i = (atom_i + 1) * orbitals_per_atom
        start_j = atom_j * orbitals_per_atom
        end_j = (atom_j + 1) * orbitals_per_atom

        bond_order = 0.0
        for mu in range(start_i, end_i):
            for nu in range(start_j, end_j):
                bond_order += density_matrix[mu, nu] * self.S[mu, nu]

        return abs(bond_order)

    def get_overlap_matrix(self) -> np.ndarray:
        """Get overlap matrix."""
        return self.S.copy()

    def get_mo_energies(self) -> np.ndarray:
        """
        Get molecular orbital energies.

        Returns:
            Array of MO energies (sorted)
        """
        energies, _ = self.compute_molecular_orbitals()
        return energies

    def get_homo_lumo_gap(self) -> float:
        """
        Compute HOMO-LUMO gap.

        Returns:
            Gap in Hartree
        """
        energies = self.get_mo_energies()
        n_occ = self.n_electrons // 2

        if n_occ < len(energies):
            gap = energies[n_occ] - energies[n_occ - 1]
            return gap
        else:
            return 0.0

    def compute_overlap_population(
        self,
        density_matrix: np.ndarray,
        mu: int,
        nu: int
    ) -> float:
        """
        Compute overlap population between orbitals μ and ν.

        OP_μν = P_μν S_μν

        Args:
            density_matrix: Density matrix
            mu: Orbital index
            nu: Orbital index

        Returns:
            Overlap population
        """
        return density_matrix[mu, nu] * self.S[mu, nu]

    def analyze_bonding(self, density_matrix: np.ndarray) -> Dict:
        """
        Comprehensive bonding analysis.

        Args:
            density_matrix: Density matrix

        Returns:
            Dictionary with bonding analysis
        """
        analysis = {}

        # MO energies
        energies, coeffs = self.compute_molecular_orbitals()
        analysis['mo_energies'] = energies
        analysis['mo_coefficients'] = coeffs

        # HOMO-LUMO gap
        analysis['homo_lumo_gap'] = self.get_homo_lumo_gap()

        # Bonding/antibonding splitting
        analysis['bonding_analysis'] = self.get_bonding_antibonding_split()

        # Bond orders (for all atom pairs)
        bond_orders = np.zeros((len(self.atoms), len(self.atoms)))
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                bo = self.compute_bond_order(density_matrix, i, j)
                bond_orders[i, j] = bo
                bond_orders[j, i] = bo

        analysis['bond_orders'] = bond_orders

        return analysis
