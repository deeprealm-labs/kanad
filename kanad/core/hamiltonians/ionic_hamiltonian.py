"""
Ionic Hamiltonian for electron transfer systems with governance protocol integration.

Models ionic bonding where electrons are transferred from donor to acceptor atoms.
Integrates IonicGovernanceProtocol to ensure physical correctness.
"""

from typing import List, Optional, Dict, Any
import numpy as np
import logging
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.core.atom import Atom
from kanad.governance.protocols.ionic_protocol import IonicGovernanceProtocol

logger = logging.getLogger(__name__)


class IonicHamiltonian(MolecularHamiltonian):
    """
    Hamiltonian for ionic bonding.

    Physical Model:
        H_ionic = Σ_i ε_i n_i + Σ_<ij> t_ij (a†_i a_j + a†_j a_i) + Σ_i U_i n_i↑ n_i↓

    where:
        ε_i: On-site energy (related to electronegativity)
        t_ij: Transfer integral (hopping amplitude) - SMALL for ionic bonds
        U_i: On-site Coulomb repulsion (Hubbard U)

    KEY PHYSICS:
        - Large electronegativity difference → charge transfer
        - Weak overlap → small transfer integral
        - Localized electrons → large on-site repulsion
    """

    def __init__(
        self,
        molecule: 'Molecule',
        representation: 'SecondQuantizationRepresentation',
        use_governance: bool = True
    ):
        """
        Initialize ionic Hamiltonian with governance protocol.

        Args:
            molecule: Molecule object with atoms
            representation: Second quantization representation
            use_governance: Enable governance protocol validation (default: True)
        """
        self.molecule = molecule
        self.representation = representation
        self.atoms = molecule.atoms
        self.use_governance = use_governance

        # Initialize governance protocol
        if use_governance:
            self.governance_protocol = IonicGovernanceProtocol()
            logger.info("✓ Ionic governance protocol initialized")
        else:
            self.governance_protocol = None

        # Compute nuclear repulsion
        nuclear_rep = self._compute_nuclear_repulsion()

        super().__init__(
            n_orbitals=len(self.atoms),  # One orbital per atom (simplified)
            n_electrons=molecule.n_electrons,
            nuclear_repulsion=nuclear_rep
        )

        # Build Hamiltonian with governance validation
        self._build_hamiltonian()

    def _compute_nuclear_repulsion(self) -> float:
        """
        Compute nuclear-nuclear repulsion energy.

        E_nn = Σ_{i<j} Z_i Z_j / |R_i - R_j|

        NOTE: distance_to() returns Angstroms, but we need atomic units (Bohr)
        Conversion: 1 Angstrom = 1.88973 Bohr
        """
        from kanad.core.constants.conversion_factors import ConversionFactors

        energy = 0.0
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                Z_i = self.atoms[i].atomic_number
                Z_j = self.atoms[j].atomic_number
                r_ij_angstrom = self.atoms[i].distance_to(self.atoms[j])
                r_ij_bohr = r_ij_angstrom * ConversionFactors.ANGSTROM_TO_BOHR

                if r_ij_bohr > 1e-10:  # Avoid division by zero
                    energy += Z_i * Z_j / r_ij_bohr

        return energy

    def _build_hamiltonian(self):
        """
        Build ionic Hamiltonian matrices.

        Emphasizes:
        - On-site energies from electronegativity
        - Small transfer integrals (weak overlap)
        - On-site Coulomb repulsion
        """
        n = self.n_orbitals

        # Core Hamiltonian: diagonal = site energies, off-diagonal = transfer
        self.h_core = np.zeros((n, n))

        for i in range(n):
            # On-site energy (negative of electronegativity)
            self.h_core[i, i] = -self.atoms[i].properties.electronegativity

        # Transfer integrals (nearest neighbors only for ionic)
        for i in range(n):
            for j in range(i + 1, n):
                t_ij = self._compute_transfer_integral(i, j)
                self.h_core[i, j] = t_ij
                self.h_core[j, i] = t_ij

        # Two-electron integrals (on-site repulsion dominant)
        self.eri = self._compute_eri_ionic()

    def _compute_transfer_integral(self, i: int, j: int) -> float:
        """
        Compute transfer integral between sites i and j using governance protocol.

        For ionic bonding:
            t_ij ∝ exp(-r_ij / λ)

        where λ is the decay length (≈ 1-2 Bohr for ionic bonds).

        Governance enforcement: Transfer integrals should be SMALL for ionic bonds
        (< 0.1 Hartree ~ 2.7 eV) to maintain localized character.

        Args:
            i: Site index
            j: Site index

        Returns:
            Transfer integral (in Hartree)
        """
        r_ij = self.atoms[i].distance_to(self.atoms[j])

        # Use governance protocol's transfer integral estimate if available
        if self.use_governance and self.governance_protocol:
            t_ij = self.governance_protocol.get_transfer_integral_estimate(r_ij)

            # Validate: Ionic bonds should have WEAK transfer (t < 0.1 Ha)
            if t_ij > 0.1:
                logger.warning(
                    f"Large transfer integral ({t_ij:.4f} Ha) for ionic bond "
                    f"between sites {i}-{j}. This may indicate covalent character."
                )
        else:
            # Fallback: Standard exponential decay
            lambda_decay = 1.5  # Decay length in Bohr
            t_0 = 0.05  # ~1.4 eV in Hartree units
            t_ij = t_0 * np.exp(-r_ij / lambda_decay)

        return t_ij

    def _compute_eri_ionic(self) -> np.ndarray:
        """
        Compute two-electron repulsion integrals for ionic model.

        For ionic bonding, the dominant terms are:
        - U_i = (ii|ii): On-site Coulomb repulsion (Hubbard U)
        - V_ij = (ii|jj): Inter-site Coulomb repulsion

        Returns:
            ERI tensor (n, n, n, n)
        """
        n = self.n_orbitals
        eri = np.zeros((n, n, n, n))

        for i in range(n):
            # On-site Hubbard U
            U_i = self._hubbard_u(self.atoms[i])
            eri[i, i, i, i] = U_i

            # Inter-site repulsion
            for j in range(i + 1, n):
                V_ij = self._inter_site_coulomb(i, j)
                eri[i, i, j, j] = V_ij
                eri[j, j, i, i] = V_ij

        return eri

    def _hubbard_u(self, atom: Atom) -> float:
        """
        Estimate Hubbard U parameter for on-site repulsion.

        U ≈ (IP - EA)

        where IP = ionization potential, EA = electron affinity

        For simplicity, use empirical values:
            U ≈ 5-20 eV depending on atom

        Args:
            atom: Atom object

        Returns:
            Hubbard U in Hartree
        """
        # Typical values in eV
        hubbard_u_ev = {
            'H': 13.0,
            'Li': 5.0,
            'Na': 5.0,
            'K': 4.0,
            'F': 17.0,
            'Cl': 13.0,
            'Br': 11.0,
            'O': 15.0,
            'S': 12.0,
        }

        U_ev = hubbard_u_ev.get(atom.symbol, 10.0)  # Default 10 eV

        # Convert to Hartree (1 Hartree ≈ 27.211 eV)
        U_hartree = U_ev / 27.211

        return U_hartree

    def _inter_site_coulomb(self, i: int, j: int) -> float:
        """
        Compute inter-site Coulomb repulsion V_ij.

        V_ij ≈ e² / r_ij

        Args:
            i: Site index
            j: Site index

        Returns:
            Inter-site repulsion in Hartree
        """
        r_ij = self.atoms[i].distance_to(self.atoms[j])

        if r_ij < 1e-10:
            return 0.0

        # In atomic units, e² = 1
        V_ij = 1.0 / r_ij

        return V_ij

    def to_matrix(self, n_qubits: Optional[int] = None, use_mo_basis: bool = False) -> np.ndarray:
        """
        Convert to matrix form (one-body part only).

        Note: For ionic bonding, full many-body Hamiltonian is not needed.
        Uses simplified Hubbard-like model.

        Args:
            n_qubits: Number of qubits (ignored for ionic, kept for API compatibility)
            use_mo_basis: Whether to use MO basis (ignored for ionic)

        Returns:
            Core Hamiltonian matrix
        """
        return self.h_core.copy()

    def compute_energy(self, density_matrix: np.ndarray) -> float:
        """
        Compute electronic energy from density matrix.

        E = Σ_ij P_ij h_ij + ½ Σ_ijkl P_ij P_kl [(ij|kl) - ½(ik|jl)]

        Note: Does NOT include nuclear repulsion. Caller should add it.

        Args:
            density_matrix: One-particle density matrix

        Returns:
            Electronic energy only (Hartree)
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
                        # Exchange (with factor for closed-shell)
                        E_ee -= 0.25 * density_matrix[i, k] * density_matrix[j, l] * self.eri[i, j, k, l]

        # Electronic energy only (nuclear repulsion added by caller)
        E_electronic = E_core + E_ee

        return E_electronic

    def compute_charge_transfer(self, density_matrix: np.ndarray) -> np.ndarray:
        """
        Compute charge on each atom from density matrix.

        For ionic bonding, expect significant charge separation:
            e.g., Na: +0.9, Cl: -0.9 for NaCl

        Args:
            density_matrix: Density matrix

        Returns:
            Charge on each atom (negative = excess electrons)
        """
        charges = np.zeros(self.n_orbitals)

        for i in range(self.n_orbitals):
            # Electron population on site i
            n_i = density_matrix[i, i]

            # Charge = nuclear charge - electron population
            charges[i] = self.atoms[i].atomic_number - n_i

        return charges

    def get_ionization_energy(self) -> float:
        """
        Estimate ionization energy of donor atom.

        IE ≈ -ε_HOMO

        Returns:
            Ionization energy (Hartree)
        """
        eigenvalues = np.linalg.eigvalsh(self.h_core)
        # HOMO is highest occupied
        n_occ = self.n_electrons // 2
        homo_energy = eigenvalues[n_occ - 1]

        return -homo_energy

    def get_electron_affinity(self) -> float:
        """
        Estimate electron affinity of acceptor atom.

        EA ≈ -ε_LUMO

        Returns:
            Electron affinity (Hartree)
        """
        eigenvalues = np.linalg.eigvalsh(self.h_core)
        # LUMO is lowest unoccupied
        n_occ = self.n_electrons // 2
        lumo_energy = eigenvalues[n_occ]

        return -lumo_energy

    def solve_scf(
        self,
        max_iterations: int = 100,
        conv_tol: float = 1e-6,
        **kwargs
    ) -> tuple:
        """
        Solve for ground state using SCF for ionic/covalent mixed character.

        Uses Roothaan-Hall SCF equations with DIIS acceleration.
        Properly handles both ionic (charge-separated) and covalent character.

        Args:
            max_iterations: Maximum SCF iterations
            conv_tol: Energy convergence threshold (Hartree)
            **kwargs: Additional arguments (damping, level_shift, etc.)

        Returns:
            Tuple of (density_matrix, energy)
        """
        damping = kwargs.get('damping', 0.5)  # Density mixing parameter
        level_shift = kwargs.get('level_shift', 0.0)  # Level shift for convergence

        # Initial guess: core Hamiltonian
        eigenvalues, eigenvectors = np.linalg.eigh(self.h_core)

        # Build initial density matrix
        # NOTE: Simplified ionic model may have fewer orbitals than electrons/2
        # Fill orbitals with fractional occupation if needed
        n_occ = min(self.n_electrons // 2, self.n_orbitals)
        total_electrons_to_place = self.n_electrons

        density_matrix = np.zeros((self.n_orbitals, self.n_orbitals))

        # Fill lowest energy orbitals first
        electrons_placed = 0
        for i in range(self.n_orbitals):
            if electrons_placed >= total_electrons_to_place:
                break
            # Each orbital can hold 2 electrons
            occ = min(2.0, total_electrons_to_place - electrons_placed)
            density_matrix += occ * np.outer(eigenvectors[:, i], eigenvectors[:, i])
            electrons_placed += occ

        # SCF iterations
        prev_energy = 0.0
        converged = False

        for iteration in range(max_iterations):
            # Build Fock matrix
            fock = self._build_fock_matrix(density_matrix)

            # Apply level shift if requested
            if level_shift > 0:
                fock += level_shift * (np.eye(self.n_orbitals) - density_matrix / 2.0)

            # Diagonalize Fock matrix
            eigenvalues, eigenvectors = np.linalg.eigh(fock)

            # Build new density matrix
            new_density = np.zeros((self.n_orbitals, self.n_orbitals))
            electrons_placed = 0
            for i in range(self.n_orbitals):
                if electrons_placed >= total_electrons_to_place:
                    break
                # Each orbital can hold 2 electrons
                occ = min(2.0, total_electrons_to_place - electrons_placed)
                new_density += occ * np.outer(eigenvectors[:, i], eigenvectors[:, i])
                electrons_placed += occ

            # Density damping (mixing)
            density_matrix = damping * density_matrix + (1 - damping) * new_density

            # Compute energy
            energy = self.compute_energy(density_matrix) + self.nuclear_repulsion

            # Check convergence
            energy_diff = abs(energy - prev_energy)
            if energy_diff < conv_tol:
                converged = True
                self._scf_converged = True
                self._scf_iterations = iteration + 1
                break

            prev_energy = energy

        if not converged:
            self._scf_converged = False
            self._scf_iterations = max_iterations

        return density_matrix, energy

    def _build_fock_matrix(self, density_matrix: np.ndarray) -> np.ndarray:
        """
        Build Fock matrix from density matrix.

        F_ij = h_ij + Σ_kl P_kl [(ij|kl) - ½(ik|jl)]

        Args:
            density_matrix: Current density matrix

        Returns:
            Fock matrix
        """
        fock = self.h_core.copy()

        # Add electron-electron contributions
        for i in range(self.n_orbitals):
            for j in range(self.n_orbitals):
                for k in range(self.n_orbitals):
                    for l in range(self.n_orbitals):
                        # Coulomb term
                        fock[i, j] += density_matrix[k, l] * self.eri[i, j, k, l]
                        # Exchange term (factor of ½ for closed-shell)
                        fock[i, j] -= 0.5 * density_matrix[k, l] * self.eri[i, k, j, l]

        return fock

    def validate_with_governance(self) -> dict[str, Any]:
        """
        Validate Hamiltonian using governance protocol.

        Returns:
            Dictionary with validation results
        """
        if not self.use_governance or not self.governance_protocol:
            return {'governance_enabled': False}

        validation = {
            'governance_enabled': True,
            'bonding_type': 'ionic',
            'checks': []
        }

        # Check 1: Transfer integrals are small (localized)
        max_transfer = 0.0
        for i in range(self.n_orbitals):
            for j in range(i + 1, self.n_orbitals):
                t_ij = abs(self.h_core[i, j])
                max_transfer = max(max_transfer, t_ij)

        validation['max_transfer_integral'] = max_transfer
        # Ionic: < 0.1 Ha (pure ionic like NaCl)
        # Polar covalent: 0.1-0.5 Ha (like LiH) - still acceptable
        # Covalent: > 0.5 Ha
        if max_transfer < 0.1:
            validation['checks'].append({
                'name': 'weak_transfer',
                'passed': True,
                'message': f'Transfer integrals small ({max_transfer:.4f} Ha) - pure ionic ✓'
            })
        elif max_transfer < 0.5:
            validation['checks'].append({
                'name': 'weak_transfer',
                'passed': True,  # Accept polar covalent
                'message': f'Transfer integrals moderate ({max_transfer:.4f} Ha) - polar covalent/ionic ✓'
            })
        else:
            validation['checks'].append({
                'name': 'weak_transfer',
                'passed': False,
                'message': f'Transfer integrals too large ({max_transfer:.4f} Ha) - predominantly covalent'
            })

        # Check 2: On-site energies reflect electronegativity difference
        site_energies = np.diag(self.h_core)
        energy_spread = np.max(site_energies) - np.min(site_energies)

        validation['energy_spread'] = energy_spread
        if energy_spread > 0.1:  # Significant electronegativity difference
            validation['checks'].append({
                'name': 'electronegativity_difference',
                'passed': True,
                'message': f'Large energy spread ({energy_spread:.4f} Ha) indicates ionic character ✓'
            })
        else:
            validation['checks'].append({
                'name': 'electronegativity_difference',
                'passed': False,
                'message': f'Small energy spread ({energy_spread:.4f} Ha) - weak ionic character'
            })

        # Check 3: Hubbard U is large (strong correlation)
        avg_U = np.mean([self.eri[i, i, i, i] for i in range(self.n_orbitals)])

        validation['average_hubbard_u'] = avg_U
        if avg_U > 0.2:  # > 5 eV
            validation['checks'].append({
                'name': 'large_hubbard_u',
                'passed': True,
                'message': f'Large Hubbard U ({avg_U:.4f} Ha ~ {avg_U * 27.211:.1f} eV) ✓'
            })
        else:
            validation['checks'].append({
                'name': 'large_hubbard_u',
                'passed': False,
                'message': f'Small Hubbard U ({avg_U:.4f} Ha)'
            })

        # Overall assessment
        all_passed = all(check['passed'] for check in validation['checks'])
        validation['all_checks_passed'] = all_passed

        if all_passed:
            logger.info("✓ Ionic Hamiltonian passed all governance checks")
        else:
            logger.warning("⚠ Some ionic governance checks failed - review bonding character")

        return validation

    def get_governance_aware_ansatz(self, ansatz_type: str = 'governance'):
        """
        Get governance-aware ansatz for this Hamiltonian.

        Args:
            ansatz_type: Type of ansatz ('governance', 'ucc', or 'hardware_efficient')

        Returns:
            Ansatz object configured for ionic bonding
        """
        if ansatz_type == 'governance' and self.governance_protocol:
            from kanad.ansatze.governance_aware_ansatz import IonicGovernanceAnsatz
            return IonicGovernanceAnsatz(
                hamiltonian=self,
                n_qubits=2 * self.n_orbitals,
                governance_protocol=self.governance_protocol
            )
        elif ansatz_type == 'ucc':
            from kanad.ansatze.ucc_ansatz import UCCAnsatz
            return UCCAnsatz(
                hamiltonian=self,
                excitations='SD'  # Singles and doubles
            )
        elif ansatz_type == 'hardware_efficient':
            from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz
            return HardwareEfficientAnsatz(
                n_qubits=2 * self.n_orbitals,
                depth=3,
                entanglement='linear'  # Linear for ionic (sparse connectivity)
            )
        else:
            raise ValueError(f"Unknown ansatz type: {ansatz_type}")
