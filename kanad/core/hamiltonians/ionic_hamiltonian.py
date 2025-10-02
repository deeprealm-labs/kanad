"""
Ionic Hamiltonian for electron transfer systems.

Models ionic bonding where electrons are transferred from donor to acceptor atoms.
"""

from typing import List, Optional
import numpy as np
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.core.atom import Atom


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
        representation: 'SecondQuantizationRepresentation'
    ):
        """
        Initialize ionic Hamiltonian.

        Args:
            molecule: Molecule object with atoms
            representation: Second quantization representation
        """
        self.molecule = molecule
        self.representation = representation
        self.atoms = molecule.atoms

        # Compute nuclear repulsion
        nuclear_rep = self._compute_nuclear_repulsion()

        super().__init__(
            n_orbitals=len(self.atoms),  # One orbital per atom (simplified)
            n_electrons=molecule.n_electrons,
            nuclear_repulsion=nuclear_rep
        )

        # Build Hamiltonian
        self._build_hamiltonian()

    def _compute_nuclear_repulsion(self) -> float:
        """
        Compute nuclear-nuclear repulsion energy.

        E_nn = Σ_{i<j} Z_i Z_j / |R_i - R_j|
        """
        energy = 0.0
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                Z_i = self.atoms[i].atomic_number
                Z_j = self.atoms[j].atomic_number
                r_ij = self.atoms[i].distance_to(self.atoms[j])
                if r_ij > 1e-10:  # Avoid division by zero
                    energy += Z_i * Z_j / r_ij

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
        Compute transfer integral between sites i and j.

        For ionic bonding:
            t_ij ∝ exp(-r_ij / λ)

        where λ is the decay length (≈ 1-2 Bohr for ionic bonds).

        Args:
            i: Site index
            j: Site index

        Returns:
            Transfer integral (in Hartree)
        """
        r_ij = self.atoms[i].distance_to(self.atoms[j])

        # Decay length for ionic bonding (in Bohr)
        lambda_decay = 1.5

        # Prefactor (typical ionic transfer integral scale)
        t_0 = 0.05  # ~1.4 eV in Hartree units

        # Exponential decay
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

        E = Σ_ij P_ij h_ij + ½ Σ_ijkl P_ij P_kl [(ij|kl) - ½(ik|jl)] + E_nn

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
                        # Exchange (with factor for closed-shell)
                        E_ee -= 0.25 * density_matrix[i, k] * density_matrix[j, l] * self.eri[i, j, k, l]

        # Total energy
        E_total = E_core + E_ee + self.nuclear_repulsion

        return E_total

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
