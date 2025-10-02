"""
Second quantization representation for ionic bonding.

In ionic bonding, electrons are localized on specific atoms.
This representation uses localized atomic orbitals with minimal
entanglement between sites.

Physical picture:
- Electron transfer: a†_i a_j (electron hops from j to i)
- Localized states: minimal entanglement
- Hubbard-like model with on-site repulsion
"""

import numpy as np
from typing import Dict, List, Optional
from kanad.core.representations.base_representation import BaseRepresentation, Molecule


class SecondQuantizationRepresentation(BaseRepresentation):
    """
    Second quantization (Fock space) representation for ionic systems.

    Uses atomic orbital basis with creation/annihilation operators.
    Optimized for systems with electron transfer between localized sites.

    Example: Na+ Cl-
    - Na: loses electron (donor site)
    - Cl: gains electron (acceptor site)
    - Hamiltonian: H = ε_Na n_Na + ε_Cl n_Cl + t(a†_Na a_Cl + h.c.)
    """

    def __init__(self, molecule: Molecule, include_spin: bool = True):
        """
        Initialize second quantization representation.

        Args:
            molecule: Molecule object
            include_spin: Whether to include spin degrees of freedom
        """
        super().__init__(molecule)
        self.include_spin = include_spin

        # Set number of orbitals (one per atom for minimal basis)
        self.n_orbitals = molecule.n_atoms

        # Number of spin orbitals
        self.n_spin_orbitals = 2 * self.n_orbitals if include_spin else self.n_orbitals

        # Number of qubits (Jordan-Wigner mapping: 1 qubit per spin-orbital)
        self.n_qubits = self.n_spin_orbitals

        # Calculate total number of electrons
        self.n_electrons = sum(atom.atomic_number for atom in molecule.atoms)

    def build_hamiltonian(self) -> 'IonicHamiltonian':
        """
        Build ionic Hamiltonian in second quantized form.

        H = Σ_i ε_i n_i + Σ_ij t_ij a†_i a_j + Σ_i U_i n_i↑ n_i↓

        where:
            ε_i: on-site energy (electronegativity)
            t_ij: hopping/transfer integral
            U_i: on-site Coulomb repulsion (Hubbard U)

        Returns:
            IonicHamiltonian object
        """
        from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian

        # Build ionic Hamiltonian with charge transfer and site energies
        return IonicHamiltonian(
            molecule=self.molecule,
            representation=self
        )

    def get_reference_state(self) -> np.ndarray:
        """
        Get reference state (Hartree-Fock or simple product state).

        For ionic systems, this is typically the charge-separated state.
        E.g., for NaCl: |Na+ Cl-⟩ = |0⟩_Na ⊗ |↑↓⟩_Cl

        Returns:
            Reference state vector
        """
        # For simplified one-orbital-per-atom model, we fill spin orbitals
        # based on the minimal representation (not full electron count)
        state_dim = 2 ** self.n_qubits
        ref_state = np.zeros(state_dim)

        # Construct reference determinant
        # For ionic bonds, typically occupy the more electronegative atom's orbitals
        # For now, fill lowest energy spin orbitals up to n_qubits
        # This gives the proper charge-separated state
        hf_occupation = 0
        n_occ = min(self.n_electrons, self.n_qubits)  # Can't exceed number of qubits
        for i in range(n_occ):
            hf_occupation |= (1 << i)

        ref_state[hf_occupation] = 1.0

        return ref_state

    def compute_observables(self, state: np.ndarray) -> Dict[str, float]:
        """
        Compute observables for ionic system.

        Args:
            state: Quantum state vector

        Returns:
            Dictionary of observables:
                - 'charge_transfer': Amount of charge transferred
                - 'site_occupations': Occupation numbers per site
                - 'energy': Expectation value of energy
        """
        observables = {}

        # Compute site occupations from quantum state
        site_occupations = np.zeros(self.n_orbitals)

        # Parse state vector to get occupations for each site
        for i in range(self.n_orbitals):
            # Count occupation of spin-up and spin-down on site i
            site_occupations[i] = self._compute_site_occupation(state, i)

        observables['site_occupations'] = site_occupations

        # Charge transfer (difference from neutral atoms)
        neutral_occupations = np.array([atom.n_valence for atom in self.molecule.atoms])
        charge_transfer = site_occupations - neutral_occupations
        observables['charge_transfer'] = charge_transfer

        # Total charge transferred (magnitude)
        observables['total_charge_transfer'] = np.sum(np.abs(charge_transfer)) / 2

        return observables

    def _compute_site_occupation(self, state: np.ndarray, site: int) -> float:
        """
        Compute occupation number for a site.

        n_i = ⟨a†_i↑ a_i↑ + a†_i↓ a_i↓⟩

        Args:
            state: Quantum state
            site: Site index

        Returns:
            Occupation number (0 to 2)
        """
        # Compute expectation value ⟨ψ|n_i|ψ⟩ where n_i = n_{i↑} + n_{i↓}
        # by summing |⟨basis|ψ⟩|² weighted by occupation in each basis state
        occupation = 0.0
        state_dim = len(state)

        for basis_state in range(state_dim):
            amplitude = state[basis_state]
            if abs(amplitude) > 1e-10:
                # Count occupation in this basis state using bit operations
                # Spin-up orbital for site i
                if basis_state & (1 << (2 * site)):
                    occupation += abs(amplitude) ** 2

                # Spin-down orbital for site i
                if basis_state & (1 << (2 * site + 1)):
                    occupation += abs(amplitude) ** 2

        return occupation

    def to_qubit_operator(self) -> Dict[str, complex]:
        """
        Map Hamiltonian to qubit operators using Jordan-Wigner transformation.

        a†_j → (∏_{k<j} Z_k) (X_j - iY_j) / 2
        a_j → (∏_{k<j} Z_k) (X_j + iY_j) / 2

        Returns:
            Dictionary mapping Pauli strings to complex coefficients
        """
        from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper

        # Build Hamiltonian if needed to get nuclear repulsion
        if not hasattr(self, 'hamiltonian') or self.hamiltonian is None:
            self.hamiltonian = self.build_hamiltonian()

        mapper = JordanWignerMapper()
        pauli_hamiltonian = {}

        # Get Hamiltonian parameters
        on_site = self.get_on_site_energies()
        hopping = self.get_hopping_matrix()
        n_sites = len(on_site)

        # Map on-site terms: ε_i n_i
        for i in range(n_sites):
            if abs(on_site[i]) > 1e-10:
                pauli_term = mapper.map_hamiltonian_term((i, i), on_site[i], n_sites)
                for pauli_string, coeff in pauli_term.items():
                    if pauli_string in pauli_hamiltonian:
                        pauli_hamiltonian[pauli_string] += coeff
                    else:
                        pauli_hamiltonian[pauli_string] = coeff

        # Map hopping terms: t_ij a†_i a_j
        for i in range(n_sites):
            for j in range(n_sites):
                if i != j and abs(hopping[i, j]) > 1e-10:
                    pauli_term = mapper.map_hamiltonian_term((i, j), hopping[i, j], n_sites)
                    for pauli_string, coeff in pauli_term.items():
                        if pauli_string in pauli_hamiltonian:
                            pauli_hamiltonian[pauli_string] += coeff
                        else:
                            pauli_hamiltonian[pauli_string] = coeff

        # Add nuclear repulsion as constant (identity)
        identity = 'I' * self.n_qubits
        nuclear_rep = self.hamiltonian.nuclear_repulsion if hasattr(self.hamiltonian, 'nuclear_repulsion') else 0.0
        if identity in pauli_hamiltonian:
            pauli_hamiltonian[identity] += nuclear_rep
        else:
            pauli_hamiltonian[identity] = nuclear_rep

        # Clean up near-zero terms
        pauli_hamiltonian = {k: v for k, v in pauli_hamiltonian.items() if abs(v) > 1e-12}

        return pauli_hamiltonian

    def get_num_qubits(self) -> int:
        """Get number of qubits (one per spin-orbital)."""
        return self.n_qubits

    def get_hopping_matrix(self) -> np.ndarray:
        """
        Compute hopping matrix t_ij.

        For ionic systems, hopping is typically nearest-neighbor only
        and exponentially decreases with distance.

        Returns:
            Hopping matrix (n_orbitals, n_orbitals)
        """
        n = self.n_orbitals
        t = np.zeros((n, n))

        # Compute distances
        D = self.molecule.distance_matrix()

        # Hopping decreases exponentially with distance
        # t_ij = t_0 exp(-d_ij / d_0)
        t_0 = 1.0  # Base hopping strength (eV or Hartree)
        d_0 = 2.0  # Decay length (Angstrom)

        for i in range(n):
            for j in range(i + 1, n):
                t[i, j] = t_0 * np.exp(-D[i, j] / d_0)
                t[j, i] = t[i, j]  # Symmetric

        return t

    def get_on_site_energies(self) -> np.ndarray:
        """
        Get on-site energies from electronegativities.

        ε_i = -χ_i (electronegativity)

        Returns:
            On-site energies (n_orbitals,)
        """
        energies = np.zeros(self.n_orbitals)

        for i, atom in enumerate(self.molecule.atoms):
            # On-site energy related to electronegativity
            # More electronegative → lower energy
            energies[i] = -atom.properties.electronegativity

        return energies

    def __repr__(self) -> str:
        """String representation."""
        return f"SecondQuantization(n_qubits={self.n_qubits}, n_electrons={self.n_electrons})"
