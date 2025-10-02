"""
LCAO (Linear Combination of Atomic Orbitals) representation for covalent bonding.

In covalent bonding, atomic orbitals hybridize to form molecular orbitals.
This representation uses hybrid orbitals (sp, sp², sp³) that form
bonding/antibonding pairs.

Physical picture:
- Orbital hybridization: sp³ → 4 tetrahedral orbitals
- Bonding/antibonding MOs: |ψ_±⟩ = (|φ_A⟩ ± |φ_B⟩)/√2
- Paired entanglement: Bell-like states for electron pairs
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from enum import Enum
from kanad.core.representations.base_representation import BaseRepresentation, Molecule


class HybridizationType(Enum):
    """Types of orbital hybridization."""
    SP = "sp"  # Linear (180°)
    SP2 = "sp2"  # Trigonal planar (120°)
    SP3 = "sp3"  # Tetrahedral (109.5°)
    NONE = "none"  # No hybridization (e.g., H)


class LCAORepresentation(BaseRepresentation):
    """
    LCAO representation for covalent bonding systems.

    Uses hybrid atomic orbitals that combine to form molecular orbitals.
    Optimized for systems with shared electron pairs and orbital overlap.

    Example: H₂O
    - O: sp³ hybridization (2 bonding + 2 lone pairs)
    - H: 1s orbitals
    - Bonding MOs: |ψ_OH⟩ = (|sp³_O⟩ + |1s_H⟩)/√2
    """

    def __init__(
        self,
        molecule: Molecule,
        hybridization: Optional[Dict[int, HybridizationType]] = None
    ):
        """
        Initialize LCAO representation.

        Args:
            molecule: Molecule object
            hybridization: Dict mapping atom index to hybridization type
                          If None, automatically determine
        """
        super().__init__(molecule)

        # Determine hybridization for each atom
        if hybridization is None:
            self.hybridization = self._auto_detect_hybridization()
        else:
            self.hybridization = hybridization

        # Build hybrid orbital basis
        self.hybrid_orbitals = self._construct_hybrid_orbitals()

        # Number of spatial orbitals = number of hybrid orbitals
        self.n_orbitals = len(self.hybrid_orbitals)

        # Number of spin orbitals
        self.n_spin_orbitals = 2 * self.n_orbitals

        # For covalent systems, use paired mapping
        # Each bonding/antibonding pair → 2 qubits
        self.n_qubits = self.n_spin_orbitals

        # Build molecular orbital pairs
        self.mo_pairs = self._construct_mo_pairs()

    def _auto_detect_hybridization(self) -> Dict[int, HybridizationType]:
        """
        Automatically detect hybridization based on atom type and bonding.

        Returns:
            Dictionary mapping atom index to hybridization type
        """
        hybridization = {}

        for i, atom in enumerate(self.molecule.atoms):
            if atom.symbol == 'H':
                # Hydrogen: no hybridization
                hybridization[i] = HybridizationType.NONE

            elif atom.symbol == 'C':
                # Carbon: default to sp³ (can be refined based on bonding)
                hybridization[i] = HybridizationType.SP3

            elif atom.symbol in ['N', 'O']:
                # Nitrogen, Oxygen: sp³
                hybridization[i] = HybridizationType.SP3

            else:
                # Default: no hybridization
                hybridization[i] = HybridizationType.NONE

        return hybridization

    def _construct_hybrid_orbitals(self) -> List[Dict]:
        """
        Construct hybrid orbitals for each atom.

        Returns:
            List of hybrid orbital dictionaries
        """
        hybrids = []

        for i, atom in enumerate(self.molecule.atoms):
            hybrid_type = self.hybridization[i]

            if hybrid_type == HybridizationType.SP3:
                # sp³: 4 tetrahedral orbitals
                for j in range(4):
                    hybrids.append({
                        'atom_index': i,
                        'type': 'sp3',
                        'index': j,
                        'direction': self._sp3_direction(j)
                    })

            elif hybrid_type == HybridizationType.SP2:
                # sp²: 3 trigonal planar + 1 p orbital
                for j in range(3):
                    hybrids.append({
                        'atom_index': i,
                        'type': 'sp2',
                        'index': j,
                        'direction': self._sp2_direction(j)
                    })
                # Pure p orbital
                hybrids.append({
                    'atom_index': i,
                    'type': 'p',
                    'index': 3,
                    'direction': np.array([0, 0, 1])
                })

            elif hybrid_type == HybridizationType.SP:
                # sp: 2 linear + 2 p orbitals
                hybrids.append({
                    'atom_index': i,
                    'type': 'sp',
                    'index': 0,
                    'direction': np.array([1, 0, 0])
                })
                hybrids.append({
                    'atom_index': i,
                    'type': 'sp',
                    'index': 1,
                    'direction': np.array([-1, 0, 0])
                })
                # Two p orbitals
                hybrids.append({
                    'atom_index': i,
                    'type': 'p',
                    'index': 2,
                    'direction': np.array([0, 1, 0])
                })
                hybrids.append({
                    'atom_index': i,
                    'type': 'p',
                    'index': 3,
                    'direction': np.array([0, 0, 1])
                })

            else:  # NONE (e.g., Hydrogen)
                # Single s orbital
                hybrids.append({
                    'atom_index': i,
                    'type': 's',
                    'index': 0,
                    'direction': np.array([0, 0, 0])  # Spherical
                })

        return hybrids

    @staticmethod
    def _sp3_direction(index: int) -> np.ndarray:
        """Get direction vector for sp³ hybrid orbital."""
        # Tetrahedral geometry
        directions = [
            np.array([1, 1, 1]),    # sp³_1
            np.array([1, -1, -1]),  # sp³_2
            np.array([-1, 1, -1]),  # sp³_3
            np.array([-1, -1, 1])   # sp³_4
        ]
        return directions[index] / np.linalg.norm(directions[index])

    @staticmethod
    def _sp2_direction(index: int) -> np.ndarray:
        """Get direction vector for sp² hybrid orbital."""
        # Trigonal planar geometry (120° apart)
        angle = index * 2 * np.pi / 3
        return np.array([np.cos(angle), np.sin(angle), 0])

    def _construct_mo_pairs(self) -> List[Tuple[int, int]]:
        """
        Construct bonding/antibonding molecular orbital pairs.

        For each bond, create a pair of indices (bonding, antibonding).

        Returns:
            List of (bonding_idx, antibonding_idx) tuples
        """
        mo_pairs = []

        # Simplified: Assume each pair of atoms forms one bond
        # Full implementation would use bond connectivity
        n_atoms = self.molecule.n_atoms

        if n_atoms == 2:
            # Diatomic: single bond
            mo_pairs.append((0, 1))  # Bonding and antibonding

        return mo_pairs

    def build_hamiltonian(self) -> 'CovalentHamiltonian':
        """
        Build covalent Hamiltonian in hybrid orbital basis.

        H = Σ_μν h_μν c†_μ c_ν + ½ Σ_μνλσ (μν|λσ) c†_μ c†_ν c_σ c_λ

        where μ,ν run over hybrid orbitals.

        Returns:
            CovalentHamiltonian object
        """
        from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

        return CovalentHamiltonian(
            molecule=self.molecule,
            representation=self
        )

    def get_reference_state(self) -> np.ndarray:
        """
        Get reference state (Hartree-Fock in MO basis).

        For covalent systems, this is the closed-shell state with
        bonding orbitals occupied.

        Returns:
            Reference state vector
        """
        state_dim = 2 ** self.n_qubits
        ref_state = np.zeros(state_dim)

        # Occupy bonding orbitals
        hf_occupation = 0
        n_occ = self.n_electrons // 2  # Pairs

        for i in range(n_occ):
            # Occupy bonding orbital (spin up and down)
            hf_occupation |= (1 << (2 * i))      # Spin up
            hf_occupation |= (1 << (2 * i + 1))  # Spin down

        ref_state[hf_occupation] = 1.0

        return ref_state

    def compute_observables(self, state: np.ndarray) -> Dict[str, float]:
        """
        Compute observables for covalent system.

        Args:
            state: Quantum state vector

        Returns:
            Dictionary of observables:
                - 'bond_orders': Bond order for each bond
                - 'overlap_populations': Mulliken overlap populations
                - 'energy': Expectation value of energy
        """
        observables = {}

        # Compute bond orders (simplified)
        bond_orders = []
        for bonding_idx, antibonding_idx in self.mo_pairs:
            # Bond order = (n_bonding - n_antibonding) / 2
            n_bonding = self._compute_orbital_occupation(state, bonding_idx)
            n_antibonding = self._compute_orbital_occupation(state, antibonding_idx)
            bond_order = (n_bonding - n_antibonding) / 2
            bond_orders.append(bond_order)

        observables['bond_orders'] = np.array(bond_orders)

        return observables

    def _compute_orbital_occupation(self, state: np.ndarray, orbital: int) -> float:
        """Compute occupation of a molecular orbital."""
        # Simplified calculation
        occupation = 0.0

        for basis_state in range(len(state)):
            amplitude = state[basis_state]
            if abs(amplitude) > 1e-10:
                # Check if orbital is occupied
                if basis_state & (1 << (2 * orbital)):  # Spin up
                    occupation += abs(amplitude) ** 2
                if basis_state & (1 << (2 * orbital + 1)):  # Spin down
                    occupation += abs(amplitude) ** 2

        return occupation

    def to_qubit_operator(self) -> 'QubitOperator':
        """
        Map to qubit operator using hybrid orbital mapper.

        Uses paired qubit encoding for bonding/antibonding orbitals.

        Returns:
            QubitOperator (placeholder)
        """
        return None

    def get_num_qubits(self) -> int:
        """Get number of qubits."""
        return self.n_qubits

    def get_bonding_antibonding_split(self, bond_idx: int) -> Dict[str, float]:
        """
        Compute bonding/antibonding energy splitting.

        Δε = E_antibonding - E_bonding

        Args:
            bond_idx: Bond index

        Returns:
            Dictionary with energies and splitting
        """
        # Placeholder - would compute from Hamiltonian eigenvalues
        return {
            'bonding_energy': -1.0,
            'antibonding_energy': 1.0,
            'splitting': 2.0
        }

    def __repr__(self) -> str:
        """String representation."""
        return f"LCAO(n_qubits={self.n_qubits}, n_orbitals={self.n_orbitals}, n_electrons={self.n_electrons})"
