"""
Active Space Selection with Governance Protocol Integration

This module provides intelligent active space selection using governance protocols
to identify core vs bonding orbitals based on bonding physics.

Key Benefits:
- Qubit reduction: H2O goes from 14 qubits to 12 qubits
- Physics-guided: Uses bonding type to determine frozen orbitals
- Governance integration: Leverages existing governance protocols
"""

import logging
from typing import List, Tuple, Optional
import numpy as np

logger = logging.getLogger(__name__)


class ActiveSpaceSelector:
    """
    Select active orbitals using governance protocol guidance.

    Different bond types require different active space strategies:
    - Covalent: Freeze core, keep valence orbitals involved in bonding
    - Ionic: Keep orbitals involved in charge transfer
    - Metallic: Keep conduction band orbitals
    """

    def __init__(self, molecule, protocol=None):
        """
        Initialize active space selector.

        Args:
            molecule: Molecule object with atoms and orbitals
            protocol: GovernanceProtocol instance (auto-detected if None)
        """
        self.molecule = molecule
        self.protocol = protocol

        # Auto-detect protocol if not provided
        if self.protocol is None:
            self.protocol = self._detect_protocol()

    def _detect_protocol(self):
        """Auto-detect appropriate governance protocol based on molecule."""
        # Import here to avoid circular dependency
        from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol

        # For now, default to covalent (most common)
        # TODO: Implement bond type detection based on electronegativity difference
        logger.info("Auto-detecting bond type: defaulting to covalent")
        return CovalentGovernanceProtocol()

    def get_active_space(self) -> Tuple[List[int], List[int]]:
        """
        Get frozen and active orbital indices using governance protocol.

        Returns:
            Tuple of (frozen_orbitals, active_orbitals)

        Example for H2O with CovalentGovernanceProtocol:
            - Frozen: [0] (O 1s core)
            - Active: [1,2,3,4,5,6] (O 2s, 2p and H 1s valence)
            - Result: 12 qubits instead of 14
        """
        bond_type = getattr(self.protocol, 'bond_type', None)

        # Convert BondingType enum to string if necessary
        if bond_type is not None and hasattr(bond_type, 'value'):
            bond_type_str = str(bond_type.value).lower()
        elif isinstance(bond_type, str):
            bond_type_str = bond_type.lower()
        else:
            bond_type_str = 'covalent'  # Default

        if bond_type_str == 'covalent':
            return self._get_covalent_active_space()
        elif bond_type_str == 'ionic':
            return self._get_ionic_active_space()
        elif bond_type_str == 'metallic':
            return self._get_metallic_active_space()
        else:
            # Fallback: no frozen orbitals (compute orbital count without triggering hamiltonian)
            logger.warning(f"Unknown bond type {bond_type_str}, using all orbitals")
            n_orbitals = sum(self._count_orbitals(atom.symbol) for atom in self.molecule.atoms)
            all_orbitals = list(range(n_orbitals))
            return [], all_orbitals

    def _get_covalent_active_space(self) -> Tuple[List[int], List[int]]:
        """
        Get active space for covalent bonding.

        Strategy:
        - Freeze core orbitals (1s for second row: O, N, C, etc.)
        - Keep valence orbitals involved in bonding

        Returns:
            (frozen_orbitals, active_orbitals)
        """
        frozen = []

        # Get total orbitals from molecule (triggers hamiltonian construction if needed)
        try:
            n_orbitals = self.molecule.n_orbitals
        except AttributeError:
            # If molecule doesn't have n_orbitals yet, compute it from atoms
            # For STO-3G: H has 1 orbital, O/N/C/F have 5 orbitals, etc.
            n_orbitals = sum(self._count_orbitals(atom.symbol) for atom in self.molecule.atoms)

        # Identify core orbitals to freeze
        orbital_idx = 0
        for atom in self.molecule.atoms:
            symbol = atom.symbol

            # First row elements (H, He) have no core
            if symbol in ['H', 'He']:
                orbital_idx += 1  # 1s orbital (valence, not frozen)

            # Li, Be have 1s core (freeze it)
            elif symbol in ['Li', 'Be']:
                frozen.append(orbital_idx)
                orbital_idx += 1  # 1s orbital (frozen)
                # 2s, 2px, 2py, 2pz are valence (not frozen)
                orbital_idx += 4

            # Second row elements (O, N, C, F, B) have 1s core
            elif symbol in ['O', 'N', 'C', 'F', 'Ne', 'B']:
                frozen.append(orbital_idx)
                orbital_idx += 1  # 1s orbital (frozen)
                # 2s, 2px, 2py, 2pz are valence (not frozen)
                orbital_idx += 4

            # Third row elements (Na-Ar) have 1s, 2s, 2p core
            elif symbol in ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']:
                # Freeze 1s, 2s, 2p (5 orbitals)
                frozen.extend([orbital_idx + i for i in range(5)])
                orbital_idx += 5
                # 3s, 3px, 3py, 3pz are valence (not frozen)
                orbital_idx += 4

            else:
                # Unknown atom, don't freeze anything
                logger.warning(f"Unknown atom {symbol}, not freezing orbitals")

        # Active orbitals are all non-frozen orbitals
        all_orbitals = list(range(n_orbitals))
        active = [i for i in all_orbitals if i not in frozen]

        logger.info(f"Covalent active space: {len(frozen)} frozen, {len(active)} active")
        logger.info(f"  Frozen orbitals: {frozen}")
        logger.info(f"  Active orbitals: {active}")
        logger.info(f"  Qubits: {len(active) * 2} (reduced from {len(all_orbitals) * 2})")

        return frozen, active

    def _count_orbitals(self, symbol: str) -> int:
        """Count number of orbitals for an atom (STO-3G basis)."""
        # STO-3G: minimal basis
        if symbol in ['H', 'He']:
            return 1  # 1s
        elif symbol in ['Li', 'Be']:
            return 5  # 1s, 2s, 2px, 2py, 2pz (first row after H/He)
        elif symbol in ['O', 'N', 'C', 'F', 'Ne', 'B']:
            return 5  # 1s, 2s, 2px, 2py, 2pz
        elif symbol in ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']:
            return 9  # 1s, 2s, 2px, 2py, 2pz, 3s, 3px, 3py, 3pz
        else:
            logger.warning(f"Unknown atom {symbol}, assuming 1 orbital")
            return 1

    def _get_ionic_active_space(self) -> Tuple[List[int], List[int]]:
        """
        Get active space for ionic bonding.

        Strategy:
        - Freeze deep core orbitals
        - Keep orbitals involved in charge transfer (valence of both atoms)

        TODO: Implement ionic-specific logic
        """
        logger.warning("Ionic active space not yet implemented, using covalent strategy")
        return self._get_covalent_active_space()

    def _get_metallic_active_space(self) -> Tuple[List[int], List[int]]:
        """
        Get active space for metallic bonding.

        Strategy:
        - Freeze deep core orbitals
        - Keep conduction band orbitals (delocalized states)

        TODO: Implement metallic-specific logic
        """
        logger.warning("Metallic active space not yet implemented, using covalent strategy")
        return self._get_covalent_active_space()

    def get_active_space_electrons(self, frozen_orbitals: List[int]) -> int:
        """
        Calculate number of electrons in active space.

        Args:
            frozen_orbitals: List of frozen orbital indices

        Returns:
            Number of electrons in active space
        """
        # Each frozen orbital contributes 2 electrons (assumed doubly occupied)
        frozen_electrons = len(frozen_orbitals) * 2
        active_electrons = self.molecule.n_electrons - frozen_electrons

        logger.info(f"Active space electrons: {active_electrons} (total: {self.molecule.n_electrons}, frozen: {frozen_electrons})")

        return active_electrons


def get_governance_active_space(molecule, protocol=None) -> Tuple[List[int], List[int], int]:
    """
    Convenience function to get governance-guided active space.

    Args:
        molecule: Molecule object
        protocol: GovernanceProtocol instance (auto-detected if None)

    Returns:
        Tuple of (frozen_orbitals, active_orbitals, n_active_electrons)

    Example for H2O:
        >>> frozen, active, n_electrons = get_governance_active_space(h2o_molecule)
        >>> # frozen = [0], active = [1,2,3,4,5,6], n_electrons = 8
        >>> # Result: 12 qubits instead of 14!
    """
    selector = ActiveSpaceSelector(molecule, protocol)
    frozen, active = selector.get_active_space()
    n_active_electrons = selector.get_active_space_electrons(frozen)

    return frozen, active, n_active_electrons
