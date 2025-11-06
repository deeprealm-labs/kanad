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
        """
        Auto-detect appropriate governance protocol based on molecule.

        Uses electronegativity differences to classify bonding:
        - ΔEN > 1.7: Ionic
        - ΔEN < 0.4 and contains metals: Metallic
        - Otherwise: Covalent
        """
        # Import here to avoid circular dependency
        from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol
        from kanad.governance.protocols.ionic_protocol import IonicGovernanceProtocol
        from kanad.governance.protocols.metallic_protocol import MetallicGovernanceProtocol

        # Pauling electronegativities
        ELECTRONEGATIVITY = {
            'H': 2.20, 'He': 0.0,
            'Li': 0.98, 'Be': 1.57, 'B': 2.04, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
            'Na': 0.93, 'Mg': 1.31, 'Al': 1.61, 'Si': 1.90, 'P': 2.19, 'S': 2.58, 'Cl': 3.16,
            'K': 0.82, 'Ca': 1.00, 'Fe': 1.83, 'Cu': 1.90, 'Zn': 1.65,
            'Ag': 1.93, 'Au': 2.54
        }

        METALS = {'Li', 'Na', 'K', 'Mg', 'Ca', 'Al', 'Fe', 'Cu', 'Zn', 'Ag', 'Au'}

        # Calculate average electronegativity difference
        en_differences = []
        has_metal = False

        for atom in self.molecule.atoms:
            if atom.symbol in METALS:
                has_metal = True

        # Calculate EN differences between bonded atoms (use molecular structure if available)
        # For simplicity, calculate for all atom pairs (approximation)
        for i, atom1 in enumerate(self.molecule.atoms):
            for atom2 in self.molecule.atoms[i+1:]:
                en1 = ELECTRONEGATIVITY.get(atom1.symbol, 2.0)  # Default to C if unknown
                en2 = ELECTRONEGATIVITY.get(atom2.symbol, 2.0)
                en_diff = abs(en1 - en2)
                en_differences.append(en_diff)

        if not en_differences:
            # Single atom - use covalent
            logger.info("Auto-detecting bond type: single atom, using covalent")
            return CovalentGovernanceProtocol()

        max_en_diff = max(en_differences)
        avg_en_diff = sum(en_differences) / len(en_differences)

        # Classification rules
        if max_en_diff > 1.7:
            # Ionic bonding
            logger.info(f"Auto-detecting bond type: max ΔEN = {max_en_diff:.2f} > 1.7 → ionic")
            return IonicGovernanceProtocol()
        elif avg_en_diff < 0.4 and has_metal:
            # Metallic bonding
            logger.info(f"Auto-detecting bond type: avg ΔEN = {avg_en_diff:.2f} < 0.4, has metal → metallic")
            return MetallicGovernanceProtocol()
        else:
            # Covalent bonding (default)
            logger.info(f"Auto-detecting bond type: ΔEN = {avg_en_diff:.2f} → covalent")
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
        - Freeze deep core orbitals (like covalent)
        - Keep HOMO/LUMO orbitals (involved in charge transfer)
        - Keep all valence orbitals on both atoms

        For ionic bonds (e.g., LiF, NaCl), the key orbitals are:
        - Donor atom's HOMO (gives electron)
        - Acceptor atom's LUMO (receives electron)
        - Neighboring orbitals for correlation

        Returns:
            (frozen_orbitals, active_orbitals)

        Example for LiF:
            - Freeze: Li 1s, F 1s (2 orbitals)
            - Active: Li 2s, F 2s, 2p (valence involved in charge transfer)
            - Result: Fewer qubits, focus on charge transfer
        """
        # For ionic bonds, use same freezing strategy as covalent
        # (freeze core, keep valence)
        # The difference is in circuit construction (handled by governance protocol)
        frozen = []

        try:
            n_orbitals = self.molecule.n_orbitals
        except AttributeError:
            n_orbitals = sum(self._count_orbitals(atom.symbol) for atom in self.molecule.atoms)

        # Freeze core orbitals (same as covalent)
        orbital_idx = 0
        for atom in self.molecule.atoms:
            symbol = atom.symbol

            if symbol in ['H', 'He']:
                orbital_idx += 1  # No core to freeze

            elif symbol in ['Li', 'Be']:
                frozen.append(orbital_idx)
                orbital_idx += 1  # 1s frozen
                orbital_idx += 4  # 2s, 2p valence (active)

            elif symbol in ['O', 'N', 'C', 'F', 'Ne', 'B']:
                frozen.append(orbital_idx)
                orbital_idx += 1  # 1s frozen
                orbital_idx += 4  # 2s, 2p valence (active)

            elif symbol in ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']:
                # Freeze 1s, 2s, 2p (5 orbitals)
                frozen.extend([orbital_idx + i for i in range(5)])
                orbital_idx += 5
                orbital_idx += 4  # 3s, 3p valence (active)

            else:
                logger.warning(f"Unknown atom {symbol}, not freezing orbitals")

        all_orbitals = list(range(n_orbitals))
        active = [i for i in all_orbitals if i not in frozen]

        logger.info(f"Ionic active space: {len(frozen)} frozen, {len(active)} active")
        logger.info(f"  Frozen orbitals: {frozen}")
        logger.info(f"  Active orbitals: {active}")
        logger.info(f"  Qubits: {len(active) * 2} (reduced from {len(all_orbitals) * 2})")
        logger.info(f"  Strategy: Keep HOMO/LUMO for charge transfer")

        return frozen, active

    def _get_metallic_active_space(self) -> Tuple[List[int], List[int]]:
        """
        Get active space for metallic bonding.

        Strategy:
        - Freeze deep core orbitals
        - Keep conduction band orbitals (delocalized states)
        - For metallic systems, valence electrons are delocalized

        For metallic bonds (e.g., metal clusters), the key orbitals are:
        - Conduction band (valence s, p orbitals)
        - Partially occupied d orbitals (for transition metals)

        Returns:
            (frozen_orbitals, active_orbitals)

        Example for Na2:
            - Freeze: 1s, 2s, 2p on each Na (10 orbitals total)
            - Active: 3s on each Na (conduction band)
            - Result: Focus on delocalized metallic bonding
        """
        frozen = []

        try:
            n_orbitals = self.molecule.n_orbitals
        except AttributeError:
            n_orbitals = sum(self._count_orbitals(atom.symbol) for atom in self.molecule.atoms)

        # For metallic bonds, freeze deep core, keep outer valence
        orbital_idx = 0
        for atom in self.molecule.atoms:
            symbol = atom.symbol

            if symbol in ['H', 'He']:
                orbital_idx += 1  # No core, all active

            elif symbol in ['Li', 'Be']:
                frozen.append(orbital_idx)
                orbital_idx += 1  # 1s frozen
                orbital_idx += 4  # 2s, 2p active (conduction band)

            elif symbol in ['O', 'N', 'C', 'F', 'Ne', 'B']:
                frozen.append(orbital_idx)
                orbital_idx += 1  # 1s frozen
                orbital_idx += 4  # 2s, 2p active

            elif symbol in ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']:
                # For metals, freeze inner shells (1s, 2s, 2p)
                frozen.extend([orbital_idx + i for i in range(5)])
                orbital_idx += 5
                orbital_idx += 4  # 3s, 3p active (conduction band)

            else:
                logger.warning(f"Unknown atom {symbol}, not freezing orbitals")

        all_orbitals = list(range(n_orbitals))
        active = [i for i in all_orbitals if i not in frozen]

        logger.info(f"Metallic active space: {len(frozen)} frozen, {len(active)} active")
        logger.info(f"  Frozen orbitals: {frozen}")
        logger.info(f"  Active orbitals: {active}")
        logger.info(f"  Qubits: {len(active) * 2} (reduced from {len(all_orbitals) * 2})")
        logger.info(f"  Strategy: Keep conduction band orbitals (delocalized)")

        return frozen, active

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
