"""
Basis Set Registry - Comprehensive basis set support using PySCF.

Provides access to all PySCF basis sets plus custom implementations.
"""

from typing import List, Set, Optional
import logging

logger = logging.getLogger(__name__)

# Try to import PySCF for comprehensive basis set support
try:
    from pyscf import gto
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False
    logger.warning("PySCF not available - only built-in basis sets (sto-3g, 6-31g) supported")


class BasisSetRegistry:
    """
    Registry of available basis sets.

    Provides unified access to both built-in and PySCF basis sets.
    """

    # Built-in basis sets (implemented in basis_sets.py)
    BUILTIN_BASIS_SETS = {
        'sto-3g': 'Minimal basis (3 Gaussians per STO)',
        '6-31g': 'Split-valence basis',
    }

    # Common PySCF basis sets (if available)
    COMMON_PYSCF_BASIS_SETS = {
        # Minimal basis sets
        'sto-3g': 'STO-3G minimal basis',
        'sto-6g': 'STO-6G minimal basis',
        'sto3g': 'STO-3G (alternate name)',
        'sto6g': 'STO-6G (alternate name)',
        'minao': 'Minimal atomic orbital basis',

        # Split-valence basis sets
        '3-21g': 'Split-valence 3-21G',
        '6-31g': 'Split-valence 6-31G',
        '6-311g': 'Triple-zeta split-valence',

        # Polarization basis sets
        '6-31g*': '6-31G with d-polarization on heavy atoms',
        '6-31g**': '6-31G with d on heavy, p on H',
        '6-31g(d)': '6-31G with d-polarization (alias for 6-31G*)',
        '6-31g(d,p)': '6-31G with d on heavy, p on H (alias for 6-31G**)',
        '6-311g*': '6-311G with d-polarization',
        '6-311g**': '6-311G with d on heavy, p on H',
        '6-311g(d)': '6-311G with d-polarization',
        '6-311g(d,p)': '6-311G with d on heavy, p on H',

        # Diffuse function basis sets
        '6-31+g': '6-31G with diffuse functions on heavy atoms',
        '6-31++g': '6-31G with diffuse on all atoms',
        '6-31+g*': '6-31+G with polarization',
        '6-31++g**': '6-31++G with full polarization',

        # Dunning correlation-consistent basis sets
        'cc-pvdz': 'Correlation-consistent double-zeta',
        'cc-pvtz': 'Correlation-consistent triple-zeta',
        'cc-pvqz': 'Correlation-consistent quadruple-zeta',
        'cc-pv5z': 'Correlation-consistent quintuple-zeta',
        'aug-cc-pvdz': 'Augmented cc-pVDZ (with diffuse)',
        'aug-cc-pvtz': 'Augmented cc-pVTZ',
        'aug-cc-pvqz': 'Augmented cc-pVQZ',

        # Ahlrichs basis sets
        'def2-svp': 'Def2 split-valence + polarization',
        'def2-svpd': 'Def2-SVP with diffuse',
        'def2-tzvp': 'Def2 triple-zeta + polarization',
        'def2-tzvpd': 'Def2-TZVP with diffuse',
        'def2-tzvpp': 'Def2 triple-zeta + double polarization',
        'def2-qzvp': 'Def2 quadruple-zeta + polarization',

        # Pople basis sets with polarization
        '4-31g': 'Pople 4-31G',
        '6-21g': 'Pople 6-21G',

        # ANO basis sets
        'ano-rcc-vdz': 'ANO-RCC double-zeta',
        'ano-rcc-vtz': 'ANO-RCC triple-zeta',

        # Periodic system basis sets
        'gth-dzvp': 'GTH double-zeta + polarization (for PBC)',
        'gth-tzvp': 'GTH triple-zeta + polarization (for PBC)',
    }

    @classmethod
    def list_available_basis_sets(cls) -> List[str]:
        """Get list of all available basis sets."""
        available = list(cls.BUILTIN_BASIS_SETS.keys())

        if PYSCF_AVAILABLE:
            # Add all PySCF basis sets
            available.extend(cls.COMMON_PYSCF_BASIS_SETS.keys())

        return sorted(set(available))

    @classmethod
    def get_basis_description(cls, basis_name: str) -> Optional[str]:
        """Get description of a basis set."""
        basis_lower = basis_name.lower()

        if basis_lower in cls.BUILTIN_BASIS_SETS:
            return cls.BUILTIN_BASIS_SETS[basis_lower]

        if basis_lower in cls.COMMON_PYSCF_BASIS_SETS:
            return cls.COMMON_PYSCF_BASIS_SETS[basis_lower]

        return None

    @classmethod
    def is_available(cls, basis_name: str) -> bool:
        """Check if a basis set is available."""
        basis_lower = basis_name.lower()

        # Check built-in
        if basis_lower in cls.BUILTIN_BASIS_SETS:
            return True

        # Check PySCF
        if PYSCF_AVAILABLE and basis_lower in cls.COMMON_PYSCF_BASIS_SETS:
            return True

        # Try PySCF anyway (might be a basis set we don't have in our list)
        if PYSCF_AVAILABLE:
            try:
                # Try to load basis for hydrogen as a test
                gto.basis.load(basis_name, 'H')
                return True
            except:
                pass

        return False

    @classmethod
    def validate_basis(cls, basis_name: str) -> str:
        """
        Validate and normalize basis set name.

        Args:
            basis_name: Requested basis set name

        Returns:
            Normalized basis set name

        Raises:
            ValueError: If basis set is not available
        """
        basis_lower = basis_name.lower()

        if not cls.is_available(basis_lower):
            available = cls.list_available_basis_sets()
            raise ValueError(
                f"Basis set '{basis_name}' not available.\n"
                f"Available basis sets: {', '.join(available[:10])}...\n"
                f"Use BasisSetRegistry.list_available_basis_sets() for full list"
            )

        return basis_lower

    @classmethod
    def recommend_basis(cls, purpose: str = 'general') -> str:
        """
        Recommend a basis set for a given purpose.

        Args:
            purpose: Purpose of calculation
                - 'minimal': Fastest, least accurate
                - 'general': Good balance of speed and accuracy
                - 'accurate': High accuracy
                - 'correlation': For correlated methods (MP2, CCSD)
                - 'diffuse': For anions, excited states
                - 'periodic': For periodic systems

        Returns:
            Recommended basis set name
        """
        recommendations = {
            'minimal': 'sto-3g',
            'general': '6-31g*' if PYSCF_AVAILABLE else '6-31g',
            'accurate': 'cc-pvtz' if PYSCF_AVAILABLE else '6-31g',
            'correlation': 'cc-pvdz' if PYSCF_AVAILABLE else '6-31g',
            'diffuse': 'aug-cc-pvdz' if PYSCF_AVAILABLE else '6-31++g**',
            'periodic': 'gth-dzvp' if PYSCF_AVAILABLE else 'sto-3g',
        }

        return recommendations.get(purpose, '6-31g')

    @classmethod
    def get_info(cls, basis_name: str) -> dict:
        """
        Get detailed information about a basis set.

        Returns:
            Dictionary with basis set information
        """
        basis_lower = basis_name.lower()

        info = {
            'name': basis_lower,
            'available': cls.is_available(basis_lower),
            'description': cls.get_basis_description(basis_lower),
            'source': None
        }

        if basis_lower in cls.BUILTIN_BASIS_SETS:
            info['source'] = 'built-in'
        elif PYSCF_AVAILABLE and basis_lower in cls.COMMON_PYSCF_BASIS_SETS:
            info['source'] = 'PySCF'

        return info


def print_available_basis_sets():
    """Print all available basis sets with descriptions."""
    print("\n" + "=" * 80)
    print("AVAILABLE BASIS SETS")
    print("=" * 80)

    registry = BasisSetRegistry()

    print("\nBuilt-in Basis Sets:")
    print("-" * 80)
    for name, desc in registry.BUILTIN_BASIS_SETS.items():
        print(f"  {name:<20} - {desc}")

    if PYSCF_AVAILABLE:
        print("\nPySCF Basis Sets:")
        print("-" * 80)

        # Group by type
        groups = {
            'Minimal': ['sto-3g', 'sto-6g', 'minao'],
            'Split-Valence': ['3-21g', '6-31g', '6-311g'],
            'Polarization': ['6-31g*', '6-31g**', '6-311g*', '6-311g**'],
            'Diffuse': ['6-31+g', '6-31++g', '6-31+g*', '6-31++g**'],
            'Correlation-Consistent': ['cc-pvdz', 'cc-pvtz', 'cc-pvqz', 'aug-cc-pvdz', 'aug-cc-pvtz'],
            'Def2': ['def2-svp', 'def2-tzvp', 'def2-tzvpp', 'def2-qzvp'],
            'Periodic': ['gth-dzvp', 'gth-tzvp'],
        }

        for group_name, basis_list in groups.items():
            print(f"\n  {group_name}:")
            for basis in basis_list:
                if basis in registry.COMMON_PYSCF_BASIS_SETS:
                    desc = registry.COMMON_PYSCF_BASIS_SETS[basis]
                    print(f"    {basis:<20} - {desc}")

    print("\n" + "=" * 80)
    print(f"Total available: {len(registry.list_available_basis_sets())} basis sets")
    print("=" * 80 + "\n")
