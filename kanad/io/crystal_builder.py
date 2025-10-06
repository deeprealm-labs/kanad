"""
Crystal Structure Builder.

Convenience functions to generate common crystal structures:
- Simple cubic, BCC, FCC, HCP
- Diamond, zinc blende, wurtzite
- Rocksalt (NaCl), CsCl, fluorite
"""

import numpy as np
from typing import Tuple, Optional, List, Dict, Any
import logging

from kanad.core.atom import Atom
from kanad.core.lattice import Lattice
from kanad.core.molecule import Molecule

logger = logging.getLogger(__name__)


def build_crystal(element: str,
                  lattice_type: str,
                  lattice_constant: float,
                  size: Tuple[int, int, int] = (1, 1, 1),
                  **kwargs) -> Molecule:
    """
    Build common crystal structures.

    Args:
        element: Element symbol (e.g., 'Cu', 'Si', 'C')
        lattice_type: Crystal structure type:
            - 'sc': Simple cubic
            - 'bcc': Body-centered cubic
            - 'fcc': Face-centered cubic
            - 'hcp': Hexagonal close-packed
            - 'diamond': Diamond structure
        lattice_constant: Lattice parameter in Angstrom
        size: Supercell size (nx, ny, nz)
        **kwargs: Additional parameters (k_points, basis, pseudo, etc.)

    Returns:
        crystal: Molecule object with periodic lattice

    Examples:
        >>> # Copper FCC crystal
        >>> cu = build_crystal('Cu', 'fcc', lattice_constant=3.61)

        >>> # Silicon diamond structure with 2x2x2 supercell
        >>> si = build_crystal('Si', 'diamond', lattice_constant=5.43, size=(2,2,2))

        >>> # Iron BCC with k-points
        >>> fe = build_crystal('Fe', 'bcc', lattice_constant=2.87, k_points=(4,4,4))
    """
    lattice_type = lattice_type.lower()

    if lattice_type == 'sc':
        crystal = _build_simple_cubic(element, lattice_constant)
    elif lattice_type == 'bcc':
        crystal = _build_bcc(element, lattice_constant)
    elif lattice_type == 'fcc':
        crystal = _build_fcc(element, lattice_constant)
    elif lattice_type == 'hcp':
        c_a_ratio = kwargs.pop('c_a_ratio', 1.633)  # Ideal HCP
        crystal = _build_hcp(element, lattice_constant, c_a_ratio)
    elif lattice_type == 'diamond':
        crystal = _build_diamond(element, lattice_constant)
    else:
        raise ValueError(f"Unknown lattice type: {lattice_type}")

    # Apply supercell expansion if needed
    if size != (1, 1, 1):
        crystal = crystal.make_supercell(size)

    # Apply k_points if specified
    if 'k_points' in kwargs:
        crystal.k_points = kwargs['k_points']

    logger.info(f"Built {lattice_type.upper()} crystal: {element}, a={lattice_constant:.3f} Å")

    return crystal


def build_binary_crystal(element_a: str,
                         element_b: str,
                         lattice_type: str,
                         lattice_constant: float,
                         size: Tuple[int, int, int] = (1, 1, 1),
                         **kwargs) -> Molecule:
    """
    Build binary compound crystals.

    Args:
        element_a: First element (cation typically)
        element_b: Second element (anion typically)
        lattice_type: Structure type:
            - 'rocksalt': NaCl structure
            - 'zincblende': ZnS (sphalerite) structure
            - 'wurtzite': ZnS (wurtzite) structure
            - 'cscl': CsCl structure
            - 'fluorite': CaF2 structure
        lattice_constant: Lattice parameter in Angstrom
        size: Supercell size
        **kwargs: Additional parameters

    Returns:
        crystal: Binary crystal Molecule

    Examples:
        >>> # Sodium chloride
        >>> nacl = build_binary_crystal('Na', 'Cl', 'rocksalt', lattice_constant=5.64)

        >>> # Gallium arsenide
        >>> gaas = build_binary_crystal('Ga', 'As', 'zincblende', lattice_constant=5.65)
    """
    lattice_type = lattice_type.lower()

    if lattice_type == 'rocksalt':
        crystal = _build_rocksalt(element_a, element_b, lattice_constant)
    elif lattice_type == 'zincblende':
        crystal = _build_zincblende(element_a, element_b, lattice_constant)
    elif lattice_type == 'wurtzite':
        u_param = kwargs.pop('u_param', 0.375)
        c_a_ratio = kwargs.pop('c_a_ratio', 1.633)
        crystal = _build_wurtzite(element_a, element_b, lattice_constant, u_param, c_a_ratio)
    elif lattice_type == 'cscl':
        crystal = _build_cscl(element_a, element_b, lattice_constant)
    elif lattice_type == 'fluorite':
        crystal = _build_fluorite(element_a, element_b, lattice_constant)
    else:
        raise ValueError(f"Unknown binary lattice type: {lattice_type}")

    if size != (1, 1, 1):
        crystal = crystal.make_supercell(size)

    logger.info(f"Built {lattice_type.upper()} crystal: {element_a}{element_b}, a={lattice_constant:.3f} Å")

    return crystal


# ===== Elemental Crystal Structures =====

def _build_simple_cubic(element: str, a: float) -> Molecule:
    """Simple cubic lattice (Po)."""
    lattice_vectors = np.array([
        [a, 0, 0],
        [0, a, 0],
        [0, 0, a]
    ])
    lattice = Lattice(lattice_vectors)

    atoms = [Atom(element, position=[0, 0, 0])]

    return Molecule(atoms, lattice=lattice, basis='gth-dzvp')


def _build_bcc(element: str, a: float) -> Molecule:
    """Body-centered cubic (Fe, Cr, W, Mo)."""
    lattice_vectors = np.array([
        [a, 0, 0],
        [0, a, 0],
        [0, 0, a]
    ])
    lattice = Lattice(lattice_vectors)

    # Atom at (0,0,0) and body center (a/2, a/2, a/2)
    atoms = [
        Atom(element, position=[0, 0, 0]),
        Atom(element, position=[a/2, a/2, a/2])
    ]

    return Molecule(atoms, lattice=lattice, basis='gth-dzvp')


def _build_fcc(element: str, a: float) -> Molecule:
    """Face-centered cubic (Cu, Al, Ag, Au, Ni, Pb)."""
    lattice_vectors = np.array([
        [0, a/2, a/2],
        [a/2, 0, a/2],
        [a/2, a/2, 0]
    ])
    lattice = Lattice(lattice_vectors)

    # Single atom at origin (FCC primitive cell)
    atoms = [Atom(element, position=[0, 0, 0])]

    return Molecule(atoms, lattice=lattice, basis='gth-dzvp')


def _build_hcp(element: str, a: float, c_a_ratio: float = 1.633) -> Molecule:
    """
    Hexagonal close-packed (Mg, Zn, Ti, Cd).

    Args:
        a: In-plane lattice constant
        c_a_ratio: c/a ratio (ideal = sqrt(8/3) ≈ 1.633)
    """
    c = a * c_a_ratio

    # Hexagonal lattice vectors
    lattice_vectors = np.array([
        [a, 0, 0],
        [-a/2, a * np.sqrt(3)/2, 0],
        [0, 0, c]
    ])
    lattice = Lattice(lattice_vectors, pbc=(True, True, True))

    # Two atoms per unit cell
    atoms = [
        Atom(element, position=[0, 0, 0]),
        Atom(element, position=[a/2, a/(2*np.sqrt(3)), c/2])
    ]

    return Molecule(atoms, lattice=lattice, basis='gth-dzvp')


def _build_diamond(element: str, a: float) -> Molecule:
    """
    Diamond structure (C, Si, Ge).

    FCC lattice with 2-atom basis.
    """
    # FCC lattice vectors
    lattice_vectors = np.array([
        [0, a/2, a/2],
        [a/2, 0, a/2],
        [a/2, a/2, 0]
    ])
    lattice = Lattice(lattice_vectors)

    # Two atoms: (0,0,0) and (a/4, a/4, a/4)
    atoms = [
        Atom(element, position=[0, 0, 0]),
        Atom(element, position=[a/4, a/4, a/4])
    ]

    return Molecule(atoms, lattice=lattice, basis='gth-dzvp')


# ===== Binary Compound Structures =====

def _build_rocksalt(element_a: str, element_b: str, a: float) -> Molecule:
    """
    Rocksalt structure (NaCl, LiF, MgO).

    Two interpenetrating FCC lattices.
    """
    lattice_vectors = np.array([
        [a, 0, 0],
        [0, a, 0],
        [0, 0, a]
    ])
    lattice = Lattice(lattice_vectors)

    # A at (0,0,0), B at (a/2, 0, 0)
    atoms = [
        Atom(element_a, position=[0, 0, 0]),
        Atom(element_b, position=[a/2, 0, 0])
    ]

    return Molecule(atoms, lattice=lattice, basis='gth-dzvp')


def _build_zincblende(element_a: str, element_b: str, a: float) -> Molecule:
    """
    Zinc blende (sphalerite) structure (GaAs, ZnS, InP).

    Diamond structure with two different atoms.
    """
    lattice_vectors = np.array([
        [0, a/2, a/2],
        [a/2, 0, a/2],
        [a/2, a/2, 0]
    ])
    lattice = Lattice(lattice_vectors)

    atoms = [
        Atom(element_a, position=[0, 0, 0]),
        Atom(element_b, position=[a/4, a/4, a/4])
    ]

    return Molecule(atoms, lattice=lattice, basis='gth-dzvp')


def _build_wurtzite(element_a: str, element_b: str, a: float,
                    u_param: float = 0.375, c_a_ratio: float = 1.633) -> Molecule:
    """
    Wurtzite structure (ZnO, GaN, CdS).

    Hexagonal lattice with 4-atom basis.

    Args:
        u_param: Internal parameter (ideal = 3/8 = 0.375)
        c_a_ratio: c/a ratio
    """
    c = a * c_a_ratio

    lattice_vectors = np.array([
        [a, 0, 0],
        [-a/2, a * np.sqrt(3)/2, 0],
        [0, 0, c]
    ])
    lattice = Lattice(lattice_vectors)

    # 4 atoms per cell
    atoms = [
        Atom(element_a, position=[0, 0, 0]),
        Atom(element_a, position=[a/2, a/(2*np.sqrt(3)), c/2]),
        Atom(element_b, position=[0, 0, u_param * c]),
        Atom(element_b, position=[a/2, a/(2*np.sqrt(3)), c/2 + u_param * c])
    ]

    return Molecule(atoms, lattice=lattice, basis='gth-dzvp')


def _build_cscl(element_a: str, element_b: str, a: float) -> Molecule:
    """
    CsCl structure.

    Simple cubic with 2-atom basis.
    """
    lattice_vectors = np.array([
        [a, 0, 0],
        [0, a, 0],
        [0, 0, a]
    ])
    lattice = Lattice(lattice_vectors)

    atoms = [
        Atom(element_a, position=[0, 0, 0]),
        Atom(element_b, position=[a/2, a/2, a/2])
    ]

    return Molecule(atoms, lattice=lattice, basis='gth-dzvp')


def _build_fluorite(element_a: str, element_b: str, a: float) -> Molecule:
    """
    Fluorite structure (CaF2).

    FCC lattice of A, with B at tetrahedral sites.
    """
    lattice_vectors = np.array([
        [a, 0, 0],
        [0, a, 0],
        [0, 0, a]
    ])
    lattice = Lattice(lattice_vectors)

    # A at FCC positions, B at (a/4, a/4, a/4) and (3a/4, 3a/4, 3a/4)
    atoms = [
        Atom(element_a, position=[0, 0, 0]),
        Atom(element_b, position=[a/4, a/4, a/4]),
        Atom(element_b, position=[3*a/4, 3*a/4, 3*a/4])
    ]

    return Molecule(atoms, lattice=lattice, basis='gth-dzvp')


# ===== High-Symmetry k-Paths =====

def get_kpath(lattice_type: str, n_points: int = 50) -> Tuple[np.ndarray, List[str], List[float]]:
    """
    Generate high-symmetry k-point path for band structure.

    Args:
        lattice_type: 'fcc', 'bcc', 'sc', 'hexagonal', 'diamond'
        n_points: Points per segment

    Returns:
        k_points: (N, 3) array of k-points
        labels: List of high-symmetry point names
        label_positions: Cumulative distances for labels

    Examples:
        >>> k_path, labels, positions = get_kpath('fcc', n_points=50)
        >>> # Use for band structure calculation
        >>> bands = crystal.compute_band_structure(k_path, n_bands=8)
    """
    lattice_type = lattice_type.lower()

    if lattice_type in ['fcc', 'diamond']:
        # FCC Brillouin zone: Γ → X → W → K → Γ → L → U → W → L → K | U → X
        return _kpath_fcc(n_points)
    elif lattice_type == 'bcc':
        # BCC: Γ → H → N → Γ → P → H | P → N
        return _kpath_bcc(n_points)
    elif lattice_type == 'sc':
        # Simple cubic: Γ → X → M → Γ → R → X | M → R
        return _kpath_sc(n_points)
    elif lattice_type == 'hexagonal':
        # Hexagonal: Γ → M → K → Γ → A → L → H → A | L → M | K → H
        return _kpath_hexagonal(n_points)
    else:
        raise ValueError(f"Unknown lattice type for k-path: {lattice_type}")


def _kpath_fcc(n_points: int) -> Tuple[np.ndarray, List[str], List[float]]:
    """FCC k-path: Γ → X → W → K → Γ → L → W"""
    # High-symmetry points (in units of 2π/a)
    Gamma = np.array([0.0, 0.0, 0.0])
    X = np.array([0.5, 0.0, 0.5])
    W = np.array([0.5, 0.25, 0.75])
    K = np.array([0.375, 0.375, 0.75])
    L = np.array([0.5, 0.5, 0.5])

    # Path segments
    segments = [
        (Gamma, X),
        (X, W),
        (W, K),
        (K, Gamma),
        (Gamma, L),
        (L, W)
    ]

    k_points, labels, positions = _build_kpath_from_segments(segments, n_points)

    return k_points, labels, positions


def _kpath_bcc(n_points: int) -> Tuple[np.ndarray, List[str], List[float]]:
    """BCC k-path: Γ → H → N → Γ → P → H"""
    Gamma = np.array([0.0, 0.0, 0.0])
    H = np.array([0.5, -0.5, 0.5])
    N = np.array([0.0, 0.0, 0.5])
    P = np.array([0.25, 0.25, 0.25])

    segments = [
        (Gamma, H),
        (H, N),
        (N, Gamma),
        (Gamma, P),
        (P, H)
    ]

    return _build_kpath_from_segments(segments, n_points)


def _kpath_sc(n_points: int) -> Tuple[np.ndarray, List[str], List[float]]:
    """Simple cubic: Γ → X → M → Γ → R → X"""
    Gamma = np.array([0.0, 0.0, 0.0])
    X = np.array([0.5, 0.0, 0.0])
    M = np.array([0.5, 0.5, 0.0])
    R = np.array([0.5, 0.5, 0.5])

    segments = [
        (Gamma, X),
        (X, M),
        (M, Gamma),
        (Gamma, R),
        (R, X)
    ]

    return _build_kpath_from_segments(segments, n_points)


def _kpath_hexagonal(n_points: int) -> Tuple[np.ndarray, List[str], List[float]]:
    """Hexagonal: Γ → M → K → Γ → A"""
    Gamma = np.array([0.0, 0.0, 0.0])
    M = np.array([0.5, 0.0, 0.0])
    K = np.array([1./3., 1./3., 0.0])
    A = np.array([0.0, 0.0, 0.5])

    segments = [
        (Gamma, M),
        (M, K),
        (K, Gamma),
        (Gamma, A)
    ]

    return _build_kpath_from_segments(segments, n_points)


def _build_kpath_from_segments(segments: List[Tuple[np.ndarray, np.ndarray]],
                                n_points: int) -> Tuple[np.ndarray, List[str], List[float]]:
    """
    Build k-path from segments.

    Args:
        segments: List of (k_start, k_end) tuples
        n_points: Points per segment

    Returns:
        k_points, labels, label_positions
    """
    all_k_points = []
    label_positions = [0.0]
    cumulative_distance = 0.0

    # Symbol mapping
    symbol_map = {
        tuple([0, 0, 0]): 'Γ',
        tuple([0.5, 0, 0.5]): 'X',
        tuple([0.5, 0.25, 0.75]): 'W',
        tuple([0.375, 0.375, 0.75]): 'K',
        tuple([0.5, 0.5, 0.5]): 'L',
        tuple([0.5, -0.5, 0.5]): 'H',
        tuple([0.0, 0.0, 0.5]): 'N',
        tuple([0.25, 0.25, 0.25]): 'P',
        tuple([0.5, 0.5, 0.0]): 'M',
        tuple([0.5, 0.0, 0.0]): 'X',
        tuple([1./3., 1./3., 0.0]): 'K',
        tuple([0.0, 0.0, 0.5]): 'A',
    }

    def get_label(k):
        k_tuple = tuple(np.round(k, 3))
        return symbol_map.get(k_tuple, '?')

    labels = [get_label(segments[0][0])]

    for i, (k_start, k_end) in enumerate(segments):
        # Generate points along segment
        segment_k = np.linspace(k_start, k_end, n_points)

        if i == 0:
            all_k_points.extend(segment_k)
        else:
            all_k_points.extend(segment_k[1:])  # Skip duplicate point

        # Update cumulative distance
        segment_distance = np.linalg.norm(k_end - k_start)
        cumulative_distance += segment_distance
        label_positions.append(cumulative_distance)
        labels.append(get_label(k_end))

    k_points = np.array(all_k_points)

    return k_points, labels, label_positions


# ===== Utility Functions =====

def get_lattice_info(lattice_type: str) -> Dict[str, Any]:
    """
    Get information about a crystal structure.

    Args:
        lattice_type: Structure type

    Returns:
        info: Dictionary with:
            - coordination: Coordination number
            - description: Text description
            - examples: List of example materials
    """
    info_dict = {
        'sc': {
            'coordination': 6,
            'description': 'Simple cubic',
            'examples': ['Po (polonium)']
        },
        'bcc': {
            'coordination': 8,
            'description': 'Body-centered cubic',
            'examples': ['Fe', 'Cr', 'W', 'Mo', 'V']
        },
        'fcc': {
            'coordination': 12,
            'description': 'Face-centered cubic',
            'examples': ['Cu', 'Al', 'Ag', 'Au', 'Ni', 'Pb']
        },
        'hcp': {
            'coordination': 12,
            'description': 'Hexagonal close-packed',
            'examples': ['Mg', 'Zn', 'Ti', 'Cd']
        },
        'diamond': {
            'coordination': 4,
            'description': 'Diamond structure (FCC with 2-atom basis)',
            'examples': ['C (diamond)', 'Si', 'Ge']
        },
        'rocksalt': {
            'coordination': 6,
            'description': 'Rocksalt (NaCl) structure',
            'examples': ['NaCl', 'LiF', 'MgO', 'TiN']
        },
        'zincblende': {
            'coordination': 4,
            'description': 'Zinc blende (sphalerite)',
            'examples': ['GaAs', 'ZnS', 'InP', 'GaN']
        },
        'wurtzite': {
            'coordination': 4,
            'description': 'Wurtzite (hexagonal)',
            'examples': ['ZnO', 'GaN', 'CdS']
        },
        'cscl': {
            'coordination': 8,
            'description': 'CsCl structure',
            'examples': ['CsCl', 'CsBr', 'CsI']
        }
    }

    return info_dict.get(lattice_type.lower(), {'description': 'Unknown'})
