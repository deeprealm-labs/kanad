"""
Crystal Lattice for Periodic Boundary Conditions.

Provides lattice vectors, reciprocal space, minimum image convention,
and supercell generation for periodic systems.
"""

import numpy as np
from typing import Tuple, Optional, List
import logging

logger = logging.getLogger(__name__)


class Lattice:
    """
    Crystal lattice with periodic boundary conditions.

    Supports 1D, 2D, and 3D periodicity for solid-state calculations.
    """

    def __init__(self,
                 lattice_vectors: np.ndarray,
                 pbc: Tuple[bool, bool, bool] = (True, True, True)):
        """
        Initialize crystal lattice.

        Args:
            lattice_vectors: (3, 3) array where rows are lattice vectors a, b, c
                            Units: Angstrom
            pbc: (bool, bool, bool) for periodicity in x, y, z directions
                 Examples:
                 - (True, True, True): 3D crystal
                 - (True, True, False): 2D sheet/slab
                 - (True, False, False): 1D chain/wire

        Examples:
            >>> # Simple cubic lattice, a = 3.0 Angstrom
            >>> lattice = Lattice(np.eye(3) * 3.0)

            >>> # FCC lattice (copper, a = 3.61 Angstrom)
            >>> a = 3.61
            >>> lattice_vectors = np.array([
            ...     [0, a/2, a/2],
            ...     [a/2, 0, a/2],
            ...     [a/2, a/2, 0]
            ... ])
            >>> lattice = Lattice(lattice_vectors)
        """
        self.lattice_vectors = np.array(lattice_vectors, dtype=float)
        assert self.lattice_vectors.shape == (3, 3), "Lattice vectors must be 3x3 array"

        self.pbc = tuple(pbc)
        assert len(self.pbc) == 3, "PBC must be 3-tuple of bools"

        # Compute derived quantities
        self._compute_reciprocal_vectors()
        self._compute_lattice_parameters()

        logger.debug(f"Created lattice: PBC={self.pbc}, volume={self.volume:.3f} Å³")

    def _compute_reciprocal_vectors(self):
        """
        Compute reciprocal lattice vectors.

        Definition:
            b_i = 2π (a_j × a_k) / V

        where V = a_1 · (a_2 × a_3) is the unit cell volume.

        Property: a_i · b_j = 2π δ_ij
        """
        a1, a2, a3 = self.lattice_vectors[0], self.lattice_vectors[1], self.lattice_vectors[2]

        # Volume (scalar triple product)
        self.volume = np.abs(np.dot(a1, np.cross(a2, a3)))

        if self.volume < 1e-10:
            raise ValueError("Lattice vectors are linearly dependent (zero volume)")

        # Reciprocal vectors
        b1 = 2 * np.pi * np.cross(a2, a3) / self.volume
        b2 = 2 * np.pi * np.cross(a3, a1) / self.volume
        b3 = 2 * np.pi * np.cross(a1, a2) / self.volume

        self.reciprocal_vectors = np.array([b1, b2, b3])

        # Verify orthogonality: a_i · b_j = 2π δ_ij
        dot_matrix = self.lattice_vectors @ self.reciprocal_vectors.T
        expected = 2 * np.pi * np.eye(3)
        if not np.allclose(dot_matrix, expected, atol=1e-8):
            logger.warning("Reciprocal vectors may be inaccurate")

    def _compute_lattice_parameters(self):
        """
        Compute lattice parameters: a, b, c, alpha, beta, gamma.

        a, b, c: lengths of lattice vectors
        alpha: angle between b and c
        beta:  angle between a and c
        gamma: angle between a and b
        """
        a_vec, b_vec, c_vec = self.lattice_vectors

        # Lengths
        self.a = np.linalg.norm(a_vec)
        self.b = np.linalg.norm(b_vec)
        self.c = np.linalg.norm(c_vec)

        # Angles (in degrees)
        self.alpha = np.degrees(np.arccos(np.dot(b_vec, c_vec) / (self.b * self.c)))
        self.beta = np.degrees(np.arccos(np.dot(a_vec, c_vec) / (self.a * self.c)))
        self.gamma = np.degrees(np.arccos(np.dot(a_vec, b_vec) / (self.a * self.b)))

    def get_reciprocal_vectors(self) -> np.ndarray:
        """
        Get reciprocal lattice vectors.

        Returns:
            reciprocal_vectors: (3, 3) array in units of 2π/Angstrom
        """
        return self.reciprocal_vectors.copy()

    def fractional_to_cartesian(self, fractional_coords: np.ndarray) -> np.ndarray:
        """
        Convert fractional coordinates to Cartesian.

        Args:
            fractional_coords: (N, 3) or (3,) array of fractional coordinates
                              (0 ≤ f_i < 1 for atoms inside unit cell)

        Returns:
            cartesian_coords: (N, 3) or (3,) array in Angstrom

        Formula:
            r = f_1 a + f_2 b + f_3 c
        """
        fractional = np.atleast_2d(fractional_coords)
        cartesian = fractional @ self.lattice_vectors

        if fractional_coords.ndim == 1:
            return cartesian[0]
        return cartesian

    def cartesian_to_fractional(self, cartesian_coords: np.ndarray) -> np.ndarray:
        """
        Convert Cartesian coordinates to fractional.

        Args:
            cartesian_coords: (N, 3) or (3,) array in Angstrom

        Returns:
            fractional_coords: (N, 3) or (3,) array

        Formula:
            f = r · (inverse of lattice matrix)
        """
        cartesian = np.atleast_2d(cartesian_coords)
        fractional = cartesian @ np.linalg.inv(self.lattice_vectors.T)

        if cartesian_coords.ndim == 1:
            return fractional[0]
        return fractional

    def wrap_to_unit_cell(self, positions: np.ndarray) -> np.ndarray:
        """
        Wrap positions into unit cell (0 ≤ f_i < 1).

        Args:
            positions: (N, 3) Cartesian coordinates in Angstrom

        Returns:
            wrapped_positions: (N, 3) Cartesian coordinates
        """
        fractional = self.cartesian_to_fractional(positions)

        # Apply PBC: wrap periodic dimensions to [0, 1)
        for i, periodic in enumerate(self.pbc):
            if periodic:
                fractional[:, i] = fractional[:, i] % 1.0

        return self.fractional_to_cartesian(fractional)

    def minimum_image_distance(self, r1: np.ndarray, r2: np.ndarray) -> float:
        """
        Compute minimum image convention distance.

        Accounts for periodic boundary conditions by considering all
        periodic images and returning the shortest distance.

        Args:
            r1, r2: Position vectors (3,) in Angstrom

        Returns:
            distance: Minimum distance in Angstrom
        """
        # Displacement vector
        dr = r2 - r1

        # Convert to fractional coordinates
        frac_dr = self.cartesian_to_fractional(dr)

        # Apply minimum image convention
        for i, periodic in enumerate(self.pbc):
            if periodic:
                # Wrap to [-0.5, 0.5)
                frac_dr[i] = frac_dr[i] - np.round(frac_dr[i])

        # Convert back to Cartesian
        min_dr = self.fractional_to_cartesian(frac_dr)

        return np.linalg.norm(min_dr)

    def minimum_image_vector(self, r1: np.ndarray, r2: np.ndarray) -> np.ndarray:
        """
        Compute minimum image displacement vector.

        Args:
            r1, r2: Position vectors (3,) in Angstrom

        Returns:
            dr_min: Minimum displacement vector (3,) in Angstrom
        """
        dr = r2 - r1
        frac_dr = self.cartesian_to_fractional(dr)

        for i, periodic in enumerate(self.pbc):
            if periodic:
                frac_dr[i] = frac_dr[i] - np.round(frac_dr[i])

        return self.fractional_to_cartesian(frac_dr)

    def make_supercell(self, size: Tuple[int, int, int]) -> 'Lattice':
        """
        Create supercell expansion.

        Args:
            size: (nx, ny, nz) number of unit cells in each direction

        Returns:
            supercell_lattice: New Lattice with expanded vectors

        Example:
            >>> # 2x2x2 supercell
            >>> supercell = lattice.make_supercell((2, 2, 2))
            >>> # Volume is 8× larger
            >>> assert np.isclose(supercell.volume, 8 * lattice.volume)
        """
        nx, ny, nz = size

        supercell_vectors = self.lattice_vectors.copy()
        supercell_vectors[0] *= nx
        supercell_vectors[1] *= ny
        supercell_vectors[2] *= nz

        return Lattice(supercell_vectors, pbc=self.pbc)

    def get_lattice_points(self, size: Tuple[int, int, int]) -> np.ndarray:
        """
        Generate lattice points for supercell.

        Args:
            size: (nx, ny, nz) supercell dimensions

        Returns:
            lattice_points: (nx*ny*nz, 3) array of lattice vectors

        Example:
            >>> # For 2x2x1 supercell of cubic lattice
            >>> points = lattice.get_lattice_points((2, 2, 1))
            >>> # Returns 4 points: (0,0,0), (a,0,0), (0,a,0), (a,a,0)
        """
        nx, ny, nz = size

        points = []
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    R = i * self.lattice_vectors[0] + \
                        j * self.lattice_vectors[1] + \
                        k * self.lattice_vectors[2]
                    points.append(R)

        return np.array(points)

    def get_nearest_neighbors(self,
                             position: np.ndarray,
                             cutoff: float = 5.0,
                             exclude_self: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """
        Find nearest neighbor lattice points within cutoff.

        Args:
            position: Reference position (3,) in Angstrom
            cutoff: Maximum distance in Angstrom
            exclude_self: Exclude (0,0,0) lattice point

        Returns:
            distances: Array of distances
            lattice_points: (N, 3) array of neighbor lattice vectors
        """
        # Generate sufficient lattice points
        max_n = int(np.ceil(cutoff / min(self.a, self.b, self.c))) + 1

        neighbors = []
        distances_list = []

        for i in range(-max_n, max_n + 1):
            for j in range(-max_n, max_n + 1):
                for k in range(-max_n, max_n + 1):
                    if exclude_self and i == 0 and j == 0 and k == 0:
                        continue

                    R = i * self.lattice_vectors[0] + \
                        j * self.lattice_vectors[1] + \
                        k * self.lattice_vectors[2]

                    dist = self.minimum_image_distance(position, position + R)

                    if dist < cutoff:
                        neighbors.append(R)
                        distances_list.append(dist)

        if len(neighbors) == 0:
            return np.array([]), np.array([]).reshape(0, 3)

        # Sort by distance
        distances = np.array(distances_list)
        lattice_points = np.array(neighbors)

        sort_idx = np.argsort(distances)

        return distances[sort_idx], lattice_points[sort_idx]

    def get_brillouin_zone_volume(self) -> float:
        """
        Get volume of first Brillouin zone.

        Returns:
            BZ_volume: Volume in (2π/Angstrom)³
        """
        # BZ volume = (2π)³ / real_space_volume
        return (2 * np.pi)**3 / self.volume

    def __repr__(self) -> str:
        """String representation."""
        return (f"Lattice(a={self.a:.3f}, b={self.b:.3f}, c={self.c:.3f}, "
                f"α={self.alpha:.1f}°, β={self.beta:.1f}°, γ={self.gamma:.1f}°, "
                f"V={self.volume:.3f} Å³, PBC={self.pbc})")

    def __str__(self) -> str:
        """Detailed string representation."""
        s = "Crystal Lattice\n"
        s += "=" * 50 + "\n"
        s += f"Lattice vectors (Angstrom):\n"
        s += f"  a = [{self.lattice_vectors[0, 0]:.4f}, {self.lattice_vectors[0, 1]:.4f}, {self.lattice_vectors[0, 2]:.4f}]\n"
        s += f"  b = [{self.lattice_vectors[1, 0]:.4f}, {self.lattice_vectors[1, 1]:.4f}, {self.lattice_vectors[1, 2]:.4f}]\n"
        s += f"  c = [{self.lattice_vectors[2, 0]:.4f}, {self.lattice_vectors[2, 1]:.4f}, {self.lattice_vectors[2, 2]:.4f}]\n"
        s += f"\nLattice parameters:\n"
        s += f"  a = {self.a:.4f} Å\n"
        s += f"  b = {self.b:.4f} Å\n"
        s += f"  c = {self.c:.4f} Å\n"
        s += f"  α = {self.alpha:.2f}°\n"
        s += f"  β = {self.beta:.2f}°\n"
        s += f"  γ = {self.gamma:.2f}°\n"
        s += f"\nUnit cell volume: {self.volume:.4f} Å³\n"
        s += f"Periodic boundary conditions: {self.pbc}\n"
        return s
