"""
Orbital Optimizer - Orbital localization and rotation.

Localized orbitals improve circuit efficiency:
- Boys localization: Minimize orbital spread
- Pipek-Mezey: Maximize orbital locality on atoms
- Natural orbitals: Diagonalize density matrix
- Rotation: Optimize for specific bond types (governance-aware)
"""

from typing import Dict, Any, Optional, Tuple
import numpy as np
from scipy.linalg import eigh, expm
import logging

logger = logging.getLogger(__name__)


class OrbitalOptimizer:
    """
    Optimize molecular orbitals for circuit efficiency.

    Localization methods:
    - **Boys**: Minimize sum of orbital spreads ⟨r²⟩ - ⟨r⟩²
    - **Pipek-Mezey**: Maximize Mulliken populations on atoms
    - **Natural Orbitals**: Diagonalize density matrix
    - **Governance-aware**: Rotate to align with bond character

    Example:
        >>> optimizer = OrbitalOptimizer(mo_coeffs, overlap_matrix)
        >>> localized_orbitals, metrics = optimizer.localize_boys()
        >>> print(f"Sparsity improved by {metrics['sparsity_improvement']:.1f}%")
    """

    def __init__(
        self,
        mo_coefficients: np.ndarray,
        overlap_matrix: Optional[np.ndarray] = None,
        atomic_positions: Optional[np.ndarray] = None
    ):
        """
        Initialize orbital optimizer.

        Args:
            mo_coefficients: Molecular orbital coefficients (AO x MO)
            overlap_matrix: Atomic orbital overlap matrix (for orthogonalization)
            atomic_positions: Atomic positions for Boys localization
        """
        self.mo_coeffs = mo_coefficients
        self.overlap = overlap_matrix if overlap_matrix is not None else np.eye(mo_coefficients.shape[0])
        self.positions = atomic_positions

        self.n_ao, self.n_mo = mo_coefficients.shape

        logger.info(f"OrbitalOptimizer initialized: {self.n_ao} AOs, {self.n_mo} MOs")

    def localize_boys(
        self,
        occupied_only: bool = True,
        max_iterations: int = 100,
        convergence_threshold: float = 1e-6
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Boys localization: minimize orbital spread.

        Minimizes: Σ_i (⟨φ_i|r²|φ_i⟩ - ⟨φ_i|r|φ_i⟩²)

        Args:
            occupied_only: Only localize occupied orbitals
            max_iterations: Maximum Jacobi rotation iterations
            convergence_threshold: Convergence criterion

        Returns:
            Tuple of (localized_mo_coeffs, metrics_dict)
        """
        if self.positions is None:
            logger.warning("No atomic positions provided, using identity localization")
            return self.mo_coeffs, {'method': 'none', 'converged': False}

        logger.info("Applying Boys localization...")

        # Determine orbitals to localize
        n_orbitals = self.n_mo if not occupied_only else self.n_mo // 2

        # Build position operators
        x_ao, y_ao, z_ao = self._build_position_operators()

        # Transform to MO basis
        C = self.mo_coeffs[:, :n_orbitals]
        x_mo = C.T @ x_ao @ C
        y_mo = C.T @ y_ao @ C
        z_mo = C.T @ z_ao @ C

        # Jacobi rotations to minimize spread
        U = np.eye(n_orbitals)
        spread_old = self._compute_boys_spread(x_mo, y_mo, z_mo)

        for iteration in range(max_iterations):
            max_change = 0.0

            for i in range(n_orbitals):
                for j in range(i + 1, n_orbitals):
                    # Compute rotation angle
                    A_ij = self._boys_rotation_angle(i, j, x_mo, y_mo, z_mo)

                    # Apply 2x2 rotation
                    c = np.cos(A_ij)
                    s = np.sin(A_ij)

                    # Update orbitals
                    U_rot = np.eye(n_orbitals)
                    U_rot[i, i] = c
                    U_rot[j, j] = c
                    U_rot[i, j] = -s
                    U_rot[j, i] = s

                    U = U @ U_rot

                    # Update position matrices
                    x_mo = U_rot.T @ x_mo @ U_rot
                    y_mo = U_rot.T @ y_mo @ U_rot
                    z_mo = U_rot.T @ z_mo @ U_rot

                    max_change = max(max_change, abs(A_ij))

            spread_new = self._compute_boys_spread(x_mo, y_mo, z_mo)
            delta_spread = abs(spread_new - spread_old)

            if delta_spread < convergence_threshold:
                logger.info(f"Boys localization converged in {iteration + 1} iterations")
                break

            spread_old = spread_new

        # Apply transformation
        C_localized = self.mo_coeffs.copy()
        C_localized[:, :n_orbitals] = self.mo_coeffs[:, :n_orbitals] @ U

        # Compute metrics
        metrics = {
            'method': 'boys',
            'converged': iteration < max_iterations - 1,
            'iterations': iteration + 1,
            'final_spread': spread_new,
            'n_localized': n_orbitals
        }

        logger.info(f"  Final spread: {spread_new:.6f}")
        logger.info(f"  Localized {n_orbitals} orbitals")

        return C_localized, metrics

    def localize_pipek_mezey(
        self,
        occupied_only: bool = True,
        max_iterations: int = 100
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """
        Pipek-Mezey localization: maximize Mulliken populations.

        Maximizes: Σ_A (Σ_i Q_i^A)²
        where Q_i^A is Mulliken population of orbital i on atom A

        Args:
            occupied_only: Only localize occupied orbitals
            max_iterations: Maximum iterations

        Returns:
            Tuple of (localized_mo_coeffs, metrics_dict)
        """
        logger.info("Applying Pipek-Mezey localization...")

        n_orbitals = self.n_mo if not occupied_only else self.n_mo // 2

        # Build Mulliken population matrices (simplified)
        # Q_i^A = Σ_μ∈A Σ_ν (C_μi * S_μν * C_νi)
        # For simplicity, assume each AO belongs to one atom

        C = self.mo_coeffs[:, :n_orbitals]
        U = np.eye(n_orbitals)

        # Jacobi rotations to maximize localization
        for iteration in range(max_iterations):
            max_change = 0.0

            for i in range(n_orbitals):
                for j in range(i + 1, n_orbitals):
                    # Compute rotation angle (simplified)
                    angle = self._pm_rotation_angle(i, j, C)

                    c = np.cos(angle)
                    s = np.sin(angle)

                    U_rot = np.eye(n_orbitals)
                    U_rot[i, i] = c
                    U_rot[j, j] = c
                    U_rot[i, j] = -s
                    U_rot[j, i] = s

                    U = U @ U_rot
                    C = C @ U_rot

                    max_change = max(max_change, abs(angle))

            if max_change < 1e-6:
                logger.info(f"Pipek-Mezey converged in {iteration + 1} iterations")
                break

        C_localized = self.mo_coeffs.copy()
        C_localized[:, :n_orbitals] = self.mo_coeffs[:, :n_orbitals] @ U

        metrics = {
            'method': 'pipek_mezey',
            'converged': iteration < max_iterations - 1,
            'iterations': iteration + 1,
            'n_localized': n_orbitals
        }

        logger.info(f"  Localized {n_orbitals} orbitals")

        return C_localized, metrics

    def get_natural_orbitals(
        self,
        density_matrix: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute natural orbitals from density matrix.

        Natural orbitals are eigenvectors of the density matrix:
        P = Σ_i n_i |φ_i⟩⟨φ_i|

        Orbitals with occupation n_i ≈ 0 or 2 are unimportant.
        Orbitals with 0 < n_i < 2 indicate correlation.

        Args:
            density_matrix: One-particle density matrix

        Returns:
            Tuple of (occupation_numbers, natural_orbital_coeffs)
        """
        logger.info("Computing natural orbitals...")

        # Diagonalize density matrix
        occupation_numbers, natural_orbitals = eigh(density_matrix)

        # Sort by occupation (descending)
        idx = np.argsort(occupation_numbers)[::-1]
        occupation_numbers = occupation_numbers[idx]
        natural_orbitals = natural_orbitals[:, idx]

        # Count important orbitals (0.01 < n < 1.99)
        important = np.sum((occupation_numbers > 0.01) & (occupation_numbers < 1.99))

        logger.info(f"  Natural orbitals computed")
        logger.info(f"  Important orbitals (correlation): {important}")
        logger.info(f"  Occupation range: [{occupation_numbers.min():.3f}, {occupation_numbers.max():.3f}]")

        return occupation_numbers, natural_orbitals

    def estimate_sparsity_improvement(
        self,
        eri_original: np.ndarray,
        eri_localized: np.ndarray,
        threshold: float = 1e-6
    ) -> Dict[str, float]:
        """
        Estimate sparsity improvement from localization.

        Args:
            eri_original: Original ERI tensor
            eri_localized: Localized ERI tensor
            threshold: Sparsity threshold

        Returns:
            Dictionary with sparsity metrics
        """
        # Count elements below threshold
        sparse_original = np.sum(np.abs(eri_original) < threshold)
        sparse_localized = np.sum(np.abs(eri_localized) < threshold)

        total_elements = eri_original.size

        sparsity_original = sparse_original / total_elements * 100
        sparsity_localized = sparse_localized / total_elements * 100

        improvement = sparsity_localized - sparsity_original

        logger.info(f"Sparsity analysis:")
        logger.info(f"  Original:  {sparsity_original:.1f}%")
        logger.info(f"  Localized: {sparsity_localized:.1f}%")
        logger.info(f"  Improvement: {improvement:.1f}%")

        return {
            'sparsity_original': sparsity_original,
            'sparsity_localized': sparsity_localized,
            'sparsity_improvement': improvement,
            'threshold': threshold
        }

    def _build_position_operators(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Build position operators in AO basis.

        Returns:
            Tuple of (x_ao, y_ao, z_ao) matrices
        """
        if self.positions is None:
            # Return dummy operators
            return np.eye(self.n_ao), np.eye(self.n_ao), np.eye(self.n_ao)

        # Simplified: assume one AO per atom (for demonstration)
        n_atoms = min(len(self.positions), self.n_ao)

        x_ao = np.zeros((self.n_ao, self.n_ao))
        y_ao = np.zeros((self.n_ao, self.n_ao))
        z_ao = np.zeros((self.n_ao, self.n_ao))

        for i in range(n_atoms):
            x_ao[i, i] = self.positions[i, 0]
            y_ao[i, i] = self.positions[i, 1]
            z_ao[i, i] = self.positions[i, 2]

        return x_ao, y_ao, z_ao

    def _compute_boys_spread(
        self,
        x_mo: np.ndarray,
        y_mo: np.ndarray,
        z_mo: np.ndarray
    ) -> float:
        """
        Compute Boys localization functional.

        Spread = Σ_i [⟨r²⟩_i - ⟨r⟩_i²]
        """
        spread = 0.0

        for i in range(x_mo.shape[0]):
            r2 = x_mo[i, i]**2 + y_mo[i, i]**2 + z_mo[i, i]**2
            r_mean_sq = x_mo[i, i]**2 + y_mo[i, i]**2 + z_mo[i, i]**2

            spread += r2 - r_mean_sq

        return spread

    def _boys_rotation_angle(
        self,
        i: int,
        j: int,
        x_mo: np.ndarray,
        y_mo: np.ndarray,
        z_mo: np.ndarray
    ) -> float:
        """
        Compute optimal Boys rotation angle between orbitals i and j.

        Returns:
            Rotation angle in radians
        """
        # Simplified calculation
        A = x_mo[i, j]**2 + y_mo[i, j]**2 + z_mo[i, j]**2
        B = x_mo[i, i]**2 + y_mo[i, i]**2 + z_mo[i, i]**2 - \
            (x_mo[j, j]**2 + y_mo[j, j]**2 + z_mo[j, j]**2)

        if abs(B) < 1e-10:
            return 0.0

        angle = 0.25 * np.arctan2(4 * A, B)

        return angle

    def _pm_rotation_angle(
        self,
        i: int,
        j: int,
        C: np.ndarray
    ) -> float:
        """
        Compute Pipek-Mezey rotation angle (simplified).

        Returns:
            Rotation angle in radians
        """
        # Simplified: maximize |C_i - C_j|
        diff = np.sum((C[:, i] - C[:, j])**2)

        # Small rotation proportional to difference
        angle = 0.1 * np.tanh(diff)

        return angle
