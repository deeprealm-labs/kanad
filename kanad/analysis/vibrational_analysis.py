"""
Vibrational Frequency Calculator

Computes vibrational frequencies and normal modes from the Hessian
(second derivatives of energy with respect to nuclear positions).

Theory:
    1. Compute Hessian: H_ij = ∂²E/∂R_i∂R_j
    2. Mass-weight: H̃_ij = H_ij / √(m_i m_j)
    3. Diagonalize: H̃ v = λ v
    4. Frequencies: ν = √λ / (2πc)
"""

import numpy as np
from typing import Dict, List, Optional, Any, Tuple
import logging

logger = logging.getLogger(__name__)


class FrequencyCalculator:
    """
    Calculate vibrational frequencies and normal modes from Hessian.

    Uses numerical differentiation of gradients to compute Hessian,
    then diagonalizes mass-weighted Hessian to obtain frequencies.

    Example:
        >>> from kanad.io import from_smiles
        >>> from kanad.analysis import FrequencyCalculator
        >>>
        >>> water = from_smiles("O")
        >>> freq_calc = FrequencyCalculator(water)
        >>> result = freq_calc.compute_frequencies(method='HF')
        >>>
        >>> print(f"Frequencies: {result['frequencies']} cm⁻¹")
        >>> print(f"ZPE: {result['zpe']:.6f} Ha")
    """

    # Physical constants (CODATA 2018)
    Ha_to_J = 4.3597447222071e-18       # Hartree to Joules
    Bohr_to_m = 5.29177210903e-11       # Bohr to meters
    Bohr_to_A = 0.529177210903          # Bohr to Angstrom
    amu_to_kg = 1.66053906660e-27       # Atomic mass unit to kg
    c_SI = 2.99792458e8                 # Speed of light (m/s)
    c_cm_s = 2.99792458e10              # Speed of light (cm/s)
    h = 6.62607015e-34                  # Planck constant (J·s)
    N_A = 6.02214076e23                 # Avogadro's number

    def __init__(self, molecule: 'Molecule'):
        """
        Initialize frequency calculator.

        Args:
            molecule: Molecule object (should be at equilibrium geometry)

        Raises:
            ValueError: If molecule has no atoms
        """
        self.molecule = molecule

        if len(molecule.atoms) == 0:
            raise ValueError("Molecule has no atoms")

        self.n_atoms = len(molecule.atoms)
        self.n_coords = 3 * self.n_atoms

        logger.info(f"FrequencyCalculator initialized for {molecule.formula}")
        logger.info(f"  {self.n_atoms} atoms, {self.n_coords} coordinates")

    def _is_linear(self) -> bool:
        """
        Check if molecule is linear.

        Returns True if all atoms are collinear.
        """
        if self.n_atoms <= 2:
            return True

        # Check if all atoms lie on a line
        # Compute inertia tensor and check if smallest moment is negligible
        positions = np.array([atom.position for atom in self.molecule.atoms])
        masses = np.array([atom.atomic_mass for atom in self.molecule.atoms])

        # Center of mass
        com = np.sum(masses[:, np.newaxis] * positions, axis=0) / np.sum(masses)

        # Positions relative to COM
        positions_rel = positions - com

        # Inertia tensor
        I = np.zeros((3, 3))
        for i, (pos, mass) in enumerate(zip(positions_rel, masses)):
            I[0, 0] += mass * (pos[1]**2 + pos[2]**2)
            I[1, 1] += mass * (pos[0]**2 + pos[2]**2)
            I[2, 2] += mass * (pos[0]**2 + pos[1]**2)
            I[0, 1] -= mass * pos[0] * pos[1]
            I[0, 2] -= mass * pos[0] * pos[2]
            I[1, 2] -= mass * pos[1] * pos[2]

        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]

        # Principal moments
        I_principal = np.linalg.eigvalsh(I)

        # Linear if smallest moment << largest
        return I_principal[0] < 1e-3 * I_principal[2]

    def compute_hessian(
        self,
        method: str = 'HF',
        step_size: float = 0.01,
        verbose: bool = True
    ) -> np.ndarray:
        """
        Compute Hessian matrix via finite differences of gradients.

        H_ij = ∂²E/∂R_i∂R_j ≈ (∂G_i/∂R_j)
             ≈ (G_i(R_j + δ) - G_i(R_j - δ)) / (2δ)

        Args:
            method: Electronic structure method ('HF', 'MP2')
            step_size: Displacement for finite difference (Bohr), default 0.01
            verbose: Print progress (default True)

        Returns:
            Hessian matrix (3N × 3N) in atomic units (Ha/Bohr²)

        Note:
            Requires 2×3N gradient evaluations (can be slow for large molecules)
        """
        from kanad.core.gradients import GradientCalculator

        # Store original positions
        orig_positions = [atom.position.copy() for atom in self.molecule.atoms]

        # Initialize Hessian
        hessian = np.zeros((self.n_coords, self.n_coords))

        if verbose:
            print(f"\nComputing Hessian ({self.n_coords}×{self.n_coords})...")
            print(f"  Method: {method}")
            print(f"  Step size: {step_size:.4f} Bohr = {step_size * self.Bohr_to_A:.4f} Å")
            print(f"  Gradient evaluations needed: {2 * self.n_coords}")
            print("-" * 70)

        # For each coordinate
        for i in range(self.n_coords):
            atom_idx = i // 3
            coord_idx = i % 3

            # Forward displacement (+δ)
            self.molecule.atoms[atom_idx].position[coord_idx] += step_size * self.Bohr_to_A
            self.molecule._hamiltonian = None  # Force rebuild
            grad_calc = GradientCalculator(self.molecule, method=method)
            result_plus = grad_calc.compute_gradient()
            grad_plus = result_plus['gradient'].flatten()

            # Backward displacement (-δ)
            self.molecule.atoms[atom_idx].position[coord_idx] -= 2 * step_size * self.Bohr_to_A
            self.molecule._hamiltonian = None
            grad_calc = GradientCalculator(self.molecule, method=method)
            result_minus = grad_calc.compute_gradient()
            grad_minus = result_minus['gradient'].flatten()

            # Restore position
            self.molecule.atoms[atom_idx].position = orig_positions[atom_idx].copy()

            # Hessian column via central difference
            hessian[:, i] = (grad_plus - grad_minus) / (2.0 * step_size)

            if verbose:
                print(f"  Progress: {i+1}/{self.n_coords} ({100*(i+1)/self.n_coords:.1f}%)", end='\r')

        if verbose:
            print()
            print("-" * 70)

        # Force rebuild with original geometry
        self.molecule._hamiltonian = None

        # Symmetrize (Hessian should be symmetric)
        hessian = 0.5 * (hessian + hessian.T)

        # Check symmetry
        asymmetry = np.max(np.abs(hessian - hessian.T))
        if verbose:
            print(f"Hessian computed. Asymmetry: {asymmetry:.2e} (should be ~0)")

        return hessian

    def _mass_weight_hessian(self, hessian: np.ndarray) -> np.ndarray:
        """
        Convert Cartesian Hessian to mass-weighted Hessian.

        H̃_ij = H_ij / √(m_i m_j)

        Args:
            hessian: Cartesian Hessian (Ha/Bohr²)

        Returns:
            Mass-weighted Hessian (Ha·amu/Bohr²)
        """
        # Get masses (repeat each 3 times for x, y, z)
        masses = np.repeat([atom.atomic_mass for atom in self.molecule.atoms], 3)

        # Mass-weighting factor: 1/√(m_i m_j)
        mass_factor = 1.0 / np.sqrt(np.outer(masses, masses))

        # Mass-weighted Hessian
        hessian_mw = hessian * mass_factor

        return hessian_mw

    def _project_translations_rotations(
        self,
        eigenvalues: np.ndarray,
        eigenvectors: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Separate vibrational modes from translations and rotations.

        Translation and rotation have zero (or near-zero) eigenvalues.

        Args:
            eigenvalues: All eigenvalues
            eigenvectors: All eigenvectors

        Returns:
            vib_eigenvalues: Vibrational eigenvalues only
            vib_eigenvectors: Vibrational eigenvectors only
        """
        # Count expected vibrational modes
        if self.n_atoms == 1:
            n_vib = 0  # Atom - no vibrations
        elif self._is_linear():
            n_vib = 3 * self.n_atoms - 5  # Linear molecule
        else:
            n_vib = 3 * self.n_atoms - 6  # Nonlinear molecule

        if n_vib == 0:
            return np.array([]), np.array([]).reshape(self.n_coords, 0)

        # Sort eigenvalues in descending order
        # Vibrational modes have largest (most positive) eigenvalues
        # Translation/rotation have smallest (near zero) eigenvalues
        indices = np.argsort(eigenvalues)[::-1]

        # Take the n_vib largest eigenvalues
        vib_indices = indices[:n_vib]

        vib_eigenvalues = eigenvalues[vib_indices]
        vib_eigenvectors = eigenvectors[:, vib_indices]

        return vib_eigenvalues, vib_eigenvectors

    def _eigenvalues_to_frequencies(self, eigenvalues: np.ndarray) -> np.ndarray:
        """
        Convert eigenvalues to frequencies in cm⁻¹.

        ν (cm⁻¹) = √λ / (2πc)

        where λ is in (Ha/Bohr²)/amu

        Args:
            eigenvalues: Eigenvalues from mass-weighted Hessian

        Returns:
            Frequencies in cm⁻¹ (negative for imaginary frequencies)
        """
        # Convert eigenvalues to SI units (s⁻²)
        # λ in (Ha/Bohr²)/amu → (J/m²)/kg = s⁻²
        eigenvalues_SI = eigenvalues * (self.Ha_to_J / self.Bohr_to_m**2) / self.amu_to_kg

        # Compute frequencies
        frequencies = np.zeros_like(eigenvalues)

        for i, λ in enumerate(eigenvalues_SI):
            if λ > 0:
                # Real frequency
                frequencies[i] = np.sqrt(λ) / (2 * np.pi * self.c_cm_s)
            else:
                # Imaginary frequency (negative indicates saddle point or numerical error)
                frequencies[i] = -np.sqrt(-λ) / (2 * np.pi * self.c_cm_s)

        return frequencies

    def _compute_zpe(self, frequencies: np.ndarray) -> float:
        """
        Compute zero-point energy from frequencies.

        ZPE = Σ_i (h·ν_i / 2)

        Args:
            frequencies: Vibrational frequencies (cm⁻¹)

        Returns:
            Zero-point energy (Ha)
        """
        zpe = 0.0

        for ν_cm in frequencies:
            if ν_cm > 0:  # Only positive (real) frequencies contribute
                ν_Hz = ν_cm * self.c_cm_s  # Convert cm⁻¹ to Hz
                zpe += 0.5 * self.h * ν_Hz  # J

        # Convert to Hartree
        zpe_Ha = zpe / self.Ha_to_J

        return zpe_Ha

    def compute_frequencies(
        self,
        method: str = 'HF',
        hessian: Optional[np.ndarray] = None,
        step_size: float = 0.01,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Compute vibrational frequencies and normal modes.

        Args:
            method: Electronic structure method ('HF', 'MP2')
            hessian: Pre-computed Hessian (if None, will compute)
            step_size: Step size for Hessian calculation (Bohr)
            verbose: Print progress and results

        Returns:
            Dictionary with:
                frequencies: Vibrational frequencies (cm⁻¹), sorted
                normal_modes: Normal mode eigenvectors (3N × n_vib)
                reduced_masses: Reduced masses for each mode (amu)
                force_constants: Force constants for each mode (mdyn/Å)
                zpe: Zero-point energy (Ha)
                hessian: Hessian matrix (Ha/Bohr²)
                n_imaginary: Number of imaginary frequencies
        """
        # Compute Hessian if not provided
        if hessian is None:
            hessian = self.compute_hessian(method=method, step_size=step_size, verbose=verbose)
        else:
            if verbose:
                print("\nUsing provided Hessian matrix")

        # Mass-weight Hessian
        if verbose:
            print("\nMass-weighting Hessian...")

        hessian_mw = self._mass_weight_hessian(hessian)

        # Diagonalize
        if verbose:
            print("Diagonalizing mass-weighted Hessian...")

        eigenvalues, eigenvectors = np.linalg.eigh(hessian_mw)

        # Project out translations and rotations
        if verbose:
            print("Projecting out translations and rotations...")

        vib_eigenvalues, vib_eigenvectors = self._project_translations_rotations(
            eigenvalues, eigenvectors
        )

        if len(vib_eigenvalues) == 0:
            if verbose:
                print("\nNo vibrational modes (atom)")
            return {
                'frequencies': np.array([]),
                'normal_modes': np.array([]).reshape(self.n_coords, 0),
                'reduced_masses': np.array([]),
                'force_constants': np.array([]),
                'zpe': 0.0,
                'hessian': hessian,
                'n_imaginary': 0
            }

        # Convert eigenvalues to frequencies
        frequencies = self._eigenvalues_to_frequencies(vib_eigenvalues)

        # Sort by frequency (ascending)
        sort_indices = np.argsort(np.abs(frequencies))
        frequencies = frequencies[sort_indices]
        vib_eigenvectors = vib_eigenvectors[:, sort_indices]
        vib_eigenvalues = vib_eigenvalues[sort_indices]

        # Compute reduced masses and force constants
        # μ_i = 1 (normalized eigenvector of mass-weighted Hessian)
        # k_i = λ_i (eigenvalue in atomic units, convert to mdyn/Å)
        reduced_masses = np.ones(len(frequencies))  # amu

        # Force constant: k = 4π²ν²μ (convert to mdyn/Å)
        # Or from eigenvalue: k = λ
        force_constants = vib_eigenvalues * (self.Ha_to_J / self.Bohr_to_m**2) # J/m² = N/m
        force_constants *= 1e-15  # N/m to mdyn/Å (1 mdyn = 10^-8 N, 1 Å = 10^-10 m)

        # Compute ZPE
        zpe = self._compute_zpe(frequencies)

        # Count imaginary frequencies
        n_imaginary = np.sum(frequencies < 0)

        if verbose:
            print("\n" + "=" * 70)
            print("VIBRATIONAL ANALYSIS RESULTS")
            print("=" * 70)
            print(f"Molecule: {self.molecule.formula}")
            print(f"Method: {method}")
            print(f"Number of atoms: {self.n_atoms}")
            print(f"Linear: {self._is_linear()}")
            print(f"Vibrational modes: {len(frequencies)}")
            print(f"Imaginary frequencies: {n_imaginary}")
            print(f"\nZero-point energy: {zpe:.8f} Ha = {zpe * 627.509:.2f} kcal/mol")
            print("\nVibrational frequencies (cm⁻¹):")
            print("-" * 70)
            for i, (ν, k) in enumerate(zip(frequencies, force_constants)):
                if ν < 0:
                    print(f"  Mode {i+1:2d}:  {ν:10.2f}i cm⁻¹  (k = {abs(k):8.2f} mdyn/Å) [IMAGINARY]")
                else:
                    print(f"  Mode {i+1:2d}:  {ν:10.2f}  cm⁻¹  (k = {k:8.2f} mdyn/Å)")
            print("=" * 70)

        return {
            'frequencies': frequencies,
            'normal_modes': vib_eigenvectors,
            'reduced_masses': reduced_masses,
            'force_constants': force_constants,
            'zpe': zpe,
            'hessian': hessian,
            'n_imaginary': n_imaginary,
            'method': method
        }
