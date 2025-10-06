"""
Molecular Property Calculator

Computes molecular properties from quantum chemistry calculations:
- Dipole moment
- Polarizability (future)
- Quadrupole moment (future)
- Molecular orbitals analysis
"""

import numpy as np
from typing import Dict, Optional, List, Any
import logging

logger = logging.getLogger(__name__)


class PropertyCalculator:
    """
    Calculate molecular properties from Hamiltonian and density matrix.

    Uses PySCF backend for integral calculations.
    """

    # Unit conversion constants
    AU_TO_DEBYE = 2.541746  # 1 a.u. = 2.541746 Debye
    DEBYE_TO_AU = 1.0 / AU_TO_DEBYE

    def __init__(self, hamiltonian: 'MolecularHamiltonian'):
        """
        Initialize property calculator.

        Args:
            hamiltonian: MolecularHamiltonian object with PySCF backend
        """
        self.hamiltonian = hamiltonian
        self.mol = hamiltonian.mol  # PySCF molecule
        self.atoms = hamiltonian.atoms

    def compute_dipole_moment(
        self,
        density_matrix: Optional[np.ndarray] = None,
        origin: Optional[np.ndarray] = None
    ) -> Dict[str, Any]:
        """
        Compute electric dipole moment.

        Theory:
            μ = -∫ ρ(r)·r dr + Σ_A Z_A R_A
              = -Σ_μν P_μν ⟨φ_μ|r|φ_ν⟩ + Σ_A Z_A R_A

        where:
            P_μν = density matrix
            ⟨φ_μ|r|φ_ν⟩ = dipole integrals
            Z_A = nuclear charge
            R_A = nuclear position

        Args:
            density_matrix: Density matrix (uses HF if None)
            origin: Origin for dipole calculation (default: [0,0,0])

        Returns:
            dict:
                dipole_vector: np.ndarray (3,) in Debye [x, y, z]
                dipole_magnitude: float in Debye
                dipole_au: np.ndarray (3,) in atomic units
                components: dict with x, y, z in Debye
                origin: np.ndarray (3,) origin used

        Examples:
            >>> calc = PropertyCalculator(water.hamiltonian)
            >>> result = calc.compute_dipole_moment()
            >>> print(f"Dipole: {result['dipole_magnitude']:.3f} D")
            Dipole: 2.123 D
        """
        # Get density matrix
        if density_matrix is None:
            density_matrix = self.hamiltonian.mf.make_rdm1()

        # Set origin (default: center of mass)
        if origin is None:
            origin = np.zeros(3)

        # Get dipole integrals from PySCF
        # int1e_r returns (3, n_ao, n_ao) for x, y, z components
        dip_ints = self.mol.intor('int1e_r')  # Shape: (3, n_ao, n_ao)

        # Shift origin if needed
        if not np.allclose(origin, 0.0):
            # Add -origin × overlap to each component
            overlap = self.mol.intor('int1e_ovlp')
            for i in range(3):
                dip_ints[i] -= origin[i] * overlap

        # Electronic contribution: -Tr(P · r)
        mu_elec = np.zeros(3)
        for i in range(3):
            mu_elec[i] = -np.einsum('ij,ji->', density_matrix, dip_ints[i])

        # Nuclear contribution: Σ Z_A (R_A - origin)
        mu_nuc = np.zeros(3)
        for atom in self.atoms:
            # Convert position from Angstroms to Bohr
            from kanad.core.constants.conversion_factors import ConversionFactors
            pos_bohr = atom.position * ConversionFactors.ANGSTROM_TO_BOHR
            mu_nuc += atom.atomic_number * (pos_bohr - origin)

        # Total dipole (atomic units)
        mu_au = mu_elec + mu_nuc

        # Convert to Debye
        mu_debye = mu_au * self.AU_TO_DEBYE
        magnitude = np.linalg.norm(mu_debye)

        logger.info(f"Dipole moment: {magnitude:.4f} D")
        logger.debug(f"  Electronic: {mu_elec * self.AU_TO_DEBYE}")
        logger.debug(f"  Nuclear: {mu_nuc * self.AU_TO_DEBYE}")

        return {
            'dipole_vector': mu_debye,
            'dipole_magnitude': magnitude,
            'dipole_au': mu_au,
            'components': {
                'x': mu_debye[0],
                'y': mu_debye[1],
                'z': mu_debye[2]
            },
            'origin': origin,
            'electronic_contribution': mu_elec * self.AU_TO_DEBYE,
            'nuclear_contribution': mu_nuc * self.AU_TO_DEBYE
        }

    def compute_center_of_mass(self) -> np.ndarray:
        """
        Compute molecular center of mass.

        Returns:
            np.ndarray: Center of mass position (Angstroms)
        """
        total_mass = 0.0
        com = np.zeros(3)

        for atom in self.atoms:
            mass = atom.atomic_mass
            com += mass * atom.position
            total_mass += mass

        return com / total_mass

    def compute_center_of_charge(self) -> np.ndarray:
        """
        Compute center of nuclear charge.

        Returns:
            np.ndarray: Center of charge (Angstroms)
        """
        total_charge = 0.0
        center = np.zeros(3)

        for atom in self.atoms:
            charge = atom.atomic_number
            center += charge * atom.position
            total_charge += charge

        return center / total_charge

    def verify_dipole_with_pyscf(
        self,
        density_matrix: Optional[np.ndarray] = None
    ) -> Dict[str, Any]:
        """
        Verify dipole calculation against PySCF's built-in method.

        Useful for debugging and validation.

        Args:
            density_matrix: Density matrix (uses HF if None)

        Returns:
            dict:
                kanad_dipole: Kanad calculation
                pyscf_dipole: PySCF calculation
                difference: Absolute difference
                agree: bool (True if difference < 0.01 D)
        """
        # Kanad calculation
        kanad_result = self.compute_dipole_moment(density_matrix)

        # PySCF calculation
        if density_matrix is None:
            density_matrix = self.hamiltonian.mf.make_rdm1()

        # PySCF dipole (returns in a.u.)
        pyscf_dipole_au = self.hamiltonian.mf.dip_moment(
            mol=self.mol,
            dm=density_matrix,
            unit='AU'
        )
        pyscf_dipole_debye = pyscf_dipole_au * self.AU_TO_DEBYE
        pyscf_magnitude = np.linalg.norm(pyscf_dipole_debye)

        # Compare
        difference = abs(kanad_result['dipole_magnitude'] - pyscf_magnitude)
        agree = difference < 0.01  # 0.01 D tolerance

        return {
            'kanad_dipole': kanad_result['dipole_magnitude'],
            'pyscf_dipole': pyscf_magnitude,
            'kanad_vector': kanad_result['dipole_vector'],
            'pyscf_vector': pyscf_dipole_debye,
            'difference': difference,
            'agree': agree
        }

    def compute_polarizability(
        self,
        method: str = 'finite_field',
        field_strength: float = 0.001
    ) -> Dict[str, Any]:
        """
        Compute static polarizability tensor.

        Polarizability (α) describes how easily a molecule's electron cloud
        is distorted by an external electric field:
            μ_induced = α · E

        Theory:
            α_ij = -∂μ_i/∂E_j |_{E=0}

        Args:
            method: 'finite_field' (numerical) or 'analytical' (future)
            field_strength: Electric field strength in a.u. (default: 0.001)
                - Too large: non-linear effects
                - Too small: numerical noise
                - Optimal: 0.001 a.u. ≈ 5×10⁸ V/m

        Returns:
            dict:
                alpha_tensor: 3×3 polarizability tensor (a.u.)
                alpha_mean: Mean polarizability (a.u. and Å³)
                alpha_anisotropy: Polarizability anisotropy (a.u.)
                eigenvalues: Principal polarizabilities (a.u.)
                alpha_xx, alpha_yy, alpha_zz: Diagonal elements (a.u.)

        Examples:
            >>> calc = PropertyCalculator(water.hamiltonian)
            >>> result = calc.compute_polarizability()
            >>> print(f"Polarizability: {result['alpha_mean']:.2f} a.u.")
            Polarizability: 9.85 a.u.

        Warning:
            Polarizability calculations require adequate basis sets:
            - STO-3G: Severely underestimates (~20-35% of experimental)
            - Minimum recommended: 6-311G(d,p)
            - For quantitative accuracy: aug-cc-pVDZ or larger + correlation (MP2)

            See POLARIZABILITY_ACCURACY_ANALYSIS.md for details.
        """
        # Check basis set and warn if minimal
        basis_name = self.hamiltonian.mol.basis.lower() if hasattr(self.hamiltonian.mol, 'basis') else 'unknown'
        minimal_basis = ['sto-3g', 'sto-6g', '3-21g', 'sto3g', 'sto6g']

        if any(basis in basis_name for basis in minimal_basis):
            logger.warning(
                f"Computing polarizability with minimal basis set '{self.hamiltonian.mol.basis}'. "
                f"Results will severely underestimate experimental values (typically 20-40% accuracy). "
                f"For quantitative results, use basis='6-311g(d,p)' or larger. "
                f"See POLARIZABILITY_ACCURACY_ANALYSIS.md for details."
            )

        if method == 'finite_field':
            alpha_tensor = self._compute_polarizability_finite_field(field_strength)
        elif method == 'analytical':
            raise NotImplementedError(
                "Analytical polarizability not yet implemented. "
                "Use method='finite_field'"
            )
        else:
            raise ValueError(f"Unknown method: {method}. Use 'finite_field'")

        # Mean polarizability: ᾱ = Tr(α)/3
        alpha_mean = np.trace(alpha_tensor) / 3.0

        # Polarizability anisotropy
        # Δα = √(3/2 ||α - ᾱI||_F)
        alpha_iso = alpha_mean * np.eye(3)
        alpha_aniso_tensor = alpha_tensor - alpha_iso
        anisotropy = np.sqrt(1.5 * np.sum(alpha_aniso_tensor**2))

        # Principal polarizabilities (eigenvalues)
        eigenvalues = np.linalg.eigvalsh(alpha_tensor)

        # Unit conversion: 1 a.u. = 0.1482 Å³
        AU_TO_ANGSTROM3 = 0.1482
        alpha_mean_angstrom = alpha_mean * AU_TO_ANGSTROM3

        logger.info(f"Polarizability: {alpha_mean:.3f} a.u. = {alpha_mean_angstrom:.3f} Å³")
        logger.debug(f"  Anisotropy: {anisotropy:.3f} a.u.")
        logger.debug(f"  Eigenvalues: {eigenvalues}")

        return {
            'alpha_tensor': alpha_tensor,  # 3×3 matrix (a.u.)
            'alpha_mean': alpha_mean,  # scalar (a.u.)
            'alpha_mean_angstrom3': alpha_mean_angstrom,  # scalar (Å³)
            'alpha_anisotropy': anisotropy,  # scalar (a.u.)
            'eigenvalues': eigenvalues,  # (3,) array (a.u.)
            'alpha_xx': alpha_tensor[0, 0],
            'alpha_yy': alpha_tensor[1, 1],
            'alpha_zz': alpha_tensor[2, 2],
            'method': method,
            'field_strength': field_strength
        }

    def _compute_polarizability_finite_field(
        self,
        field_strength: float
    ) -> np.ndarray:
        """
        Compute polarizability via finite field method.

        Applies small electric fields in ±x, ±y, ±z directions,
        computes induced dipole moments, and uses finite differences:
            α_ij ≈ -[μ_i(+E_j) - μ_i(-E_j)] / (2E_j)

        Args:
            field_strength: Electric field magnitude (a.u.)

        Returns:
            np.ndarray: 3×3 polarizability tensor (a.u.)
        """
        alpha = np.zeros((3, 3))

        logger.debug(f"Computing polarizability with field strength {field_strength} a.u.")

        # Apply field in each direction (x, y, z)
        for direction in range(3):
            # Positive field
            field_vec_plus = np.zeros(3)
            field_vec_plus[direction] = field_strength
            dipole_plus = self._compute_dipole_with_field(field_vec_plus)

            # Negative field
            field_vec_minus = np.zeros(3)
            field_vec_minus[direction] = -field_strength
            dipole_minus = self._compute_dipole_with_field(field_vec_minus)

            # Finite difference: α_ij = -dμ_i/dE_j
            for component in range(3):
                alpha[component, direction] = -(
                    dipole_plus[component] - dipole_minus[component]
                ) / (2.0 * field_strength)

        # Symmetrize tensor (α should be symmetric)
        alpha_sym = 0.5 * (alpha + alpha.T)

        # Check symmetry
        asymmetry = np.max(np.abs(alpha - alpha.T))
        if asymmetry > 0.1:
            logger.warning(f"Polarizability tensor asymmetry: {asymmetry:.4f} a.u.")

        return alpha_sym

    def _compute_dipole_with_field(self, field_vector: np.ndarray) -> np.ndarray:
        """
        Compute dipole moment with external electric field applied.

        Modifies the core Hamiltonian to include field interaction:
            H' = H₀ - μ·E = H₀ - r·E

        Then runs SCF to self-consistency and computes dipole.

        Args:
            field_vector: Electric field [Ex, Ey, Ez] in atomic units

        Returns:
            np.ndarray: Dipole moment vector in atomic units
        """
        from pyscf import scf

        # Build modified core Hamiltonian
        # H' = T + V_ne - r·E
        h1e = self.mol.intor('int1e_kin') + self.mol.intor('int1e_nuc')

        # Add field term: -r·E
        dip_ints = self.mol.intor('int1e_r')  # (3, n_ao, n_ao)
        for i in range(3):
            h1e -= field_vector[i] * dip_ints[i]

        # Create new SCF object with modified Hamiltonian
        if self.hamiltonian.spin == 0:
            mf_field = scf.RHF(self.mol)
        else:
            mf_field = scf.ROHF(self.mol)

        # Override get_hcore to use our modified H
        mf_field.get_hcore = lambda *args: h1e

        # Run SCF (suppress output)
        mf_field.verbose = 0
        mf_field.kernel()

        # Check convergence
        if not mf_field.converged:
            logger.warning(f"SCF with field {field_vector} did not converge")

        # Get density matrix
        dm_field = mf_field.make_rdm1()

        # Compute dipole with this density
        # Important: Use original (unperturbed) dipole calculation
        # The field is already accounted for in the density
        result = self.compute_dipole_moment(density_matrix=dm_field)

        return result['dipole_au']

    # =========================================================================
    # MP2 Polarizability (with electron correlation)
    # =========================================================================

    def compute_polarizability_mp2(
        self,
        field_strength: float = 0.001
    ) -> Dict[str, Any]:
        """
        Compute polarizability with MP2 electron correlation.

        Significantly more accurate than HF (~20-30% improvement over HF).
        Uses finite field method with MP2-correlated density matrices.

        Expected accuracy:
            - MP2/6-311G(d,p):  ~70-80% of experimental
            - MP2/aug-cc-pVDZ:  ~85-95% of experimental
            - MP2/aug-cc-pVTZ:  ~90-98% of experimental

        Computational cost: ~6-10x more expensive than HF polarizability
        (6 field directions × MP2 iterations)

        Args:
            field_strength: Electric field strength in a.u. (default: 0.001)

        Returns:
            dict: Same format as compute_polarizability(), plus:
                method: 'mp2_finite_field'

        Example:
            >>> calc = PropertyCalculator(water.hamiltonian)
            >>> result = calc.compute_polarizability_mp2()
            >>> print(f"MP2 polarizability: {result['alpha_mean']:.2f} a.u.")
            MP2 polarizability: 9.15 a.u.

        Note:
            Requires converged HF reference. Will run MP2 calculation
            at each field point (6 total: ±x, ±y, ±z).
        """
        from pyscf import mp

        logger.info("Computing MP2 polarizability (this will take longer than HF)...")

        # Check that HF has converged
        if not self.hamiltonian.mf.converged:
            raise ValueError("HF must converge before MP2 polarizability")

        alpha_tensor = np.zeros((3, 3))

        # Apply field in each direction (x, y, z)
        for direction in range(3):
            # Positive field
            field_vec_plus = np.zeros(3)
            field_vec_plus[direction] = field_strength
            dipole_plus = self._compute_dipole_with_field_mp2(field_vec_plus)

            # Negative field
            field_vec_minus = np.zeros(3)
            field_vec_minus[direction] = -field_strength
            dipole_minus = self._compute_dipole_with_field_mp2(field_vec_minus)

            # Finite difference: α_ij = -dμ_i/dE_j
            for component in range(3):
                alpha_tensor[component, direction] = -(
                    dipole_plus[component] - dipole_minus[component]
                ) / (2.0 * field_strength)

            logger.debug(f"Direction {direction} ({'xyz'[direction]}): complete")

        # Symmetrize tensor (α should be symmetric)
        alpha_sym = 0.5 * (alpha_tensor + alpha_tensor.T)

        # Check symmetry
        asymmetry = np.max(np.abs(alpha_tensor - alpha_tensor.T))
        if asymmetry > 0.1:
            logger.warning(f"MP2 polarizability tensor asymmetry: {asymmetry:.4f} a.u.")

        # Mean polarizability: ᾱ = Tr(α)/3
        alpha_mean = np.trace(alpha_sym) / 3.0

        # Convert to Angstrom^3 (1 a.u. = 0.1482 Å³)
        AU_TO_ANGSTROM3 = 0.1482
        alpha_mean_angstrom = alpha_mean * AU_TO_ANGSTROM3

        # Polarizability anisotropy: Δα = √(3/2 ||α - ᾱI||_F)
        alpha_iso = alpha_mean * np.eye(3)
        alpha_aniso_tensor = alpha_sym - alpha_iso
        alpha_anisotropy = np.sqrt(1.5 * np.sum(alpha_aniso_tensor**2))

        # Principal polarizabilities (eigenvalues)
        eigenvalues = np.linalg.eigvalsh(alpha_sym)

        logger.info(f"MP2 mean polarizability: {alpha_mean:.4f} a.u. = {alpha_mean_angstrom:.4f} Å³")

        return {
            'alpha_tensor': alpha_sym,
            'alpha_mean': alpha_mean,
            'alpha_mean_angstrom3': alpha_mean_angstrom,
            'alpha_anisotropy': alpha_anisotropy,
            'eigenvalues': eigenvalues,
            'alpha_xx': alpha_sym[0, 0],
            'alpha_yy': alpha_sym[1, 1],
            'alpha_zz': alpha_sym[2, 2],
            'method': 'mp2_finite_field',
            'field_strength': field_strength
        }

    def _compute_dipole_with_field_mp2(self, field_vector: np.ndarray) -> np.ndarray:
        """
        Compute MP2 dipole moment with external electric field applied.

        Runs HF with field, then MP2 on top, and computes dipole from MP2 density.

        Args:
            field_vector: Electric field [Ex, Ey, Ez] in atomic units

        Returns:
            np.ndarray: MP2 dipole moment vector in atomic units
        """
        from pyscf import scf, mp

        # Build modified core Hamiltonian
        # H' = T + V_ne - r·E
        h1e = self.mol.intor('int1e_kin') + self.mol.intor('int1e_nuc')

        # Add field term: -r·E
        dip_ints = self.mol.intor('int1e_r')  # (3, n_ao, n_ao)
        for i in range(3):
            h1e -= field_vector[i] * dip_ints[i]

        # Create new SCF object with modified Hamiltonian
        if self.hamiltonian.spin == 0:
            mf_field = scf.RHF(self.mol)
        else:
            mf_field = scf.ROHF(self.mol)

        # Override get_hcore to use our modified H
        mf_field.get_hcore = lambda *args: h1e

        # Run SCF (suppress output)
        mf_field.verbose = 0
        mf_field.kernel()

        # Check convergence
        if not mf_field.converged:
            logger.warning(f"SCF with field {field_vector} did not converge for MP2")

        # Run MP2 on top of field-perturbed HF
        mp2_solver = mp.MP2(mf_field)
        mp2_solver.verbose = 0
        e_corr, t2 = mp2_solver.kernel()

        # Get MP2 density matrix
        dm_mp2 = mp2_solver.make_rdm1(ao_repr=True)

        # Compute dipole with MP2 density
        # Use the unperturbed dipole calculation (field already in density)
        result = self.compute_dipole_moment(density_matrix=dm_mp2)

        return result['dipole_au']
