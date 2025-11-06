"""
Energy gradients (forces) for molecular systems.

Provides analytical and numerical gradients for:
- Geometry optimization
- Molecular dynamics
- Transition state searches
- Vibrational frequency calculations
"""

import numpy as np
from typing import Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)


class GradientCalculator:
    """
    Calculate energy gradients (forces) for molecules.

    Computes ∂E/∂R_A (gradient) and F_A = -∂E/∂R_A (force) for each atom.

    Provides both:
    - Analytical gradients (from PySCF, fast and accurate)
    - Numerical gradients (finite difference, for validation)

    At equilibrium geometry, all forces should be zero.

    Attributes:
        hamiltonian: Hamiltonian object with molecular system
        method: Electronic structure method ('HF', 'MP2')

    Example:
        >>> from kanad.io import from_smiles
        >>> from kanad.core import GradientCalculator
        >>>
        >>> mol = from_smiles("O")
        >>> grad_calc = GradientCalculator(mol)
        >>>
        >>> # Compute gradient
        >>> result = grad_calc.compute_gradient()
        >>> print(f"Max force: {result['max_force']:.6f} Ha/Bohr")
        >>> print(f"RMS force: {result['rms_force']:.6f} Ha/Bohr")
        >>>
        >>> # Verify against numerical gradient
        >>> verification = grad_calc.verify_analytical_gradient()
        >>> print(f"Agreement: {verification['agree']}")
    """

    def __init__(self, molecule, method: str = 'HF'):
        """
        Initialize gradient calculator.

        Args:
            molecule: Molecule object with converged electronic structure
            method: Electronic structure method ('HF', 'MP2')

        Raises:
            ValueError: If method is not supported
        """
        self.molecule = molecule
        if hasattr(molecule, 'hamiltonian'):
            self.hamiltonian = molecule.hamiltonian
        else:
            self.hamiltonian = None
        self.method = method.upper()
        self._mf_cache = None  # Cache for mean-field object

        if self.method not in ['HF', 'MP2']:
            raise ValueError(f"Method '{method}' not supported. Use 'HF' or 'MP2'.")

        logger.debug(f"GradientCalculator initialized: {self.method}")

    def _get_pyscf_mf(self):
        """
        Get or create PySCF mean-field object.

        Returns:
            PySCF mean-field object (RHF or RKS)
        """
        # Try to get from hamiltonian if available
        if self.hamiltonian is not None and hasattr(self.hamiltonian, 'mf'):
            return self.hamiltonian.mf

        # Create fresh PySCF calculation
        from pyscf import gto, scf

        # Build PySCF molecule
        atoms = []
        if hasattr(self.molecule, 'atoms'):
            # Full molecule
            for atom in self.molecule.atoms:
                pos = atom.position  # Already in Bohr
                atoms.append([atom.symbol, pos])
        elif hasattr(self.molecule, 'atom_1') and hasattr(self.molecule, 'atom_2'):
            # Bond (diatomic)
            atoms.append([self.molecule.atom_1.symbol, self.molecule.atom_1.position])
            atoms.append([self.molecule.atom_2.symbol, self.molecule.atom_2.position])
        else:
            raise ValueError("Cannot extract atoms from molecule/bond object")

        # Determine charge and spin
        charge = getattr(self.molecule, 'charge', 0)
        multiplicity = getattr(self.molecule, 'multiplicity', 1)

        # Create PySCF molecule
        mol = gto.Mole()
        mol.atom = atoms
        mol.basis = getattr(self.molecule, 'basis_set', 'sto-3g')
        mol.charge = charge
        mol.spin = multiplicity - 1
        mol.unit = 'Bohr'
        mol.verbose = 0
        mol.build()

        # Run HF
        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.kernel()

        # Cache it
        self._mf_cache = mf

        return mf

    def compute_gradient(self) -> Dict[str, Any]:
        """
        Compute analytical energy gradient using PySCF.

        Returns gradient ∂E/∂R_A and forces F_A = -∂E/∂R_A for each atom.

        Returns:
            dict:
                gradient (np.ndarray): Energy gradient (n_atoms, 3) in Ha/Bohr
                forces (np.ndarray): Forces on atoms = -gradient (Ha/Bohr)
                max_force (float): Maximum force component magnitude (Ha/Bohr)
                rms_force (float): RMS force (Ha/Bohr)
                max_force_atom (int): Index of atom with maximum force
                method (str): 'analytical'

        Example:
            >>> result = grad_calc.compute_gradient()
            >>> if result['max_force'] < 1e-3:
            ...     print("Geometry is at minimum (forces near zero)")
        """
        logger.info(f"Computing analytical gradient ({self.method})...")

        # Get PySCF mean-field object
        mf = self._get_pyscf_mf()

        # Get PySCF gradient object
        if self.method == 'HF':
            grad_obj = mf.Gradients()
        elif self.method == 'MP2':
            from pyscf import mp
            # Run MP2 if not already done
            mp2_solver = mp.MP2(mf)
            mp2_solver.verbose = 0
            mp2_solver.kernel()
            grad_obj = mp2_solver.Gradients()
        else:
            raise ValueError(f"Unsupported method: {self.method}")

        grad_obj.verbose = 0

        # Compute gradient (n_atoms, 3) in Ha/Bohr
        gradient = grad_obj.kernel()

        # Forces = -gradient
        forces = -gradient

        # Compute statistics
        force_magnitudes = np.linalg.norm(gradient, axis=1)
        max_force_atom = np.argmax(force_magnitudes)
        max_force = np.max(np.abs(gradient))
        rms_force = np.sqrt(np.mean(gradient**2))

        logger.info(f"Max force: {max_force:.6f} Ha/Bohr (atom {max_force_atom})")
        logger.info(f"RMS force: {rms_force:.6f} Ha/Bohr")

        # Get energy from mean-field object
        if self.method == 'HF':
            energy = mf.e_tot
        elif self.method == 'MP2':
            energy = mp2_solver.e_tot
        else:
            energy = mf.e_tot

        logger.info(f"Energy: {energy:.6f} Ha")

        return {
            'gradient': gradient,
            'forces': forces,
            'energy': energy,
            'max_force': max_force,
            'rms_force': rms_force,
            'max_force_atom': max_force_atom,
            'method': 'analytical'
        }

    def compute_numerical_gradient(
        self,
        step: float = 0.001,
        method: str = 'central'
    ) -> Dict[str, Any]:
        """
        Compute numerical gradient using finite differences.

        Useful for validating analytical gradients.

        Args:
            step: Finite difference step size in Bohr (default: 0.001)
            method: Finite difference method:
                - 'central': (E(R+h) - E(R-h)) / 2h (more accurate)
                - 'forward': (E(R+h) - E(R)) / h (faster)

        Returns:
            dict: Same format as compute_gradient()

        Warning:
            Numerical gradients are expensive (6*N energy evaluations for
            central differences, where N = number of atoms).
        """
        logger.info(f"Computing numerical gradient (method={method}, step={step} Bohr)...")

        from pyscf import gto, scf

        n_atoms = len(self.molecule.atoms)
        gradient = np.zeros((n_atoms, 3))

        # Get current geometry in Bohr
        coords = np.array([atom.position for atom in self.molecule.atoms])
        coords_bohr = coords / 0.529177  # Angstrom → Bohr

        # Finite difference for each atom and coordinate
        for i in range(n_atoms):
            for j in range(3):  # x, y, z
                if method == 'central':
                    # E(R + h)
                    coords_plus = coords_bohr.copy()
                    coords_plus[i, j] += step
                    e_plus = self._compute_energy_at_geometry(coords_plus)

                    # E(R - h)
                    coords_minus = coords_bohr.copy()
                    coords_minus[i, j] -= step
                    e_minus = self._compute_energy_at_geometry(coords_minus)

                    # Central difference: (E(R+h) - E(R-h)) / 2h
                    gradient[i, j] = (e_plus - e_minus) / (2 * step)

                elif method == 'forward':
                    # E(R)
                    e_0 = self._compute_energy_at_geometry(coords_bohr)

                    # E(R + h)
                    coords_plus = coords_bohr.copy()
                    coords_plus[i, j] += step
                    e_plus = self._compute_energy_at_geometry(coords_plus)

                    # Forward difference: (E(R+h) - E(R)) / h
                    gradient[i, j] = (e_plus - e_0) / step

                logger.debug(f"Numerical gradient: atom {i}, coord {j}: {gradient[i,j]:.6f}")

        forces = -gradient
        max_force = np.max(np.abs(gradient))
        rms_force = np.sqrt(np.mean(gradient**2))

        logger.info(f"Numerical gradient computed: max={max_force:.6f}, rms={rms_force:.6f}")

        return {
            'gradient': gradient,
            'forces': forces,
            'max_force': max_force,
            'rms_force': rms_force,
            'method': f'numerical_{method}'
        }

    def _compute_energy_at_geometry(self, coords_bohr: np.ndarray) -> float:
        """
        Compute energy at given geometry.

        Args:
            coords_bohr: Atomic coordinates in Bohr (n_atoms, 3)

        Returns:
            float: Total energy in Hartree
        """
        from pyscf import gto, scf

        # Convert Bohr to Angstrom for PySCF
        coords_ang = coords_bohr * 0.529177

        # Build atom string
        atom_strs = []
        for i, atom in enumerate(self.molecule.atoms):
            x, y, z = coords_ang[i]
            atom_strs.append(f"{atom.symbol} {x:.10f} {y:.10f} {z:.10f}")
        atom_str = "; ".join(atom_strs)

        # Create PySCF molecule
        mol = gto.M(
            atom=atom_str,
            basis=self.molecule.basis,
            charge=self.molecule.charge,
            spin=self.molecule.spin,
            unit='Angstrom'
        )

        # Run electronic structure calculation
        if self.molecule.spin == 0:
            mf = scf.RHF(mol)
        else:
            mf = scf.ROHF(mol)

        mf.verbose = 0
        energy = mf.kernel()

        return float(energy)

    def verify_analytical_gradient(
        self,
        numerical_step: float = 0.001,
        tolerance: float = 1e-5
    ) -> Dict[str, Any]:
        """
        Verify analytical gradient against numerical gradient.

        Computes both analytical and numerical gradients and compares them.
        Useful for debugging and validation.

        Args:
            numerical_step: Step size for numerical gradient (Bohr)
            tolerance: Acceptable difference (Ha/Bohr)

        Returns:
            dict:
                agree (bool): Whether gradients agree within tolerance
                max_difference (float): Maximum absolute difference (Ha/Bohr)
                rms_difference (float): RMS difference (Ha/Bohr)
                analytical (np.ndarray): Analytical gradient
                numerical (np.ndarray): Numerical gradient
                difference (np.ndarray): Absolute difference

        Example:
            >>> verification = grad_calc.verify_analytical_gradient()
            >>> if verification['agree']:
            ...     print("Analytical gradient is correct!")
            >>> else:
            ...     print(f"Difference: {verification['max_difference']:.2e}")
        """
        logger.info("Verifying analytical gradient against numerical gradient...")

        # Compute both gradients
        result_analytical = self.compute_gradient()
        result_numerical = self.compute_numerical_gradient(step=numerical_step, method='central')

        grad_analytical = result_analytical['gradient']
        grad_numerical = result_numerical['gradient']

        # Compute differences
        difference = np.abs(grad_analytical - grad_numerical)
        max_diff = np.max(difference)
        rms_diff = np.sqrt(np.mean(difference**2))
        agree = max_diff < tolerance

        logger.info(f"Max difference: {max_diff:.2e} Ha/Bohr")
        logger.info(f"RMS difference: {rms_diff:.2e} Ha/Bohr")
        logger.info(f"Agreement (<{tolerance:.2e}): {'YES ✓' if agree else 'NO ✗'}")

        return {
            'agree': agree,
            'max_difference': max_diff,
            'rms_difference': rms_diff,
            'analytical': grad_analytical,
            'numerical': grad_numerical,
            'difference': difference,
            'tolerance': tolerance
        }
