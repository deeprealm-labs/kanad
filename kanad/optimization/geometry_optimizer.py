"""
Geometry optimization for molecular structures.

Finds minimum energy geometries using analytical gradients from PySCF
combined with scipy's optimization algorithms.
"""

import numpy as np
from scipy.optimize import minimize
from typing import Dict, Any, List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


class GeometryOptimizer:
    """
    Geometry optimizer using analytical gradients.

    Finds the minimum energy molecular structure by optimizing
    atomic positions using gradient-based methods (BFGS, L-BFGS, CG).

    Uses PySCF for analytical gradients (fast, accurate) and scipy
    for optimization algorithms (well-tested, robust).

    Attributes:
        molecule: Molecule object to optimize
        method: Electronic structure method ('HF', 'MP2')
        basis: Basis set for calculations
        energies: Energy history during optimization
        geometries: Geometry history during optimization
        gradients: Gradient history during optimization

    Example:
        >>> from kanad.io import from_smiles
        >>> from kanad.optimization import GeometryOptimizer
        >>>
        >>> # Create molecule with initial geometry
        >>> mol = from_smiles("O")
        >>>
        >>> # Optimize geometry
        >>> opt = GeometryOptimizer(mol, method='HF', basis='6-31g(d)')
        >>> result = opt.optimize()
        >>>
        >>> print(f"Optimized energy: {result['energy']:.6f} Ha")
        >>> print(f"Converged: {result['success']}")
        >>> print(f"Steps: {result['n_steps']}")
    """

    def __init__(
        self,
        molecule,
        method: str = 'HF',
        basis: Optional[str] = None
    ):
        """
        Initialize geometry optimizer.

        Args:
            molecule: Molecule object to optimize
            method: Electronic structure method ('HF', 'MP2')
            basis: Basis set (if None, use molecule's current basis)
        """
        self.molecule = molecule
        self.method = method.upper()
        self.basis = basis if basis is not None else getattr(molecule, 'basis', 'sto-3g')

        # Optimization history
        self.energies = []
        self.geometries = []
        self.gradients = []
        self.step_count = 0

        # Convergence criteria (standard values)
        self.energy_tol = 1e-6  # Ha (very tight)
        self.gradient_tol = 3e-4  # Ha/Bohr (standard in quantum chemistry)
        self.displacement_tol = 1.2e-3  # Bohr

        logger.info(f"GeometryOptimizer initialized: {self.method}/{self.basis}")

    def optimize(
        self,
        algorithm: str = 'BFGS',
        max_steps: int = 100,
        gradient_tol: Optional[float] = None,
        **optimizer_kwargs
    ) -> Dict[str, Any]:
        """
        Optimize molecular geometry to minimum energy.

        Args:
            algorithm: Optimization algorithm
                - 'BFGS': Quasi-Newton (default, recommended)
                - 'L-BFGS-B': Limited-memory BFGS (for large molecules)
                - 'CG': Conjugate gradient (robust but slower)
            max_steps: Maximum optimization steps (default: 100)
            gradient_tol: Gradient convergence threshold (default: 3e-4 Ha/Bohr)
            **optimizer_kwargs: Additional arguments for scipy.optimize.minimize

        Returns:
            dict:
                success (bool): Whether optimization converged
                energy (float): Final energy (Ha)
                geometry (np.ndarray): Optimized atomic positions (Nx3, Angstrom)
                max_gradient (float): Maximum gradient component (Ha/Bohr)
                rms_gradient (float): RMS gradient (Ha/Bohr)
                n_steps (int): Number of optimization steps
                trajectory (list): List of geometries during optimization
                energy_history (list): Energy at each step
                message (str): Optimization status message

        Example:
            >>> result = opt.optimize(algorithm='BFGS', max_steps=50)
            >>> if result['success']:
            ...     print("Optimization converged!")
            ...     print(f"Final energy: {result['energy']:.6f} Ha")
        """
        if gradient_tol is not None:
            self.gradient_tol = gradient_tol

        logger.info(f"Starting geometry optimization with {algorithm}")
        logger.info(f"Max steps: {max_steps}, Gradient tol: {self.gradient_tol:.2e} Ha/Bohr")

        # Get initial geometry as flat array
        initial_coords = self._get_flat_coords()

        # Store initial geometry
        self.geometries.append(initial_coords.copy())

        # Optimize using scipy
        result = minimize(
            fun=self._compute_energy_and_gradient,
            x0=initial_coords,
            jac=True,  # Gradient provided by function
            method=algorithm,
            options={
                'maxiter': max_steps,
                'gtol': self.gradient_tol,
                'disp': False,  # Suppress scipy output
                **optimizer_kwargs
            }
        )

        # Extract final results
        final_coords = result.x
        final_energy = result.fun
        converged = result.success

        # Update molecule with optimized geometry
        self._set_geometry(final_coords)

        # Compute final gradient statistics
        if len(self.gradients) > 0:
            final_gradient = self.gradients[-1]
            max_grad = np.max(np.abs(final_gradient))
            rms_grad = np.sqrt(np.mean(final_gradient**2))
        else:
            max_grad = 0.0
            rms_grad = 0.0

        logger.info(f"Optimization {'converged' if converged else 'did not converge'}")
        logger.info(f"Final energy: {final_energy:.8f} Ha")
        logger.info(f"Max gradient: {max_grad:.6f} Ha/Bohr")
        logger.info(f"RMS gradient: {rms_grad:.6f} Ha/Bohr")
        logger.info(f"Steps taken: {self.step_count}")

        return {
            'success': converged,
            'energy': final_energy,
            'geometry': self._get_atomic_positions(),  # Nx3 array in Angstrom
            'max_gradient': max_grad,
            'rms_gradient': rms_grad,
            'n_steps': self.step_count,
            'trajectory': [g.reshape(-1, 3) for g in self.geometries],  # List of Nx3 arrays
            'energy_history': self.energies.copy(),
            'gradient_history': [g.reshape(-1, 3) for g in self.gradients],
            'message': result.message
        }

    def _compute_energy_and_gradient(
        self,
        coords: np.ndarray
    ) -> Tuple[float, np.ndarray]:
        """
        Compute energy and gradient at given geometry.

        This is called by scipy.optimize at each optimization step.

        Args:
            coords: Flat array of atomic coordinates [x1,y1,z1,x2,y2,z2,...] in Bohr

        Returns:
            tuple: (energy in Ha, gradient in Ha/Bohr as flat array)
        """
        from pyscf import gto, scf

        self.step_count += 1

        # Convert Bohr to Angstrom for internal use
        coords_angstrom = coords * 0.529177  # Bohr → Angstrom

        # Build PySCF molecule (PySCF uses Angstrom by default)
        atom_str = self._make_pyscf_atom_string(coords_angstrom)

        mol = gto.M(
            atom=atom_str,
            basis=self.basis,
            charge=self.molecule.charge,
            spin=self.molecule.spin,
            unit='Angstrom'
        )

        # Run electronic structure calculation
        if self.method == 'HF':
            if self.molecule.spin == 0:
                mf = scf.RHF(mol)
            else:
                mf = scf.ROHF(mol)
        elif self.method == 'MP2':
            from pyscf import mp
            mf = scf.RHF(mol) if self.molecule.spin == 0 else scf.ROHF(mol)
            mf.verbose = 0
            mf.kernel()
            mf = mp.MP2(mf)
        else:
            raise ValueError(f"Unsupported method: {self.method}. Use 'HF' or 'MP2'.")

        mf.verbose = 0
        energy = mf.kernel()

        # Compute analytical gradient
        grad_obj = mf.Gradients()
        grad_obj.verbose = 0
        gradient = grad_obj.kernel()  # Shape: (n_atoms, 3) in Ha/Bohr

        # Flatten gradient
        grad_flat = gradient.flatten()

        # Store in history
        self.energies.append(float(energy))
        self.geometries.append(coords.copy())
        self.gradients.append(grad_flat.copy())

        max_grad = np.max(np.abs(grad_flat))
        logger.debug(f"Step {self.step_count}: E={energy:.8f} Ha, max|g|={max_grad:.6f} Ha/Bohr")

        return float(energy), grad_flat

    def _get_flat_coords(self) -> np.ndarray:
        """
        Get atomic coordinates as flat array in Bohr.

        Returns:
            np.ndarray: Flat array [x1,y1,z1,x2,y2,z2,...] in Bohr
        """
        coords_angstrom = []
        for atom in self.molecule.atoms:
            coords_angstrom.extend(atom.position)
        coords_angstrom = np.array(coords_angstrom)

        # Convert Angstrom to Bohr for scipy optimizer
        coords_bohr = coords_angstrom / 0.529177

        return coords_bohr

    def _set_geometry(self, coords: np.ndarray):
        """
        Update molecule geometry from flat coordinate array.

        Args:
            coords: Flat array in Bohr
        """
        coords_bohr = coords.reshape(-1, 3)
        coords_angstrom = coords_bohr * 0.529177  # Bohr → Angstrom

        for i, atom in enumerate(self.molecule.atoms):
            atom.position = coords_angstrom[i]

    def _make_pyscf_atom_string(self, coords: np.ndarray) -> str:
        """
        Create PySCF atom string from coordinates.

        Args:
            coords: Flat array of coordinates in Angstrom

        Returns:
            str: PySCF atom specification string
        """
        coords_3d = coords.reshape(-1, 3)
        atom_strs = []
        for i, atom in enumerate(self.molecule.atoms):
            x, y, z = coords_3d[i]
            atom_strs.append(f"{atom.symbol} {x:.10f} {y:.10f} {z:.10f}")
        return "; ".join(atom_strs)

    def _get_atomic_positions(self) -> np.ndarray:
        """
        Get current atomic positions.

        Returns:
            np.ndarray: Atomic positions as Nx3 array in Angstrom
        """
        return np.array([atom.position for atom in self.molecule.atoms])
