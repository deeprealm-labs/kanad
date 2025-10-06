"""
Electron correlation methods.

Provides post-Hartree-Fock methods for improved accuracy:
- MP2 (Møller-Plesset 2nd order perturbation theory)
- Future: MP3, CCSD, CCSD(T)

These methods add electron correlation on top of the mean-field HF solution,
significantly improving accuracy for energies and molecular properties.
"""

import numpy as np
from typing import Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)


class MP2Solver:
    """
    MP2 (Møller-Plesset 2nd order perturbation theory) solver.

    Adds electron correlation correction to Hartree-Fock:
        E_MP2 = E_HF + E_corr

    where the correlation energy is:
        E_corr = -1/4 Σ_ijab |⟨ij||ab⟩|² / (ε_i + ε_j - ε_a - ε_b)

    Notation:
        i,j: occupied molecular orbitals
        a,b: virtual (unoccupied) molecular orbitals
        ⟨ij||ab⟩: antisymmetrized two-electron integral
        ε: orbital energies

    MP2 typically recovers 80-90% of the correlation energy and improves
    molecular properties significantly over HF.

    Uses PySCF backend for efficient, well-tested implementation.

    Attributes:
        hamiltonian: Hamiltonian object with converged HF reference
        mf: PySCF mean-field object
        mol: PySCF molecule object

    Example:
        >>> from kanad.io import from_smiles
        >>> mol = from_smiles("O", basis='6-311g(d,p)')
        >>> mp2 = MP2Solver(mol.hamiltonian)
        >>> result = mp2.compute_energy()
        >>> print(f"MP2 energy: {result['e_mp2']:.6f} Ha")
        >>> print(f"Correlation: {result['e_corr']:.6f} Ha")
    """

    def __init__(self, hamiltonian):
        """
        Initialize MP2 solver.

        Args:
            hamiltonian: Hamiltonian object with converged HF reference

        Raises:
            ValueError: If HF calculation has not converged
        """
        self.hamiltonian = hamiltonian
        self.mf = hamiltonian.mf
        self.mol = hamiltonian.mol

        if not self.mf.converged:
            raise ValueError(
                "HF calculation must converge before MP2. "
                "Check that hamiltonian.mf.converged is True."
            )

        logger.debug(f"MP2Solver initialized for {self.mol.atom}")

    def compute_energy(self) -> Dict[str, Any]:
        """
        Compute MP2 correlation energy.

        Runs MP2 calculation on top of the converged HF reference
        and returns energy components.

        Returns:
            dict:
                e_hf (float): Hartree-Fock energy (Ha)
                e_corr (float): MP2 correlation energy (Ha, negative)
                e_mp2 (float): Total MP2 energy = e_hf + e_corr (Ha)
                t2_amplitudes (np.ndarray): T2 cluster amplitudes (for advanced use)
                mp2_solver (pyscf.mp.MP2): PySCF MP2 solver object

        Example:
            >>> result = mp2.compute_energy()
            >>> print(f"HF:        {result['e_hf']:.8f} Ha")
            >>> print(f"MP2 corr:  {result['e_corr']:.8f} Ha")
            >>> print(f"MP2 total: {result['e_mp2']:.8f} Ha")
        """
        from pyscf import mp

        logger.info("Computing MP2 correlation energy...")
        logger.debug(f"HF reference energy: {self.mf.e_tot:.8f} Ha")

        # Create MP2 solver
        # For RHF reference, use standard MP2
        # For ROHF reference, use ROMP2
        if isinstance(self.mf, (type(self.mf).__bases__[0])):
            mp2_solver = mp.MP2(self.mf)
        else:
            mp2_solver = mp.MP2(self.mf)

        # Compute correlation energy
        # kernel() returns (e_corr, t2_amplitudes)
        e_corr, t2 = mp2_solver.kernel()
        e_mp2 = self.mf.e_tot + e_corr

        logger.info(f"MP2 correlation energy: {e_corr:.8f} Ha")
        logger.info(f"MP2 total energy:       {e_mp2:.8f} Ha")
        logger.debug(f"Correlation recovery:   {abs(e_corr):.6f} Ha")

        return {
            'e_hf': self.mf.e_tot,
            'e_corr': e_corr,
            'e_mp2': e_mp2,
            't2_amplitudes': t2,
            'mp2_solver': mp2_solver
        }

    def make_rdm1(self, ao_repr: bool = True) -> np.ndarray:
        """
        Compute MP2 1-electron reduced density matrix.

        The MP2 density matrix includes correlation effects and can be used
        for computing correlated molecular properties (dipole moment, etc.).

        The density matrix satisfies:
            Tr(P) = N_electrons
            E_MP2 = Tr(P·H) + V_nn + E_corr

        Args:
            ao_repr (bool): If True, return density in AO basis (default).
                           If False, return in MO basis.

        Returns:
            np.ndarray: MP2 density matrix
                - Shape: (n_ao, n_ao) if ao_repr=True
                - Shape: (n_mo, n_mo) if ao_repr=False

        Example:
            >>> dm_mp2 = mp2.make_rdm1()
            >>> n_elec = np.trace(dm_mp2)
            >>> print(f"Electron count: {n_elec:.1f}")
        """
        from pyscf import mp

        logger.debug("Computing MP2 density matrix...")

        # Create MP2 solver and run
        mp2_solver = mp.MP2(self.mf)
        mp2_solver.kernel()

        # Get MP2 density matrix
        # By default, PySCF returns in AO basis
        if ao_repr:
            dm_mp2 = mp2_solver.make_rdm1(ao_repr=True)
        else:
            dm_mp2 = mp2_solver.make_rdm1(ao_repr=False)

        logger.debug(f"MP2 density matrix shape: {dm_mp2.shape}")
        logger.debug(f"Trace (electron count): {np.trace(dm_mp2):.6f}")

        return dm_mp2

    def make_rdm2(self, ao_repr: bool = False) -> np.ndarray:
        """
        Compute MP2 2-electron reduced density matrix.

        This is the 2-particle density matrix including correlation.
        Primarily for advanced calculations.

        Args:
            ao_repr (bool): If True, return in AO basis (very large!).
                           If False (default), return in MO basis.

        Returns:
            np.ndarray: MP2 2-electron density matrix
                - Shape: (n_mo, n_mo, n_mo, n_mo)

        Warning:
            The 2-RDM is very large (N⁴) and expensive to compute.
            Only use for small molecules.
        """
        from pyscf import mp

        logger.warning("Computing MP2 2-RDM - this is expensive for large molecules!")

        mp2_solver = mp.MP2(self.mf)
        mp2_solver.kernel()

        dm2_mp2 = mp2_solver.make_rdm2(ao_repr=ao_repr)

        logger.debug(f"MP2 2-RDM shape: {dm2_mp2.shape}")

        return dm2_mp2


class MP3Solver:
    """
    MP3 (Møller-Plesset 3rd order) solver.

    Placeholder for future implementation.
    MP3 adds 3rd-order correlation correction.
    """

    def __init__(self, hamiltonian):
        raise NotImplementedError(
            "MP3 not yet implemented. Use MP2Solver for now."
        )


class CCSDSolver:
    """
    CCSD (Coupled Cluster Singles and Doubles) solver.

    Placeholder for future implementation.
    CCSD is more accurate than MP2 (typically 95%+ of correlation energy)
    but more expensive.
    """

    def __init__(self, hamiltonian):
        raise NotImplementedError(
            "CCSD not yet implemented. Use MP2Solver for now."
        )
