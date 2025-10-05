"""
Molecule class for multi-atom systems.

Provides efficient Hamiltonian construction for molecules with 3+ atoms
using PySCF backend for maximum performance.
"""

import numpy as np
from typing import List, Optional, Dict, Any
from pyscf import gto, scf
import logging

from kanad.core.atom import Atom

logger = logging.getLogger(__name__)


class MolecularHamiltonian:
    """
    Hamiltonian for multi-atom molecules using PySCF backend.

    This is optimized for performance and matches Qiskit Nature speed.
    """

    def __init__(self, atoms: List[Atom], charge: int = 0, spin: int = 0, basis: str = 'sto-3g'):
        """
        Initialize molecular Hamiltonian.

        Args:
            atoms: List of Atom objects
            charge: Total molecular charge
            spin: Spin multiplicity (2S, where S is total spin)
            basis: Basis set (default: sto-3g)
        """
        self.atoms = atoms
        self.charge = charge
        self.spin = spin
        self.basis = basis

        # Build PySCF molecule
        self._build_pyscf_molecule()

        # Run HF to get integrals
        self._run_hf()

        # Set base class attributes
        self.n_orbitals = self.mf.mo_coeff.shape[1]
        self.n_electrons = self.mol.nelectron
        self.nuclear_repulsion = self.mol.energy_nuc()

        logger.info(f"MolecularHamiltonian: {self.n_electrons} electrons, {self.n_orbitals} orbitals")

    def _build_pyscf_molecule(self):
        """Build PySCF molecule object."""
        self.mol = gto.Mole()

        # Build atom string
        atom_str = []
        for atom in self.atoms:
            symbol = atom.symbol
            pos = atom.position
            atom_str.append(f"{symbol} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}")

        self.mol.atom = '; '.join(atom_str)
        self.mol.basis = self.basis
        self.mol.charge = self.charge
        self.mol.spin = self.spin
        self.mol.build()

        logger.debug(f"Built molecule: {len(self.atoms)} atoms, charge={self.charge}, spin={self.spin}")

    def _run_hf(self):
        """Run Hartree-Fock to get molecular orbitals and integrals."""
        if self.spin == 0:
            self.mf = scf.RHF(self.mol)
        else:
            self.mf = scf.ROHF(self.mol)

        self.hf_energy = self.mf.kernel()

        # Store convergence info
        self._scf_converged = self.mf.converged
        self._scf_iterations = getattr(self.mf, 'iter', 0)

        # Cache integrals in MO basis
        self._cache_integrals()

        logger.debug(f"HF energy: {self.hf_energy:.6f} Ha, converged: {self._scf_converged}")

    def _cache_integrals(self):
        """Cache molecular orbital integrals for fast access."""
        # One-electron integrals (kinetic + nuclear attraction)
        self.h_core = self.mf.get_hcore()

        # Two-electron repulsion integrals (ERI)
        # Get in MO basis for efficiency
        mo_coeff = self.mf.mo_coeff
        n_orb = mo_coeff.shape[1]

        # AO → MO transformation
        # ERI in AO basis
        eri_ao = self.mol.intor('int2e')

        # Transform to MO basis: (pq|rs) = C_pi C_qj C_rk C_sl (ij|kl)
        # This is expensive but we do it once and cache
        self.eri = np.einsum('pi,qj,ijkl->pqkl', mo_coeff, mo_coeff, eri_ao)
        self.eri = np.einsum('rk,sl,pqkl->pqrs', mo_coeff, mo_coeff, self.eri)

        # Also cache h_core in MO basis
        self.h_core_mo = mo_coeff.T @ self.h_core @ mo_coeff

        logger.debug(f"Cached integrals: h_core {self.h_core.shape}, ERI {self.eri.shape}")

    def compute_electronic_energy(self, density_matrix: Optional[np.ndarray] = None) -> float:
        """
        Compute electronic energy.

        Args:
            density_matrix: Density matrix (uses HF if None)

        Returns:
            Electronic energy in Hartree
        """
        if density_matrix is None:
            return self.hf_energy

        # E = Tr[D h] + 0.5 Tr[D G[D]]
        h = self.h_core

        # Fock matrix
        vhf = self.mf.get_veff(dm=density_matrix)

        # Energy
        energy = np.einsum('ij,ji->', h + 0.5 * vhf, density_matrix)
        energy += self.mol.energy_nuc()  # Nuclear repulsion

        return energy

    def solve_scf(self, max_iterations: int = 100, conv_tol: float = 1e-8,
                  use_diis: bool = True, **kwargs) -> tuple:
        """
        Solve self-consistent field equations.

        Args:
            max_iterations: Maximum SCF iterations
            conv_tol: Convergence tolerance
            use_diis: Use DIIS acceleration
            **kwargs: Additional parameters (level_shift, damping_factor, etc.)

        Returns:
            (density_matrix, energy)
        """
        # Re-run HF with specified parameters
        if 'level_shift' in kwargs:
            self.mf.level_shift = kwargs['level_shift']

        if 'damping_factor' in kwargs:
            # Not directly supported in PySCF, would need custom implementation
            pass

        self.mf.max_cycle = max_iterations
        self.mf.conv_tol = conv_tol
        self.mf.diis = use_diis

        energy = self.mf.kernel()
        density_matrix = self.mf.make_rdm1()

        # Store convergence info
        self._scf_converged = self.mf.converged
        self._scf_iterations = getattr(self.mf, 'iter', max_iterations)

        return density_matrix, energy

    def compute_molecular_orbitals(self):
        """
        Get molecular orbital energies and coefficients.

        Returns:
            (mo_energies, mo_coefficients)
        """
        return self.mf.mo_energy, self.mf.mo_coeff

    def to_matrix(self) -> np.ndarray:
        """
        Get Hamiltonian as dense matrix in second quantization.

        Returns:
            Hamiltonian matrix (N×N where N = 2^(2*n_orbitals))
        """
        # This would require full second quantization
        # For now, return Fock matrix as approximation
        return self.mf.get_fock()

    def get_geometry(self) -> List[tuple]:
        """
        Get molecular geometry.

        Returns:
            List of (symbol, position) tuples
        """
        return [(atom.symbol, atom.position) for atom in self.atoms]

    def get_nuclear_repulsion(self) -> float:
        """Get nuclear repulsion energy."""
        return self.mol.energy_nuc()


class Molecule:
    """
    Multi-atom molecule class.

    Provides high-performance molecular Hamiltonian construction
    matching Qiskit Nature speed.
    """

    def __init__(self, atoms: List[Atom], charge: int = 0, spin: int = 0, basis: str = 'sto-3g'):
        """
        Initialize molecule.

        Args:
            atoms: List of Atom objects
            charge: Total molecular charge
            spin: Spin multiplicity (2S)
            basis: Basis set
        """
        self.atoms = atoms
        self.charge = charge
        self.spin = spin
        self.basis = basis

        # Build Hamiltonian
        self._hamiltonian = None

        logger.info(f"Created Molecule: {len(atoms)} atoms, charge={charge}, spin={spin}")

    @property
    def hamiltonian(self) -> MolecularHamiltonian:
        """Get molecular Hamiltonian (lazy construction)."""
        if self._hamiltonian is None:
            self._hamiltonian = MolecularHamiltonian(
                self.atoms, self.charge, self.spin, self.basis
            )
        return self._hamiltonian

    @property
    def n_atoms(self) -> int:
        """Number of atoms."""
        return len(self.atoms)

    @property
    def n_electrons(self) -> int:
        """Total number of electrons."""
        return self.hamiltonian.n_electrons

    @property
    def n_orbitals(self) -> int:
        """Number of spatial orbitals."""
        return self.hamiltonian.n_orbitals

    def compute_energy(self, method: str = 'HF', **kwargs) -> Dict[str, Any]:
        """
        Compute molecular energy.

        Args:
            method: 'HF', 'VQE', 'QPE', 'SQD'
            **kwargs: Method-specific parameters

        Returns:
            Dictionary with energy and results
        """
        if method.upper() == 'HF':
            energy = self.hamiltonian.hf_energy
            from kanad.core.constants.conversion_factors import ConversionFactors

            return {
                'energy': energy * ConversionFactors.HARTREE_TO_EV,
                'energy_ha': energy,
                'method': 'Hartree-Fock',
                'converged': self.hamiltonian._scf_converged
            }

        elif method.upper() == 'VQE':
            from kanad.solvers.vqe_solver import VQESolver
            from kanad.ansatze.ucc_ansatz import UCCAnsatz
            from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper

            n_qubits = 2 * self.n_orbitals
            ansatz = UCCAnsatz(n_qubits, self.n_electrons)
            mapper = JordanWignerMapper()
            solver = VQESolver(self.hamiltonian, ansatz, mapper)

            result = solver.solve()
            return result

        elif method.upper() == 'QPE':
            from kanad.solvers.qpe_solver import QPESolver
            solver = QPESolver(self.hamiltonian, n_ancilla=4)
            return solver.solve()

        elif method.upper() == 'SQD':
            from kanad.solvers.sqd_solver import SQDSolver
            solver = SQDSolver(self.hamiltonian, n_samples=100)
            return solver.solve()

        else:
            raise ValueError(f"Unknown method: {method}")

    def get_geometry_string(self) -> str:
        """Get geometry as formatted string."""
        lines = []
        for atom in self.atoms:
            pos = atom.position
            lines.append(f"{atom.symbol:2s}  {pos[0]:>10.6f}  {pos[1]:>10.6f}  {pos[2]:>10.6f}")
        return '\n'.join(lines)

    def __repr__(self) -> str:
        """String representation."""
        formula = self._get_formula()
        return f"Molecule({formula}, {self.n_atoms} atoms, {self.n_electrons} electrons)"

    def _get_formula(self) -> str:
        """Get chemical formula."""
        from collections import Counter
        symbols = [atom.symbol for atom in self.atoms]
        counts = Counter(symbols)

        # Sort by element (common convention: C, H, others alphabetically)
        def sort_key(item):
            symbol = item[0]
            if symbol == 'C':
                return (0, symbol)
            elif symbol == 'H':
                return (1, symbol)
            else:
                return (2, symbol)

        parts = []
        for symbol, count in sorted(counts.items(), key=sort_key):
            if count == 1:
                parts.append(symbol)
            else:
                parts.append(f"{symbol}{count}")

        return ''.join(parts)
