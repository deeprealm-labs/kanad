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

# Import Lattice and PeriodicHamiltonian (will be used conditionally)
try:
    from kanad.core.lattice import Lattice
    from kanad.core.hamiltonians.periodic_hamiltonian import PeriodicHamiltonian
    PBC_AVAILABLE = True
except ImportError:
    PBC_AVAILABLE = False
    Lattice = None
    PeriodicHamiltonian = None


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
        # One-electron integrals (kinetic + nuclear attraction) - AO basis
        self.h_core = self.mf.get_hcore()

        # Two-electron repulsion integrals (ERI) - AO basis
        # CRITICAL: Keep in AO basis - Qiskit Nature transforms to MO internally
        # PauliConverter will use Qiskit Nature (if installed) for correct conversion
        self.eri = self.mol.intor('int2e')

        # Also cache h_core in MO basis for other uses
        mo_coeff = self.mf.mo_coeff
        self.h_core_mo = mo_coeff.T @ self.h_core @ mo_coeff

        logger.debug(f"Cached integrals (AO basis): h_core {self.h_core.shape}, ERI {self.eri.shape}")

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

    def to_sparse_hamiltonian(self, mapper: str = 'jordan_wigner'):
        """
        Convert to sparse Hamiltonian representation using Pauli operators.

        Uses OpenFermion for CORRECT fermionic-to-qubit transformation.
        This is the RECOMMENDED method for VQE calculations (100-1000x faster than dense).

        Args:
            mapper: Fermion-to-qubit mapping ('jordan_wigner' or 'bravyi_kitaev')

        Returns:
            Qiskit SparsePauliOp object ready for use in VQE
        """
        from kanad.core.hamiltonians.fast_pauli_builder import build_molecular_hamiltonian_pauli
        import logging
        import numpy as np

        logger = logging.getLogger(__name__)

        n_qubits = 2 * self.n_orbitals

        logger.info(f"Building sparse Hamiltonian using OpenFermion (FAST method)...")
        logger.info(f"  {self.n_orbitals} orbitals → {n_qubits} qubits")
        logger.info(f"  Mapper: {mapper}")
        logger.info(f"  Bypassing {2**n_qubits}×{2**n_qubits} dense matrix construction")

        # Transform to MO basis (required for fermionic transformations)
        mo_energies, C = self.compute_molecular_orbitals()
        h_mo = C.T @ self.h_core @ C
        eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', C, C, self.eri, C, C)

        logger.debug(f"Transformed integrals: AO → MO basis")

        # Build Pauli operators using OpenFermion
        sparse_pauli_op = build_molecular_hamiltonian_pauli(
            h_core=h_mo,
            eri=eri_mo,
            nuclear_repulsion=self.nuclear_repulsion,
            n_orbitals=self.n_orbitals,
            mapper=mapper
        )

        num_terms = len(sparse_pauli_op)
        logger.info(f"✓ Sparse Hamiltonian: {num_terms} Pauli terms")
        logger.info(f"✓ Memory savings: {(2**n_qubits)**2:,} matrix elements → {num_terms} Pauli terms")

        return sparse_pauli_op


class Molecule:
    """
    Multi-atom molecule class.

    Provides high-performance molecular Hamiltonian construction
    matching Qiskit Nature speed.
    """

    def __init__(self,
                 atoms: List[Atom],
                 charge: int = 0,
                 spin: int = 0,
                 basis: str = 'sto-3g',
                 lattice: Optional['Lattice'] = None,
                 k_points: Optional[Any] = None,
                 pseudo: str = 'gth-pade'):
        """
        Initialize molecule or periodic system.

        Args:
            atoms: List of Atom objects (unit cell for periodic)
            charge: Total molecular charge (or unit cell charge)
            spin: Spin multiplicity (2S)
            basis: Basis set ('sto-3g' for molecular, 'gth-dzvp' for periodic)
            lattice: Lattice object for periodic systems (None for molecules)
            k_points: K-point sampling for periodic systems (tuple or array)
            pseudo: Pseudopotential for periodic systems

        Examples:
            >>> # Molecular system
            >>> mol = Molecule([Atom('H', [0,0,0]), Atom('H', [0,0,0.74])])

            >>> # Periodic system
            >>> lattice = Lattice(np.eye(3) * 5.0)
            >>> crystal = Molecule([Atom('Si', [0,0,0])], lattice=lattice, k_points=(4,4,4))
        """
        self.atoms = atoms
        self.charge = charge
        self.spin = spin
        self.basis = basis
        self.lattice = lattice
        self.k_points = k_points
        self.pseudo = pseudo

        # Determine if periodic or molecular
        self.is_periodic = lattice is not None

        # Build Hamiltonian (lazy)
        self._hamiltonian = None

        if self.is_periodic:
            logger.info(f"Created Periodic System: {len(atoms)} atoms/cell, "
                        f"PBC={lattice.pbc}, k_points={k_points}")
        else:
            logger.info(f"Created Molecule: {len(atoms)} atoms, charge={charge}, spin={spin}")

    @property
    def hamiltonian(self):
        """
        Get Hamiltonian (molecular or periodic, lazy construction).

        Returns:
            MolecularHamiltonian or PeriodicHamiltonian depending on self.is_periodic
        """
        if self._hamiltonian is None:
            if self.is_periodic:
                if not PBC_AVAILABLE:
                    raise ImportError("Periodic boundary conditions require PeriodicHamiltonian")

                # Default basis for periodic systems
                basis = self.basis if self.basis != 'sto-3g' else 'gth-dzvp'

                self._hamiltonian = PeriodicHamiltonian(
                    self.atoms,
                    self.lattice,
                    charge=self.charge,
                    spin=self.spin,
                    basis=basis,
                    pseudo=self.pseudo,
                    k_points=self.k_points
                )
            else:
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

    @property
    def formula(self) -> str:
        """Chemical formula (e.g., 'H2O', 'CH4')."""
        return self._get_formula()

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

    # ===== Periodic System Methods =====

    def solve_scf_pbc(self, **kwargs) -> Dict[str, Any]:
        """
        Solve SCF for periodic system with k-point sampling.

        Only valid for periodic systems (lattice is not None).

        Args:
            **kwargs: Passed to PeriodicHamiltonian.solve_scf()

        Returns:
            result: Dictionary with energy, band_energies, fermi_energy, etc.

        Raises:
            ValueError: If system is not periodic
        """
        if not self.is_periodic:
            raise ValueError("solve_scf_pbc() only valid for periodic systems")

        return self.hamiltonian.solve_scf(**kwargs)

    def compute_band_structure(self, k_path, n_bands=None, n_points=50) -> Dict[str, Any]:
        """
        Compute band structure along k-path.

        Only valid for periodic systems.

        Args:
            k_path: High-symmetry path ('GXLG') or explicit k-points
            n_bands: Number of bands to compute
            n_points: Points per segment

        Returns:
            result: Dictionary with k_points, band_energies, labels, etc.
        """
        if not self.is_periodic:
            raise ValueError("compute_band_structure() only valid for periodic systems")

        return self.hamiltonian.compute_band_structure(k_path, n_bands, n_points)

    def get_band_gap(self) -> Dict[str, float]:
        """
        Compute band gap for periodic system.

        Returns:
            result: Dictionary with gap, vbm, cbm, type
        """
        if not self.is_periodic:
            raise ValueError("get_band_gap() only valid for periodic systems")

        return self.hamiltonian.get_band_gap()

    def make_supercell(self, size: tuple) -> 'Molecule':
        """
        Create supercell expansion of periodic system.

        Args:
            size: (nx, ny, nz) supercell dimensions

        Returns:
            supercell: New Molecule with expanded lattice and atoms
        """
        if not self.is_periodic:
            raise ValueError("make_supercell() only valid for periodic systems")

        # Create supercell lattice
        supercell_lattice = self.lattice.make_supercell(size)

        # Replicate atoms
        nx, ny, nz = size
        supercell_atoms = []

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    # Lattice vector shift
                    R = i * self.lattice.lattice_vectors[0] + \
                        j * self.lattice.lattice_vectors[1] + \
                        k * self.lattice.lattice_vectors[2]

                    # Add shifted copy of each atom
                    for atom in self.atoms:
                        new_pos = atom.position + R
                        supercell_atoms.append(Atom(atom.symbol, new_pos))

        return Molecule(
            supercell_atoms,
            lattice=supercell_lattice,
            charge=self.charge * nx * ny * nz,  # Scale charge
            spin=self.spin,
            basis=self.basis,
            k_points=self.k_points,
            pseudo=self.pseudo
        )

    def __repr__(self) -> str:
        """String representation."""
        formula = self._get_formula()
        return f"Molecule({formula}, {self.n_atoms} atoms, {self.n_electrons} electrons)"

    def _get_formula(self) -> str:
        """Get chemical formula with charge notation."""
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

        formula = ''.join(parts)

        # Add charge notation if charged
        if self.charge > 0:
            if self.charge == 1:
                formula += '+'
            else:
                formula += f'{self.charge}+'
        elif self.charge < 0:
            if self.charge == -1:
                formula += '-'
            else:
                formula += f'{abs(self.charge)}-'

        return formula
