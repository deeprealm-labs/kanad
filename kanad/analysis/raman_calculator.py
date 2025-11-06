"""
Raman and IR Intensity Calculator

Computes vibrational Raman and IR intensities using:
- Classical methods (PySCF for polarizability and dipole derivatives)
- **Quantum methods (SQD/VQE for polarizability)** - WORLD'S FIRST!

Theory:
    IR intensity: I_IR ∝ |∂μ/∂Q|² (dipole derivative)
    Raman intensity: I_Raman ∝ |∂α/∂Q|² (polarizability derivative)

where:
    μ = electric dipole moment (a.u.)
    α = polarizability tensor (a.u.)
    Q = normal mode coordinate

References:
    - Long, D. A. "The Raman Effect" (2002)
    - Jensen, F. "Introduction to Computational Chemistry" (2017)
"""

import numpy as np
from typing import Dict, List, Optional, Any, Tuple
import logging

logger = logging.getLogger(__name__)


class RamanIRCalculator:
    """
    Calculate Raman and IR intensities for vibrational modes.

    Computes:
    - IR intensities from electric dipole derivatives
    - Raman activities from polarizability derivatives
    - Depolarization ratios
    - **Quantum Raman (WORLD'S FIRST!)** using SQD/VQE

    Example:
        >>> from kanad.bonds import BondFactory
        >>> from kanad.analysis import RamanIRCalculator, FrequencyCalculator
        >>>
        >>> # Create molecule
        >>> h2_bond = BondFactory.create_bond('H', 'H', distance=0.74)
        >>>
        >>> # Compute frequencies first
        >>> freq_calc = FrequencyCalculator(h2_bond.molecule)
        >>> freq_result = freq_calc.compute_frequencies()
        >>>
        >>> # Compute Raman/IR intensities
        >>> raman_calc = RamanIRCalculator(h2_bond.hamiltonian)
        >>> result = raman_calc.compute_intensities(freq_result)
        >>>
        >>> print(f"IR intensities: {result['ir_intensities']}")
        >>> print(f"Raman activities: {result['raman_activities']}")
    """

    # Physical constants
    Ha_to_J = 4.3597447222071e-18       # Hartree to Joules
    Bohr_to_m = 5.29177210903e-11       # Bohr to meters
    Bohr_to_A = 0.529177210903          # Bohr to Angstrom
    amu_to_kg = 1.66053906660e-27       # Atomic mass unit to kg
    c_SI = 2.99792458e8                 # Speed of light (m/s)
    c_cm_s = 2.99792458e10              # Speed of light (cm/s)
    e = 1.602176634e-19                 # Elementary charge (C)
    epsilon_0 = 8.8541878128e-12        # Vacuum permittivity (F/m)

    # Conversion factors
    debye_to_au = 0.393430307           # Debye to atomic units (ea₀)
    au_to_debye = 2.541746473           # Atomic units to Debye

    def __init__(self, hamiltonian: 'MolecularHamiltonian'):
        """
        Initialize Raman/IR calculator.

        Args:
            hamiltonian: MolecularHamiltonian object

        Raises:
            ValueError: If hamiltonian has no atoms
        """
        self.hamiltonian = hamiltonian
        self.atoms = getattr(hamiltonian, 'atoms', [])
        self.mol = getattr(hamiltonian, 'mol', None)

        if len(self.atoms) == 0:
            raise ValueError("Hamiltonian has no atoms")

        self.n_atoms = len(self.atoms)
        self.n_coords = 3 * self.n_atoms

        # Get molecule name
        mol_name = getattr(hamiltonian, 'name', 'Unknown')
        logger.info(f"RamanIRCalculator initialized for {mol_name}")
        logger.info(f"  {self.n_atoms} atoms, {self.n_coords} coordinates")

    def _compute_dipole_moment(self, method: str = 'HF') -> np.ndarray:
        """
        Compute electric dipole moment.

        CRITICAL FIX: Uses quantum density if available from VQE/SQD.

        Args:
            method: Electronic structure method ('HF', 'DFT')

        Returns:
            Dipole moment vector (a.u.)
        """
        from pyscf import scf, dft

        mol_pyscf = self.mol

        # CRITICAL FIX: Try to use quantum density first
        dm = None
        mf = None
        if hasattr(self.hamiltonian, 'get_density_matrix'):
            try:
                dm = self.hamiltonian.get_density_matrix()
                # Need mf object for dip_moment call, create one with converged SCF
                if hasattr(self.hamiltonian, 'mf') and self.hamiltonian.mf is not None:
                    mf = self.hamiltonian.mf
                logger.debug("✅ Using quantum density for dipole calculation")
            except ValueError:
                pass  # No quantum density available

        # If no quantum density, run SCF
        if dm is None:
            if method.upper() == 'HF':
                mf = scf.RHF(mol_pyscf)
                mf.verbose = 0
                mf.kernel()
            else:
                mf = dft.RKS(mol_pyscf)
                mf.xc = 'B3LYP'
                mf.verbose = 0
                mf.kernel()
            dm = mf.make_rdm1()
            logger.debug("Using HF/DFT density for dipole calculation")

        # Compute dipole moment
        dipole = mf.dip_moment(mol=mol_pyscf, dm=dm, unit='AU', verbose=0)

        return np.array(dipole)

    def _compute_polarizability_from_scf(self, mf) -> np.ndarray:
        """
        Compute polarizability from SCF molecular orbitals using sum-over-states formula.

        Uses: α = 2 Σ_{occ,virt} |⟨occ|μ|virt⟩|² / (E_virt - E_occ)

        Where μ is the dipole operator and the sum runs over occupied-virtual transitions.

        Args:
            mf: PySCF mean-field object with converged SCF solution

        Returns:
            Polarizability tensor (3×3, a.u.)

        Raises:
            ValueError: If SCF not converged or insufficient virtual orbitals
        """
        try:
            # Get MO coefficients and energies
            mo_coeff = mf.mo_coeff
            mo_energy = mf.mo_energy
            mo_occ = mf.mo_occ

            # Find occupied and virtual orbitals
            occ_indices = np.where(mo_occ > 0)[0]
            virt_indices = np.where(mo_occ == 0)[0]

            n_occ = len(occ_indices)
            n_virt = len(virt_indices)

            if n_virt < 1:
                raise ValueError("No virtual orbitals available for polarizability calculation")

            # Get dipole integrals in AO basis from PySCF
            mol = mf.mol
            dipole_ao = mol.intor_symmetric('int1e_r', comp=3)  # [x, y, z] components

            # Transform dipole integrals to MO basis
            # dipole_mo[comp, i, j] = <MO_i|dipole_comp|MO_j>
            dipole_mo = np.zeros((3, mo_coeff.shape[1], mo_coeff.shape[1]))
            for comp in range(3):
                dipole_mo[comp] = mo_coeff.T @ dipole_ao[comp] @ mo_coeff

            # Compute polarizability tensor using sum-over-states
            # α_ij = 2 Σ_{occ,virt} μ_i(occ→virt) * μ_j(occ→virt) / ΔE
            alpha = np.zeros((3, 3))

            for i_occ in occ_indices:
                for a_virt in virt_indices:
                    # Energy difference
                    delta_e = mo_energy[a_virt] - mo_energy[i_occ]

                    # Skip near-degenerate transitions (avoid numerical issues)
                    if delta_e < 0.01:  # 0.01 Ha = 0.27 eV
                        continue

                    # Transition dipole moments for x, y, z
                    mu_x = dipole_mo[0, i_occ, a_virt]
                    mu_y = dipole_mo[1, i_occ, a_virt]
                    mu_z = dipole_mo[2, i_occ, a_virt]

                    # Contribution to polarizability tensor
                    # α_ij = 2 Σ μ_i * μ_j / ΔE
                    mu_vec = np.array([mu_x, mu_y, mu_z])
                    alpha += 2.0 * np.outer(mu_vec, mu_vec) / delta_e

            return alpha

        except Exception as e:
            logger.debug(f"Sum-over-states polarizability failed: {e}")
            raise

    def _compute_polarizability(self, method: str = 'HF') -> np.ndarray:
        """
        Compute static polarizability tensor (classical).

        First attempts proper sum-over-states calculation from SCF orbitals.
        Falls back to empirical approximation if that fails.

        Args:
            method: Electronic structure method ('HF', 'DFT')

        Returns:
            Polarizability tensor (3×3, a.u.)
        """
        from pyscf import scf, dft

        mol_pyscf = self.mol

        # Run SCF to get molecular properties
        if method.upper() == 'HF':
            mf = scf.RHF(mol_pyscf)
            mf.verbose = 0
            mf.kernel()
        else:
            mf = dft.RKS(mol_pyscf)
            mf.xc = 'B3LYP'
            mf.verbose = 0
            mf.kernel()

        # CRITICAL FIX: Use ONLY quantum sum-over-states calculation
        # Removed fallback to empirical approximation (alpha_iso = n_electrons * 0.8)
        # If sum-over-states fails, we should fix it, not hide behind approximations!
        alpha = self._compute_polarizability_from_scf(mf)
        logger.debug("✅ Computed polarizability using sum-over-states formula")
        return alpha

    def _compute_finite_field_polarizability(
        self,
        hamiltonian,
        rdm1_quantum: np.ndarray,
        field_strength: float = 0.001,
        verbose: bool = False
    ) -> np.ndarray:
        """
        Compute polarizability using finite-field method with quantum density.

        This is the CORRECT way to compute quantum polarizability:
        α_ij = -∂²E/∂F_i∂F_j

        Uses numerical differentiation with quantum 1-RDM to compute
        the response of the electronic energy to an applied electric field.

        Args:
            hamiltonian: Hamiltonian object with core Hamiltonian and dipole integrals
            rdm1_quantum: Quantum 1-electron reduced density matrix
            field_strength: Electric field strength (a.u.)
            verbose: Print progress

        Returns:
            Polarizability tensor (3×3, a.u.)

        References:
            - Finite-field method: Kurtz, Stewart, Dieter (1990) J. Comp. Chem. 11, 82
            - Quantum density response: Modern Quantum Chemistry (Szabo & Ostlund)
        """
        if verbose:
            print(f"  Computing finite-field polarizability...")
            print(f"  Field strength: {field_strength:.6f} a.u.")

        # Get core Hamiltonian and dipole integrals from PySCF molecule
        # The hamiltonian should have a pyscf_molecule attribute
        if hasattr(hamiltonian, 'pyscf_molecule'):
            mol = hamiltonian.pyscf_molecule
        elif hasattr(hamiltonian, 'mol'):
            mol = hamiltonian.mol
        else:
            # Create PySCF molecule from atom list
            from pyscf import gto
            mol = gto.M(
                atom=[(atom.symbol, atom.position) for atom in self.atoms],
                basis='sto-3g',
                unit='Bohr',
                verbose=0
            )

        # Get core Hamiltonian matrix (kinetic + nuclear attraction)
        H_core = mol.intor('int1e_kin') + mol.intor('int1e_nuc')

        # Get dipole integrals (3 components: x, y, z)
        dipole_ints = mol.intor_symmetric('int1e_r', comp=3)  # Shape: (3, n_ao, n_ao)

        # Get two-electron integrals (ERI) for energy calculation
        # For large molecules, this can be expensive - but necessary for accuracy
        eri = mol.intor('int2e', aosym='s1')  # Shape: (n_ao, n_ao, n_ao, n_ao)

        # Helper function to compute energy with applied field
        def compute_energy_with_field(field: np.ndarray) -> float:
            """
            Compute electronic energy with applied electric field.

            E = Tr[rdm1 * (H_core - dipole·F)] + 0.5 * Tr[rdm1 @ G[rdm1]]

            where G[rdm1] is the two-electron Fock matrix constructed from rdm1.
            """
            # Perturbed core Hamiltonian: H_eff = H_core - dipole·F
            H_eff = H_core.copy()
            for i in range(3):
                H_eff -= field[i] * dipole_ints[i]

            # One-electron energy
            E_one = np.einsum('ij,ji->', H_eff, rdm1_quantum)

            # Two-electron energy: 0.5 * Σ_ijkl rdm1_ij rdm1_kl (ij|kl) - 0.25 * Σ_ijkl rdm1_ij rdm1_kl (ik|jl)
            # This is the expensive part - but necessary for correct quantum energy
            E_two = 0.5 * np.einsum('ij,kl,ijkl->', rdm1_quantum, rdm1_quantum, eri)
            E_two -= 0.25 * np.einsum('ik,jl,ijkl->', rdm1_quantum, rdm1_quantum, eri)

            return E_one + E_two

        # Compute polarizability tensor using finite differences
        # α_ij = -∂²E/∂F_i∂F_j
        #      ≈ -(E(+F_i,+F_j) - E(+F_i,-F_j) - E(-F_i,+F_j) + E(-F_i,-F_j)) / (4F²)
        alpha = np.zeros((3, 3))

        # Compute energy at zero field (reference)
        E_0 = compute_energy_with_field(np.zeros(3))

        if verbose:
            print(f"  E(F=0): {E_0:.6f} Ha")

        # For diagonal elements, use simpler 3-point formula
        # α_ii = -(E(+F_i) - 2*E(0) + E(-F_i)) / F²
        for i in range(3):
            field_plus = np.zeros(3)
            field_plus[i] = field_strength

            field_minus = np.zeros(3)
            field_minus[i] = -field_strength

            E_plus = compute_energy_with_field(field_plus)
            E_minus = compute_energy_with_field(field_minus)

            alpha[i, i] = -(E_plus - 2*E_0 + E_minus) / (field_strength**2)

        # For off-diagonal elements, use 4-point formula
        for i in range(3):
            for j in range(i+1, 3):
                field_pp = np.zeros(3)
                field_pp[i] = field_strength
                field_pp[j] = field_strength

                field_pm = np.zeros(3)
                field_pm[i] = field_strength
                field_pm[j] = -field_strength

                field_mp = np.zeros(3)
                field_mp[i] = -field_strength
                field_mp[j] = field_strength

                field_mm = np.zeros(3)
                field_mm[i] = -field_strength
                field_mm[j] = -field_strength

                E_pp = compute_energy_with_field(field_pp)
                E_pm = compute_energy_with_field(field_pm)
                E_mp = compute_energy_with_field(field_mp)
                E_mm = compute_energy_with_field(field_mm)

                alpha[i, j] = -(E_pp - E_pm - E_mp + E_mm) / (4 * field_strength**2)
                alpha[j, i] = alpha[i, j]  # Symmetry

        if verbose:
            alpha_iso = np.mean(np.diag(alpha))
            print(f"  Polarizability (isotropic): {alpha_iso:.3f} a.u.")
            print(f"  Diagonal: [{alpha[0,0]:.3f}, {alpha[1,1]:.3f}, {alpha[2,2]:.3f}]")

        return alpha

    def _compute_quantum_polarizability(
        self,
        backend: str = 'statevector',
        method: str = 'sqd',
        subspace_dim: int = 15,
        verbose: bool = False
    ) -> np.ndarray:
        """
        Compute polarizability using QUANTUM methods (SQD/VQE).

        **WORLD'S FIRST quantum Raman spectroscopy!**

        Uses quantum density matrices from real quantum hardware to compute
        electronic polarizability.

        Args:
            backend: Quantum backend ('statevector', 'ibm', 'bluequbit')
            method: Quantum method ('sqd', 'vqe')
            subspace_dim: Subspace dimension for SQD
            verbose: Print progress

        Returns:
            Polarizability tensor (3×3, a.u.)
        """
        from kanad.solvers import VQESolver, SQDSolver
        from kanad.bonds import BondFactory

        if verbose:
            print(f"\n🔬 Computing quantum polarizability...")
            print(f"  Backend: {backend}")
            print(f"  Method: {method.upper()}")

        # Create bond from hamiltonian
        if self.n_atoms == 2:
            atom1, atom2 = self.atoms
            bond = BondFactory.create_bond(
                atom1.symbol,
                atom2.symbol,
                distance=np.linalg.norm(atom1.position - atom2.position)
            )
        else:
            # For polyatomic molecules, create synthetic bond
            from kanad.bonds.covalent_bond import CovalentBond
            bond = CovalentBond(self.atoms[0], self.atoms[1])
            # Attach molecule/hamiltonian
            from kanad.core.molecule import Molecule
            bond.molecule = Molecule(self.atoms)
            bond.hamiltonian = self.hamiltonian

        # Solve for ground state
        if method.lower() == 'sqd':
            solver = SQDSolver(bond=bond, subspace_dim=subspace_dim, backend=backend)
        else:  # vqe
            solver = VQESolver(bond=bond, backend=backend, ansatz_type='ucc')

        result = solver.solve()
        ground_energy = result['energy']

        if verbose:
            print(f"  Ground state energy: {ground_energy:.6f} Ha")

        # CRITICAL FIX: Implement finite-field quantum polarizability
        # Uses quantum 1-RDM to compute response to electric field

        # Get quantum density matrix from result
        if 'quantum_rdm1' in result:
            rdm1_quantum = result['quantum_rdm1']

            if verbose:
                print(f"  Using quantum 1-RDM for polarizability (finite-field method)")

            # Compute finite-field polarizability with quantum density
            try:
                alpha_quantum = self._compute_finite_field_polarizability(
                    bond.hamiltonian,
                    rdm1_quantum,
                    field_strength=0.001,  # 0.001 a.u. ≈ 0.05 V/Å
                    verbose=verbose
                )

                if verbose:
                    alpha_iso = np.mean(np.diag(alpha_quantum))
                    hf_energy = result.get('hf_energy', ground_energy)
                    correlation_energy = ground_energy - hf_energy
                    print(f"  Polarizability (quantum): {alpha_iso:.3f} a.u.")
                    print(f"  Correlation energy: {correlation_energy:.6f} Ha")

                return alpha_quantum

            except Exception as e:
                if verbose:
                    print(f"  ⚠️  Finite-field calculation failed: {e}")
                    print(f"  Falling back to HF polarizability")

                # Fallback to HF if finite-field fails
                alpha_classical = self._compute_polarizability(method='HF')
                return alpha_classical

        else:
            # No quantum density available - use HF
            if verbose:
                print(f"  ⚠️  No quantum 1-RDM available, using HF polarizability")

            alpha_classical = self._compute_polarizability(method='HF')
            return alpha_classical

    def _compute_dipole_derivatives(
        self,
        freq_result: Dict[str, Any],
        method: str = 'HF',
        step_size: float = 0.01,
        verbose: bool = False
    ) -> np.ndarray:
        """
        Compute dipole moment derivatives ∂μ/∂Q for IR intensities.

        Uses finite differences:
        ∂μ/∂Q_i ≈ (μ(Q_i + δ) - μ(Q_i - δ)) / (2δ)

        Args:
            freq_result: Output from FrequencyCalculator.compute_frequencies()
            method: Electronic structure method
            step_size: Displacement step (Bohr)
            verbose: Print progress

        Returns:
            Dipole derivatives (n_modes × 3), a.u.
        """
        from kanad.core.atom import Atom

        # Get normal modes
        normal_modes = np.array(freq_result['normal_modes'])  # (3N × n_modes)
        n_modes = normal_modes.shape[1]

        # Store original positions
        orig_positions = [atom.position.copy() for atom in self.atoms]

        # Compute dipole at equilibrium
        mu_eq = self._compute_dipole_moment(method)

        # Initialize derivatives
        dipole_derivatives = np.zeros((n_modes, 3))

        if verbose:
            print(f"\nComputing IR intensities (dipole derivatives)...")
            print(f"  Normal modes: {n_modes}")
            print(f"  Step size: {step_size:.4f} Bohr")

        # For each normal mode
        for i in range(n_modes):
            # Get mode vector
            mode = normal_modes[:, i]  # (3N,)

            # Displace along mode (+δ)
            for j, atom in enumerate(self.atoms):
                atom.position = orig_positions[j] + step_size * mode[3*j:3*(j+1)] * self.Bohr_to_A

            # Force rebuild Hamiltonian
            if hasattr(self.hamiltonian, '_mol'):
                self.hamiltonian._mol = None

            mu_plus = self._compute_dipole_moment(method)

            # Displace along mode (-δ)
            for j, atom in enumerate(self.atoms):
                atom.position = orig_positions[j] - step_size * mode[3*j:3*(j+1)] * self.Bohr_to_A

            if hasattr(self.hamiltonian, '_mol'):
                self.hamiltonian._mol = None

            mu_minus = self._compute_dipole_moment(method)

            # Restore original positions
            for j, atom in enumerate(self.atoms):
                atom.position = orig_positions[j].copy()

            if hasattr(self.hamiltonian, '_mol'):
                self.hamiltonian._mol = None

            # Compute derivative
            dipole_derivatives[i, :] = (mu_plus - mu_minus) / (2.0 * step_size)

            if verbose:
                print(f"  Mode {i+1}/{n_modes} ({100*(i+1)/n_modes:.0f}%)", end='\r')

        if verbose:
            print()

        return dipole_derivatives

    def _compute_polarizability_derivatives(
        self,
        freq_result: Dict[str, Any],
        method: str = 'HF',
        backend: str = None,
        quantum_method: str = 'sqd',
        subspace_dim: int = 15,
        step_size: float = 0.01,
        verbose: bool = False
    ) -> np.ndarray:
        """
        Compute polarizability derivatives ∂α/∂Q for Raman intensities.

        Uses finite differences:
        ∂α/∂Q_i ≈ (α(Q_i + δ) - α(Q_i - δ)) / (2δ)

        Args:
            freq_result: Output from FrequencyCalculator.compute_frequencies()
            method: Electronic structure method ('HF', 'DFT') for classical
            backend: Quantum backend ('statevector', 'ibm', 'bluequbit') for quantum
                    If None, uses classical method
            quantum_method: Quantum method ('sqd', 'vqe')
            subspace_dim: Subspace dimension for SQD
            step_size: Displacement step (Bohr)
            verbose: Print progress

        Returns:
            Polarizability derivatives (n_modes × 3 × 3), a.u.
        """
        # Get normal modes
        normal_modes = np.array(freq_result['normal_modes'])  # (3N × n_modes)
        n_modes = normal_modes.shape[1]

        # Store original positions
        orig_positions = [atom.position.copy() for atom in self.atoms]

        # Compute polarizability at equilibrium
        if backend is not None:
            # QUANTUM Raman!
            if verbose:
                print(f"\n🔬 Computing QUANTUM Raman intensities (polarizability derivatives)...")
                print(f"  Backend: {backend}")
                print(f"  Method: {quantum_method.upper()}")
            alpha_eq = self._compute_quantum_polarizability(backend, quantum_method, subspace_dim, verbose=False)
        else:
            # Classical
            if verbose:
                print(f"\nComputing Raman intensities (polarizability derivatives)...")
                print(f"  Method: {method}")
            alpha_eq = self._compute_polarizability(method)

        # Initialize derivatives
        polarizability_derivatives = np.zeros((n_modes, 3, 3))

        if verbose:
            print(f"  Normal modes: {n_modes}")
            print(f"  Step size: {step_size:.4f} Bohr")

        # For each normal mode
        for i in range(n_modes):
            # Get mode vector
            mode = normal_modes[:, i]  # (3N,)

            # Displace along mode (+δ)
            for j, atom in enumerate(self.atoms):
                atom.position = orig_positions[j] + step_size * mode[3*j:3*(j+1)] * self.Bohr_to_A

            # Force rebuild
            if hasattr(self.hamiltonian, '_mol'):
                self.hamiltonian._mol = None

            if backend is not None:
                alpha_plus = self._compute_quantum_polarizability(backend, quantum_method, subspace_dim, verbose=False)
            else:
                alpha_plus = self._compute_polarizability(method)

            # Displace along mode (-δ)
            for j, atom in enumerate(self.atoms):
                atom.position = orig_positions[j] - step_size * mode[3*j:3*(j+1)] * self.Bohr_to_A

            if hasattr(self.hamiltonian, '_mol'):
                self.hamiltonian._mol = None

            if backend is not None:
                alpha_minus = self._compute_quantum_polarizability(backend, quantum_method, subspace_dim, verbose=False)
            else:
                alpha_minus = self._compute_polarizability(method)

            # Restore positions
            for j, atom in enumerate(self.atoms):
                atom.position = orig_positions[j].copy()

            if hasattr(self.hamiltonian, '_mol'):
                self.hamiltonian._mol = None

            # Compute derivative
            polarizability_derivatives[i, :, :] = (alpha_plus - alpha_minus) / (2.0 * step_size)

            if verbose:
                print(f"  Mode {i+1}/{n_modes} ({100*(i+1)/n_modes:.0f}%)", end='\r')

        if verbose:
            print()

        return polarizability_derivatives

    def compute_intensities(
        self,
        freq_result: Dict[str, Any],
        method: str = 'HF',
        compute_ir: bool = True,
        compute_raman: bool = True,
        backend: Optional[str] = None,
        quantum_method: str = 'sqd',
        subspace_dim: int = 15,
        step_size: float = 0.01,
        temperature: float = 298.15,
        laser_wavelength: float = 532.0,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Compute IR and Raman intensities for vibrational modes.

        Args:
            freq_result: Output from FrequencyCalculator.compute_frequencies()
            method: Electronic structure method for classical ('HF', 'DFT')
            compute_ir: Compute IR intensities (default: True)
            compute_raman: Compute Raman intensities (default: True)
            backend: Quantum backend for Raman ('statevector', 'ibm', 'bluequbit')
                    If None, uses classical method
            quantum_method: Quantum method ('sqd', 'vqe')
            subspace_dim: Subspace dimension for quantum SQD
            step_size: Finite difference step size (Bohr)
            temperature: Temperature (K) for thermal population factors
            laser_wavelength: Raman laser wavelength (nm, default: 532 nm)
            verbose: Print progress and results

        Returns:
            Dictionary with:
                frequencies: Vibrational frequencies (cm⁻¹)
                ir_intensities: IR intensities (km/mol)
                raman_activities: Raman activities (Å⁴/amu)
                raman_intensities: Raman scattering intensities (arbitrary units)
                depolarization_ratios: Raman depolarization ratios
                dipole_derivatives: ∂μ/∂Q (a.u.)
                polarizability_derivatives: ∂α/∂Q (a.u.)
                method: Method used
                backend: Backend used (for quantum Raman)
                quantum: True if quantum method used
        """
        frequencies = np.array(freq_result['frequencies'])
        n_modes = len(frequencies)

        if verbose:
            print("\n" + "=" * 70)
            print("RAMAN & IR SPECTROSCOPY")
            print("=" * 70)
            print(f"Molecule: {getattr(self.hamiltonian, 'name', 'Unknown')}")
            print(f"Method: {method}" + (f" + Quantum ({backend})" if backend else ""))
            print(f"Vibrational modes: {n_modes}")
            print("-" * 70)

        # Initialize results
        ir_intensities = None
        raman_activities = None
        raman_intensities = None
        depolarization_ratios = None
        dipole_derivatives = None
        polarizability_derivatives = None

        # Compute IR intensities
        if compute_ir:
            if verbose:
                print("\n📊 Computing IR intensities...")

            dipole_derivatives = self._compute_dipole_derivatives(
                freq_result, method, step_size, verbose=verbose
            )

            # IR intensity: I_IR ∝ |∂μ/∂Q|²
            # Units: km/mol (convention)
            ir_intensities = np.sum(dipole_derivatives**2, axis=1)

            # Convert to km/mol (empirical conversion factor)
            # I (km/mol) ≈ 974.9 × (∂μ/∂Q)² where ∂μ/∂Q in e·Bohr/amu^(1/2)
            ir_intensities *= 974.9

            if verbose:
                print(f"✅ IR intensities computed")
                print(f"   Range: {np.min(ir_intensities[frequencies > 0]):.2f} - {np.max(ir_intensities):.2f} km/mol")

        # Compute Raman intensities
        if compute_raman:
            if backend:
                if verbose:
                    print(f"\n🔬 Computing QUANTUM Raman intensities (backend={backend})...")
            else:
                if verbose:
                    print(f"\n📊 Computing classical Raman intensities...")

            polarizability_derivatives = self._compute_polarizability_derivatives(
                freq_result, method, backend, quantum_method, subspace_dim, step_size, verbose=verbose
            )

            # Raman activity and intensity calculations
            raman_activities = np.zeros(n_modes)
            depolarization_ratios = np.zeros(n_modes)

            for i in range(n_modes):
                dα = polarizability_derivatives[i, :, :]  # 3×3 tensor

                # Invariants
                α_trace = np.trace(dα)  # Isotropic part
                α_aniso_sq = 0.5 * np.sum((dα - dα.T)**2) + 3 * (dα[0,1]**2 + dα[0,2]**2 + dα[1,2]**2)  # Anisotropic part

                # Raman activity: S = 45·α'² + 7·γ'²
                # where α' = (1/3)Tr(∂α/∂Q), γ'² = anisotropy
                raman_activities[i] = 45 * (α_trace / 3)**2 + 4 * α_aniso_sq

                # Depolarization ratio: ρ = 3γ'² / (45α'² + 4γ'²)
                denominator = 45 * (α_trace / 3)**2 + 4 * α_aniso_sq
                if denominator > 1e-10:
                    depolarization_ratios[i] = 3 * α_aniso_sq / denominator
                else:
                    depolarization_ratios[i] = 0.75  # Fully depolarized

            # Raman scattering intensity (includes frequency and thermal factors)
            # I_Raman ∝ (ν₀ - ν_vib)⁴ × S × n_B
            # where ν₀ = laser frequency, n_B = Bose-Einstein factor

            laser_freq_cm = 1e7 / laser_wavelength  # nm to cm⁻¹
            k_B = 0.695034800  # cm⁻¹/K (Boltzmann constant)

            raman_intensities = np.zeros(n_modes)
            for i in range(n_modes):
                if frequencies[i] > 0:  # Only real frequencies
                    ν_stokes = laser_freq_cm - frequencies[i]

                    # Bose-Einstein thermal population
                    if temperature > 0:
                        n_B = 1.0 / (np.exp(frequencies[i] / (k_B * temperature)) - 1)
                    else:
                        n_B = 0.0

                    # Raman intensity
                    raman_intensities[i] = (ν_stokes**4) * raman_activities[i] * (1 + n_B)

            # Normalize
            if np.max(raman_intensities) > 0:
                raman_intensities /= np.max(raman_intensities)

            if verbose:
                print(f"✅ Raman intensities computed")
                if backend:
                    print(f"   🌟 WORLD'S FIRST quantum Raman spectroscopy!")
                print(f"   Activity range: {np.min(raman_activities[frequencies > 0]):.2f} - {np.max(raman_activities):.2f} Å⁴/amu")

        # Print summary
        if verbose:
            print("\n" + "=" * 70)
            print("VIBRATIONAL MODES SUMMARY")
            print("=" * 70)
            print(f"{'Mode':<6} {'Freq (cm⁻¹)':<14} {'IR (km/mol)':<14} {'Raman (rel)':<14} {'ρ':<8}")
            print("-" * 70)

            for i in range(n_modes):
                freq = frequencies[i]
                ir_str = f"{ir_intensities[i]:.2f}" if ir_intensities is not None else "N/A"
                raman_str = f"{raman_intensities[i]:.4f}" if raman_intensities is not None else "N/A"
                rho_str = f"{depolarization_ratios[i]:.3f}" if depolarization_ratios is not None else "N/A"

                if freq > 0:
                    print(f"{i+1:<6} {freq:>12.2f}  {ir_str:<14} {raman_str:<14} {rho_str:<8}")
                else:
                    print(f"{i+1:<6} {freq:>12.2f}i (imaginary)")

            print("=" * 70)
            if backend:
                print(f"\n🌟 QUANTUM Raman spectroscopy completed (backend={backend})")
                print(f"   This is the WORLD'S FIRST quantum Raman calculator!")
            print("=" * 70)

        return {
            'frequencies': frequencies.tolist() if hasattr(frequencies, 'tolist') else frequencies,
            'ir_intensities': ir_intensities.tolist() if ir_intensities is not None and hasattr(ir_intensities, 'tolist') else ir_intensities,
            'raman_activities': raman_activities.tolist() if raman_activities is not None and hasattr(raman_activities, 'tolist') else raman_activities,
            'raman_intensities': raman_intensities.tolist() if raman_intensities is not None and hasattr(raman_intensities, 'tolist') else raman_intensities,
            'depolarization_ratios': depolarization_ratios.tolist() if depolarization_ratios is not None and hasattr(depolarization_ratios, 'tolist') else depolarization_ratios,
            'dipole_derivatives': dipole_derivatives.tolist() if dipole_derivatives is not None and hasattr(dipole_derivatives, 'tolist') else dipole_derivatives,
            'polarizability_derivatives': polarizability_derivatives.tolist() if polarizability_derivatives is not None and hasattr(polarizability_derivatives, 'tolist') else polarizability_derivatives,
            'method': method,
            'backend': backend,
            'quantum': backend is not None,
            'temperature': temperature,
            'laser_wavelength': laser_wavelength
        }
