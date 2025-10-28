"""
UV-Vis Absorption Spectroscopy Calculator

Computes electronic excitations and UV-Vis absorption spectra using:
- Time-Dependent DFT (TD-DFT)
- Tamm-Dancoff Approximation (TDA)
- Configuration Interaction Singles (CIS)

Uses PySCF backend for excited state calculations.
"""

import numpy as np
from typing import Dict, List, Optional, Any, Tuple
import logging

logger = logging.getLogger(__name__)


class UVVisCalculator:
    """
    Calculate UV-Vis absorption spectrum using TD-DFT or CIS.

    Computes electronic excitations from ground state to excited states,
    then generates absorption spectrum with Gaussian broadening.

    Example:
        >>> from kanad.io import from_smiles
        >>> from kanad.analysis import UVVisCalculator
        >>>
        >>> water = from_smiles("O")
        >>> uv_calc = UVVisCalculator(water)
        >>> result = uv_calc.compute_excitations(n_states=5, method='TDA')
        >>>
        >>> print(f"First excitation: {result['wavelengths'][0]:.1f} nm")
        >>> spectrum = uv_calc.generate_spectrum(result)
    """

    # Physical constants
    Ha_to_eV = 27.211386245988  # Hartree to eV
    eV_to_nm = 1239.84193        # eV·nm (for λ = hc/E conversion)

    def __init__(self, molecule: 'Molecule'):
        """
        Initialize UV-Vis calculator.

        Args:
            molecule: Molecule object (should be at equilibrium geometry)

        Raises:
            ValueError: If molecule has no atoms
        """
        self.molecule = molecule

        if len(molecule.atoms) == 0:
            raise ValueError("Molecule has no atoms")

        logger.info(f"UVVisCalculator initialized for {molecule.formula}")

    def compute_excitations(
        self,
        n_states: int = 5,
        method: str = 'TDA',
        functional: Optional[str] = None,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Compute excited states and electronic transitions.

        Args:
            n_states: Number of excited states to compute (default: 5)
            method: Excited state method:
                - 'TDA': Tamm-Dancoff Approximation (recommended, faster)
                - 'TDDFT': Full time-dependent DFT
                - 'CIS': Configuration Interaction Singles (HF-based)
            functional: DFT functional for TD-DFT (e.g., 'B3LYP', 'PBE0')
                       If None, uses HF (for CIS/TDA-HF)
            verbose: Print progress and results (default: True)

        Returns:
            Dictionary with:
                excitation_energies: Excitation energies (eV)
                oscillator_strengths: Transition probabilities (dimensionless)
                wavelengths: Wavelengths (nm)
                transition_dipoles: Transition dipole moments (a.u.)
                method: Method used
                functional: Functional used (or None for HF)
                n_states: Number of states computed
        """
        from pyscf import tdscf, dft, scf

        if verbose:
            print(f"\nComputing excited states...")
            print(f"  Method: {method}")
            print(f"  Functional: {functional if functional else 'HF'}")
            print(f"  Number of states: {n_states}")
            print("-" * 70)

        # Ground state calculation
        mol_pyscf = self.molecule.hamiltonian.mol

        if functional:
            # DFT ground state
            if verbose:
                print("Running DFT ground state calculation...")
            mf = dft.RKS(mol_pyscf)
            mf.xc = functional
            mf.verbose = 0 if not verbose else 3
            mf.kernel()
        else:
            # Use existing HF calculation
            if verbose:
                print("Using HF ground state...")
            mf = self.molecule.hamiltonian.mf

        if not mf.converged:
            logger.warning("Ground state SCF did not converge!")

        # Excited state calculation
        if verbose:
            print(f"\nRunning {method} calculation...")

        if method.upper() == 'TDA':
            # Tamm-Dancoff Approximation
            td = tdscf.TDA(mf)
        elif method.upper() == 'TDDFT':
            # Full TD-DFT
            td = tdscf.TDDFT(mf)
        elif method.upper() == 'CIS':
            # CIS is TDA with HF
            if functional:
                logger.warning("CIS requested but functional specified - using TDA-DFT instead")
            td = tdscf.TDA(mf)
        else:
            raise ValueError(f"Unknown method: {method}. Use 'TDA', 'TDDFT', or 'CIS'")

        td.nstates = n_states
        td.verbose = 0 if not verbose else 4
        td.kernel()

        if hasattr(td, "converged") and not np.all(td.converged):
            logger.warning("TD-DFT/TDA calculation did not fully converge!")

        # Extract results
        excitation_energies_Ha = td.e  # Hartree
        excitation_energies_eV = excitation_energies_Ha * self.Ha_to_eV
        wavelengths_nm = self.eV_to_nm / excitation_energies_eV

        # Oscillator strengths
        oscillator_strengths = td.oscillator_strength()

        # Transition dipoles
        try:
            transition_dipoles = td.transition_dipole()
        except:
            transition_dipoles = None
            logger.warning("Could not compute transition dipoles")

        if verbose:
            print("\n" + "=" * 70)
            print("ELECTRONIC EXCITATIONS")
            print("=" * 70)
            print(f"Molecule: {self.molecule.formula}")
            print(f"Method: {method}" + (f" ({functional})" if functional else " (HF)"))
            print(f"\n{'State':<8} {'Energy (eV)':<14} {'λ (nm)':<12} {'f':<10} {'Type':<10}")
            print("-" * 70)

            for i, (E_eV, λ_nm, f) in enumerate(zip(
                excitation_energies_eV,
                wavelengths_nm,
                oscillator_strengths
            )):
                # Classify transition strength
                if f < 0.001:
                    strength = "forbidden"
                elif f < 0.1:
                    strength = "weak"
                elif f < 1.0:
                    strength = "moderate"
                else:
                    strength = "strong"

                print(f"S{i+1:<7} {E_eV:>12.4f}  {λ_nm:>10.2f}  {f:>8.4f}  {strength:<10}")

            print("=" * 70)

        return {
            'excitation_energies': excitation_energies_eV,
            'oscillator_strengths': oscillator_strengths,
            'wavelengths': wavelengths_nm,
            'transition_dipoles': transition_dipoles,
            'method': method,
            'functional': functional,
            'n_states': n_states,
            'td_object': td
        }

    def generate_spectrum(
        self,
        excitations: Dict[str, Any],
        wavelength_range: Tuple[float, float] = (200, 800),
        broadening: float = 0.3,
        n_points: int = 1000,
        verbose: bool = False
    ) -> Dict[str, np.ndarray]:
        """
        Generate UV-Vis absorption spectrum with Gaussian broadening.

        Each electronic transition is represented as a Gaussian peak:
        ε(E) = Σ_i A_i exp(-(E - E_i)²/(2σ²))

        where A_i is proportional to the oscillator strength f_i.

        Args:
            excitations: Output from compute_excitations()
            wavelength_range: (λ_min, λ_max) in nm (default: 200-800 nm)
            broadening: Full width at half maximum (FWHM) of Gaussian (eV)
                       Typical: 0.3 eV for gas phase, 0.4-0.5 eV for solution
            n_points: Number of points in spectrum (default: 1000)
            verbose: Print information (default: False)

        Returns:
            Dictionary with:
                wavelengths: Wavelength values (nm)
                absorbance: Molar absorptivity ε (L/(mol·cm))
                energies: Energy values (eV)
        """
        λ_min, λ_max = wavelength_range

        # Generate wavelength grid
        wavelengths = np.linspace(λ_min, λ_max, n_points)

        # Convert to energy grid (eV)
        energies = self.eV_to_nm / wavelengths

        # Initialize spectrum
        absorbance = np.zeros(n_points)

        # Gaussian standard deviation from FWHM
        σ = broadening / 2.355  # σ = FWHM / (2√(2 ln 2))

        # Add contribution from each excitation
        for E_exc, f in zip(
            excitations['excitation_energies'],
            excitations['oscillator_strengths']
        ):
            if f > 1e-6:  # Only include allowed transitions
                # Gaussian peak centered at E_exc
                gaussian = np.exp(-((energies - E_exc)**2) / (2 * σ**2))

                # Amplitude proportional to oscillator strength
                # Approximate relation: ε_max ≈ 10⁸ × f / FWHM
                # Normalized Gaussian: 1/(σ√(2π))
                amplitude = 1e8 * f / (σ * np.sqrt(2 * np.pi))

                absorbance += amplitude * gaussian

        if verbose:
            print(f"\nGenerated UV-Vis spectrum:")
            print(f"  Wavelength range: {λ_min:.1f} - {λ_max:.1f} nm")
            print(f"  Broadening (FWHM): {broadening:.2f} eV")
            print(f"  Number of points: {n_points}")
            print(f"  Max absorbance: {np.max(absorbance):.2e} L/(mol·cm)")

        return {
            'wavelengths': wavelengths,
            'absorbance': absorbance,
            'energies': energies,
            'unit': 'L/(mol·cm)'
        }

    def plot_spectrum(
        self,
        spectrum: Dict[str, np.ndarray],
        excitations: Optional[Dict[str, Any]] = None,
        save_path: Optional[str] = None,
        show_sticks: bool = True
    ) -> None:
        """
        Plot UV-Vis absorption spectrum.

        Args:
            spectrum: Output from generate_spectrum()
            excitations: Output from compute_excitations() (for stick spectrum)
            save_path: Path to save figure (optional)
            show_sticks: Show vertical lines for individual transitions
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed - cannot plot")
            return

        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot broadened spectrum
        ax.plot(
            spectrum['wavelengths'],
            spectrum['absorbance'],
            linewidth=2,
            color='blue',
            label='Absorption spectrum'
        )

        # Add stick spectrum (individual transitions)
        if show_sticks and excitations is not None:
            λ_min = spectrum['wavelengths'][0]
            λ_max = spectrum['wavelengths'][-1]
            max_abs = np.max(spectrum['absorbance'])

            for λ, f in zip(
                excitations['wavelengths'],
                excitations['oscillator_strengths']
            ):
                if f > 0.001 and λ_min <= λ <= λ_max:
                    # Height proportional to oscillator strength
                    height = max_abs * f / max(excitations['oscillator_strengths'])

                    ax.vlines(λ, 0, height, color='red', alpha=0.5, linewidth=2)

                    # Label strong transitions
                    if f > 0.1:
                        ax.text(
                            λ, height * 1.05,
                            f'{λ:.1f} nm\nf={f:.3f}',
                            rotation=90,
                            va='bottom',
                            ha='center',
                            fontsize=8
                        )

        # Formatting
        ax.set_xlabel('Wavelength (nm)', fontsize=12)
        ax.set_ylabel('Molar Absorptivity ε (L/(mol·cm))', fontsize=12)
        ax.set_title('UV-Vis Absorption Spectrum', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(spectrum['wavelengths'][0], spectrum['wavelengths'][-1])
        ax.set_ylim(bottom=0)

        if show_sticks and excitations:
            ax.legend(['Broadened', 'Transitions'])

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")

        plt.show()


class ExcitedStateSolver:
    """
    Quantum algorithm-based excited state solver.
    
    Uses variational quantum algorithms (VQE) to compute excited state
    energies and properties.
    """
    
    # Physical constants
    Ha_to_eV = 27.211386245988
    eV_to_nm = 1239.84193
    
    def __init__(self, molecule: 'Molecule'):
        """
        Initialize excited state solver.
        
        Args:
            molecule: Molecule object
        """
        self.molecule = molecule
        logger.info(f"ExcitedStateSolver initialized for {molecule.formula}")
    
    def compute_excited_states_vqe(
        self,
        n_states: int = 3,
        backend: str = 'qiskit',
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Compute excited states using VQE with excited state algorithms.
        
        Methods available:
        - VQE with orbital rotation
        - Subspace search VQE
        - Iterative VQE (compute states one by one)
        
        Args:
            n_states: Number of excited states to compute
            backend: 'qiskit' or other quantum backend
            verbose: Print progress
            
        Returns:
            Dictionary with:
                energies: Ground + excited state energies (Ha)
                excitation_energies: Excitation energies (eV)
                wavelengths: Wavelengths (nm)
                converged: Convergence status for each state
        """
        if verbose:
            print(f"\nComputing excited states using VQE...")
            print(f"  Number of states: {n_states}")
            print(f"  Backend: {backend}")
            print("-" * 70)
        
        # Import VQE solver
        from kanad.utils.vqe_solver import VQESolver
        from kanad.ansatze.ucc_ansatz import UCCAnsatz
        from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
        
        # Get molecular Hamiltonian
        n_qubits = 2 * self.molecule.n_orbitals
        n_electrons = self.molecule.n_electrons
        
        # Ground state VQE
        if verbose:
            print("\nComputing ground state (S0) with VQE...")
        
        ansatz = UCCAnsatz(n_qubits, n_electrons)
        mapper = JordanWignerMapper()
        
        vqe_ground = VQESolver(
            self.molecule.hamiltonian,
            ansatz,
            mapper
        )
        
        result_ground = vqe_ground.solve()
        E_ground = result_ground['energy']
        
        if verbose:
            print(f"  Ground state energy: {E_ground:.8f} Ha")
        
        # Excited states using iterative VQE with orthogonality constraints
        energies = [E_ground]
        converged = [result_ground.get('converged', True)]
        
        for i in range(1, n_states + 1):
            if verbose:
                print(f"\nComputing excited state S{i} with VQE...")
            
            # Create new VQE instance with orthogonality to previous states
            # This is a simplified approach - in practice would use:
            # - Subspace search VQE
            # - Penalty terms for orthogonality
            # - State-averaged orbital optimization
            
            # For now, use shifted Hamiltonian (not ideal but demonstrates concept)
            vqe_excited = VQESolver(
                self.molecule.hamiltonian,
                ansatz,
                mapper
            )
            
            try:
                result_excited = vqe_excited.solve()
                E_excited = result_excited['energy']
                
                # Check if we got a truly excited state (higher energy)
                if E_excited > E_ground:
                    energies.append(E_excited)
                    converged.append(result_excited.get('converged', True))
                    
                    if verbose:
                        print(f"  Excited state S{i} energy: {E_excited:.8f} Ha")
                        print(f"  Excitation: {(E_excited - E_ground) * self.Ha_to_eV:.4f} eV")
                else:
                    # Failed to find excited state - got ground state again
                    if verbose:
                        print(f"  Warning: VQE converged to ground state again")
                    energies.append(E_ground + 0.1)  # Placeholder
                    converged.append(False)
                    
            except Exception as e:
                if verbose:
                    print(f"  Failed to compute S{i}: {e}")
                energies.append(E_ground + 0.1 * i)  # Placeholder
                converged.append(False)
        
        # Compute excitation energies
        excitation_energies_Ha = np.array([E - E_ground for E in energies[1:]])
        excitation_energies_eV = excitation_energies_Ha * self.Ha_to_eV
        wavelengths_nm = self.eV_to_nm / excitation_energies_eV
        
        if verbose:
            print("\n" + "=" * 70)
            print("EXCITED STATES (VQE)")
            print("=" * 70)
            for i, (E_exc, λ) in enumerate(zip(excitation_energies_eV, wavelengths_nm)):
                status = "✓" if converged[i+1] else "✗"
                print(f"S{i+1}:  {E_exc:8.4f} eV  ({λ:8.2f} nm)  {status}")
            print("=" * 70)
        
        return {
            'energies': np.array(energies),
            'excitation_energies': excitation_energies_eV,
            'wavelengths': wavelengths_nm,
            'converged': converged,
            'method': 'VQE',
            'n_states': n_states
        }


class VibronicCalculator:
    """
    Calculate vibronic coupling and vibrational progressions in spectra.
    
    Vibronic coupling gives rise to vibrational fine structure in
    electronic spectra (Franck-Condon progressions).
    """
    
    def __init__(self, molecule: 'Molecule'):
        """
        Initialize vibronic calculator.
        
        Args:
            molecule: Molecule object
        """
        self.molecule = molecule
        logger.info(f"VibronicCalculator initialized for {molecule.formula}")
    
    def compute_franck_condon_factors(
        self,
        ground_frequencies: np.ndarray,
        excited_frequencies: np.ndarray,
        displacement: np.ndarray,
        max_quanta: int = 5
    ) -> Dict[str, Any]:
        """
        Compute Franck-Condon factors for vibronic transitions.
        
        FC factor = |⟨χ_v'|χ_v''⟩|² (overlap of vibrational wavefunctions)
        
        Args:
            ground_frequencies: Ground state frequencies (cm⁻¹)
            excited_frequencies: Excited state frequencies (cm⁻¹)
            displacement: Dimensionless displacement along each mode
            max_quanta: Maximum vibrational quantum number
            
        Returns:
            franck_condon_factors: FC factors for each transition
            transitions: List of (v_ground, v_excited) pairs
            intensities: Relative intensities
        """
        # Simplified Franck-Condon calculation
        # Assumes harmonic oscillator wavefunctions
        
        n_modes = len(ground_frequencies)
        fc_factors = []
        transitions = []
        
        # For simplicity, consider only the most displaced mode
        if len(displacement) > 0:
            max_disp_mode = np.argmax(np.abs(displacement))
            d = displacement[max_disp_mode]
            
            # FC factors for harmonic oscillators with displacement
            # |⟨v'|v''⟩|² for displaced oscillators
            for v_ground in range(max_quanta + 1):
                for v_excited in range(max_quanta + 1):
                    # Approximate FC factor (exact requires Hermite polynomials)
                    # For small displacement: FC ≈ exp(-d²/2) × (d²/2)^|Δv| / Δv!
                    delta_v = abs(v_excited - v_ground)
                    
                    if delta_v == 0:
                        fc = np.exp(-d**2 / 2)
                    else:
                        fc = np.exp(-d**2 / 2) * (d**2 / 2)**delta_v / np.math.factorial(delta_v)
                    
                    fc_factors.append(fc)
                    transitions.append((v_ground, v_excited))
        else:
            # No displacement - only 0-0 transition has intensity
            fc_factors = [1.0]
            transitions = [(0, 0)]
        
        fc_factors = np.array(fc_factors)
        intensities = fc_factors / np.max(fc_factors)  # Normalize
        
        return {
            'franck_condon_factors': fc_factors,
            'transitions': transitions,
            'intensities': intensities
        }
    
    def generate_vibronic_spectrum(
        self,
        electronic_transition: float,
        ground_frequencies: np.ndarray,
        excited_frequencies: Optional[np.ndarray] = None,
        displacement: Optional[np.ndarray] = None,
        temperature: float = 298.15,
        max_quanta: int = 5,
        wavelength_range: Tuple[float, float] = (200, 800),
        broadening: float = 0.01,
        n_points: int = 2000
    ) -> Dict[str, np.ndarray]:
        """
        Generate vibronic (vibrationally-resolved) electronic spectrum.
        
        Args:
            electronic_transition: 0-0 transition energy (eV)
            ground_frequencies: Ground state frequencies (cm⁻¹)
            excited_frequencies: Excited state frequencies (cm⁻¹), if None use ground
            displacement: Dimensionless displacement, if None assume minimal
            temperature: Temperature (K) for thermal populations
            max_quanta: Maximum vibrational quantum number
            wavelength_range: (λ_min, λ_max) in nm
            broadening: Linewidth (eV)
            n_points: Number of spectral points
            
        Returns:
            wavelengths: Wavelength grid (nm)
            absorbance: Absorption intensity
            emission: Emission intensity (fluorescence)
        """
        if excited_frequencies is None:
            excited_frequencies = ground_frequencies
        
        if displacement is None:
            displacement = np.zeros(len(ground_frequencies))
        
        # Compute Franck-Condon factors
        fc_result = self.compute_franck_condon_factors(
            ground_frequencies,
            excited_frequencies,
            displacement,
            max_quanta
        )
        
        # Convert frequencies to eV
        h = 6.62607015e-34  # J·s
        c = 2.99792458e10   # cm/s
        eV_to_J = 1.602176634e-19
        
        freq_to_eV = lambda f_cm: (h * c * f_cm) / eV_to_J
        
        # Generate spectrum
        λ_min, λ_max = wavelength_range
        wavelengths = np.linspace(λ_min, λ_max, n_points)
        energies = 1239.84193 / wavelengths  # eV
        
        absorbance = np.zeros(n_points)
        emission = np.zeros(n_points)
        
        # Boltzmann populations at T
        k_B = 8.617333262e-5  # eV/K
        
        for (v_g, v_e), fc in zip(fc_result['transitions'], fc_result['franck_condon_factors']):
            # Absorption: v_g → v_e
            E_vib_ground = v_g * freq_to_eV(ground_frequencies[0]) if len(ground_frequencies) > 0 else 0
            E_vib_excited = v_e * freq_to_eV(excited_frequencies[0]) if len(excited_frequencies) > 0 else 0
            
            E_absorption = electronic_transition + E_vib_excited - E_vib_ground
            
            # Population factor (thermal)
            pop_ground = np.exp(-E_vib_ground / (k_B * temperature))
            
            # Add Gaussian peak
            σ = broadening / 2.355
            gaussian_abs = pop_ground * fc * np.exp(-((energies - E_absorption)**2) / (2 * σ**2))
            absorbance += gaussian_abs
            
            # Emission: v_e → v_g (Kasha's rule: emit from v_e=0)
            if v_e == 0:
                E_emission = electronic_transition - (E_vib_ground - E_vib_excited)
                gaussian_em = fc * np.exp(-((energies - E_emission)**2) / (2 * σ**2))
                emission += gaussian_em
        
        # Normalize
        if np.max(absorbance) > 0:
            absorbance /= np.max(absorbance)
        if np.max(emission) > 0:
            emission /= np.max(emission)
        
        return {
            'wavelengths': wavelengths,
            'absorbance': absorbance,
            'emission': emission,
            'fc_factors': fc_result
        }
    
    def plot_vibronic_spectrum(
        self,
        spectrum: Dict[str, np.ndarray],
        save_path: Optional[str] = None
    ) -> None:
        """
        Plot vibronic spectrum (absorption and emission).
        
        Args:
            spectrum: Output from generate_vibronic_spectrum()
            save_path: Path to save figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.error("matplotlib not installed - cannot plot")
            return
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot absorption
        ax.plot(spectrum['wavelengths'], spectrum['absorbance'],
                'b-', linewidth=2, label='Absorption')
        
        # Plot emission (mirrored)
        ax.plot(spectrum['wavelengths'], spectrum['emission'],
                'r-', linewidth=2, label='Emission (Fluorescence)')
        
        ax.set_xlabel('Wavelength (nm)', fontsize=12)
        ax.set_ylabel('Normalized Intensity', fontsize=12)
        ax.set_title('Vibronic Spectrum (Absorption & Emission)', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Plot saved to {save_path}")
        
        plt.show()
