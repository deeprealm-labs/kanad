"""
AnalysisService - Decoupled Analysis System for Kanad

This service provides on-demand analysis capabilities that work independently
of the solver execution. It enables:
- Domain-specific analysis profiles (Chemistry, Spectroscopy, Materials, etc.)
- Custom analysis configuration
- Async analysis execution with progress tracking
- Analysis result caching and management

Architecture:
- Loads experiment results from database
- Runs analysis modules independently
- Stores results in analysis_results table
- Supports WebSocket progress streaming
"""

import json
import time
import logging
from typing import Dict, List, Optional, Any, Callable
from datetime import datetime
import numpy as np

# Import all analysis modules
from kanad.analysis import (
    EnergyAnalyzer,
    BondingAnalyzer,
    CorrelationAnalyzer,
    PropertyCalculator,
    BondLengthScanner,
    ThermochemistryCalculator,
    FrequencyCalculator,
    UVVisCalculator,
    ExcitedStateSolver,
    VibronicCalculator,
    DOSCalculator,
    UncertaintyAnalyzer,
)
from kanad.analysis.adme_calculator import ADMECalculator

# Import profiles (Chemistry removed - all analyses already done inline)
from kanad.services.profiles.spectroscopy import SPECTROSCOPY_PROFILE
from kanad.services.profiles.materials import MATERIALS_PROFILE
from kanad.services.profiles.drug_discovery import DRUG_DISCOVERY_PROFILE
from kanad.services.profiles.catalysis import CATALYSIS_PROFILE

logger = logging.getLogger(__name__)


def make_json_serializable(obj):
    """
    Recursively convert numpy arrays and other non-JSON types to JSON-serializable formats.

    Args:
        obj: Object to convert (dict, list, numpy array, etc.)

    Returns:
        JSON-serializable version of the object
    """
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {key: make_json_serializable(value) for key, value in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [make_json_serializable(item) for item in obj]
    elif isinstance(obj, (np.integer, np.floating)):
        return float(obj)
    elif isinstance(obj, (int, float, str, bool, type(None))):
        return obj
    else:
        # For complex objects that can't be serialized, convert to string representation
        return str(obj)


class AnalysisService:
    """
    Decoupled analysis service that operates on completed experiment results.

    This service enables on-demand analysis without re-running expensive
    quantum computations. It supports:
    - Domain-specific profiles
    - Custom analysis selection
    - Async execution with progress tracking
    - Result caching
    """

    # Registry of available analysis modules
    ANALYSIS_MODULES = {
        'energy': EnergyAnalyzer,
        'bonding': BondingAnalyzer,
        'correlation': CorrelationAnalyzer,
        'properties': PropertyCalculator,
        'bond_scan': BondLengthScanner,
        'thermochemistry': ThermochemistryCalculator,
        'adme': ADMECalculator,  # ADME properties for drug discovery
        'nmr': None,  # NMR with Google ECHOES (to be implemented)
        'frequencies': FrequencyCalculator,
        'uv_vis': UVVisCalculator,
        'excited_states': ExcitedStateSolver,
        'vibronic': VibronicCalculator,
        'dos': DOSCalculator,
        'uncertainty': UncertaintyAnalyzer,
    }

    # Registry of domain profiles (only UNIQUE analyses, not already in inline analysis)
    DOMAIN_PROFILES = {
        'drug_discovery': DRUG_DISCOVERY_PROFILE,     # ADME properties
        'spectroscopy': SPECTROSCOPY_PROFILE,         # UV-Vis, frequencies, excited states
        'materials': MATERIALS_PROFILE,               # Density of states
        'catalysis': CATALYSIS_PROFILE,               # Transition state frequencies
    }

    def __init__(
        self,
        experiment_results: Dict[str, Any],
        molecule_data: Optional[Dict[str, Any]] = None,
        progress_callback: Optional[Callable] = None
    ):
        """
        Initialize AnalysisService with experiment results.

        Args:
            experiment_results: Dictionary containing solver results
                Required keys: 'energies', 'final_energy', 'optimal_params'
                Optional: 'rdm1', 'rdm2', 'hamiltonian', 'wavefunction'
            molecule_data: Dictionary with molecule info (atoms, geometry, charge, spin)
            progress_callback: Optional callback for progress updates
        """
        self.results = experiment_results
        self.molecule = molecule_data or {}
        self.progress_callback = progress_callback

        # Extract key data from results
        # Support both 'final_energy' (old format) and 'energy' (new format)
        self.final_energy = experiment_results.get('final_energy') or experiment_results.get('energy')
        self.energies = experiment_results.get('energies', [])
        self.optimal_params = experiment_results.get('optimal_params', [])

        # Optional advanced data
        self.rdm1 = experiment_results.get('rdm1')
        self.rdm2 = experiment_results.get('rdm2')
        self.hamiltonian = experiment_results.get('hamiltonian')
        self.wavefunction = experiment_results.get('wavefunction')

        logger.info(f"AnalysisService initialized with energy={self.final_energy}")

    def run_analysis_profile(
        self,
        profile_name: str,
        parameters: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Run a predefined domain-specific analysis profile.

        Args:
            profile_name: Name of profile ('chemistry', 'spectroscopy', etc.)
            parameters: Optional profile-specific parameters

        Returns:
            Dictionary with analysis results for all modules in the profile
        """
        if profile_name not in self.DOMAIN_PROFILES:
            available = ', '.join(self.DOMAIN_PROFILES.keys())
            raise ValueError(f"Unknown profile '{profile_name}'. Available: {available}")

        profile = self.DOMAIN_PROFILES[profile_name]
        logger.info(f"Running {profile_name} analysis profile ({len(profile['analyses'])} analyses)")

        # Merge default profile params with user params
        params = profile.get('default_parameters', {})
        if parameters:
            params.update(parameters)

        # Run each analysis in the profile
        results = {
            'profile': profile_name,
            'description': profile['description'],
            'timestamp': datetime.utcnow().isoformat(),
            'analyses': {}
        }

        analyses = profile['analyses']
        for i, analysis_name in enumerate(analyses):
            if self.progress_callback:
                progress = (i / len(analyses)) * 100
                self.progress_callback(f"Running {analysis_name}", progress)

            try:
                # Get analysis-specific parameters
                analysis_params = params.get(analysis_name, {})
                result = self.run_analysis(analysis_name, analysis_params)
                results['analyses'][analysis_name] = result

                # Store result in self.results so subsequent analyses can access it
                # (e.g., thermochemistry needs frequencies)
                self.results[analysis_name] = result

                logger.info(f"✓ {analysis_name} completed")
            except Exception as e:
                logger.error(f"✗ {analysis_name} failed: {e}")
                results['analyses'][analysis_name] = {
                    'error': str(e),
                    'status': 'failed'
                }

        if self.progress_callback:
            self.progress_callback("Analysis profile complete", 100)

        # Convert numpy arrays and other non-JSON types to JSON-serializable formats
        return make_json_serializable(results)

    def run_analysis(
        self,
        analysis_name: str,
        parameters: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Run a single analysis module.

        Args:
            analysis_name: Name of analysis module
            parameters: Analysis-specific parameters

        Returns:
            Dictionary with analysis results
        """
        if analysis_name not in self.ANALYSIS_MODULES:
            available = ', '.join(self.ANALYSIS_MODULES.keys())
            raise ValueError(f"Unknown analysis '{analysis_name}'. Available: {available}")

        params = parameters or {}
        analyzer_class = self.ANALYSIS_MODULES[analysis_name]

        start_time = time.time()

        try:
            # Route to appropriate analysis method
            if analysis_name == 'energy':
                result = self._run_energy_analysis(params)
            elif analysis_name == 'bonding':
                result = self._run_bonding_analysis(params)
            elif analysis_name == 'correlation':
                result = self._run_correlation_analysis(params)
            elif analysis_name == 'properties':
                result = self._run_properties_analysis(params)
            elif analysis_name == 'adme':
                result = self._run_adme_analysis(params)
            elif analysis_name == 'thermochemistry':
                result = self._run_thermochemistry_analysis(params)
            elif analysis_name == 'frequencies':
                result = self._run_frequencies_analysis(params)
            elif analysis_name == 'uv_vis':
                result = self._run_uvvis_analysis(params)
            elif analysis_name == 'excited_states':
                result = self._run_excited_states_analysis(params)
            elif analysis_name == 'vibronic':
                result = self._run_vibronic_analysis(params)
            elif analysis_name == 'dos':
                result = self._run_dos_analysis(params)
            elif analysis_name == 'uncertainty':
                result = self._run_uncertainty_analysis(params)
            else:
                raise NotImplementedError(f"Analysis '{analysis_name}' not yet implemented")

            computation_time = time.time() - start_time
            result['computation_time'] = round(computation_time, 3)
            result['status'] = 'success'

            return result

        except Exception as e:
            logger.error(f"Analysis {analysis_name} failed: {e}")
            return {
                'status': 'failed',
                'error': str(e),
                'computation_time': time.time() - start_time
            }

    def _run_energy_analysis(self, params: Dict) -> Dict:
        """Run energy decomposition analysis."""
        # Simple energy analysis without requiring hamiltonian object
        nuclear_repulsion = self.results.get('nuclear_repulsion', 0.0)
        total_energy = self.final_energy
        electronic_energy = total_energy - nuclear_repulsion
        n_electrons = self.molecule.get('n_electrons', 2)

        result = {
            'total_energy': total_energy,
            'nuclear_repulsion': nuclear_repulsion,
            'electronic_energy': electronic_energy,
            'energy_per_electron': electronic_energy / n_electrons if n_electrons > 0 else None,
        }

        # Add convergence analysis if energy history available
        if self.energies and len(self.energies) > 1:
            result['convergence'] = {
                'initial_energy': float(self.energies[0]),
                'final_energy': float(self.energies[-1]),
                'energy_change': float(self.energies[-1] - self.energies[0]),
                'iterations': len(self.energies),
            }

        return result

    def _run_bonding_analysis(self, params: Dict) -> Dict:
        """Run bonding and molecular orbital analysis."""
        if not self.rdm1:
            return {'error': 'RDM1 (density matrix) required for bonding analysis'}

        # TODO: Implement bonding analysis with RDM1
        # This requires molecular orbital data
        return {
            'message': 'Bonding analysis requires MO coefficients',
            'available_data': 'rdm1' if self.rdm1 else None
        }

    def _run_correlation_analysis(self, params: Dict) -> Dict:
        """Run electron correlation analysis."""
        hf_energy = self.results.get('hf_energy')
        if not hf_energy:
            return {'error': 'Hartree-Fock reference energy required'}

        # Simple correlation analysis without requiring hamiltonian object
        vqe_energy = self.final_energy
        correlation_energy = vqe_energy - hf_energy

        result = {
            'vqe_energy': vqe_energy,
            'hf_energy': hf_energy,
            'correlation_energy': correlation_energy,
        }

        # Calculate correlation percentage
        if abs(hf_energy) > 1e-10:
            result['correlation_percentage'] = abs(correlation_energy / hf_energy) * 100

        # Compare with FCI if available
        fci_energy = self.results.get('fci_energy')
        if fci_energy:
            result['fci_energy'] = fci_energy
            result['error_vs_fci'] = vqe_energy - fci_energy
            total_correlation = fci_energy - hf_energy
            if abs(total_correlation) > 1e-10:
                result['correlation_recovered'] = abs(correlation_energy / total_correlation) * 100

        return result

    def _run_properties_analysis(self, params: Dict) -> Dict:
        """Calculate molecular properties (dipole, polarizability)."""
        if not self.rdm1:
            return {'error': 'RDM1 required for property calculation'}

        # TODO: Implement properties calculation with RDM1
        return {'message': 'Properties calculation requires implementation'}

    def _run_adme_analysis(self, params: Dict) -> Dict:
        """Calculate ADME properties for drug discovery."""
        if not self.molecule.get('geometry'):
            return {'error': 'Molecular geometry required for ADME calculation'}

        try:
            # Create ADME calculator
            calc = ADMECalculator(
                geometry=self.molecule['geometry'],
                smiles=self.molecule.get('smiles'),
                charge=self.molecule.get('charge', 0),
                multiplicity=self.molecule.get('multiplicity', 1)
            )

            # Calculate molecular descriptors
            descriptors = calc.calculate_descriptors(
                rdm1=np.array(self.results.get('rdm1')) if self.results.get('rdm1') else None,
                orbital_energies=self.results.get('orbital_energies'),
                dipole=np.array(self.results.get('dipole')) if self.results.get('dipole') else None,
                polarizability=self.results.get('polarizability')
            )

            # Predict ADME properties
            adme_props = calc.predict_adme(descriptors)

            # Convert to dict using dataclass fields
            from dataclasses import asdict

            result = {
                'descriptors': asdict(descriptors),
                'adme_properties': asdict(adme_props),
            }

            return result

        except Exception as e:
            import traceback
            logger.error(f"ADME calculation failed: {e}")
            logger.error(traceback.format_exc())
            return {'error': f'ADME calculation failed: {str(e)}'}

    def _run_thermochemistry_analysis(self, params: Dict) -> Dict:
        """Calculate thermochemical properties."""
        temperature = params.get('temperature', 298.15)
        pressure = params.get('pressure', 101325.0)

        if not self.molecule.get('geometry'):
            return {'error': 'Molecular geometry required for thermochemistry'}

        # Check if we have vibrational frequencies
        if not self.results.get('frequencies'):
            return {'error': 'Vibrational frequencies required. Run frequencies analysis first.'}

        try:
            # Reconstruct molecule object
            from kanad.core.molecule import Molecule

            # Get atoms from molecule_data
            atoms = self.molecule.get('atoms', [])
            if not atoms:
                return {'error': 'Atom data not available'}

            # Create molecule with existing Atom objects
            mol = Molecule(
                atoms=atoms,
                charge=self.molecule.get('charge', 0),
                spin=self.molecule.get('multiplicity', 1) - 1,  # spin = multiplicity - 1
                basis=self.molecule.get('basis', 'sto-3g')
            )

            # Set the hamiltonian directly if available
            if self.molecule.get('hamiltonian'):
                mol._hamiltonian = self.molecule.get('hamiltonian')

            # Extract frequencies from the frequencies analysis result
            freq_result = self.results['frequencies']
            frequencies_cm = freq_result.get('frequencies', [])

            # Create thermochemistry calculator
            calc = ThermochemistryCalculator(
                molecule=mol,
                frequencies=frequencies_cm
            )

            # Compute thermochemistry at specified temperature and pressure
            thermo = calc.compute_thermochemistry(
                temperature=temperature,
                pressure=pressure
            )
            return thermo
        except Exception as e:
            logger.error(f"Thermochemistry calculation failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return {'error': f'Thermochemistry calculation failed: {str(e)}'}


    def _run_nmr_analysis(self, params: Dict) -> Dict:
        """Calculate NMR chemical shifts using Google ECHOES algorithm."""
        if not self.molecule.get('geometry'):
            return {'error': 'Molecular geometry required for NMR'}

        # TODO: Implement Google ECHOES NMR spectroscopy
        return {
            'message': 'NMR spectroscopy with Google ECHOES algorithm',
            'future_features': [
                '¹H NMR chemical shifts',
                '¹³C NMR chemical shifts',
                'Coupling constants (J-values)',
                'ECHOES quantum algorithm',
                'Shielding tensors'
            ]
        }

    def _run_frequencies_analysis(self, params: Dict) -> Dict:
        """Calculate vibrational frequencies."""
        if not self.molecule.get('geometry'):
            return {'error': 'Molecular geometry required'}

        if not self.molecule.get('hamiltonian'):
            return {'error': 'Hamiltonian object required for frequency calculation. This requires running with molecular solver.'}

        try:
            # Reconstruct molecule-like object for FrequencyCalculator
            from kanad.core.molecule import Molecule

            # Get atoms from molecule_data
            atoms = self.molecule.get('atoms', [])
            if not atoms:
                return {'error': 'Atom data not available'}

            # Create molecule with existing Atom objects
            mol = Molecule(
                atoms=atoms,
                charge=self.molecule.get('charge', 0),
                spin=self.molecule.get('multiplicity', 1) - 1,  # spin = multiplicity - 1
                basis=self.molecule.get('basis', 'sto-3g')
            )

            # Set the hamiltonian directly (bypass lazy construction)
            mol._hamiltonian = self.molecule.get('hamiltonian')

            calc = FrequencyCalculator(mol)
            result = calc.compute_frequencies(method='HF')  # Use HF for now
            return result
        except Exception as e:
            logger.error(f"Frequency calculation failed: {e}")
            return {'error': f'Frequency calculation failed: {str(e)}. This feature requires Hessian computation which is computationally expensive.'}

    def _run_uvvis_analysis(self, params: Dict) -> Dict:
        """Calculate UV-Vis absorption spectrum."""
        n_states = params.get('n_states', 5)
        broadening = params.get('broadening', 0.3)

        if not self.molecule.get('geometry'):
            return {'error': 'Molecular geometry required'}

        if not self.molecule.get('hamiltonian'):
            return {'error': 'Hamiltonian object required for UV-Vis calculation. This requires running with molecular solver.'}

        try:
            # Reconstruct molecule-like object
            from kanad.core.molecule import Molecule

            # Get atoms from molecule_data
            atoms = self.molecule.get('atoms', [])
            if not atoms:
                return {'error': 'Atom data not available'}

            # Create molecule with existing Atom objects
            mol = Molecule(
                atoms=atoms,
                charge=self.molecule.get('charge', 0),
                spin=self.molecule.get('multiplicity', 1) - 1,  # spin = multiplicity - 1
                basis=self.molecule.get('basis', 'sto-3g')
            )

            # Set the hamiltonian directly (bypass lazy construction)
            mol._hamiltonian = self.molecule.get('hamiltonian')

            calc = UVVisCalculator(mol)
            result = calc.compute_excitations(n_states=n_states, method='TDA')

            # Generate spectrum
            spectrum = calc.generate_spectrum(result, broadening=broadening)
            result['spectrum'] = spectrum

            return result
        except Exception as e:
            logger.error(f"UV-Vis calculation failed: {e}")
            return {'error': f'UV-Vis calculation failed: {str(e)}. This feature requires excited state calculations.'}

    def _run_excited_states_analysis(self, params: Dict) -> Dict:
        """Calculate excited states."""
        n_states = params.get('n_states', 3)
        method = params.get('method', 'tda')

        if not self.molecule.get('geometry'):
            return {'error': 'Molecular geometry required'}

        if not self.molecule.get('hamiltonian'):
            return {'error': 'Hamiltonian object required for excited state calculation. This requires running with molecular solver.'}

        try:
            # Reconstruct molecule-like object
            from kanad.core.molecule import Molecule

            # Get atoms from molecule_data
            atoms = self.molecule.get('atoms', [])
            if not atoms:
                return {'error': 'Atom data not available'}

            # Create molecule with existing Atom objects
            mol = Molecule(
                atoms=atoms,
                charge=self.molecule.get('charge', 0),
                spin=self.molecule.get('multiplicity', 1) - 1,  # spin = multiplicity - 1
                basis=self.molecule.get('basis', 'sto-3g')
            )

            # Set the hamiltonian directly (bypass lazy construction)
            mol._hamiltonian = self.molecule.get('hamiltonian')

            # Use TDA method (Tamm-Dancoff Approximation) for excited states
            # This is the standard method and works with HF/DFT ground states
            calc = UVVisCalculator(mol)
            result = calc.compute_excitations(n_states=n_states, method=method.upper())

            return result
        except Exception as e:
            logger.error(f"Excited state calculation failed: {e}")
            return {'error': f'Excited state calculation failed: {str(e)}. This feature requires excited state calculations.'}

    def _run_vibronic_analysis(self, params: Dict) -> Dict:
        """Calculate vibronic spectrum (vibrationally-resolved electronic spectrum)."""
        if not self.molecule.get('geometry'):
            return {'error': 'Molecular geometry required'}

        if not self.molecule.get('hamiltonian'):
            return {'error': 'Hamiltonian object required for vibronic calculation. This requires running with molecular solver.'}

        try:
            # Reconstruct molecule object
            from kanad.core.molecule import Molecule

            atoms = self.molecule.get('atoms', [])
            if not atoms:
                return {'error': 'Atom data not available'}

            mol = Molecule(
                atoms=atoms,
                charge=self.molecule.get('charge', 0),
                spin=self.molecule.get('multiplicity', 1) - 1,
                basis=self.molecule.get('basis', 'sto-3g')
            )
            mol._hamiltonian = self.molecule.get('hamiltonian')

            # First compute ground state frequencies
            freq_calc = FrequencyCalculator(mol)
            freq_result = freq_calc.compute_frequencies(method='HF')
            ground_frequencies = np.array(freq_result['frequencies'])

            # Then compute excited states for electronic transition energy
            uv_calc = UVVisCalculator(mol)
            excited_result = uv_calc.compute_excitations(n_states=1, method='TDA')

            if not excited_result.get('excitation_energies'):
                return {'error': 'Could not compute excited states for vibronic analysis'}

            electronic_transition = excited_result['excitation_energies'][0]  # First excited state (eV)

            # Generate vibronic spectrum
            vibronic_calc = VibronicCalculator(mol)
            spectrum = vibronic_calc.generate_vibronic_spectrum(
                electronic_transition=electronic_transition,
                ground_frequencies=ground_frequencies,
                excited_frequencies=None,  # Use ground state frequencies
                displacement=None,  # Assume minimal displacement
                temperature=params.get('temperature', 298.15),
                max_quanta=params.get('max_quanta', 5),
                wavelength_range=params.get('wavelength_range', (200, 800)),
                broadening=params.get('broadening', 0.01),
                n_points=params.get('n_points', 2000)
            )

            return {
                'electronic_transition': electronic_transition,
                'ground_frequencies': ground_frequencies.tolist(),
                'wavelengths': spectrum['wavelengths'].tolist(),
                'absorbance': spectrum['absorbance'].tolist(),
                'emission': spectrum['emission'].tolist() if 'emission' in spectrum else None,
                'method': 'Franck-Condon'
            }

        except Exception as e:
            logger.error(f"Vibronic calculation failed: {e}")
            return {'error': f'Vibronic calculation failed: {str(e)}'}

    def _run_dos_analysis(self, params: Dict) -> Dict:
        """Calculate density of states from molecular orbitals."""

        # Try to get orbital energies from results or compute from hamiltonian
        orbital_energies = None

        if self.results.get('orbital_energies') is not None:
            orbital_energies = np.array(self.results['orbital_energies'])
        elif self.molecule.get('hamiltonian'):
            # Compute orbital energies from hamiltonian
            try:
                hamiltonian = self.molecule.get('hamiltonian')
                if hasattr(hamiltonian, 'mo_energy'):
                    # Convert from Hartree to eV
                    Ha_to_eV = 27.211386245988
                    orbital_energies = np.array(hamiltonian.mo_energy) * Ha_to_eV
                elif hasattr(hamiltonian, 'orbital_energies'):
                    orbital_energies = np.array(hamiltonian.orbital_energies)
            except Exception as e:
                logger.error(f"Failed to get orbital energies from hamiltonian: {e}")

        if orbital_energies is None or len(orbital_energies) == 0:
            return {'error': 'Orbital energies required for DOS. Run experiment with a molecular solver (HF, VQE, SQD) first.'}

        broadening = params.get('broadening', 0.1)  # eV
        n_points = params.get('n_points', 500)

        # For molecules, compute DOS from discrete orbital energies
        orbital_energies = np.array(orbital_energies)  # in eV

        # Create energy range centered around the orbital energies
        e_min = np.min(orbital_energies) - 5  # eV
        e_max = np.max(orbital_energies) + 5  # eV

        energies = np.linspace(e_min, e_max, n_points)
        dos = np.zeros(n_points)

        # Gaussian broadening: each orbital contributes a Gaussian peak
        sigma = broadening
        for E_orbital in orbital_energies:
            dos += np.exp(-((energies - E_orbital)**2) / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi))

        # Find HOMO and LUMO
        n_electrons = self.molecule.get('n_electrons', len(orbital_energies))
        n_occupied = n_electrons // 2  # Assuming closed-shell

        homo_energy = orbital_energies[n_occupied - 1] if n_occupied > 0 else None
        lumo_energy = orbital_energies[n_occupied] if n_occupied < len(orbital_energies) else None
        band_gap = (lumo_energy - homo_energy) if (homo_energy is not None and lumo_energy is not None) else None

        return {
            'energies': energies.tolist(),
            'dos': dos.tolist(),
            'orbital_energies': orbital_energies.tolist(),
            'homo_energy': float(homo_energy) if homo_energy is not None else None,
            'lumo_energy': float(lumo_energy) if lumo_energy is not None else None,
            'band_gap': float(band_gap) if band_gap is not None else None,
            'n_orbitals': len(orbital_energies),
            'broadening': broadening,
            'units': 'eV'
        }

    def _run_uncertainty_analysis(self, params: Dict) -> Dict:
        """Analyze statistical uncertainty."""
        if not self.results.get('shots'):
            return {'error': 'Shot-based data required for uncertainty analysis'}

        analyzer = UncertaintyAnalyzer(
            energies=self.energies,
            shots=self.results['shots']
        )

        return {
            'uncertainty': analyzer.estimate_uncertainty(),
            'confidence_interval': analyzer.confidence_interval(),
            'required_shots': analyzer.required_shots_for_accuracy(target_error=0.001)
        }

    def get_available_analyses(self) -> Dict[str, Any]:
        """
        Get list of available analysis modules with descriptions.

        Returns:
            Dictionary mapping analysis names to metadata
        """
        return {
            'energy': {
                'name': 'Energy Analysis',
                'description': 'Energy decomposition, convergence, per-electron energy',
                'requires': ['final_energy'],
                'optional': ['nuclear_repulsion', 'energies']
            },
            'bonding': {
                'name': 'Bonding Analysis',
                'description': 'Bond orders, HOMO-LUMO gap, charge analysis',
                'requires': ['rdm1'],
                'optional': ['mo_coefficients']
            },
            'correlation': {
                'name': 'Correlation Energy',
                'description': 'Electron correlation energy and recovery percentage',
                'requires': ['final_energy', 'hf_energy'],
                'optional': ['fci_energy']
            },
            'properties': {
                'name': 'Molecular Properties',
                'description': 'Dipole moment, polarizability',
                'requires': ['rdm1', 'geometry'],
                'optional': []
            },
            'thermochemistry': {
                'name': 'Thermochemistry',
                'description': 'Enthalpy, entropy, Gibbs free energy',
                'requires': ['geometry', 'frequencies'],
                'optional': ['temperature', 'pressure']
            },
            'frequencies': {
                'name': 'Vibrational Frequencies',
                'description': 'Normal modes, ZPE, IR intensities',
                'requires': ['geometry'],
                'optional': ['hessian']
            },
            'uv_vis': {
                'name': 'UV-Vis Spectrum',
                'description': 'Electronic absorption spectrum',
                'requires': ['geometry'],
                'optional': ['n_states', 'broadening']
            },
            'excited_states': {
                'name': 'Excited States',
                'description': 'Electronic excited states energies and properties',
                'requires': ['geometry'],
                'optional': ['n_states', 'method']
            },
            'dos': {
                'name': 'Density of States',
                'description': 'Total and projected DOS',
                'requires': ['orbital_energies'],
                'optional': ['kpoints', 'broadening']
            },
            'uncertainty': {
                'name': 'Uncertainty Analysis',
                'description': 'Statistical uncertainty and confidence intervals',
                'requires': ['energies', 'shots'],
                'optional': []
            }
        }

    def get_available_profiles(self) -> Dict[str, Any]:
        """
        Get list of available domain profiles with descriptions.

        Returns:
            Dictionary mapping profile names to metadata
        """
        return {
            name: {
                'name': profile['name'],
                'description': profile['description'],
                'analyses': profile['analyses'],
                'recommended_for': profile.get('recommended_for', []),
                'typical_use_cases': profile.get('use_cases', [])
            }
            for name, profile in self.DOMAIN_PROFILES.items()
        }

    @classmethod
    def from_experiment_id(cls, experiment_id: str, db=None):
        """
        Create AnalysisService from an experiment ID by loading from database.

        Args:
            experiment_id: ID of completed experiment
            db: Database connection (ExperimentDB)

        Returns:
            AnalysisService instance
        """
        if db is None:
            from api.core.database import ExperimentDB
            db = ExperimentDB

        experiment = db.get(experiment_id)
        if not experiment:
            raise ValueError(f"Experiment {experiment_id} not found")

        if experiment.get('status') != 'completed':
            logger.warning(f"Experiment {experiment_id} is not completed (status={experiment.get('status')})")

        results = experiment.get('results', {})
        molecule_data = experiment.get('molecule', {})

        return cls(
            experiment_results=results,
            molecule_data=molecule_data
        )


def create_analysis_service(
    experiment_id: Optional[str] = None,
    experiment_results: Optional[Dict] = None,
    molecule_data: Optional[Dict] = None,
    progress_callback: Optional[Callable] = None
) -> AnalysisService:
    """
    Factory function to create AnalysisService from various sources.

    Args:
        experiment_id: Load from database by ID
        experiment_results: Direct results dictionary
        molecule_data: Molecule information
        progress_callback: Optional progress callback

    Returns:
        AnalysisService instance
    """
    if experiment_id:
        return AnalysisService.from_experiment_id(experiment_id)
    elif experiment_results:
        return AnalysisService(experiment_results, molecule_data, progress_callback)
    else:
        raise ValueError("Must provide either experiment_id or experiment_results")
