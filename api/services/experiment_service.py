"""
Experiment execution service - integrates Kanad quantum chemistry framework.

This service wraps all Kanad framework capabilities and exposes them through a
clean API interface. Handles VQE, QPE, SQD solvers with support for:
- Multiple ansatze (UCC, Hardware-Efficient, Governance-aware)
- Multiple mappers (Jordan-Wigner, Bravyi-Kitaev, Hybrid Orbital)
- Multiple backends (Classical, IBM Quantum, BlueQubit)
- SMILES parsing and validation
- Progress tracking and convergence monitoring
"""

import logging
import numpy as np
from typing import Dict, Any, Optional, Callable
from datetime import datetime
import traceback

# Kanad imports
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.io.smiles_parser import from_smiles, validate_smiles
from kanad.solvers.vqe_solver import VQESolver
from kanad.solvers.sqd_solver import SQDSolver
from kanad.solvers.excited_states_solver import ExcitedStatesSolver
from kanad.bonds import BondFactory

logger = logging.getLogger(__name__)


class ExperimentService:
    """
    Service for executing quantum chemistry experiments using Kanad framework.

    This class provides a high-level interface for running VQE, QPE, and SQD calculations
    with various configurations, backends, and analysis options.
    """

    def __init__(self):
        """Initialize experiment service."""
        self.logger = logger

    def validate_molecule_config(self, molecule_data: Dict[str, Any]) -> tuple[bool, str]:
        """
        Validate molecule configuration.

        Args:
            molecule_data: Dictionary with molecule information

        Returns:
            (is_valid, error_message)
        """
        # Check if SMILES or atoms are provided
        if not molecule_data.get('smiles') and not molecule_data.get('atoms'):
            return False, "Either 'smiles' or 'atoms' must be provided"

        # Validate SMILES if provided
        if molecule_data.get('smiles'):
            is_valid, error = validate_smiles(molecule_data['smiles'])
            if not is_valid:
                return False, f"Invalid SMILES: {error}"

        # Validate basis set
        valid_bases = ['sto-3g', '6-31g', '6-31g*', '6-31g**', 'cc-pvdz', 'cc-pvtz']
        basis = molecule_data.get('basis', 'sto-3g').lower()
        if basis not in valid_bases:
            return False, f"Invalid basis set: {basis}. Valid options: {', '.join(valid_bases)}"

        return True, ""

    def create_molecule(self, molecule_data: Dict[str, Any]) -> Molecule:
        """
        Create Kanad Molecule from configuration.

        Args:
            molecule_data: Dictionary containing:
                - smiles: SMILES string (optional)
                - atoms: List of atom dictionaries (optional)
                - basis: Basis set (default: sto-3g)
                - charge: Molecular charge (default: 0)
                - multiplicity: Spin multiplicity (default: 1)

        Returns:
            Kanad Molecule object

        Raises:
            ValueError: If molecule data is invalid
        """
        # Extract parameters
        smiles = molecule_data.get('smiles')
        atoms_data = molecule_data.get('atoms', [])
        basis = molecule_data.get('basis', 'sto-3g').lower()
        charge = molecule_data.get('charge', 0)
        multiplicity = molecule_data.get('multiplicity', 1)

        # Convert multiplicity to spin (2S where S is total spin)
        spin = multiplicity - 1

        # Create molecule from SMILES or atoms
        if smiles:
            logger.info(f"Creating molecule from SMILES: {smiles}")
            molecule = from_smiles(
                smiles,
                basis=basis,
                optimize_geometry=True
            )
            # Override charge/spin if provided
            if charge != 0:
                molecule.charge = charge
            if spin != 0:
                molecule.spin = spin
        else:
            # Create from atom list
            logger.info(f"Creating molecule from {len(atoms_data)} atoms")
            atoms = []
            for atom_data in atoms_data:
                # Support both 'element' and 'symbol' keys
                element = atom_data.get('element') or atom_data.get('symbol')
                if not element:
                    raise ValueError(f"Atom data missing 'element' or 'symbol' field: {atom_data}")

                atom = Atom(
                    symbol=element,
                    position=np.array([
                        atom_data.get('x', 0.0),
                        atom_data.get('y', 0.0),
                        atom_data.get('z', 0.0)
                    ])
                )
                atoms.append(atom)

            molecule = Molecule(
                atoms=atoms,
                charge=charge,
                spin=spin,
                basis=basis
            )

        logger.info(f"Molecule created: {molecule.formula}, {molecule.n_electrons} electrons")
        return molecule

    def create_bond(self, molecule: Molecule, config: Dict[str, Any]) -> Any:
        """
        Create bond from molecule for bond-based API.

        For simple diatomic molecules, creates a proper bond object.
        For larger molecules, returns the molecule's Hamiltonian directly.

        Args:
            molecule: Kanad Molecule
            config: Configuration dictionary

        Returns:
            Bond or Hamiltonian object
        """
        # For diatomic molecules, create a bond
        if molecule.n_atoms == 2:
            atom1, atom2 = molecule.atoms
            distance = np.linalg.norm(atom1.position - atom2.position)

            bond = BondFactory.create_bond(
                atom1.symbol,
                atom2.symbol,
                distance=distance,
                basis=molecule.basis
            )
            logger.info(f"Created bond: {atom1.symbol}-{atom2.symbol}, distance={distance:.3f} Å")
            return bond
        else:
            # For polyatomic molecules, return Hamiltonian directly
            logger.info(f"Using molecular Hamiltonian for {molecule.n_atoms}-atom system")
            return molecule.hamiltonian

    def execute_vqe(
        self,
        molecule: Molecule,
        config: Dict[str, Any],
        progress_callback: Optional[Callable] = None,
        db_session=None
    ) -> Dict[str, Any]:
        """
        Execute VQE calculation.

        Args:
            molecule: Kanad Molecule object
            config: Configuration dictionary with:
                - ansatz: Ansatz type (ucc, uccsd, hardware_efficient, governance_aware, two_local)
                - mapper: Qubit mapper (jordan_wigner, bravyi_kitaev, hybrid_orbital, parity)
                - optimizer: Classical optimizer (SLSQP, COBYLA, L-BFGS-B, ADAM, etc.)
                - backend: Backend type (classical, ibm_quantum, bluequbit_gpu)
                - backend_name: Specific backend name (e.g., ibm_torino)
                - max_iterations: Max optimization iterations
                - conv_threshold: Convergence threshold
            progress_callback: Optional callback function(iteration, energy, params)
            db_session: Database session for retrieving cloud credentials

        Returns:
            Dictionary with results
        """
        # Extract configuration
        ansatz_type = config.get('ansatz', 'ucc')
        # Normalize ansatz names
        if ansatz_type == 'governance':
            ansatz_type = 'governance_aware'

        mapper_type = config.get('mapper', 'jordan_wigner')
        optimizer = config.get('optimizer', 'SLSQP')
        backend = config.get('backend', 'classical')
        backend_name = config.get('backend_name')
        max_iterations = config.get('max_iterations', 1000)
        conv_threshold = config.get('conv_threshold', 1e-6)

        logger.info(f"Executing VQE: ansatz={ansatz_type}, mapper={mapper_type}, backend={backend}")

        # Create bond or use Hamiltonian directly
        try:
            bond = self.create_bond(molecule, config)
        except Exception as e:
            logger.warning(f"Could not create bond, using Hamiltonian directly: {e}")
            bond = molecule.hamiltonian

        # Normalize backend name
        if backend == 'classical':
            backend = 'statevector'
        elif backend == 'ibm_quantum':
            backend = 'ibm'

        # Initialize VQE solver
        try:
            if hasattr(bond, 'bond_type'):
                # Bond-based API (preferred for diatomic molecules)
                solver = VQESolver(
                    bond=bond,
                    ansatz_type=ansatz_type,
                    mapper_type=mapper_type,
                    optimizer=optimizer,
                    max_iterations=max_iterations,
                    conv_threshold=conv_threshold,
                    backend=backend,
                    enable_analysis=True,
                    enable_optimization=True
                )
            else:
                # Hamiltonian-based API (for polyatomic molecules)
                solver = VQESolver(
                    hamiltonian=bond,
                    molecule=molecule,  # Pass molecule explicitly
                    ansatz_type=ansatz_type,
                    mapper_type=mapper_type,
                    optimizer=optimizer,
                    max_iterations=max_iterations,
                    conv_threshold=conv_threshold,
                    backend=backend,
                    enable_analysis=True,
                    enable_optimization=True
                )
        except Exception as e:
            logger.error(f"Failed to initialize VQE solver: {e}")
            raise

        # Run VQE with progress tracking
        result = solver.solve(callback=progress_callback)

        # Format results for API response
        formatted_result = {
            'energy': float(result.get('energy', 0.0)),
            'hf_energy': float(result.get('hf_energy', 0.0)),
            'correlation_energy': float(result.get('correlation_energy', 0.0)),
            'converged': result.get('converged', False),
            'iterations': result.get('iterations', 0),
            'method': 'VQE',
            'ansatz': ansatz_type,
            'mapper': mapper_type,
            'optimizer': optimizer,
            'backend': backend,
            'convergence_history': self._format_convergence_history(result.get('energy_history')),
            'properties': self._extract_properties(result),
            'analysis': result.get('analysis', {}),
            'validation': result.get('validation', {}),
        }

        logger.info(f"VQE completed: E={formatted_result['energy']:.6f} Ha, converged={formatted_result['converged']}")

        # Convert all numpy types to Python native types for JSON serialization
        return self._convert_numpy_to_python(formatted_result)

    def execute_sqd(
        self,
        molecule: Molecule,
        config: Dict[str, Any],
        progress_callback: Optional[Callable] = None,
        db_session=None
    ) -> Dict[str, Any]:
        """
        Execute SQD (Subspace Quantum Diagonalization) calculation.

        Args:
            molecule: Kanad Molecule object
            config: Configuration dictionary
            progress_callback: Optional callback function
            db_session: Database session for retrieving cloud credentials

        Returns:
            Dictionary with results
        """
        logger.info("Executing SQD calculation")

        # Create bond or Hamiltonian
        try:
            bond = self.create_bond(molecule, config)
        except Exception as e:
            logger.warning(f"Could not create bond, using Hamiltonian directly: {e}")
            bond = molecule.hamiltonian

        # Initialize SQD solver
        n_states = config.get('n_states', 3)
        try:
            solver = SQDSolver(
                bond=bond,
                subspace_dim=config.get('subspace_dim', 10),
                backend=config.get('backend', 'statevector'),
                enable_analysis=True,
                enable_optimization=True
            )
        except Exception as e:
            logger.error(f"Failed to initialize SQD solver: {e}")
            raise

        # Run SQD with requested number of states
        result = solver.solve(n_states=n_states)

        # Format results
        formatted_result = {
            'energies': [float(e) for e in result.get('energies', [])],
            'ground_state_energy': float(result.get('ground_state_energy', 0.0)),
            'excited_state_energies': [float(e) for e in result.get('excited_state_energies', [])],
            'method': 'SQD',
            'n_states': n_states,
            'subspace_dim': config.get('subspace_dim', 10),
            'converged': result.get('converged', False),
            'hf_energy': float(result.get('hf_energy', 0.0)) if result.get('hf_energy') else None,
            'correlation_energy': float(result.get('correlation_energy', 0.0)) if result.get('correlation_energy') else None,
            'properties': self._extract_properties(result),
            'analysis': result.get('analysis', {}),
        }

        logger.info(f"SQD completed: Ground state E={formatted_result['ground_state_energy']:.6f} Ha")
        return self._convert_numpy_to_python(formatted_result)

    def execute_excited_states(
        self,
        molecule: Molecule,
        config: Dict[str, Any],
        progress_callback: Optional[Callable] = None,
        db_session=None
    ) -> Dict[str, Any]:
        """
        Execute excited states calculation.

        Args:
            molecule: Kanad Molecule object
            config: Configuration dictionary with:
                - excited_method: CIS, TDDFT, etc.
                - n_states: Number of excited states
            progress_callback: Optional callback function
            db_session: Database session for retrieving cloud credentials

        Returns:
            Dictionary with results
        """
        logger.info("Executing excited states calculation")

        # Create bond or Hamiltonian
        try:
            bond = self.create_bond(molecule, config)
        except Exception as e:
            logger.warning(f"Could not create bond, using Hamiltonian directly: {e}")
            bond = molecule.hamiltonian

        # Get excited states parameters
        excited_method = config.get('excited_method', 'cis').lower()
        n_states = config.get('n_states', 5)

        # Initialize solver
        try:
            solver = ExcitedStatesSolver(
                bond=bond,
                method=excited_method,
                n_states=n_states,
                enable_analysis=True,
                enable_optimization=False
            )
        except Exception as e:
            logger.error(f"Failed to initialize excited states solver: {e}")
            raise

        # Run calculation
        result = solver.solve()

        # Format results
        formatted_result = {
            'method': f'ExcitedStates-{excited_method.upper()}',
            'excited_method': excited_method,
            'n_states': n_states,
            'ground_state_energy': float(result.get('ground_state_energy', 0.0)),
            'excited_state_energies': [float(e) for e in result.get('excited_state_energies', [])],
            'excitation_energies_ev': [float(e) for e in result.get('excitation_energies_ev', [])],
            'oscillator_strengths': [float(f) for f in result.get('oscillator_strengths', [])] if result.get('oscillator_strengths') is not None else [],
            'dominant_transitions': result.get('dominant_transitions', []),
            'converged': result.get('converged', False),
            'analysis': result.get('analysis', {}),
        }

        logger.info(f"Excited states completed: {len(formatted_result['excitation_energies_ev'])} excitations found")
        return self._convert_numpy_to_python(formatted_result)

    def execute_experiment(
        self,
        molecule_data: Dict[str, Any],
        config: Dict[str, Any],
        progress_callback: Optional[Callable] = None,
        experiment_id: Optional[int] = None,
        db_session=None
    ) -> Dict[str, Any]:
        """
        Main entry point for executing experiments.

        Args:
            molecule_data: Molecule configuration
            config: Computation configuration
            progress_callback: Optional progress callback
            experiment_id: Optional experiment ID for cancellation checks
            db_session: Optional database session for cloud credentials

        Returns:
            Dictionary with results

        Raises:
            ValueError: If configuration is invalid
            Exception: If execution fails
        """
        start_time = datetime.now()

        try:
            # Validate molecule configuration
            is_valid, error = self.validate_molecule_config(molecule_data)
            if not is_valid:
                raise ValueError(error)

            # Create molecule
            molecule = self.create_molecule(molecule_data)

            # Get computation method
            method = config.get('method', 'VQE').upper()

            # Execute based on method
            if method == 'VQE':
                result = self.execute_vqe(molecule, config, progress_callback, db_session)
            elif method == 'SQD':
                result = self.execute_sqd(molecule, config, progress_callback, db_session)
            elif method == 'EXCITED' or method == 'EXCITED_STATES' or method == 'EXCITEDSTATES':
                result = self.execute_excited_states(molecule, config, progress_callback, db_session)
            elif method == 'HF':
                # Simple Hartree-Fock calculation
                result = {
                    'energy': float(molecule.hamiltonian.hf_energy),
                    'method': 'HF',
                    'converged': molecule.hamiltonian._scf_converged,
                    'iterations': molecule.hamiltonian._scf_iterations,
                }
            else:
                raise ValueError(f"Unknown method: {method}")

            # Add execution time
            execution_time = (datetime.now() - start_time).total_seconds()
            result['execution_time'] = execution_time
            result['molecule_formula'] = molecule.formula
            result['n_electrons'] = molecule.n_electrons
            result['n_orbitals'] = molecule.n_orbitals

            logger.info(f"Experiment completed successfully in {execution_time:.2f}s")
            return result

        except Exception as e:
            logger.error(f"Experiment failed: {e}")
            logger.error(traceback.format_exc())
            raise

    def _convert_numpy_to_python(self, obj):
        """Recursively convert numpy types to Python native types for JSON serialization."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.generic):
            return obj.item()
        elif isinstance(obj, dict):
            return {key: self._convert_numpy_to_python(value) for key, value in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [self._convert_numpy_to_python(item) for item in obj]
        else:
            return obj

    def _format_convergence_history(self, energy_history) -> list:
        """Format energy history for API response."""
        if energy_history is None:
            return []

        # Convert numpy array to list of dicts
        convergence = []
        if isinstance(energy_history, np.ndarray):
            for i, energy in enumerate(energy_history):
                convergence.append({
                    'iteration': i + 1,
                    'energy': float(energy)
                })

        return convergence

    def _extract_properties(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Extract molecular properties from result."""
        properties = {}

        # Extract from analysis if available
        if 'analysis' in result and result['analysis']:
            analysis = result['analysis']
            if 'properties' in analysis:
                properties = analysis['properties']

        # Add standard properties if not present
        if 'dipole_moment' not in properties and 'dipole_moment' in result:
            properties['dipole_moment'] = result['dipole_moment']

        return properties


# Global service instance
experiment_service = ExperimentService()
