"""
Core service for executing quantum chemistry experiments
"""

import traceback
import logging
import numpy as np
from typing import Dict, Any, Optional
from datetime import datetime

logger = logging.getLogger(__name__)

from rdkit import Chem
from rdkit.Chem import AllChem

from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.bonds import BondFactory
from kanad.solvers import VQESolver, SQDSolver, ExcitedStatesSolver
from kanad.backends.ibm import IBMBackend
from kanad.backends.bluequbit import BlueQubitBackend

from api.core.database import ExperimentDB, JobDB
from api.core.config import get_settings
from api.routes.websockets import manager as ws_manager
from api.utils import (
    broadcast_log_sync as _broadcast_log_sync,
    broadcast_convergence_sync as _broadcast_convergence_sync,
    broadcast_status_sync as _broadcast_status_sync,
    set_main_event_loop
)
import asyncio
import threading


class ExperimentCancelledException(Exception):
    """Exception raised when an experiment is cancelled by the user."""
    pass


def create_molecule_from_config(molecule_config: Dict[str, Any]) -> Molecule:
    """Create Molecule object from configuration."""
    if molecule_config.get('smiles'):
        # Create from SMILES using RDKit
        smiles = molecule_config['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        # Extract atoms and positions
        conf = mol.GetConformer()
        atoms = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            atoms.append(
                Atom(
                    atom.GetSymbol(),
                    position=np.array([pos.x, pos.y, pos.z])
                )
            )

        molecule = Molecule(
            atoms,
            basis=molecule_config.get('basis', 'sto-3g'),
            charge=molecule_config.get('charge', 0),
            spin=molecule_config.get('multiplicity', 1) - 1
        )

    elif molecule_config.get('atoms'):
        # Create from atoms
        atoms = [
            Atom(
                atom['symbol'],
                position=np.array([atom['x'], atom['y'], atom['z']])
            )
            for atom in molecule_config['atoms']
        ]

        # Calculate total electrons
        atomic_numbers = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10}
        charge = molecule_config.get('charge', 0)
        n_electrons = sum(atomic_numbers.get(atom['symbol'], 0) for atom in molecule_config['atoms']) - charge

        # Auto-correct multiplicity if needed
        multiplicity = molecule_config.get('multiplicity', 1)
        spin = multiplicity - 1  # spin = 2S where multiplicity = 2S+1

        # Check if electron count is consistent with multiplicity
        # For odd electrons, multiplicity must be even (doublet=2, quartet=4, etc.)
        # For even electrons, multiplicity must be odd (singlet=1, triplet=3, etc.)
        if (n_electrons % 2 == 1 and multiplicity % 2 == 1) or (n_electrons % 2 == 0 and multiplicity % 2 == 0):
            # Inconsistent - auto-correct to lowest valid multiplicity
            if n_electrons % 2 == 1:
                # Odd electrons -> use doublet (multiplicity=2)
                multiplicity = 2
                print(f"⚠️  Auto-corrected multiplicity to {multiplicity} (doublet) for {n_electrons} electrons")
            else:
                # Even electrons -> use singlet (multiplicity=1)
                multiplicity = 1
                print(f"⚠️  Auto-corrected multiplicity to {multiplicity} (singlet) for {n_electrons} electrons")
            spin = multiplicity - 1

        molecule = Molecule(
            atoms,
            basis=molecule_config.get('basis', 'sto-3g'),
            charge=charge,
            spin=spin
        )
    else:
        raise ValueError("Must provide either SMILES or atoms")

    return molecule


def get_cloud_credentials(provider: str) -> Dict[str, Any]:
    """
    Retrieve cloud credentials from database.

    Args:
        provider: 'ibm' or 'bluequbit'

    Returns:
        Dictionary with credentials

    Raises:
        ValueError: If credentials not found
    """
    from api.core.database import get_db
    import json

    with get_db() as conn:
        cursor = conn.cursor()
        cursor.execute(
            "SELECT credentials FROM cloud_credentials WHERE provider = ?",
            (provider,)
        )
        row = cursor.fetchone()

        if not row:
            raise ValueError(
                f"{provider.upper()} credentials not configured. "
                f"Please configure them in Settings → Backend"
            )

        return json.loads(row['credentials'])


def get_backend_kwargs(backend_config: Dict[str, Any], experiment_id: str = None) -> tuple:
    """
    Get backend type and kwargs for VQE solver.

    Returns:
        tuple: (backend_type: str, backend_kwargs: dict)
    """
    backend_type = backend_config.get('backend', 'classical')
    print(f"🔧 get_backend_kwargs called with backend_type: {backend_type}")

    if backend_type == 'classical' or backend_type == 'statevector':
        msg = f"📍 Using statevector simulation"
        if experiment_id:
            _broadcast_log_sync(experiment_id, msg)
        else:
            print(msg)
        return ('statevector', {})

    elif backend_type == 'ibm_quantum':
        msg = f"🌐 Configuring IBM Quantum backend..."
        if experiment_id:
            _broadcast_log_sync(experiment_id, msg)
        else:
            print(msg)

        # Try database first, fall back to environment variables
        try:
            credentials = get_cloud_credentials('ibm')
            api_token = credentials.get('api_token')
            crn = credentials.get('crn')
            print(f"✅ IBM credentials loaded from database")
        except ValueError:
            # Fallback to environment variables
            settings = get_settings()
            api_token = settings.IBM_API_TOKEN
            crn = settings.IBM_CRN
            print(f"✅ IBM credentials loaded from environment")

            if not api_token:
                raise ValueError(
                    "IBM Quantum credentials not configured. "
                    "Please configure them in Settings → Backend or set IBM_API and IBM_CRN environment variables"
                )

        backend_name = backend_config.get('backend_name', 'ibm_torino')
        msg = f"✅ Connected to IBM Quantum: {backend_name}"
        if experiment_id:
            _broadcast_log_sync(experiment_id, msg)
        else:
            print(msg)

        msg = f"🔗 Track jobs at: https://quantum.ibm.com/jobs"
        if experiment_id:
            _broadcast_log_sync(experiment_id, msg)
        else:
            print(msg)

        return ('ibm', {
            'backend_name': backend_name,
            'api_token': api_token,
            'instance': crn
        })

    elif backend_type == 'bluequbit':
        msg = f"🌐 Configuring BlueQubit backend..."
        if experiment_id:
            _broadcast_log_sync(experiment_id, msg)
        else:
            print(msg)
        print(f"🔍 Backend config received: {backend_config}")

        # Try database first, fall back to environment variables
        try:
            credentials = get_cloud_credentials('bluequbit')
            api_token = credentials.get('api_token')
            print(f"✅ BlueQubit credentials loaded from database")
        except ValueError:
            # Fallback to environment variables
            settings = get_settings()
            api_token = settings.BLUEQUBIT_API_TOKEN
            print(f"✅ BlueQubit credentials loaded from environment")

            if not api_token:
                raise ValueError(
                    "BlueQubit credentials not configured. "
                    "Please configure them in Settings → Backend or set BLUE_TOKEN environment variable"
                )

        # Get BlueQubit device/method options
        device = backend_config.get('bluequbit_device', 'gpu')
        print(f"🔧 BlueQubit device from config: {device}")
        bluequbit_options = {}

        # Add device-specific options
        if device.startswith('mps'):
            bluequbit_options['mps_bond_dimension'] = backend_config.get('mps_bond_dimension', 100)
            print(f"📊 MPS bond dimension: {bluequbit_options['mps_bond_dimension']}")
        elif device == 'pauli-path':
            truncation = backend_config.get('pauli_path_truncation_threshold')
            if truncation is not None:
                bluequbit_options['pauli_path_truncation_threshold'] = truncation
            transpilation_level = backend_config.get('pauli_path_circuit_transpilation_level')
            if transpilation_level is not None:
                bluequbit_options['pauli_path_circuit_transpilation_level'] = transpilation_level
            print(f"📊 Pauli-path options: {bluequbit_options}")

        msg = f"✅ Connected to BlueQubit cloud: device={device}"
        if experiment_id:
            _broadcast_log_sync(experiment_id, msg)
        else:
            print(msg)

        msg = f"🔗 Track jobs at: https://app.bluequbit.io/jobs"
        if experiment_id:
            _broadcast_log_sync(experiment_id, msg)
        else:
            print(msg)
        print(f"🎯 Backend kwargs: {{'device': '{device}', **bluequbit_options}}")
        return ('bluequbit', {
            'api_token': api_token,
            'device': device,
            **bluequbit_options
        })

    else:
        raise ValueError(f"Unknown backend type: {backend_type}")


def execute_hartree_fock(molecule: Molecule) -> Dict[str, Any]:
    """Execute Hartree-Fock calculation."""
    result = molecule.compute_energy(method='HF')

    return {
        'energy': float(result['energy']),
        'hf_energy': float(result['energy']),
        'converged': result['converged'],
        'iterations': result.get('n_iterations', 0),
        'n_orbitals': molecule.n_orbitals,
        'n_electrons': molecule.n_electrons,
        'method': 'HF'
    }


def check_cancellation(experiment_id: str, job_id: str):
    """Check if experiment/job has been cancelled."""
    experiment = ExperimentDB.get(experiment_id)
    job = JobDB.get(job_id)

    if experiment and experiment['status'] == 'cancelled':
        raise ExperimentCancelledException(f"Experiment {experiment_id} was cancelled")
    if job and job['status'] == 'cancelled':
        raise ExperimentCancelledException(f"Job {job_id} was cancelled")


def convert_numpy_types(obj):
    """Recursively convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [convert_numpy_types(item) for item in obj]
    elif isinstance(obj, (int, float, str, bool, type(None))):
        return obj
    else:
        # Skip non-serializable objects (like MolecularHamiltonian, TDA, etc.)
        # Return string representation as fallback
        try:
            import json
            json.dumps(obj)
            return obj
        except (TypeError, ValueError):
            logger.warning(f"Skipping non-serializable object of type {type(obj).__name__}")
            return None


def execute_vqe(
    molecule: Molecule,
    config: Dict[str, Any],
    job_id: str,
    experiment_id: str = None
) -> Dict[str, Any]:
    """Execute VQE calculation."""
    # Check for cancellation before starting
    if experiment_id:
        check_cancellation(experiment_id, job_id)

    # Create bond (for VQE high-level API)
    if molecule.n_atoms == 2:
        # Diatomic molecule - use bond API
        atom_1, atom_2 = molecule.atoms
        distance = atom_1.distance_to(atom_2)
        bond = BondFactory.create_bond(
            atom_1,
            atom_2,
            distance=distance,
            basis=molecule.basis
        )

        # Normalize ansatz type for VQESolver
        ansatz_type = config.get('ansatz', 'hardware_efficient')
        # Map specific governance ansatz types to generic 'governance'
        if ansatz_type in ['covalent_governance', 'ionic_governance', 'adaptive_governance', 'metallic_governance']:
            ansatz_type = 'governance'

        # Get backend type and kwargs
        backend_type, backend_kwargs = get_backend_kwargs(config, experiment_id)

        # Create VQE solver
        solver = VQESolver(
            bond=bond,
            ansatz_type=ansatz_type,
            mapper_type=config.get('mapper', 'jordan_wigner'),
            optimizer=config.get('optimizer', 'SLSQP'),
            max_iterations=config.get('max_iterations', 1000),
            backend=backend_type,
            shots=config.get('shots', 1024) if backend_type != 'statevector' else None,
            experiment_id=experiment_id,  # Pass for WebSocket broadcasting
            job_id=job_id,  # Pass for cancellation checking
            **backend_kwargs
        )

    else:
        # Multi-atom molecule - use low-level API
        from kanad.ansatze import UCCAnsatz, HardwareEfficientAnsatz
        from kanad.core.mappers import JordanWignerMapper

        hamiltonian = molecule.hamiltonian
        n_qubits = molecule.n_orbitals * 2  # Spin orbitals

        # Select ansatz
        ansatz_type = config.get('ansatz', 'hardware_efficient')
        if ansatz_type == 'ucc':
            ansatz = UCCAnsatz(n_qubits=n_qubits, n_electrons=molecule.n_electrons)
        else:
            ansatz = HardwareEfficientAnsatz(
                n_qubits=n_qubits,
                n_electrons=molecule.n_electrons,
                n_layers=3
            )

        mapper = JordanWignerMapper()

        # Get backend type and kwargs
        backend_type, backend_kwargs = get_backend_kwargs(config, experiment_id)

        solver = VQESolver(
            hamiltonian=hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            molecule=molecule,  # Pass molecule for analysis
            optimizer=config.get('optimizer', 'SLSQP'),
            max_iterations=config.get('max_iterations', 1000),
            backend=backend_type,
            shots=config.get('shots', 1024) if backend_type != 'statevector' else None,
            enable_analysis=True,  # Explicitly enable analysis
            experiment_id=experiment_id,  # Pass for WebSocket broadcasting
            job_id=job_id,  # Pass for cancellation checking
            **backend_kwargs
        )

    # Progress callback with cancellation check and real-time updates
    function_eval_count = [0]  # Use list to allow modification in nested function
    last_broadcasted_iter = [0]  # Track last iteration we broadcasted

    def progress_callback(iteration: int, energy: float, parameters: np.ndarray):
        # Check for cancellation
        if experiment_id:
            check_cancellation(experiment_id, job_id)

        # Track function evaluations
        function_eval_count[0] = iteration

        # Estimate actual optimizer iteration based on optimizer type
        optimizer = config.get('optimizer', 'SLSQP')
        if optimizer in ['SLSQP', 'L-BFGS-B']:
            # Gradient-based optimizers use ~40 function evals per iteration
            estimated_iteration = max(1, iteration // 40)
        elif optimizer in ['COBYLA', 'POWELL']:
            # Direct search optimizers use ~2-3 function evals per iteration
            estimated_iteration = max(1, iteration // 2)
        else:
            # Default: assume function eval = iteration
            estimated_iteration = iteration

        # Only broadcast when iteration actually changes (reduce spam)
        if estimated_iteration > last_broadcasted_iter[0]:
            last_broadcasted_iter[0] = estimated_iteration

            max_iter = config.get('max_iterations', 1000)
            progress = min(100.0, (estimated_iteration / max_iter) * 100.0)

            # Update job database with estimated iteration
            JobDB.update_progress(
                job_id,
                progress=progress,
                current_iteration=estimated_iteration,
                current_energy=float(energy)
            )

            # Send real-time WebSocket update
            if experiment_id:
                try:
                    # Use asyncio.run() for better compatibility
                    import concurrent.futures

                    async def send_update():
                        await ws_manager.broadcast_convergence(
                            experiment_id,
                            iteration=estimated_iteration,
                            energy=float(energy),
                            parameters=parameters.tolist() if parameters is not None else None
                        )

                    # Run async code in a thread pool to avoid event loop conflicts
                    with concurrent.futures.ThreadPoolExecutor() as executor:
                        future = executor.submit(asyncio.run, send_update())
                        future.result(timeout=1.0)  # 1 second timeout

                except Exception as e:
                    # Don't fail the experiment if WebSocket fails
                    print(f"⚠️  WebSocket broadcast failed: {e}")

    # Execute VQE with callback
    result = solver.solve(callback=progress_callback)

    # Debug: Check if analysis was generated
    if 'analysis' in result and result['analysis']:
        print(f"✅ Analysis data generated: {list(result['analysis'].keys())}")
    else:
        print(f"⚠️  No analysis data in results. Keys: {list(result.keys())}")

    # Build results dictionary with all data including analysis
    results_dict = {
        'energy': float(result['energy']),
        'hf_energy': float(result.get('hf_energy', 0.0)),
        'correlation_energy': float(result.get('correlation_energy', 0.0)),
        'converged': result['converged'],
        'iterations': result['iterations'],
        'convergence_history': [
            {'iteration': i, 'energy': float(e)}
            for i, e in enumerate(result.get('energy_history', []))
        ],
        'parameters': result['parameters'].tolist() if isinstance(result.get('parameters'), np.ndarray) else [],
        'method': 'VQE',
        'ansatz': config.get('ansatz'),
        'mapper': config.get('mapper'),
        # Enhanced data for analysis service
        'geometry': result.get('geometry'),
        'atoms': result.get('atoms'),
        'n_atoms': result.get('n_atoms'),
        'n_electrons': result.get('n_electrons'),
        'charge': result.get('charge'),
        'multiplicity': result.get('multiplicity'),
        'nuclear_repulsion': result.get('nuclear_repulsion'),
        'rdm1': result.get('rdm1'),
        'orbital_energies': result.get('orbital_energies'),
        'dipole': result.get('dipole'),
    }

    # Add analysis data if present (convert numpy types for JSON serialization)
    if 'analysis' in result and result['analysis']:
        results_dict['analysis'] = convert_numpy_types(result['analysis'])
        print(f"✅ Added analysis to results_dict")
    else:
        print(f"⚠️  No analysis to add to results")

    return results_dict


def execute_sqd(
    molecule: Molecule,
    config: Dict[str, Any],
    job_id: str,
    experiment_id: str = None
) -> Dict[str, Any]:
    """Execute SQD calculation."""
    print(f"🚀 execute_sqd called for experiment {experiment_id}")
    # Check for cancellation before starting
    if experiment_id:
        check_cancellation(experiment_id, job_id)
    print(f"📊 Starting SQD setup...")
    # SQD requires bond API - only works for diatomic molecules
    if molecule.n_atoms != 2:
        raise NotImplementedError(
            f"SQD currently only supports diatomic molecules. "
            f"Your molecule has {molecule.n_atoms} atoms. "
            f"Please use VQE for multi-atom molecules."
        )

    # Diatomic molecule - use bond API
    atom_1, atom_2 = molecule.atoms
    distance = atom_1.distance_to(atom_2)
    bond = BondFactory.create_bond(
        atom_1,
        atom_2,
        distance=distance,
        basis=molecule.basis
    )
    print(f"✅ Created bond for diatomic molecule: {atom_1.symbol}-{atom_2.symbol}")

    # Get backend configuration
    backend_type, backend_kwargs = get_backend_kwargs(config, experiment_id)
    print(f"🔧 SQD backend_type: {backend_type}")
    print(f"🔧 SQD backend_kwargs: {backend_kwargs}")

    # Create SQD solver
    solver = SQDSolver(
        bond=bond,
        subspace_dim=config.get('subspace_dim', 10),
        circuit_depth=config.get('circuit_depth', 3),
        backend=backend_type,
        shots=config.get('shots', 1024) if backend_type != 'statevector' else None,
        enable_analysis=True,
        enable_optimization=True,
        experiment_id=experiment_id,  # Pass for WebSocket broadcasting
        **backend_kwargs
    )
    print(f"✅ SQD solver created with backend: {solver.backend}")

    # Update progress
    JobDB.update_progress(job_id, progress=10.0)

    # Define progress callback for real-time updates
    # For ground state SQD, we only need 1 state (the ground state)
    # n_states is only used for excited states calculations
    n_states = 1  # Ground state only
    total_stages = 4 + n_states  # 0=init, 1=basis, 2=projection, 3=diag, 4=ground state

    # Track energy history for convergence graph
    energy_history = []

    def sqd_progress_callback(stage: int, energy: float, message: str):
        """Progress callback for SQD execution."""
        # Check for cancellation
        if experiment_id:
            check_cancellation(experiment_id, job_id)

        # Calculate progress based on stage
        progress = 10.0 + (stage / total_stages) * 85.0
        progress = min(progress, 95.0)

        # Update database
        JobDB.update_progress(job_id, progress=progress, current_iteration=stage)

        # Log progress
        print(f"📊 SQD Progress: Stage {stage}/{total_stages-1} - {message}, E = {energy:.8f} Ha")

        # Store energy for history
        energy_history.append(float(energy))

        # Send WebSocket update
        if experiment_id:
            _broadcast_convergence_sync(experiment_id, stage, float(energy))

    # Execute SQD with callback
    result = solver.solve(n_states=n_states, callback=sqd_progress_callback)

    # Update progress
    JobDB.update_progress(job_id, progress=100.0)

    results_dict = {
        'energy': float(result['energy']),
        'hf_energy': float(result.get('hf_energy', 0.0)),
        'correlation_energy': float(result.get('correlation_energy', 0.0)),
        'converged': result['converged'],
        'iterations': result['iterations'],
        'method': 'SQD',
        'subspace_dim': result['subspace_dim'],
        'energies': [float(e) for e in result['energies']],
        'excited_state_energies': [float(e) for e in result.get('excited_state_energies', [])],
        'circuit_depth': result.get('circuit_depth', 3),
        'energy_history': energy_history,  # For convergence graph
        # Enhanced data for analysis service
        'geometry': result.get('geometry'),
        'atoms': result.get('atoms'),
        'n_atoms': result.get('n_atoms'),
        'n_electrons': result.get('n_electrons'),
        'charge': result.get('charge'),
        'multiplicity': result.get('multiplicity'),
        'nuclear_repulsion': result.get('nuclear_repulsion'),
        'rdm1': result.get('rdm1'),
        'orbital_energies': result.get('orbital_energies'),
        'dipole': result.get('dipole'),
    }

    print(f"📊 SQD energy_history: {len(energy_history)} stages recorded")

    # Debug: Check what's in result
    print(f"🔍 SQD result keys: {result.keys()}")
    print(f"🔍 Has analysis: {'analysis' in result}")
    if 'analysis' in result:
        print(f"🔍 Analysis keys: {result['analysis'].keys() if result['analysis'] else 'None'}")

    # Add analysis data if present
    if 'analysis' in result and result['analysis']:
        results_dict['analysis'] = convert_numpy_types(result['analysis'])
        print(f"✅ Added analysis to results_dict: {list(results_dict['analysis'].keys())}")
    else:
        print(f"⚠️  No analysis data in SQD result!")

    return results_dict


def execute_excited_states(
    molecule: Molecule,
    config: Dict[str, Any],
    job_id: str,
    experiment_id: str = None
) -> Dict[str, Any]:
    """Execute excited states calculation."""
    # Check for cancellation before starting
    if experiment_id:
        check_cancellation(experiment_id, job_id)
    # Check if using quantum backend
    backend_type = config.get('backend', 'classical')

    # Get backend configuration
    backend_kwargs = {}

    # Get user-selected method FIRST - respect user's choice!
    # Support both snake_case (excited_method) and camelCase (excitedMethod)
    user_method = config.get('excited_method') or config.get('excitedMethod', 'cis')

    print(f"🔍 User selected excited states method: {user_method}")

    # Only configure backend_kwargs if method requires quantum execution
    if user_method == 'vqe' and backend_type in ['bluequbit', 'ibm_quantum']:
        print(f"⚠️  VQE excited states with quantum backend - requires many jobs")
        # Get backend credentials for VQE
        backend_type, backend_kwargs = get_backend_kwargs(config, experiment_id)
    elif backend_type in ['bluequbit', 'ibm_quantum'] and user_method != 'vqe':
        # User selected CIS or TDDFT with quantum backend - override to classical
        print(f"📊 {user_method.upper()} is classical-only, switching backend to classical")
        backend_type = 'classical'

    # Classical excited states (CIS/TDDFT) - create bond wrapper
    if molecule.n_atoms == 2:
        # Diatomic molecule - use bond API
        atom_1, atom_2 = molecule.atoms
        distance = atom_1.distance_to(atom_2)
        bond = BondFactory.create_bond(
            atom_1,
            atom_2,
            distance=distance,
            basis=molecule.basis
        )
        print(f"✅ Created bond for diatomic molecule: {atom_1.symbol}-{atom_2.symbol}")
    else:
        # Multi-atom molecule - create a minimal bond wrapper
        # CIS/TDDFT only need hamiltonian and molecule, not actual bond properties
        print(f"⚠️  Multi-atom molecule ({molecule.n_atoms} atoms) - creating wrapper for Excited States")

        # Create a simple wrapper class that mimics a bond
        class MoleculeBondWrapper:
            def __init__(self, molecule):
                self.molecule = molecule
                self.hamiltonian = molecule.hamiltonian
                self.atoms = molecule.atoms
                self.bond_type = "molecular"  # Not a traditional bond
                self.bond_order = None
                self.distance = None

        bond = MoleculeBondWrapper(molecule)
        print(f"✅ Created molecule wrapper for classical excited states")

    # Use the user-selected method (already read above)
    method = user_method

    # Support both snake_case and camelCase for n_states
    # Try each key in order, using the first non-None, non-zero value
    _excited_n_states = config.get('excited_n_states')
    _excitedNStates = config.get('excitedNStates')
    _n_states = config.get('n_states')

    print(f"🔍 DEBUG Excited States Config:")
    print(f"   config.get('excited_method'): {config.get('excited_method')}")
    print(f"   config.get('excitedMethod'): {config.get('excitedMethod')}")
    print(f"   config.get('excited_n_states'): {_excited_n_states} (type: {type(_excited_n_states)})")
    print(f"   config.get('excitedNStates'): {_excitedNStates} (type: {type(_excitedNStates)})")
    print(f"   config.get('n_states'): {_n_states} (type: {type(_n_states)})")
    print(f"   FINAL method: {method}")

    # Use the first non-None value
    if _excited_n_states is not None:
        n_states = _excited_n_states
    elif _excitedNStates is not None:
        n_states = _excitedNStates
    elif _n_states is not None:
        n_states = _n_states
    else:
        n_states = 5

    print(f"   FINAL n_states: {n_states}")
    print(f"   Backend type: {backend_type}")

    # Broadcast initial status to frontend
    if experiment_id:
        _broadcast_log_sync(experiment_id, f"🔬 Starting {method.upper()} excited states calculation...")
        _broadcast_log_sync(experiment_id, f"📊 Computing {n_states} excited states")

    # Build solver kwargs
    solver_kwargs = {
        'bond': bond,
        'method': method,
        'n_states': n_states,
        'enable_analysis': True,
        'enable_optimization': False,
        'experiment_id': experiment_id  # Pass for WebSocket broadcasting
    }

    # Add method-specific backend settings
    if method == 'vqe':
        solver_kwargs['backend'] = backend_type
        solver_kwargs['ansatz'] = config.get('ansatz', 'ucc')  # Fixed: 'ucc' not 'uccsd'
        solver_kwargs['optimizer'] = config.get('optimizer', 'COBYLA')
        solver_kwargs['max_iterations'] = config.get('max_iterations', 100)
        solver_kwargs['penalty_weight'] = config.get('penalty_weight', 5.0)
        solver_kwargs['backend_kwargs'] = backend_kwargs  # Pass credentials

        # Create progress callback for VQE (broadcasts convergence to frontend)
        def vqe_progress_callback(iteration, energy, parameters):
            """Broadcast VQE progress updates during optimization"""
            try:
                _broadcast_convergence_sync(experiment_id, iteration, float(energy))

                # Update progress bar (20% to 80% during optimization)
                max_iters = solver_kwargs.get('max_iterations', 100)
                if max_iters > 0:
                    progress = 20.0 + (iteration / max_iters) * 60.0
                    progress = min(progress, 80.0)
                    print(f"🔍 DEBUG: Updating job_id={job_id}, progress={progress:.1f}%, iter={iteration}")
                    JobDB.update_progress(job_id, progress=progress, current_iteration=iteration)
                    print(f"📊 Live VQE update: iter={iteration}, E={energy:.8f}, progress={progress:.1f}%")
            except Exception as e:
                logger.error(f"VQE progress broadcast failed: {e}")

        solver_kwargs['vqe_callback'] = vqe_progress_callback  # Add progress callback

        print(f"🔧 VQE Excited States configuration:")
        print(f"   Backend: {backend_type}")
        print(f"   Ansatz: {solver_kwargs['ansatz']}")
        print(f"   Optimizer: {solver_kwargs['optimizer']}")
        print(f"   Max iterations: {solver_kwargs['max_iterations']}")
        print(f"   Penalty weight: {solver_kwargs['penalty_weight']}")
        print(f"   Backend kwargs keys: {list(backend_kwargs.keys())}")

    elif method == 'sqd':
        solver_kwargs['backend'] = backend_type
        solver_kwargs['subspace_dim'] = config.get('subspace_dim', config.get('subspaceDim', 10))
        solver_kwargs['circuit_depth'] = config.get('circuit_depth', config.get('circuitDepth', 3))
        solver_kwargs['backend_kwargs'] = backend_kwargs  # Pass credentials for cloud backends

        print(f"🔧 SQD Excited States configuration:")
        print(f"   Backend: {backend_type}")
        print(f"   Subspace dimension: {solver_kwargs['subspace_dim']}")
        print(f"   Circuit depth: {solver_kwargs['circuit_depth']}")
        print(f"   Number of states: {n_states}")
        print(f"   Backend kwargs keys: {list(backend_kwargs.keys())}")

    solver = ExcitedStatesSolver(**solver_kwargs)

    # Update progress
    JobDB.update_progress(job_id, progress=20.0, current_iteration=0)
    backend_desc = f"quantum ({backend_type})" if method == 'vqe' and backend_type in ['bluequbit', 'ibm_quantum'] else "classical"
    print(f"📊 Excited States: Computing states using {method.upper()} ({backend_desc})...")

    # Broadcast that calculation is running
    if experiment_id:
        _broadcast_log_sync(experiment_id, f"⚙️ Running {method.upper()} solver...")

    # Execute
    result = solver.solve()

    # Broadcast completion
    if experiment_id:
        _broadcast_log_sync(experiment_id, f"✅ {method.upper()} calculation complete")

    # DEBUG: Check result structure
    print(f"🔍 DEBUG: Excited states result keys: {result.keys()}")
    print(f"🔍 DEBUG: Has 'energies' key: {'energies' in result}")
    if 'energies' in result:
        print(f"🔍 DEBUG: energies type: {type(result['energies'])}, length: {len(result['energies'])}")
        print(f"🔍 DEBUG: energies values: {result['energies']}")

    # Broadcast final energies for all computed states
    if experiment_id and 'energies' in result:
        actual_states = len(result['energies'])
        print(f"📊 Excited States result: requested={n_states}, found={actual_states} states (ground + {actual_states-1} excited)")
        print(f"📊 Broadcasting {actual_states} state energies...")
        _broadcast_log_sync(experiment_id, f"📊 Broadcasting {actual_states} state energies to frontend...")

        for i, energy in enumerate(result['energies']):
            # Check for cancellation during broadcasting
            check_cancellation(experiment_id, job_id)

            JobDB.update_progress(job_id, progress=50.0 + (i / actual_states) * 50.0, current_iteration=i+1)
            print(f"📊 Broadcasting state {i+1}/{actual_states}: iter={i+1}, E={energy:.8f}")
            _broadcast_convergence_sync(experiment_id, i+1, float(energy))
            # Small delay to ensure WebSocket messages are sent sequentially
            import time
            time.sleep(0.1)

        print(f"✅ Finished broadcasting {actual_states} convergence points")
    else:
        print(f"⚠️  WARNING: Cannot broadcast convergence - experiment_id={experiment_id}, has_energies={'energies' in result}")

    # Update final progress
    JobDB.update_progress(job_id, progress=100.0)

    # Generate analysis if not present in result
    analysis_data = result.get('analysis')
    if not analysis_data:
        # Generate analysis using solver's built-in method
        try:
            # Get HF density matrix for analysis
            density_matrix, _ = solver.hamiltonian.solve_scf(
                max_iterations=100,
                conv_tol=1e-8,
                use_diis=True
            )
            # Use BaseSolver's analysis method
            solver._add_analysis_to_results(
                energy=result['ground_state_energy'],
                density_matrix=density_matrix
            )
            analysis_data = solver.results.get('analysis')
            if analysis_data:
                print(f"✅ Generated analysis for EXCITED_STATES")
            else:
                print(f"⚠️  No analysis generated")
        except Exception as e:
            print(f"⚠️  Analysis generation failed: {e}")
            traceback.print_exc()
            analysis_data = None

    # Build energy history from all states for graph
    energy_history = []
    if 'energies' in result and result['energies'] is not None:
        # Include ground state and excited states
        energy_history = [float(e) for e in result['energies'][:n_states]]
    elif 'excited_state_energies' in result:
        # Fallback: just ground state and first few excited states
        energy_history = [float(result['ground_state_energy'])]
        energy_history.extend([float(e) for e in result['excited_state_energies'][:n_states-1]])

    # Format results
    results_dict = {
        'energy': float(result['energy']),  # Ground state energy
        'ground_state_energy': float(result.get('ground_state_energy', result['energy'])),
        'excited_state_energies': [float(e) for e in result.get('excited_state_energies', [])],
        'excitation_energies_ev': [float(e) for e in result.get('excitation_energies_ev', [])],
        'oscillator_strengths': [float(f) for f in result.get('oscillator_strengths', [])],
        'dominant_transitions': result.get('dominant_transitions', []),
        'converged': result['converged'],
        'iterations': result['iterations'],
        'method': 'EXCITED_STATES',
        'excited_method': method,
        'n_states': n_states,
        'energy_history': energy_history,  # For graph display
    }

    # Add analysis if available
    if analysis_data:
        results_dict['analysis'] = convert_numpy_types(analysis_data)
        print(f"✅ Added analysis to EXCITED_STATES results")

    return results_dict


def execute_experiment(experiment_id: str, job_id: str):
    """
    Execute a quantum chemistry experiment.

    This is called as a background task.
    """
    import time
    start_time = time.time()  # Track execution time

    try:
        print(f"⚛️  Starting experiment {experiment_id}")

        # Update status
        ExperimentDB.update_status(experiment_id, 'running')
        JobDB.update_status(job_id, 'running')

        # Get experiment data
        experiment = ExperimentDB.get(experiment_id)
        if not experiment:
            raise ValueError(f"Experiment {experiment_id} not found")

        molecule_config = experiment['molecule']
        config = experiment['configuration']

        # Create molecule
        print(f"📊 Creating molecule: {molecule_config}")
        molecule = create_molecule_from_config(molecule_config)
        print(f"✅ Molecule created: {molecule.formula} ({molecule.n_electrons} electrons, {molecule.n_orbitals} orbitals)")

        # Execute based on method
        method = config.get('method', 'VQE').upper()

        if method == 'HF':
            print("🔬 Running Hartree-Fock calculation...")
            results = execute_hartree_fock(molecule)

        elif method == 'VQE':
            print("🔬 Running VQE calculation...")
            results = execute_vqe(molecule, config, job_id, experiment_id)

        elif method == 'SQD':
            print("🔬 Running SQD calculation...")
            results = execute_sqd(molecule, config, job_id, experiment_id)

        elif method == 'EXCITED_STATES':
            print("🔬 Running Excited States calculation...")
            results = execute_excited_states(molecule, config, job_id, experiment_id)

        else:
            raise ValueError(f"Unknown method: {method}")

        # Add molecule info to results
        results['molecule_info'] = {
            'formula': molecule.formula,
            'n_atoms': molecule.n_atoms,
            'n_electrons': molecule.n_electrons,
            'n_orbitals': molecule.n_orbitals,
        }

        print(f"✅ Experiment completed: Energy = {results['energy']:.8f} Ha")

        # Run advanced analysis if enabled
        if config.get('advancedAnalysisEnabled') and config.get('advancedAnalysisProfile'):
            try:
                print(f"🔬 Running advanced analysis: {config['advancedAnalysisProfile']}")
                _broadcast_log_sync(experiment_id, f"Running {config['advancedAnalysisProfile']} analysis...")

                from kanad.services.analysis_service import AnalysisService

                # Prepare molecule data for analysis
                # Extract orbital energies from hamiltonian for DOS analysis (if not already in results)
                # Don't include the hamiltonian object itself to avoid JSON serialization issues
                if not results.get('orbital_energies') and hasattr(molecule, 'hamiltonian'):
                    try:
                        hamiltonian = molecule.hamiltonian
                        if hasattr(hamiltonian, 'mf') and hasattr(hamiltonian.mf, 'mo_energy'):
                            Ha_to_eV = 27.211386245988
                            results['orbital_energies'] = (hamiltonian.mf.mo_energy * Ha_to_eV).tolist()
                            logger.info(f"✅ Extracted orbital energies from hamiltonian for DOS analysis")
                    except Exception as e:
                        logger.warning(f"Could not extract orbital energies from hamiltonian: {e}")

                molecule_data = {
                    'geometry': molecule.geometry.tolist() if hasattr(molecule, 'geometry') else results.get('geometry'),
                    'atoms': molecule.atoms if hasattr(molecule, 'atoms') else results.get('atoms'),
                    'n_atoms': molecule.n_atoms if hasattr(molecule, 'n_atoms') else results.get('n_atoms'),
                    'n_electrons': molecule.n_electrons if hasattr(molecule, 'n_electrons') else results.get('n_electrons'),
                    'charge': molecule.charge if hasattr(molecule, 'charge') else results.get('charge', 0),
                    'multiplicity': molecule.multiplicity if hasattr(molecule, 'multiplicity') else results.get('multiplicity', 1),
                    'formula': molecule.formula if hasattr(molecule, 'formula') else None,
                    'smiles': molecule.smiles if hasattr(molecule, 'smiles') else None,
                    'basis': molecule.basis if hasattr(molecule, 'basis') else results.get('basis', 'sto-3g'),
                    # Include hamiltonian for advanced analysis (frequencies, etc.)
                    # It will be filtered out by convert_numpy_types when storing results
                    'hamiltonian': molecule.hamiltonian if hasattr(molecule, 'hamiltonian') else None,
                }

                # Run advanced analysis - AnalysisService expects experiment_results and molecule_data
                analyzer = AnalysisService(
                    experiment_results=results,
                    molecule_data=molecule_data
                )

                analysis_results = analyzer.run_analysis_profile(config['advancedAnalysisProfile'])

                # Add advanced analysis results (convert numpy types for JSON serialization)
                results['advanced_analysis'] = {
                    'profile': config['advancedAnalysisProfile'],
                    'results': convert_numpy_types(analysis_results.get('analyses', {})),  # Convert numpy types
                    'status': 'completed',
                }

                print(f"✅ Advanced analysis completed")
                _broadcast_log_sync(experiment_id, "Advanced analysis completed")

            except Exception as e:
                print(f"⚠️  Advanced analysis failed: {e}")
                traceback.print_exc()
                _broadcast_log_sync(experiment_id, f"Advanced analysis failed: {str(e)}")
                results['advanced_analysis'] = {
                    'profile': config['advancedAnalysisProfile'],
                    'status': 'failed',
                    'error': str(e)
                }

        # Update experiment with results
        # Convert numpy arrays to lists for JSON serialization
        from kanad.services.analysis_service import make_json_serializable
        results_serializable = make_json_serializable(results)

        ExperimentDB.update_status(
            experiment_id,
            'completed',
            results=results_serializable
        )

        # Extract and save cloud provider job information
        if results.get('cloud_provider'):
            provider = results['cloud_provider']
            job_ids = results.get('cloud_job_ids', [])
            job_urls = results.get('cloud_job_urls', [])
            execution_mode = results.get('execution_mode')

            # Save the first job ID and URL to the experiment record
            provider_job_id = job_ids[0] if job_ids else None
            provider_job_url = job_urls[0] if job_urls else None

            if provider_job_id:
                ExperimentDB.update_cloud_job_info(
                    experiment_id,
                    provider_job_id=provider_job_id,
                    provider_job_url=provider_job_url,
                    execution_mode=execution_mode
                )
                print(f"✅ Saved cloud job info: {provider} job {provider_job_id}")

        # Update job
        JobDB.update_progress(job_id, progress=100.0)
        JobDB.update_status(job_id, 'completed')

        # Broadcast completion status via WebSocket
        _broadcast_status_sync(experiment_id, status='completed', progress=100.0)

    except ExperimentCancelledException as e:
        # Handle cancellation gracefully
        execution_time = time.time() - start_time
        print(f"⚠️  Experiment cancelled after {execution_time:.2f} seconds: {e}")

        # Save partial results with execution time if available
        try:
            experiment = ExperimentDB.get(experiment_id)
            if experiment and experiment.get('results'):
                results = experiment['results']
                results['execution_time'] = round(execution_time, 2)
                results['status_note'] = 'cancelled'
                ExperimentDB.update_status(experiment_id, 'cancelled', results=results)
            else:
                ExperimentDB.update_status(experiment_id, 'cancelled')
        except:
            ExperimentDB.update_status(experiment_id, 'cancelled')

        JobDB.update_status(job_id, 'cancelled')

        # CRITICAL: Broadcast cancellation status via WebSocket so frontend updates
        _broadcast_status_sync(experiment_id, status='cancelled', progress=0.0)
        _broadcast_log_sync(experiment_id, "🚫 Experiment cancelled by user")

    except Exception as e:
        execution_time = time.time() - start_time
        print(f"❌ Experiment failed after {execution_time:.2f} seconds: {e}")
        traceback.print_exc()

        # Update experiment as failed with execution time
        ExperimentDB.update_status(
            experiment_id,
            'failed',
            error_message=str(e),
            results={'execution_time': round(execution_time, 2), 'status_note': 'failed'}
        )

        # Update job as failed
        JobDB.update_status(job_id, 'failed')
