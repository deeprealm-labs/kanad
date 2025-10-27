"""
Core service for executing quantum chemistry experiments
"""

import traceback
import numpy as np
from typing import Dict, Any, Optional
from datetime import datetime

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
from api.utils import broadcast_log_sync as _broadcast_log_sync, set_main_event_loop
import asyncio
import threading

# Store reference to main event loop (set during startup)
_main_event_loop = None


class ExperimentCancelledException(Exception):
    """Exception raised when an experiment is cancelled by the user."""
    pass


def _broadcast_convergence_sync(experiment_id: str, iteration: int, energy: float):
    """
    Broadcast convergence update from sync context.

    This function safely broadcasts WebSocket updates from synchronous code
    by submitting the coroutine to the main event loop.
    """
    global _main_event_loop

    try:
        # Try to get the stored main loop first
        loop = _main_event_loop

        # Fallback: try to get running loop (works in async context)
        if loop is None:
            try:
                loop = asyncio.get_running_loop()
            except RuntimeError:
                pass

        # If still no loop, we can't broadcast
        if loop is None:
            print(f"⚠️  No event loop available for WebSocket broadcast (iteration {iteration})")
            return

        # Create the coroutine
        coro = ws_manager.broadcast_convergence(
            experiment_id,
            iteration=iteration,
            energy=energy
        )

        # Schedule it on the main loop (thread-safe)
        future = asyncio.run_coroutine_threadsafe(coro, loop)

        # Wait briefly for completion (non-blocking)
        future.result(timeout=0.5)

    except RuntimeError as e:
        if "no running event loop" in str(e).lower() or "no current event loop" in str(e).lower():
            print(f"⚠️  No event loop available for WebSocket broadcast (iteration {iteration})")
        else:
            print(f"⚠️  WebSocket broadcast failed: {e}")
    except Exception as e:
        print(f"⚠️  WebSocket broadcast failed: {e}")


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
    else:
        return obj


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
            **backend_kwargs
        )

    # Progress callback with cancellation check and real-time updates
    function_eval_count = [0]  # Use list to allow modification in nested function

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
            estimated_iteration = iteration

        max_iter = config.get('max_iterations', 1000)
        progress = min(100.0, (estimated_iteration / max_iter) * 100.0)

        # Update job database with estimated iteration
        JobDB.update_progress(
            job_id,
            progress=progress,
            current_iteration=estimated_iteration,  # Use estimated iteration, not function evals
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
                        iteration=estimated_iteration,  # Send estimated iteration
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
    n_states = config.get('n_states', 3)
    total_stages = 4 + n_states  # 0=init, 1=basis, 2=projection, 3=diag, 4+=states

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

    # For quantum backends, use SQD for diatomic or VQE for multi-atom molecules
    if backend_type in ['bluequbit', 'ibm_quantum']:
        if molecule.n_atoms == 2:
            # Use SQD for diatomic molecules (more efficient)
            print(f"⚠️  Excited States with quantum backend detected - using SQD")
            print(f"📊 Quantum excited states computed via Subspace Quantum Diagonalization")
            return execute_sqd(molecule, config, job_id, experiment_id)
        else:
            # Use VQE-based approach for multi-atom molecules
            print(f"⚠️  Multi-atom molecule with quantum backend detected - using VQE-based excited states")
            print(f"📊 Computing multiple eigenstates via repeated VQE with orthogonality constraints")

            # For now, fall through to classical excited states
            # TODO: Implement proper VQE-based excited states with orthogonality constraints
            print(f"⚠️  VQE-based excited states not yet implemented - falling back to classical")
            backend_type = 'classical'  # Override to use classical method

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

    # Classical excited states (CIS, TDDFT)
    method = config.get('excited_method', 'cis')  # CIS, TDDFT
    n_states = config.get('n_states', 5)

    # Broadcast initial status to frontend
    if experiment_id:
        _broadcast_log_sync(experiment_id, f"🔬 Starting {method.upper()} excited states calculation...")
        _broadcast_log_sync(experiment_id, f"📊 Computing {n_states} excited states")

    solver = ExcitedStatesSolver(
        bond=bond,
        method=method,
        n_states=n_states,
        enable_analysis=True,
        enable_optimization=False
    )

    # Update progress (no initial broadcast - CIS/TDDFT complete instantly)
    JobDB.update_progress(job_id, progress=20.0, current_iteration=0)
    print(f"📊 Excited States: Computing states using {method.upper()} (classical)...")

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
            import traceback
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

        # Update experiment with results
        ExperimentDB.update_status(
            experiment_id,
            'completed',
            results=results
        )

        # Update job
        JobDB.update_progress(job_id, progress=100.0)
        JobDB.update_status(job_id, 'completed')

        # Broadcast completion status via WebSocket
        global _main_event_loop
        try:
            loop = _main_event_loop
            if loop is None:
                try:
                    loop = asyncio.get_running_loop()
                except RuntimeError:
                    pass

            if loop is not None:
                coro = ws_manager.broadcast_status(
                    experiment_id,
                    status='completed',
                    progress=100.0
                )
                future = asyncio.run_coroutine_threadsafe(coro, loop)
                future.result(timeout=0.5)
                print(f"✅ Broadcasted completion status via WebSocket")
            else:
                print(f"⚠️  No event loop available for completion broadcast")
        except Exception as e:
            print(f"⚠️  WebSocket completion broadcast failed: {e}")

    except ExperimentCancelledException as e:
        # Handle cancellation gracefully
        print(f"⚠️  Experiment cancelled: {e}")
        # Ensure both experiment and job are marked as cancelled
        ExperimentDB.update_status(experiment_id, 'cancelled')
        JobDB.update_status(job_id, 'cancelled')

    except Exception as e:
        print(f"❌ Experiment failed: {e}")
        traceback.print_exc()

        # Update experiment as failed
        ExperimentDB.update_status(
            experiment_id,
            'failed',
            error_message=str(e)
        )

        # Update job as failed
        JobDB.update_status(job_id, 'failed')
