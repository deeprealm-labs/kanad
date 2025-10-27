"""
Test classical excited states execution to see if convergence broadcasting works.
"""
import sys
import asyncio

# Add api to path
sys.path.insert(0, '/home/mk/deeprealm/kanad')
sys.path.insert(0, '/home/mk/deeprealm/kanad/api')

from kanad.core.molecule import Molecule
from api.services.experiment_service import execute_excited_states
from api.core.database import init_db, JobDB
from api.utils import set_main_event_loop
import uuid

async def test_excited_states():
    """Test excited states with classical backend."""
    # Set event loop
    loop = asyncio.get_running_loop()
    set_main_event_loop(loop)

    # Initialize database
    init_db()

    # Create molecule (H2)
    molecule = Molecule.from_smiles('[H][H]', basis='sto-3g')

    # Config
    config = {
        'backend': 'classical',
        'excited_method': 'cis',
        'n_states': 3,
    }

    # Create job
    job_id = str(uuid.uuid4())
    experiment_id = str(uuid.uuid4())

    job_data = {
        'id': job_id,
        'experiment_id': experiment_id,
        'status': 'running',
        'priority': 0,
        'max_iterations': 100,
    }
    JobDB.create(job_data)

    print(f"🧪 Testing excited states: experiment_id={experiment_id}")
    print(f"🧪 Job ID: {job_id}")
    print("")

    # Execute
    result = execute_excited_states(molecule, config, job_id, experiment_id)

    print("")
    print(f"✅ Result keys: {result.keys()}")
    print(f"✅ Ground state energy: {result.get('ground_state_energy', 'N/A')}")
    print(f"✅ Excited energies: {result.get('excited_state_energies', [])}")
    print(f"✅ Excitation energies (eV): {result.get('excitation_energies_ev', [])}")
    print(f"✅ Energy history: {result.get('energy_history', [])}")

if __name__ == '__main__':
    asyncio.run(test_excited_states())
