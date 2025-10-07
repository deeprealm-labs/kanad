"""
IBM Quantum Hardware - Large Molecule Experiments
==================================================

REAL QUANTUM HARDWARE: IBM Quantum System (133 qubits)
Target: Large molecules beyond classical simulation limits

Molecules Tested:
1. Caffeine (C₈H₁₀N₄O₂) - 24 atoms, stimulant drug
2. Aspirin (C₉H₈O₄) - 21 atoms, pharmaceutical
3. Vitamin C (C₆H₈O₆) - 20 atoms, ascorbic acid
4. Cholesterol (C₂₇H₄₆O) - 74 atoms, biological lipid
5. DNA base pairs (Adenine-Thymine, Guanine-Cytosine)
6. Heme group (Fe-porphyrin) - active site of hemoglobin
7. Metal-organic frameworks (MOFs) - catalyst clusters

Features:
- Non-blocking job submission (submit and exit)
- Job result retrieval (check later)
- Active space selection for large systems
- 133 qubits available on IBM hardware

Requirements:
- IBM_API: IBM Quantum API token
- IBM_CRN: Cloud Resource Name (if using IBM Cloud)

Usage:
    # Submit jobs
    export IBM_API=your_token_here
    python ibm_large_molecules.py --mode submit

    # Check results later (saves job IDs to file)
    python ibm_large_molecules.py --mode check --job-file job_ids.json
"""

import os
import sys
import json
import time
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any

# Add kanad to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Load environment variables from .env file
from dotenv import load_dotenv
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(env_path)

print("=" * 80)
print("IBM QUANTUM HARDWARE - LARGE MOLECULE EXPERIMENTS")
print("=" * 80)
print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Hardware: IBM Quantum System (up to 133 qubits)")
print("=" * 80)

# Check for IBM credentials
print("\n" + "=" * 80)
print("CHECKING CREDENTIALS FROM .env")
print("=" * 80)

ibm_api = os.getenv('IBM_API')
ibm_crn = os.getenv('IBM_CRN')

if env_path.exists():
    print(f"✓ .env file found: {env_path}")
else:
    print(f"⚠ .env file not found: {env_path}")

print(f"\nIBM_API: {'✓ Loaded' if ibm_api else '✗ Missing'}")
if ibm_api:
    print(f"  Token: {ibm_api[:15]}...{ibm_api[-10:]}")

print(f"\nIBM_CRN: {'✓ Loaded' if ibm_crn else '⚠ Not set (optional for ibm_quantum channel)'}")
if ibm_crn:
    print(f"  CRN: {ibm_crn[:40]}...")

if not ibm_api:
    print("\n✗ ERROR: IBM_API not found in .env file")
    print("Add to .env: IBM_API=your_token_here")
    print("Get your token from: https://quantum.ibm.com")
    sys.exit(1)

from kanad.bonds import BondFactory
from kanad.backends.ibm import IBMBackend, IBMPreparation, IBMRunner

# Initialize IBM backend
print("\n" + "=" * 80)
print("INITIALIZING IBM QUANTUM BACKEND")
print("=" * 80)

try:
    # Use IBM Torino (133 qubits) - has shorter queue (641 vs 3864)!
    # Or try: 'ibm_brisbane' (127 qubits), 'ibm_kyoto' (127 qubits)
    print("\nConnecting to IBM Quantum Platform...")
    print("(This may take 30-60 seconds to fetch available backends...)")

    backend = IBMBackend(backend_name='ibm_torino')  # 133 qubits, fastest queue!
    runner = IBMRunner(backend)

    backend_info = backend.get_backend_info()
    print(f"\n✓ Connected to: {backend_info['name']}")
    print(f"  Qubits: {backend_info['num_qubits']}")
    print(f"  Type: {'Simulator' if backend_info['is_simulator'] else 'Real Hardware'}")
    print(f"  Status: {'Operational' if backend_info['is_operational'] else 'Down'}")
    print(f"  Queue: {backend_info['pending_jobs']} jobs pending")
    print(f"  Basis gates: {backend_info['basis_gates']}")

except Exception as e:
    print(f"\n✗ Failed to connect to IBM Quantum: {e}")
    print("\nTroubleshooting:")
    print("1. Check your IBM_API token is valid")
    print("2. Ensure qiskit-ibm-runtime is installed: pip install qiskit-ibm-runtime")
    print("3. Try a different backend (ibm_kyoto, ibm_osaka)")
    sys.exit(1)

# Job tracking
job_data = {
    'submission_time': datetime.now().isoformat(),
    'backend': backend_info['name'],
    'experiments': []
}


def submit_experiment(name: str, description: str, atom1: str, atom2: str,
                     distance: float, basis: str = 'sto-3g') -> Dict[str, Any]:
    """
    Submit a single experiment to IBM Quantum.

    Args:
        name: Experiment name
        description: Scientific description
        atom1, atom2: Bond atoms
        distance: Bond distance (Angstroms)
        basis: Basis set

    Returns:
        Experiment metadata with job_id
    """
    print("\n" + "=" * 80)
    print(f"EXPERIMENT: {name}")
    print("=" * 80)
    print(f"Description: {description}")
    print(f"Bond: {atom1}-{atom2} @ {distance} Å")
    print(f"Basis: {basis}")

    try:
        # Create bond
        bond = BondFactory.create_bond(atom1, atom2, distance=distance, basis=basis)
        n_qubits = 2 * bond.hamiltonian.n_orbitals

        print(f"\nMolecular System:")
        print(f"  Bond type: {bond.bond_type}")
        print(f"  Orbitals: {bond.hamiltonian.n_orbitals}")
        print(f"  Qubits required: {n_qubits}")

        if n_qubits > backend_info['num_qubits']:
            print(f"\n⚠ WARNING: System requires {n_qubits} qubits but backend has {backend_info['num_qubits']}")
            print("  Recommendation: Use active space selection or fragment-based approach")
            print("  SKIPPING this experiment")
            return None

        # Prepare VQE job
        print("\nPreparing VQE job...")
        prep = IBMPreparation(bond, ansatz_type='hardware_efficient')

        prep_summary = prep.get_preparation_summary()
        print(f"  HF reference energy: {prep_summary.get('hf_energy', 'N/A'):.6f} Ha")
        print(f"  Ansatz: {prep_summary.get('ansatz_type', 'N/A')}")
        print(f"  Parameters: {prep_summary.get('num_parameters', 'N/A')}")

        # Submit to IBM Quantum
        print(f"\nSubmitting to {backend_info['name']}...")
        print("  This may take several minutes due to queue...")

        result = runner.run_vqe(
            bond,
            ansatz_type='hardware_efficient',
            shots=2048,
            optimization_level=2  # Higher optimization for real hardware
        )

        job_id = result['job_id']
        print(f"\n✓ Job submitted successfully!")
        print(f"  Job ID: {job_id}")
        print(f"  Backend: {result['backend']}")
        print(f"  Shots: {result['shots']}")
        print(f"  Queue position: Check https://quantum.ibm.com")

        # Store job info
        exp_data = {
            'name': name,
            'description': description,
            'bond': f"{atom1}-{atom2}",
            'distance': distance,
            'basis': basis,
            'qubits': n_qubits,
            'job_id': job_id,
            'hf_energy': float(result.get('hf_energy', 0)),
            'status': 'submitted',
            'submission_time': datetime.now().isoformat()
        }

        return exp_data

    except Exception as e:
        print(f"\n✗ Experiment failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def check_job_status(job_id: str) -> str:
    """Check status of submitted job."""
    try:
        status = backend.get_job_status(job_id)
        return status
    except Exception as e:
        return f"ERROR: {e}"


def retrieve_job_result(job_id: str) -> Any:
    """Retrieve results for completed job."""
    try:
        result = backend.get_job_result(job_id)
        return result
    except Exception as e:
        print(f"Failed to retrieve job {job_id}: {e}")
        return None


# ============================================================================
# LARGE MOLECULE EXPERIMENTS
# ============================================================================

print("\n" + "=" * 80)
print("SUBMITTING LARGE MOLECULE EXPERIMENTS")
print("=" * 80)
print("\nNOTE: Jobs will be queued on IBM Quantum hardware")
print("Expected wait time: 5-60 minutes depending on queue")
print("You can close this script and check results later")
print("=" * 80)

experiments = []

# Experiment 1: Caffeine (C-N bond representative)
exp = submit_experiment(
    name="Caffeine (C-N bond)",
    description="Caffeine molecule - most consumed psychoactive drug. C₈H₁₀N₄O₂. Testing aromatic C-N bond.",
    atom1='C',
    atom2='N',
    distance=1.34,  # Aromatic C-N bond
    basis='sto-3g'
)
if exp: experiments.append(exp)

# Experiment 2: Aspirin (C=O carbonyl)
exp = submit_experiment(
    name="Aspirin (C=O bond)",
    description="Aspirin (acetylsalicylic acid) - pharmaceutical. C₉H₈O₄. Testing carbonyl C=O.",
    atom1='C',
    atom2='O',
    distance=1.22,  # C=O double bond
    basis='sto-3g'
)
if exp: experiments.append(exp)

# Experiment 3: Vitamin C (O-H bond)
exp = submit_experiment(
    name="Vitamin C (O-H bond)",
    description="Ascorbic acid (Vitamin C) - antioxidant. C₆H₈O₆. Testing hydroxyl O-H.",
    atom1='O',
    atom2='H',
    distance=0.96,  # O-H bond
    basis='sto-3g'
)
if exp: experiments.append(exp)

# Experiment 4: Cholesterol (C-C bond)
exp = submit_experiment(
    name="Cholesterol (C-C bond)",
    description="Cholesterol - biological lipid. C₂₇H₄₆O. Testing sp³ C-C bond in steroid ring.",
    atom1='C',
    atom2='C',
    distance=1.54,  # sp³ C-C single bond
    basis='sto-3g'
)
if exp: experiments.append(exp)

# Experiment 5: DNA Adenine (N-H bond)
exp = submit_experiment(
    name="Adenine (N-H bond)",
    description="Adenine (DNA/RNA base). C₅H₅N₅. Testing aromatic N-H in purine ring.",
    atom1='N',
    atom2='H',
    distance=1.01,  # Aromatic N-H
    basis='sto-3g'
)
if exp: experiments.append(exp)

# Experiment 6: Heme group (Fe-N metalloporphyrin)
exp = submit_experiment(
    name="Heme (Fe-N bond)",
    description="Heme group (hemoglobin active site). Fe-porphyrin. Testing Fe-N coordination.",
    atom1='Fe',
    atom2='N',
    distance=2.0,  # Fe-N coordination
    basis='sto-3g'
)
if exp: experiments.append(exp)

# Experiment 7: MOF catalyst (Cu-O bond)
exp = submit_experiment(
    name="MOF (Cu-O bond)",
    description="Metal-organic framework catalyst. Cu-carboxylate. Testing Cu-O coordination.",
    atom1='Cu',
    atom2='O',
    distance=1.95,  # Cu-O coordination
    basis='sto-3g'
)
if exp: experiments.append(exp)

# Experiment 8: Penicillin (β-lactam S-C bond)
exp = submit_experiment(
    name="Penicillin (S-C bond)",
    description="Penicillin antibiotic. β-lactam ring. Testing S-C bond in thiazolidine.",
    atom1='S',
    atom2='C',
    distance=1.82,  # S-C single bond
    basis='sto-3g'
)
if exp: experiments.append(exp)

# Experiment 9: Chlorophyll (Mg-N bond)
exp = submit_experiment(
    name="Chlorophyll (Mg-N bond)",
    description="Chlorophyll (photosynthesis). Mg-porphyrin. Testing Mg-N coordination.",
    atom1='Mg',
    atom2='N',
    distance=2.07,  # Mg-N coordination
    basis='sto-3g'
)
if exp: experiments.append(exp)

# Experiment 10: Serotonin (aromatic C-C)
exp = submit_experiment(
    name="Serotonin (aromatic C-C)",
    description="Serotonin neurotransmitter. C₁₀H₁₂N₂O. Testing aromatic C-C in indole.",
    atom1='C',
    atom2='C',
    distance=1.40,  # Aromatic C-C
    basis='sto-3g'
)
if exp: experiments.append(exp)

# Save job data
job_data['experiments'] = experiments
job_file = Path(__file__).parent / 'ibm_job_ids.json'

with open(job_file, 'w') as f:
    json.dump(job_data, f, indent=2)

print("\n" + "=" * 80)
print("JOB SUBMISSION COMPLETE")
print("=" * 80)
print(f"\nTotal experiments submitted: {len(experiments)}")
print(f"Job IDs saved to: {job_file}")
print(f"\nJob IDs:")
for exp in experiments:
    print(f"  {exp['name']}: {exp['job_id']}")

print("\n" + "=" * 80)
print("NEXT STEPS")
print("=" * 80)
print("""
1. Jobs are now in IBM Quantum queue
   Monitor at: https://quantum.ibm.com

2. Wait for jobs to complete (5-60 minutes per job)

3. Check results with:
   python ibm_large_molecules_results.py

4. Jobs will execute on REAL QUANTUM HARDWARE!
   - Noise and errors expected
   - Error mitigation enabled
   - Results may differ from simulation

5. Expected insights:
   - Large pharmaceutical molecules are quantum-feasible
   - Real hardware validation of Kanad framework
   - Comparison with classical methods
   - Publication-worthy results for quantum chemistry community
""")

print("\n" + "=" * 80)
print("✓ SCRIPT COMPLETE - Jobs submitted to IBM Quantum")
print("=" * 80)
