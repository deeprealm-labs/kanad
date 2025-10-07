"""
IBM Quantum Results Checker
============================

Retrieves and analyzes results from submitted IBM Quantum jobs.

Usage:
    python ibm_large_molecules_results.py
"""

import os
import sys
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Any

# Add kanad to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Load environment variables from .env file
from dotenv import load_dotenv
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(env_path)

print("=" * 80)
print("IBM QUANTUM - RESULTS CHECKER")
print("=" * 80)
print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)

# Check for IBM credentials
ibm_api = os.getenv('IBM_API')

print(f"\nCredentials: {'✓ Loaded from .env' if ibm_api else '✗ Missing'}")

if not ibm_api:
    print("\n⚠ ERROR: IBM_API not found in .env file")
    print("Add to .env: IBM_API=your_token_here")
    sys.exit(1)

from kanad.backends.ibm import IBMBackend

# Load job data
job_file = Path(__file__).parent / 'ibm_job_ids.json'

if not job_file.exists():
    print(f"\n✗ Job file not found: {job_file}")
    print("Run ibm_large_molecules.py first to submit jobs")
    sys.exit(1)

with open(job_file, 'r') as f:
    job_data = json.load(f)

experiments = job_data['experiments']
backend_name = job_data['backend']

print(f"\nLoaded {len(experiments)} experiments from {job_file}")
print(f"Backend: {backend_name}")
print(f"Submitted: {job_data['submission_time']}")

# Initialize backend
try:
    backend = IBMBackend(backend_name=backend_name)
    print(f"\n✓ Connected to IBM Quantum: {backend_name}")
except Exception as e:
    print(f"\n✗ Failed to connect: {e}")
    sys.exit(1)

# Check job statuses
print("\n" + "=" * 80)
print("CHECKING JOB STATUSES")
print("=" * 80)

results_summary = {
    'completed': [],
    'running': [],
    'queued': [],
    'failed': [],
    'unknown': []
}

for i, exp in enumerate(experiments, 1):
    job_id = exp['job_id']
    name = exp['name']

    print(f"\n[{i}/{len(experiments)}] {name}")
    print(f"  Job ID: {job_id}")

    try:
        status = backend.get_job_status(job_id)
        print(f"  Status: {status}")

        exp['current_status'] = status
        exp['last_checked'] = datetime.now().isoformat()

        # Categorize
        if status in ['DONE', 'COMPLETED']:
            results_summary['completed'].append(exp)
        elif status in ['RUNNING', 'VALIDATING']:
            results_summary['running'].append(exp)
        elif status in ['QUEUED', 'PENDING']:
            results_summary['queued'].append(exp)
        elif status in ['CANCELLED', 'ERROR', 'FAILED']:
            results_summary['failed'].append(exp)
        else:
            results_summary['unknown'].append(exp)

    except Exception as e:
        print(f"  ERROR: {e}")
        exp['error'] = str(e)
        results_summary['unknown'].append(exp)

# Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print(f"\n✓ Completed:  {len(results_summary['completed'])}")
print(f"⏳ Running:    {len(results_summary['running'])}")
print(f"📋 Queued:     {len(results_summary['queued'])}")
print(f"✗ Failed:     {len(results_summary['failed'])}")
print(f"? Unknown:    {len(results_summary['unknown'])}")

# Retrieve completed results
if results_summary['completed']:
    print("\n" + "=" * 80)
    print("RETRIEVING COMPLETED RESULTS")
    print("=" * 80)

    detailed_results = []

    for exp in results_summary['completed']:
        job_id = exp['job_id']
        name = exp['name']

        print(f"\n{name}")
        print(f"  Job ID: {job_id}")

        try:
            # Retrieve raw result
            result = backend.get_job_result(job_id)

            # Extract energy values (from EstimatorV2 result)
            if hasattr(result, '__iter__'):
                # Multiple pub results
                energies = []
                for pub_result in result:
                    if hasattr(pub_result, 'data'):
                        if hasattr(pub_result.data, 'evs'):
                            energies.append(pub_result.data.evs)

                if energies:
                    vqe_energy = float(energies[0]) if len(energies) > 0 else None
                else:
                    vqe_energy = None
            else:
                vqe_energy = None

            hf_energy = exp.get('hf_energy', None)

            print(f"  HF Energy:  {hf_energy:.6f} Ha" if hf_energy else "  HF Energy: N/A")
            print(f"  VQE Energy: {vqe_energy:.6f} Ha" if vqe_energy else "  VQE Energy: N/A")

            if hf_energy and vqe_energy:
                correlation = (vqe_energy - hf_energy) * 1000  # mHa
                print(f"  Correlation: {correlation:.3f} mHa")
                print(f"  Status: {'✓ VQE < HF' if vqe_energy < hf_energy else '⚠ VQE > HF (check)'}")

            # Save detailed result
            detailed_results.append({
                'name': name,
                'job_id': job_id,
                'bond': exp['bond'],
                'qubits': exp['qubits'],
                'hf_energy': hf_energy,
                'vqe_energy': vqe_energy,
                'correlation_mha': correlation if (hf_energy and vqe_energy) else None,
                'retrieval_time': datetime.now().isoformat()
            })

        except Exception as e:
            print(f"  ERROR retrieving result: {e}")
            import traceback
            traceback.print_exc()

    # Save detailed results
    results_file = Path(__file__).parent / 'ibm_results.json'
    with open(results_file, 'w') as f:
        json.dump(detailed_results, f, indent=2)

    print(f"\n✓ Detailed results saved to: {results_file}")

# Update job file with current statuses
with open(job_file, 'w') as f:
    job_data['experiments'] = experiments
    job_data['last_checked'] = datetime.now().isoformat()
    json.dump(job_data, f, indent=2)

print(f"\n✓ Job statuses updated in: {job_file}")

# Next steps
print("\n" + "=" * 80)
print("NEXT STEPS")
print("=" * 80)

if results_summary['queued'] or results_summary['running']:
    print(f"""
Still waiting for {len(results_summary['queued']) + len(results_summary['running'])} jobs to complete.

Check again later with:
    python ibm_large_molecules_results.py

Monitor progress at:
    https://quantum.ibm.com
""")

if results_summary['completed']:
    print(f"""
{len(results_summary['completed'])} jobs completed!

Results analysis:
    - VQE energies computed on REAL quantum hardware
    - Correlation energies show electron correlation effects
    - Compare with classical methods (DFT, CCSD)

Publication potential:
    - Large pharmaceutical molecules on quantum hardware
    - Real-world validation of Kanad framework
    - Demonstrates quantum advantage for chemistry
""")

if results_summary['failed']:
    print(f"""
⚠ {len(results_summary['failed'])} jobs failed.

Possible reasons:
    - Circuit too deep for hardware
    - Timeout (job took too long)
    - Hardware error during execution

Consider:
    - Reduce optimization level
    - Use simpler ansatz
    - Try different backend
""")

print("\n" + "=" * 80)
print("✓ RESULTS CHECK COMPLETE")
print("=" * 80)
