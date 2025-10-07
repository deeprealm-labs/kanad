"""
IBM Quantum Results Retrieval
==============================

Simple script to check job statuses and retrieve completed results.

Usage:
    # Check all jobs from job file
    python get_ibm_results.py

    # Check specific job ID
    python get_ibm_results.py --job-id d3ie3mpb641c738n3d80
"""

import os
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime
from dotenv import load_dotenv

# Load credentials
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(env_path)

sys.path.insert(0, str(Path(__file__).parent.parent))

print("=" * 80)
print("IBM QUANTUM - RESULTS RETRIEVAL")
print("=" * 80)
print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--job-id', help='Specific job ID to check')
parser.add_argument('--job-file', default='ibm_job_ids.json', help='Job IDs file')
args = parser.parse_args()

# Check credentials
ibm_api = os.getenv('IBM_API')
if not ibm_api:
    print("✗ IBM_API not found in .env")
    sys.exit(1)

print("✓ Credentials loaded\n")

# Connect to IBM Quantum
print("Connecting to IBM Quantum...")

try:
    from qiskit_ibm_runtime import QiskitRuntimeService

    service = QiskitRuntimeService(
        channel='ibm_quantum_platform',
        token=ibm_api
    )

    print("✓ Connected!\n")

    # Check specific job or all jobs from file
    if args.job_id:
        job_ids = [args.job_id]
        job_names = {args.job_id: "Manual check"}
    else:
        job_file = Path(__file__).parent / args.job_file

        if not job_file.exists():
            print(f"✗ Job file not found: {job_file}")
            print("Run ibm_large_molecules.py first or provide --job-id")
            sys.exit(1)

        with open(job_file, 'r') as f:
            job_data = json.load(f)

        experiments = job_data['experiments']
        job_ids = [exp['job_id'] for exp in experiments]
        job_names = {exp['job_id']: exp['name'] for exp in experiments}

        print(f"Loaded {len(job_ids)} jobs from {job_file}\n")

    # Check each job
    print("=" * 80)
    print("CHECKING JOB STATUSES")
    print("=" * 80)

    results = []

    for i, job_id in enumerate(job_ids, 1):
        name = job_names.get(job_id, "Unknown")

        print(f"\n[{i}/{len(job_ids)}] {name}")
        print(f"  Job ID: {job_id}")

        try:
            job = service.job(job_id)
            status = job.status()

            print(f"  Status: {status}")

            result_data = {
                'job_id': job_id,
                'name': name,
                'status': str(status),
                'checked_at': datetime.now().isoformat()
            }

            # If completed, retrieve results
            status_str = str(status)
            if 'DONE' in status_str or 'COMPLETED' in status_str:
                print(f"  ✓ Job completed! Retrieving results...")

                try:
                    result = job.result()

                    # Extract data based on primitive type
                    # For SamplerV2, results are measurement counts
                    for pub_result in result:
                        if hasattr(pub_result.data, 'meas'):
                            counts = pub_result.data.meas.get_counts()
                            print(f"  Measurement counts: {counts}")
                            result_data['counts'] = counts
                        elif hasattr(pub_result.data, 'evs'):
                            # EstimatorV2 energy values
                            energy = float(pub_result.data.evs)
                            print(f"  Energy: {energy:.6f} Ha")
                            result_data['energy'] = energy

                    # Save metrics
                    metrics = job.metrics()
                    if metrics:
                        result_data['execution_time'] = metrics.get('execution_time')
                        result_data['queue_time'] = metrics.get('estimated_start_time')

                    print(f"  ✓ Results retrieved successfully")

                except Exception as e:
                    print(f"  ⚠ Could not retrieve result data: {e}")
                    result_data['result_error'] = str(e)

            elif 'QUEUED' in status_str or 'PENDING' in status_str or 'VALIDATING' in status_str:
                print(f"  ⏳ Still in queue...")

            elif 'RUNNING' in status_str:
                print(f"  ⚡ Currently executing...")

            elif 'CANCELLED' in status_str or 'ERROR' in status_str or 'FAILED' in status_str:
                print(f"  ✗ Job failed")
                try:
                    error = job.error_message()
                except:
                    error = "No error message available"
                if error:
                    print(f"  Error: {error}")
                    result_data['error'] = error

            results.append(result_data)

        except Exception as e:
            print(f"  ✗ Error checking job: {e}")
            results.append({
                'job_id': job_id,
                'name': name,
                'error': str(e),
                'checked_at': datetime.now().isoformat()
            })

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    completed = sum(1 for r in results if r.get('status', '').startswith('JobStatus.DONE'))
    queued = sum(1 for r in results if 'QUEUED' in r.get('status', '') or 'PENDING' in r.get('status', ''))
    running = sum(1 for r in results if 'RUNNING' in r.get('status', ''))
    failed = sum(1 for r in results if 'FAILED' in r.get('status', '') or 'CANCELLED' in r.get('status', '') or 'ERROR' in r.get('status', ''))

    print(f"\n✓ Completed:  {completed}/{len(results)}")
    print(f"⏳ Queued:     {queued}/{len(results)}")
    print(f"⚡ Running:    {running}/{len(results)}")
    print(f"✗ Failed:     {failed}/{len(results)}")

    # Save results
    results_file = Path(__file__).parent / 'ibm_results.json'
    with open(results_file, 'w') as f:
        json.dump({
            'checked_at': datetime.now().isoformat(),
            'total_jobs': len(results),
            'completed': completed,
            'queued': queued,
            'running': running,
            'failed': failed,
            'results': results
        }, f, indent=2)

    print(f"\n✓ Results saved to: {results_file}")

    # Show completed energies
    if completed > 0:
        print("\n" + "=" * 80)
        print("COMPLETED JOB ENERGIES")
        print("=" * 80)

        for r in results:
            if r.get('status', '').startswith('JobStatus.DONE'):
                name = r['name']
                energy = r.get('primary_energy')
                if energy:
                    print(f"{name}: {energy:.6f} Ha")
                else:
                    print(f"{name}: No energy data")

    print("\n" + "=" * 80)
    print("✓ CHECK COMPLETE")
    print("=" * 80)

    if queued or running:
        print(f"\nℹ {queued + running} jobs still pending. Run this script again later.")

    if completed:
        print(f"\n✓ {completed} jobs completed! Results saved to {results_file}")

except Exception as e:
    print(f"\n✗ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
