"""
Check IBM Quantum Job Status

Usage:
    export IBM_API='your_token'
    export IBM_CRN='your_crn'
    python check_ibm_jobs.py [job_id]

If no job_id is provided, lists all recent jobs.
"""
import os
import sys
from kanad.backends.ibm.backend import IBMBackend

IBM_API = os.getenv('IBM_API')
IBM_CRN = os.getenv('IBM_CRN')

if not IBM_API:
    print("ERROR: IBM_API not set")
    exit(1)

print("="*80)
print("IBM QUANTUM JOB STATUS CHECKER")
print("="*80)

# Initialize backend (use torino by default, or specify in args)
backend = IBMBackend(
    backend_name='ibm_torino',
    api_token=IBM_API,
    channel='ibm_cloud' if IBM_CRN else 'ibm_quantum_platform',
    crn=IBM_CRN
)

print(f"\n✓ Connected to: {backend.backend.name}")

if len(sys.argv) > 1:
    # Check specific job
    job_id = sys.argv[1]

    print(f"\n📋 Checking job: {job_id}")

    try:
        status = backend.get_job_status(job_id)
        print(f"   Status: {status}")

        if status == 'DONE':
            print(f"\n✅ Job completed! Retrieving results...")
            results = backend.get_job_result(job_id)

            print(f"\n📊 Results:")
            print(f"   Values: {results.values}")
            print(f"\n   These are the diagonal energies for each configuration.")
            print(f"   Use these with classical diagonalization to get ground state.")

        elif status == 'ERROR':
            print(f"\n✗ Job failed!")
            try:
                job = backend.service.job(job_id)
                print(f"   Error details: {job.error_message()}")
            except:
                print(f"   Could not retrieve error details")

        elif status == 'CANCELLED':
            print(f"\n⏭️  Job was cancelled")

        else:
            print(f"\n⏳ Job is still {status}")
            print(f"   Check back later...")

    except Exception as e:
        print(f"\n✗ Error checking job: {e}")

else:
    # List all recent jobs
    print(f"\n📋 Recent jobs on {backend.backend.name}:")

    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
        service = backend.service

        # Get recent jobs
        jobs = service.jobs(backend_name='ibm_torino', limit=10)

        print(f"\n{'Job ID':<50} {'Status':<15} {'Created'}")
        print("-"*80)

        for job in jobs:
            job_id = job.job_id()
            status = str(job.status())
            created = job.creation_date.strftime('%Y-%m-%d %H:%M:%S') if hasattr(job, 'creation_date') else 'N/A'

            print(f"{job_id:<50} {status:<15} {created}")

        print(f"\n💡 To check a specific job:")
        print(f"   python check_ibm_jobs.py <job_id>")

    except Exception as e:
        print(f"\n✗ Error listing jobs: {e}")

print("\n" + "="*80)
