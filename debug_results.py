import os
from kanad.backends.ibm.backend import IBMBackend

IBM_API = os.getenv('IBM_API')
IBM_CRN = os.getenv('IBM_CRN')

backend = IBMBackend(
    backend_name='ibm_torino',
    api_token=IBM_API,
    channel='ibm_cloud' if IBM_CRN else 'ibm_quantum_platform',
    crn=IBM_CRN
)

job_id = 'd44motnlcjfs73atu3n0'
print(f"Getting results for job: {job_id}")

results = backend.get_job_result(job_id)

print(f"\nResults type: {type(results)}")
print(f"Results: {results}")

# Try different ways to access the data
try:
    print(f"\nlen(results): {len(results)}")
    for i, result in enumerate(results):
        print(f"\n  Result {i}:")
        print(f"    Type: {type(result)}")
        print(f"    Dir: {[x for x in dir(result) if not x.startswith('_')]}")
        if hasattr(result, 'data'):
            print(f"    data: {result.data}")
            if hasattr(result.data, 'evs'):
                print(f"    data.evs: {result.data.evs}")
            if hasattr(result.data, 'stds'):
                print(f"    data.stds: {result.data.stds}")
        if hasattr(result, 'metadata'):
            print(f"    metadata: {result.metadata}")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
