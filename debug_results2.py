import os
import numpy as np
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
print(f"Getting results for job: {job_id}\n")

results = backend.get_job_result(job_id)

print("Configuration Energies from IBM Torino:")
print("="*60)

for i, result in enumerate(results):
    print(f"\nConfiguration {i}:")
    print(f"  evs (final): {result.data.evs}")
    print(f"  stds (error): {result.data.stds}")
    
    if hasattr(result.data, 'evs_noise_factors'):
        print(f"  evs_noise_factors: {result.data.evs_noise_factors}")
    
    if hasattr(result.data, 'evs_extrapolated'):
        print(f"  evs_extrapolated shape: {result.data.evs_extrapolated.shape}")
        print(f"  evs_extrapolated: {result.data.evs_extrapolated}")
        
        # The extrapolated values might be in here
        if result.data.evs_extrapolated.size > 0:
            flat = result.data.evs_extrapolated.flatten()
            print(f"  Extrapolated values (flattened): {flat}")
            print(f"  Using extrapolated[0]: {flat[0]}")

# Try to get the actual energies
print("\n\n" + "="*60)
print("ATTEMPTING TO EXTRACT ACTUAL ENERGIES")
print("="*60)

energies = []
for i, result in enumerate(results):
    # Try different ways
    if not np.isnan(result.data.evs):
        energy = result.data.evs
    elif hasattr(result.data, 'evs_extrapolated') and result.data.evs_extrapolated.size > 0:
        # Use the first extrapolated value (noise factor = 0)
        energy = result.data.evs_extrapolated.flatten()[0]
    elif hasattr(result.data, 'evs_noise_factors') and not np.all(np.isnan(result.data.evs_noise_factors)):
        # Use the first noise factor measurement (noise = 1.0)
        energy = result.data.evs_noise_factors[0]
    else:
        energy = np.nan
    
    energies.append(energy)
    print(f"Config {i}: {energy:.8f} Ha")

print(f"\nAll energies: {energies}")
