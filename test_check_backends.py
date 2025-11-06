import os
from qiskit_ibm_runtime import QiskitRuntimeService

IBM_API = os.getenv('IBM_API')
IBM_CRN = os.getenv('IBM_CRN')

# Save account
QiskitRuntimeService.save_account(
    channel='ibm_cloud',
    token=IBM_API,
    instance=IBM_CRN,
    overwrite=True
)

# Initialize service
service = QiskitRuntimeService(channel='ibm_cloud')

# List all backends
print("Available backends:")
backends = service.backends()
for backend in backends:
    print(f"  - {backend.name}: {backend.num_qubits} qubits, simulator={backend.simulator}")
