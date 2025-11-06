"""
IBM Quantum Backend for Kanad Framework

Provides access to IBM Quantum hardware and simulators via Qiskit Runtime.

Features:
- Batch mode job submission (for non-premium users)
- Quantum hardware access
- Cloud simulators
- Qiskit Runtime primitives (Sampler, Estimator)

Authentication:
    Set environment variables:
    - IBM_API: Your IBM Quantum API token
    - IBM_CRN: Cloud Resource Name (for IBM Cloud channel)

    Channels: 'ibm_quantum_platform' (default) or 'ibm_cloud'

Usage:
    from kanad.backends.ibm import IBMBackend, IBMPreparation

    # Prepare job
    prep = IBMPreparation(bond, ansatz='ucc')
    circuits = prep.prepare_circuits()

    # Execute on IBM
    backend = IBMBackend(backend_name='ibm_brisbane')
    results = backend.run_batch(circuits, prep)
"""

from kanad.backends.ibm.backend import IBMBackend
from kanad.backends.ibm.preparation import IBMPreparation
from kanad.backends.ibm.runner import IBMRunner
from kanad.backends.ibm.error_mitigation import ErrorMitigationStrategy
from kanad.backends.ibm.governance_error_mitigation import GovernanceAwareErrorMitigation

__all__ = [
    'IBMBackend',
    'IBMPreparation',
    'IBMRunner',
    'ErrorMitigationStrategy',
    'GovernanceAwareErrorMitigation'
]
