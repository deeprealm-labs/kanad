"""
Quantum backends and solver modules.

Modular architecture organized by provider:
- kanad.backends.ibm: IBM Quantum Platform (VQE, QPE, SQD solvers)
- kanad.backends.bluequbit: BlueQubit Cloud
- kanad.backends.qiskit_backend: Generic Qiskit backends
- kanad.backends.cuquantum_backend: GPU-accelerated simulation

Usage:
    # New modular way (recommended):
    from kanad.backends.ibm import IBMVQESolver, IBMRuntimeBackend

    # Legacy imports (backward compatible):
    from kanad.backends import IBMRuntimeBackend, QiskitBackend
"""

from kanad.backends.qiskit_backend import QiskitBackend
from kanad.backends.ibm.backend import IBMRuntimeBackend  # Backward compatibility
from kanad.backends.cuquantum_backend import (
    CuQuantumBackend,
    check_cuquantum_available,
    get_gpu_info
)

# Import submodules for easy access
from kanad.backends import ibm
from kanad.backends import bluequbit

__all__ = [
    'QiskitBackend',
    'IBMRuntimeBackend',
    'CuQuantumBackend',
    'check_cuquantum_available',
    'get_gpu_info',
    'ibm',
    'bluequbit',
]
