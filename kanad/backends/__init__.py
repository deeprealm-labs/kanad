"""
Backend module for Qiskit integration.

Provides interfaces for executing quantum circuits on:
- Classical simulators (NumPy, Aer)
- GPU-accelerated simulators (cuQuantum)
- Real quantum hardware (IBM Quantum)
- Third-party providers
"""

from kanad.backends.qiskit_backend import QiskitBackend
from kanad.backends.cuquantum_backend import (
    CuQuantumBackend,
    check_cuquantum_available,
    get_gpu_info
)

__all__ = [
    'QiskitBackend',
    'CuQuantumBackend',
    'check_cuquantum_available',
    'get_gpu_info',
]
