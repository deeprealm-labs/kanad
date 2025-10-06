"""
BlueQubit Cloud Backend for Kanad Framework

Provides access to BlueQubit's GPU and CPU quantum simulators.

Features:
- Free GPU simulators (36 qubits)
- CPU simulators (34 qubits)
- MPS tensor network simulators (40+ qubits)
- NVIDIA cuQuantum acceleration

Authentication:
    Set environment variable BLUEQUBIT_API_TOKEN or pass token directly

Usage:
    from kanad.backends.bluequbit import BlueQubitBackend

    backend = BlueQubitBackend(device='gpu')  # or 'cpu', 'mps.gpu'
    result = backend.run_vqe(bond, ansatz='ucc')
"""

from kanad.backends.bluequbit.backend import BlueQubitBackend
from kanad.backends.bluequbit.runner import BlueQubitRunner

__all__ = ['BlueQubitBackend', 'BlueQubitRunner']
