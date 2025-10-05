"""
IBM Quantum Platform integration.

Provides complete quantum chemistry solver suite optimized for IBM Quantum:
- IBMRuntimeBackend: Core backend with Runtime primitives
- IBMVQESolver: Variational Quantum Eigensolver
- IBMQPESolver: Quantum Phase Estimation
- IBMSQDSolver: Sample-based Quantum Diagonalization

All solvers include:
- Automatic transpilation with optimization levels
- Error mitigation (resilience levels)
- Session management for efficient execution
- Support for both simulators and real hardware
"""

from kanad.backends.ibm.backend import IBMRuntimeBackend
from kanad.backends.ibm.vqe_solver import IBMVQESolver
from kanad.backends.ibm.qpe_solver import IBMQPESolver
from kanad.backends.ibm.sqd_solver import IBMSQDSolver

__all__ = [
    'IBMRuntimeBackend',
    'IBMVQESolver',
    'IBMQPESolver',
    'IBMSQDSolver',
]
