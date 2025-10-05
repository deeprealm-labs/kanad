"""
VQE Optimization Module

Techniques to speed up VQE execution:
- Circuit depth reduction
- Gate count optimization
- Faster convergence
- Evaluation caching
"""

from kanad.vqe_optimization.circuit_optimizer import CircuitOptimizer
from kanad.vqe_optimization.fast_vqe import FastVQE

__all__ = ['CircuitOptimizer', 'FastVQE']
