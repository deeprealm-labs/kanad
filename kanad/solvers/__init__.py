"""
Quantum Chemistry Solvers.

Provides quantum and hybrid quantum-classical algorithms for computing
molecular properties, ground state energies, and excited states.

Available Solvers:
- VQESolver: Variational Quantum Eigensolver (production-ready)
- SQDSolver: Sample-Based Quantum Diagonalization
- QPESolver: Quantum Phase Estimation
- FCISolver: Full Configuration Interaction (classical)
- ExcitedStatesSolver: Excited state calculations
- VibrationalSolver: Vibrational spectroscopy
- AlloyFormationSolver: Metallic alloy properties
- ProteinFoldingSolver: Protein structure prediction
"""

from kanad.solvers.vqe_solver import VQESolver
from kanad.solvers.sqd_solver import SQDSolver
from kanad.solvers.qpe_solver import QPESolver
from kanad.solvers.fci_solver import FCISolver
from kanad.solvers.excited_states_solver import ExcitedStatesSolver
from kanad.solvers.vibrational_solver import VibrationalSolver
from kanad.solvers.alloy_solver import AlloySolver
from kanad.solvers.protein_folding_solver import ProteinFoldingSolver

__all__ = [
    'VQESolver',
    'SQDSolver',
    'QPESolver',
    'FCISolver',
    'ExcitedStatesSolver',
    'VibrationalSolver',
    'AlloySolver',
    'ProteinFoldingSolver',
]
