"""
Quantum Chemistry Solvers - Rebuilt with Bonds Module Integration.

New Architecture:
- All solvers work with bonds module as primary interface
- Automatic integration with analysis tools
- Automatic integration with optimization tools
- Rich, comprehensive results
- Unified API across all solvers

Available Solvers:
- VQESolver: Variational Quantum Eigensolver (ground state)
- SQDSolver: Subspace Quantum Diagonalization (ground + excited)
- ExcitedStatesSolver: Molecular excited states (CIS, TDDFT)
- FCISolver: Full Configuration Interaction (exact, legacy)
- QPESolver: Quantum Phase Estimation (legacy)

Usage Example:
    from kanad.bonds import BondFactory
    from kanad.solvers import VQESolver

    # Create molecule
    bond = BondFactory.create_bond('H', 'H', distance=0.74)

    # Run VQE
    solver = VQESolver(bond, ansatz_type='ucc')
    result = solver.solve()

    # Print comprehensive results
    solver.print_summary()

    # Access results
    print(f"Energy: {result['energy']:.6f} Hartree")
    print(f"Correlation: {result['correlation_energy']:.6f} Hartree")
    print(f"Analysis: {result['analysis']}")
"""

# New Solvers (Bonds Module Integration)
from kanad.solvers.base_solver import BaseSolver
from kanad.solvers.vqe_solver import VQESolver
from kanad.solvers.sqd_solver import SQDSolver
from kanad.solvers.excited_states_solver import ExcitedStatesSolver

# Legacy Solvers (for backward compatibility)
# from kanad.solvers.fci_solver import FCISolver
# from kanad.solvers.qpe_solver import QPESolver

__all__ = [
    # Base Class
    'BaseSolver',

    # New Solvers
    'VQESolver',
    'SQDSolver',
    'ExcitedStatesSolver',

    # Legacy (uncomment when fixed)
    # 'FCISolver',
    # 'QPESolver',
]
