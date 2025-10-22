#!/usr/bin/env python3
"""
Test SQD solver with real-time progress callbacks
"""

from kanad.core.atom import Atom
from kanad.bonds import BondFactory
from kanad.solvers import SQDSolver


def progress_callback(stage: int, energy: float, message: str):
    """Simple progress callback for testing."""
    print(f"📊 Stage {stage}: {message} - E = {energy:.8f} Ha")


def test_sqd_with_callback():
    """Test SQD solver with progress callback."""
    print("=" * 60)
    print("Testing SQD Solver with Real-Time Progress Callbacks")
    print("=" * 60)

    # Create H2 molecule
    bond = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')

    # Create SQD solver with classical backend
    solver = SQDSolver(
        bond=bond,
        subspace_dim=10,
        circuit_depth=3,
        backend='statevector',
        enable_analysis=True,
        enable_optimization=True
    )

    print(f"\n✅ SQD Solver initialized")
    print(f"   Molecule: H2 (0.74 Å)")
    print(f"   Subspace: 10 dimensions")
    print(f"   Backend: statevector")

    # Solve with callback
    print(f"\n🚀 Starting SQD calculation with progress callbacks...\n")
    result = solver.solve(n_states=3, callback=progress_callback)

    # Print results
    print(f"\n{'='*60}")
    print("Results:")
    print(f"{'='*60}")
    print(f"Ground State Energy: {result['ground_state_energy']:.8f} Ha")
    print(f"Converged: {result['converged']}")
    print(f"Subspace Dimension: {result['subspace_dim']}")

    if 'excited_state_energies' in result:
        print(f"\nExcited States:")
        for i, E in enumerate(result['excited_state_energies'], start=1):
            excitation = (E - result['ground_state_energy']) * 27.2114  # Convert to eV
            print(f"  State {i}: {E:.8f} Ha (ΔE = {excitation:.4f} eV)")

    if 'analysis' in result:
        print(f"\n✅ Analysis data generated")

    print(f"\n{'='*60}")
    print("✅ SQD Test Complete!")
    print(f"{'='*60}")

    return result


if __name__ == "__main__":
    test_sqd_with_callback()
