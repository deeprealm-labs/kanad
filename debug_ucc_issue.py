"""
Debug UCC Ansatz Issue

The UCC ansatz is converging to HF energy instead of capturing correlation.
This script traces through the VQE optimization to find the bug.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.solvers.vqe_solver import VQESolver

# Create H2 molecule
H1 = Atom('H', position=[0.0, 0.0, 0.0])
H2 = Atom('H', position=[0.74, 0.0, 0.0])
bond = CovalentBond(H1, H2, basis='sto-3g')

print("="*80)
print("H2 MOLECULE - UCC ANSATZ DEBUG")
print("="*80)

# Reference energies
print("\nReference Energies:")
print(f"  Exact (FCI):  -1.137284 Ha")
print(f"  HF (SCF):     -1.116759 Ha")
print(f"  Correlation:   0.020525 Ha (20.525 mHa)")

# Get HF energy
hf_result = bond.compute_energy(method='HF')
print(f"\nHartree-Fock:")
print(f"  Energy: {hf_result['energy']:.8f} Ha")
print(f"  Converged: {hf_result['converged']}")

# Create VQE solver with UCC
print("\n" + "="*80)
print("VQE WITH UCC ANSATZ")
print("="*80)

solver = VQESolver(
    bond=bond,
    ansatz_type='ucc',
    mapper_type='jordan_wigner',
    optimizer='SLSQP',
    max_iterations=100
)

print(f"\nAnsatz: {solver.ansatz.__class__.__name__}")
print(f"  Parameters: {solver.ansatz.n_parameters}")
print(f"  Excitations: {solver.ansatz.excitations}")

print(f"\nMapper: {solver.mapper.__class__.__name__}")

# Run VQE
print("\nRunning VQE optimization...")
result = solver.solve()

print("\n" + "="*80)
print("RESULTS")
print("="*80)

print(f"\nVQE Energy: {result['energy']:.8f} Ha")
print(f"Converged: {result['converged']}")
print(f"Iterations: {result['iterations']}")

# Compare to references
error_vs_exact = (result['energy'] - (-1.137284)) * 1000  # mHa
error_vs_hf = (result['energy'] - hf_result['energy']) * 1000  # mHa

print(f"\nError Analysis:")
print(f"  vs Exact:  {error_vs_exact:+.3f} mHa")
print(f"  vs HF:     {error_vs_hf:+.3f} mHa")

# Check if converged to HF
if abs(error_vs_hf) < 0.01:  # Less than 0.01 mHa from HF
    print(f"\n{'⚠'*40}")
    print("WARNING: VQE converged to HF energy!")
    print("UCC ansatz is NOT capturing electron correlation")
    print(f"{'⚠'*40}")

    print("\nPossible causes:")
    print("1. Circuit not applying excitation operators correctly")
    print("2. Initial parameters cause zero excitations")
    print("3. Optimizer stuck in local minimum (HF)")
    print("4. Hamiltonian missing correlation terms")

    # Check final parameters
    if 'parameters' in result:
        params = result['parameters']
        print(f"\nFinal parameters: {params}")
        print(f"Max param magnitude: {np.max(np.abs(params)):.6f}")
        print(f"Are all params zero? {np.allclose(params, 0, atol=1e-6)}")
else:
    print("\n✓ VQE successfully captured correlation energy")

# Check energy history
if 'energy_history' in result and len(result['energy_history']) > 0:
    print(f"\nEnergy History:")
    history = result['energy_history']
    print(f"  Initial:  {history[0]:.8f} Ha")
    print(f"  Final:    {history[-1]:.8f} Ha")
    print(f"  Change:   {(history[-1] - history[0])*1000:.3f} mHa")
    print(f"  Min:      {min(history):.8f} Ha")

    if len(history) > 5:
        print(f"\n  First 5 energies:")
        for i, E in enumerate(history[:5]):
            print(f"    {i}: {E:.8f} Ha")

print("\n" + "="*80)
