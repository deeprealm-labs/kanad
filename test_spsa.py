"""Quick test: SPSA optimizer with H2"""
from kanad.bonds import BondFactory
from kanad.utils.vqe_solver import VQESolver

print("="*80)
print("SPSA OPTIMIZER TEST - H2")
print("="*80)

# Create H2 bond
bond = BondFactory.create_bond('H', 'H', distance=0.74)

print(f"\nH2 Molecule:")
print(f"  Atoms: {bond.molecule.n_atoms}")
print(f"  Electrons: {bond.molecule.n_electrons}")

# Test SPSA
print(f"\n{'─'*80}")
print(f"Testing SPSA Optimizer")
print(f"{'─'*80}")

solver = VQESolver(
    bond=bond,
    ansatz_type='hardware_efficient',
    optimizer='SPSA',
    max_iterations=10,
    backend='statevector'
)

result = solver.solve()

print(f"\n✅ SPSA Result:")
print(f"   Energy: {result['energy']:.8f} Ha")
print(f"   Function evals: {result['function_evaluations']}")
print(f"   Iterations: {result['iterations']}")
print(f"   Evals/iteration: {result['function_evaluations'] / result['iterations']:.1f}")

# Verify 2 evals per iteration
expected_evals = result['iterations'] * 2
print(f"\n📊 Efficiency Check:")
print(f"   Expected evals (2 per iter): {expected_evals}")
print(f"   Actual evals: {result['function_evaluations']}")
print(f"   {'✅ PASS: 2 evals/iteration' if result['function_evaluations'] == expected_evals else '❌ FAIL'}")
