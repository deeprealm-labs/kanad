import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom

# Create H2
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2_atom = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_bond = BondFactory.create_bond(H1, H2_atom)

hamiltonian = h2_bond.hamiltonian

print(f"n_orbitals: {hamiltonian.n_orbitals}")
print(f"h_core shape: {hamiltonian.h_core.shape}")
print(f"ERI shape: {hamiltonian.eri.shape}")

# Test building one excitation operator
print(f"\nTrying to build excitation operator:")
print(f"  2*0 = 0, 2*1 = 2")
print(f"  Expected n_qubits: {2 * hamiltonian.n_orbitals}")

# Call _jordan_wigner_excitation directly
n_qubits = 4
result = hamiltonian._jordan_wigner_excitation(0, 0, n_qubits)
print(f"  Result shape: {result.shape}")
