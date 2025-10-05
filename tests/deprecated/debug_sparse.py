"""Debug sparse Hamiltonian construction."""

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom

# Create H2
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_bond = BondFactory.create_bond(H1, H2)
hamiltonian = h2_bond.hamiltonian

print("Analyzing ERI terms...")
eri = hamiltonian.eri
n_orbitals = hamiltonian.n_orbitals

total_terms = 0
number_number_terms = 0
exchange_terms = 0
general_terms = 0

for i in range(n_orbitals):
    for j in range(n_orbitals):
        for k in range(n_orbitals):
            for l in range(n_orbitals):
                eri_val = eri[i, k, j, l]

                if abs(eri_val) > 1e-12:
                    # Check each spin combination
                    for spin_combo in [(2*i, 2*j, 2*l, 2*k),  # alpha-alpha
                                       (2*i, 2*j+1, 2*l+1, 2*k),  # alpha-beta
                                       (2*i+1, 2*j, 2*l, 2*k+1),  # beta-alpha
                                       (2*i+1, 2*j+1, 2*l+1, 2*k+1)]:  # beta-beta
                        p, q, r, s = spin_combo
                        total_terms += 1

                        if p == r and q == s:
                            number_number_terms += 1
                        elif p == s and q == r:
                            exchange_terms += 1
                        else:
                            general_terms += 1
                            print(f"  General term: a†_{p} a†_{q} a_{r} a_{s} (coeff={0.5*eri_val:.6f})")

print(f"\nTotal two-body terms: {total_terms}")
print(f"  Number-number: {number_number_terms}")
print(f"  Exchange: {exchange_terms}")
print(f"  General: {general_terms}")
print(f"  Coverage: {(number_number_terms + exchange_terms) / total_terms * 100:.1f}%")
