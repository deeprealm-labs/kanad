#!/usr/bin/env python3
"""Quick test of single VQE run to see new messages"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.utils.vqe_solver import VQESolver

atom1 = Atom('H', [0.0, 0.0, 0.0])
atom2 = Atom('H', [0.0, 0.0, 0.74])
bond = BondFactory.create_bond(atom1, atom2, distance=0.74, basis='sto-3g')

solver = VQESolver(
    bond=bond,
    ansatz_type='governance',
    mapper_type='jordan_wigner',
    optimizer='COBYLA',
    max_iterations=100,
    backend='statevector',
    shots=None
)

print("Running single VQE (no multi-start)...")
result = solver.solve()

print(f"\nResult: {result['energy']:.8f} Ha")
print(f"HF: {result.get('hf_energy', 'N/A')} Ha")
