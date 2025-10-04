88/usr/bin/env python
"""
Test QPE and SQD solvers on H2 molecule.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds import CovalentBond
from kanad.solvers import QPESolver, SQDSolver

# Create H2
h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

# Build Hamiltonian via CovalentBond
bond = CovalentBond(h1, h2)
hamiltonian = bond.hamiltonian

# HF baseline
print("=" * 60)
print("Testing QPE and SQD Solvers on H2")
print("=" * 60)

# solve_scf returns (density_matrix, energy)
_, hf_energy = hamiltonian.solve_scf()
# Convert from numpy scalar to float
if hasattr(hf_energy, 'item'):
    hf_energy = hf_energy.item()
else:
    hf_energy = float(hf_energy)
hf_ev = hf_energy * 27.2114
print(f"\nHF Energy: {hf_energy:.4f} Ha ({hf_ev:.4f} eV)")

# Test QPE
print("\n" + "=" * 60)
print("Testing QPE Solver")
print("=" * 60)

qpe = QPESolver(hamiltonian, n_ancilla=8)
qpe_result = qpe.solve()

qpe_energy_ha = qpe_result['energy']
qpe_energy_ev = qpe_energy_ha * 27.2114

print(f"QPE Energy: {qpe_energy_ha:.4f} Ha ({qpe_energy_ev:.4f} eV)")
print(f"Difference from HF: {abs(qpe_energy_ha - hf_energy):.4f} Ha")
print(f"Method: {qpe_result['method']}")
print(f"Phase: {qpe_result['phase']:.6f}")
print(f"Precision: {qpe_result['precision']:.2e}")

# Check if QPE is reasonable
if abs(qpe_energy_ha - hf_energy) < 0.5:
    print("\n✅ QPE PASS - Energy is reasonable")
    qpe_pass = True
else:
    print(f"\n❌ QPE FAIL - Energy difference too large: {abs(qpe_energy_ha - hf_energy):.4f} Ha")
    qpe_pass = False

# Test SQD
print("\n" + "=" * 60)
print("Testing SQD Solver")
print("=" * 60)

sqd = SQDSolver(hamiltonian, n_samples=100, max_iterations=5)  # Reduced for speed
sqd_result = sqd.solve()

sqd_energy_ha = sqd_result['energy']
sqd_energy_ev = sqd_energy_ha * 27.2114

print(f"SQD Energy: {sqd_energy_ha:.4f} Ha ({sqd_energy_ev:.4f} eV)")
print(f"Difference from HF: {abs(sqd_energy_ha - hf_energy):.4f} Ha")
print(f"Method: {sqd_result['method']}")

# Check if SQD is reasonable
if abs(sqd_energy_ha - hf_energy) < 0.5:
    print("\n✅ SQD PASS - Energy is reasonable")
    sqd_pass = True
else:
    print(f"\n❌ SQD FAIL - Energy difference too large: {abs(sqd_energy_ha - hf_energy):.4f} Ha")
    sqd_pass = False

# Summary
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"HF:  {hf_energy:.4f} Ha ({hf_ev:.4f} eV)")
print(f"QPE: {qpe_energy_ha:.4f} Ha ({qpe_energy_ev:.4f} eV) - {'✅ PASS' if qpe_pass else '❌ FAIL'}")
print(f"SQD: {sqd_energy_ha:.4f} Ha ({sqd_energy_ev:.4f} eV) - {'✅ PASS' if sqd_pass else '❌ FAIL'}")

if qpe_pass and sqd_pass:
    print("\n✅ ALL TESTS PASSED (2/2)")
else:
    failed = []
    if not qpe_pass:
        failed.append("QPE")
    if not sqd_pass:
        failed.append("SQD")
    print(f"\n❌ TESTS FAILED: {', '.join(failed)}")
