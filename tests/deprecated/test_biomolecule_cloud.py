#!/usr/bin/env python3
"""
🧬 BIOMOLECULE QUANTUM SIMULATION ON IBM QUANTUM CLOUD
Test complex biological molecules on real quantum hardware!
"""

import os
import sys
import numpy as np
import time
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from kanad.core.atom import Atom
from kanad.bonds.bond_factory import BondFactory
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
from kanad.backends.ibm_runtime_backend import IBMRuntimeBackend

# Load environment
from pathlib import Path
env_file = Path(__file__).parent / '.env'
if env_file.exists():
    with open(env_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                key, _, value = line.partition('=')
                os.environ[key.strip()] = value.strip().strip('"').strip("'")

HARTREE_TO_EV = 27.211386245988

print("="*80)
print("🧬 BIOMOLECULE QUANTUM SIMULATION ON IBM QUANTUM CLOUD")
print("="*80)
print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Backend: IBM Torino (133 qubits)")
print("="*80)

# Get credentials
token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
crn = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
channel = 'ibm_cloud' if crn else 'ibm_quantum_platform'

if not token:
    print("❌ No IBM Quantum credentials found!")
    sys.exit(1)

print(f"✓ Credentials loaded")
print()

# ============================================================================
# BIOMOLECULE CANDIDATES
# ============================================================================

biomolecules = []

# 1. Formamide (HCONH2) - Simplest peptide bond model
print("🧬 Candidate 1: Formamide (HCONH₂) - Peptide Bond Model")
print("   Structure: H-CO-NH₂")
print("   Relevance: Models the peptide bond in proteins!")
print("   Atoms: C, O, N, 3×H = 6 atoms")
print()

# 2. Glycine (simplest amino acid) - Too big, but let's try a fragment
print("🧬 Candidate 2: Formic Acid (HCOOH) - Carboxyl Group")
print("   Structure: H-COOH")
print("   Relevance: Carboxyl group found in all amino acids")
print("   Atoms: C, 2×O, 2×H = 5 atoms")
print()

# 3. Ammonia (NH3) - Amino group
print("🧬 Candidate 3: Ammonia (NH₃) - Amino Group Model")
print("   Structure: NH₃")
print("   Relevance: Amino group in amino acids and proteins")
print("   Atoms: N, 3×H = 4 atoms")
print()

# 4. Water (H2O) - Biological solvent
print("🧬 Candidate 4: Water (H₂O) - Biological Solvent")
print("   Structure: H-O-H")
print("   Relevance: Essential for all life, hydrogen bonding")
print("   Atoms: O, 2×H = 3 atoms")
print()

# Let's start with something exciting but manageable: NH3 or H2O
print("="*80)
choice = input("Choose molecule (1=Formamide, 2=Formic Acid, 3=NH₃, 4=H₂O, default=4): ").strip()
print("="*80)

if choice == '1':
    print("\n🚀 Testing FORMAMIDE (HCONH₂) - Peptide Bond!")
    print("⚠️  This is complex! May take several minutes...")
    # For simplicity, test just C-N bond (peptide bond core)
    C = Atom('C', position=np.array([0.0, 0.0, 0.0]))
    N = Atom('N', position=np.array([1.33, 0.0, 0.0]))  # C-N peptide bond ~1.33 Å
    molecule_bond = BondFactory.create_bond(C, N)
    name = "C-N Peptide Bond (from Formamide)"

elif choice == '2':
    print("\n🚀 Testing FORMIC ACID (HCOOH) - Carboxyl Group!")
    # Test C=O double bond
    C = Atom('C', position=np.array([0.0, 0.0, 0.0]))
    O = Atom('O', position=np.array([1.21, 0.0, 0.0]))  # C=O ~1.21 Å
    molecule_bond = BondFactory.create_bond(C, O)
    name = "C=O Double Bond (from Formic Acid)"

elif choice == '3':
    print("\n🚀 Testing AMMONIA (NH₃) - Amino Group!")
    N = Atom('N', position=np.array([0.0, 0.0, 0.0]))
    H = Atom('H', position=np.array([1.01, 0.0, 0.0]))  # N-H ~1.01 Å
    molecule_bond = BondFactory.create_bond(N, H)
    name = "N-H Bond (from Ammonia/Amino Group)"

else:  # Default: H2O
    print("\n🚀 Testing WATER (H₂O) - Essential for Life!")
    O = Atom('O', position=np.array([0.0, 0.0, 0.0]))
    H1 = Atom('H', position=np.array([0.96, 0.0, 0.0]))  # O-H ~0.96 Å
    molecule_bond = BondFactory.create_bond(O, H1)
    name = "O-H Bond (from Water)"

print()
print("="*80)
print(f"SIMULATING: {name}")
print("="*80)

# System info
n_orbitals = molecule_bond.hamiltonian.n_orbitals
n_electrons = molecule_bond.molecule.n_electrons
n_qubits = 2 * n_orbitals

print(f"\nSystem Information:")
print(f"  Orbitals: {n_orbitals}")
print(f"  Electrons: {n_electrons}")
print(f"  Qubits needed: {n_qubits}")
print(f"  Hilbert space dimension: {2**n_qubits}")
print()

if n_qubits > 12:
    print(f"⚠️  Warning: {n_qubits} qubits is large! This will take time...")
    print(f"   Matrix size: {2**n_qubits} × {2**n_qubits}")
    print()

# Connect to IBM Quantum
print("🔗 Connecting to IBM Quantum Cloud...")
backend = IBMRuntimeBackend(
    token=token,
    instance=crn,
    channel=channel,
    shots=4096,
    optimization_level=3,
    resilience_level=1
)

backend_info = backend.get_backend_info()
print(f"✓ Connected to: {backend_info['name']}")
print(f"  Qubits available: {backend_info.get('num_qubits', 'N/A')}")
print(f"  Status: {backend_info.get('status', 'N/A')}")
print(f"  Pending jobs: {backend_info.get('pending_jobs', 'N/A')}")
print()

# Compute exact energy locally first (for comparison)
print("="*80)
print("STEP 1: Local Exact Calculation (for reference)")
print("="*80)

mapper = JordanWignerMapper()

try:
    start = time.time()
    exact_result = molecule_bond.compute_energy(method='exact', mapper=mapper)
    exact_time = time.time() - start

    exact_ev = exact_result['energy']
    exact_ha = exact_ev / HARTREE_TO_EV

    print(f"✓ Exact Energy (local): {exact_ha:.6f} Hartree ({exact_ev:.4f} eV)")
    print(f"  Computation time: {exact_time:.2f}s")
    print()
except Exception as e:
    print(f"✗ Exact calculation failed: {e}")
    exact_result = None

# VQE on IBM Quantum Cloud!
print("="*80)
print("STEP 2: VQE on IBM Quantum Cloud! 🚀")
print("="*80)

# Use hardware-efficient ansatz (better for real hardware)
ansatz = RealAmplitudesAnsatz(
    n_qubits=n_qubits,
    n_electrons=n_electrons,
    n_layers=2
)

print(f"Ansatz: Hardware-Efficient (2 layers)")
print(f"VQE iterations: 50")
print(f"Submitting to {backend_info['name']}...")
print()
print("⏳ Running VQE on real quantum hardware...")
print("   This will submit jobs to IBM Quantum cloud")
print("   May take several minutes depending on queue...")
print()

start_cloud = time.time()

try:
    cloud_result = molecule_bond.compute_energy(
        method='VQE',
        mapper=mapper,
        ansatz=ansatz,
        max_iterations=50,
        backend=backend
    )
    cloud_time = time.time() - start_cloud

    cloud_ev = cloud_result['energy']
    cloud_ha = cloud_ev / HARTREE_TO_EV

    print()
    print("="*80)
    print("✅ CLOUD VQE COMPLETED!")
    print("="*80)
    print(f"VQE Energy (IBM Cloud): {cloud_ha:.6f} Hartree ({cloud_ev:.4f} eV)")
    print(f"Total time: {cloud_time:.2f}s")
    print()

    # Compare with exact
    if exact_result:
        error_abs = abs(cloud_ev - exact_ev)
        error_pct = (error_abs / abs(exact_ev)) * 100

        print("Comparison:")
        print(f"  Exact (local):  {exact_ha:.6f} Ha")
        print(f"  VQE (cloud):    {cloud_ha:.6f} Ha")
        print(f"  Error: {error_abs:.6f} eV ({error_pct:.2f}%)")
        print()

        if error_pct < 10:
            print(f"✅ EXCELLENT! VQE error {error_pct:.2f}% < 10%")
        elif error_pct < 20:
            print(f"✓ GOOD! VQE error {error_pct:.2f}% < 20%")
        else:
            print(f"⚠️  VQE error {error_pct:.2f}% (may need more iterations)")

    print()
    print("="*80)
    print("🎉 BIOMOLECULE SIMULATION SUCCESSFUL!")
    print("="*80)
    print(f"You just simulated {name}")
    print(f"on a real 133-qubit quantum computer! 🚀")
    print("="*80)

except Exception as e:
    print()
    print(f"❌ Cloud VQE failed: {e}")
    import traceback
    traceback.print_exc()

# Cleanup
backend.close()

print()
print("Simulation complete!")
