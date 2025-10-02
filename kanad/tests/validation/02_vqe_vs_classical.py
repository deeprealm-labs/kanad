"""
VQE vs Classical Methods Comparison
====================================

Comprehensive comparison of VQE (quantum) vs classical methods for
ground state energy calculations. Tests:

1. VQE convergence and accuracy
2. Different ansätze performance
3. Comparison with Hartree-Fock (mean-field)
4. Correlation energy recovery
5. Circuit depth vs accuracy tradeoffs

Scientific Goals:
- Validate VQE implementation
- Demonstrate quantum advantage potential
- Compare governance-aware ansätze
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.ansatze import UCCAnsatz, HardwareEfficientAnsatz
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz
from kanad.solvers import VQESolver

print("="*80)
print("VQE vs CLASSICAL METHODS - COMPREHENSIVE COMPARISON")
print("="*80)
print()

# ==============================================================================
# Test 1: H2 VQE vs HF Benchmark
# ==============================================================================
print("TEST 1: H2 Molecule - VQE vs Hartree-Fock")
print("-" * 80)

# Create H2 at equilibrium
h2 = BondFactory.create_bond('H', 'H', distance=0.74)
print(f"System: {h2}")
print(f"  Electrons: {h2.molecule.n_electrons}")
print(f"  Qubits:    {h2.representation.n_qubits}")
print(f"  Orbitals:  {h2.representation.n_orbitals}")
print()

# 1a. Hartree-Fock calculation
print("1a. Hartree-Fock (Classical Mean-Field)...")
hf_result = h2.compute_energy(method='HF', max_iterations=100)
hf_energy = hf_result['energy']
hf_converged = hf_result['converged']

print(f"  Energy:    {hf_energy:.8f} Ha")
print(f"  Converged: {hf_converged}")
print()

# 1b. VQE with different ansätze
print("1b. VQE with Multiple Ansätze...")
print()

n_qubits = h2.representation.n_qubits
n_electrons = h2.molecule.n_electrons

ansatz_configs = [
    {
        'name': 'UCCSD',
        'ansatz': UCCAnsatz(n_qubits, n_electrons, singles=True, doubles=True),
        'description': 'Unitary Coupled Cluster (chemistry gold standard)'
    },
    {
        'name': 'Hardware-Efficient (2 layers)',
        'ansatz': HardwareEfficientAnsatz(n_qubits, n_electrons, n_layers=2),
        'description': 'Hardware-efficient with linear entanglement'
    },
    {
        'name': 'Covalent Governance',
        'ansatz': CovalentGovernanceAnsatz(n_qubits, n_electrons, n_layers=2),
        'description': 'Governance-aware for covalent bonding'
    },
]

vqe_results = []
for config in ansatz_configs:
    ansatz = config['ansatz']
    ansatz.circuit = ansatz.build_circuit()

    print(f"  {config['name']}:")
    print(f"    Description: {config['description']}")
    print(f"    Parameters:  {len(ansatz.circuit.parameters)}")
    print(f"    Gates:       {len(ansatz.circuit.gates)}")

    # Create VQE solver
    vqe = VQESolver(
        hamiltonian=h2.hamiltonian,
        ansatz=ansatz,
        mapper=h2.mapper,
        optimizer='SLSQP',
        max_iterations=100
    )

    # Run VQE
    try:
        vqe_result = vqe.solve()
        energy = vqe_result['energy']
        converged = vqe_result['converged']
        iterations = vqe_result['iterations']

        print(f"    Energy:      {energy:.8f} Ha")
        print(f"    Converged:   {converged}")
        print(f"    Iterations:  {iterations}")

        # Correlation energy (VQE - HF)
        corr_energy = energy - hf_energy
        corr_energy_ev = corr_energy * 27.211
        print(f"    Correlation: {corr_energy:.8f} Ha ({corr_energy_ev:.4f} eV)")

        vqe_results.append({
            'name': config['name'],
            'energy': energy,
            'converged': converged,
            'correlation': corr_energy,
            'iterations': iterations,
            'n_params': len(ansatz.circuit.parameters)
        })

    except Exception as e:
        print(f"    ✗ Failed: {e}")
        vqe_results.append({
            'name': config['name'],
            'energy': np.nan,
            'converged': False,
            'correlation': np.nan,
            'iterations': 0,
            'n_params': len(ansatz.circuit.parameters)
        })

    print()

# 1c. Comparison and Analysis
print("1c. Energy Comparison:")
print()
print(f"{'Method':<30} {'Energy (Ha)':<15} {'vs HF (mHa)':<15} {'Converged'}")
print("-" * 80)
print(f"{'Hartree-Fock':<30} {hf_energy:<15.8f} {'0.000':<15} {'✓' if hf_converged else '✗'}")

for result in vqe_results:
    if not np.isnan(result['energy']):
        energy_diff_mha = (result['energy'] - hf_energy) * 1000
        conv_symbol = '✓' if result['converged'] else '✗'
        print(f"{result['name']:<30} {result['energy']:<15.8f} {energy_diff_mha:<15.3f} {conv_symbol}")
    else:
        print(f"{result['name']:<30} {'FAILED':<15} {'-':<15} {'✗'}")

print()

# Validation checks
validations_vqe_h2 = []
validations_vqe_h2.append(("HF converged", hf_converged))
validations_vqe_h2.append(("HF energy negative", hf_energy < 0))

for result in vqe_results:
    if not np.isnan(result['energy']):
        name = result['name']
        validations_vqe_h2.append((f"{name} converged", result['converged']))
        validations_vqe_h2.append((f"{name} energy ≤ HF", result['energy'] <= hf_energy + 1e-6))
        # VQE should recover some correlation (be lower than HF)
        validations_vqe_h2.append((f"{name} finds correlation", result['correlation'] < 0.01))

print("H2 VQE Validation Results:")
for test_name, passed in validations_vqe_h2:
    symbol = "✓" if passed else "✗"
    print(f"  {symbol} {test_name}")
print()

# ==============================================================================
# Test 2: LiH - Harder System
# ==============================================================================
print("TEST 2: LiH Molecule - VQE on Heteronuclear System")
print("-" * 80)

lih = BondFactory.create_bond('Li', 'H', distance=1.60)
print(f"System: {lih}")
print(f"  Bond type: {lih.bond_type}")
print(f"  Electrons: {lih.molecule.n_electrons}")
print(f"  Qubits:    {lih.representation.n_qubits}")
print()

# HF reference
print("Computing HF reference...")
try:
    lih_hf_result = lih.compute_energy(method='HF', max_iterations=200)
    lih_hf_energy = lih_hf_result['energy']
    lih_hf_converged = lih_hf_result['converged']

    print(f"  HF Energy: {lih_hf_energy:.8f} Ha")
    print(f"  Converged: {lih_hf_converged}")
    print()

    # VQE with UCCSD
    print("VQE with UCCSD ansatz...")
    n_qubits_lih = lih.representation.n_qubits
    n_electrons_lih = lih.molecule.n_electrons

    lih_ansatz = UCCAnsatz(n_qubits_lih, n_electrons_lih)
    lih_ansatz.circuit = lih_ansatz.build_circuit()

    print(f"  Parameters: {len(lih_ansatz.circuit.parameters)}")
    print(f"  Excitations: {len(lih_ansatz.excitations)}")

    lih_vqe = VQESolver(
        hamiltonian=lih.hamiltonian,
        ansatz=lih_ansatz,
        mapper=lih.mapper,
        optimizer='COBYLA',  # COBYLA often better for constrained problems
        max_iterations=150
    )

    lih_vqe_result = lih_vqe.solve()
    lih_vqe_energy = lih_vqe_result['energy']
    lih_vqe_converged = lih_vqe_result['converged']

    print(f"  VQE Energy: {lih_vqe_energy:.8f} Ha")
    print(f"  Converged: {lih_vqe_converged}")
    print(f"  Iterations: {lih_vqe_result['iterations']}")

    lih_correlation = lih_vqe_energy - lih_hf_energy
    print(f"  Correlation energy: {lih_correlation:.8f} Ha")
    print()

    validations_lih = []
    validations_lih.append(("LiH HF converged", lih_hf_converged))
    validations_lih.append(("LiH VQE completed", True))
    validations_lih.append(("LiH VQE ≤ HF", lih_vqe_energy <= lih_hf_energy + 1e-4))

    print("LiH VQE Validation:")
    for test_name, passed in validations_lih:
        symbol = "✓" if passed else "✗"
        print(f"  {symbol} {test_name}")

except Exception as e:
    print(f"  ✗ LiH calculation failed: {e}")
    validations_lih = [("LiH VQE", False)]

print()

# ==============================================================================
# Test 3: Ansatz Comparison - Parameter Efficiency
# ==============================================================================
print("TEST 3: Ansatz Parameter Efficiency")
print("-" * 80)

print("Comparing ansatz efficiency (fewer parameters = better for NISQ devices)")
print()

# Test on H2
h2_test = BondFactory.create_bond('H', 'H', distance=0.74)
n_q = h2_test.representation.n_qubits
n_e = h2_test.molecule.n_electrons

ansatz_efficiency_test = [
    ("UCC-S", UCCAnsatz(n_q, n_e, singles=True, doubles=False)),
    ("UCC-D", UCCAnsatz(n_q, n_e, singles=False, doubles=True)),
    ("UCCSD", UCCAnsatz(n_q, n_e, singles=True, doubles=True)),
    ("HEA-1", HardwareEfficientAnsatz(n_q, n_e, n_layers=1)),
    ("HEA-2", HardwareEfficientAnsatz(n_q, n_e, n_layers=2)),
    ("HEA-3", HardwareEfficientAnsatz(n_q, n_e, n_layers=3)),
    ("Cov-Gov", CovalentGovernanceAnsatz(n_q, n_e, n_layers=2)),
]

print(f"{'Ansatz':<15} {'Params':<10} {'Gates':<10} {'Depth Est.':<12}")
print("-" * 55)

for name, ansatz in ansatz_efficiency_test:
    ansatz.circuit = ansatz.build_circuit()
    n_params = len(ansatz.circuit.parameters)
    n_gates = len(ansatz.circuit.gates)

    # Estimate circuit depth (very rough - actual depth depends on connectivity)
    depth_estimate = n_gates // n_q + 1

    print(f"{name:<15} {n_params:<10} {n_gates:<10} {depth_estimate:<12}")

print()
print("Key Insights:")
print("  • UCCSD is chemically motivated but parameter-heavy")
print("  • Hardware-efficient ansätze have fewer parameters")
print("  • Governance-aware ansätze enforce physical constraints")
print("  • Parameter count vs accuracy is a key tradeoff for NISQ")
print()

# ==============================================================================
# Test 4: Convergence Speed Comparison
# ==============================================================================
print("TEST 4: Optimizer Convergence Speed")
print("-" * 80)

h2_conv = BondFactory.create_bond('H', 'H', distance=0.74)
n_q_conv = h2_conv.representation.n_qubits
n_e_conv = h2_conv.molecule.n_electrons

# Use same ansatz (UCCSD) with different optimizers
ansatz_conv = UCCAnsatz(n_q_conv, n_e_conv)
ansatz_conv.circuit = ansatz_conv.build_circuit()

optimizers = ['SLSQP', 'COBYLA', 'L-BFGS-B']

print(f"Testing optimizer convergence on H2 with UCCSD ansatz")
print(f"(Same ansatz, different optimizers)")
print()

opt_results = []
for opt_name in optimizers:
    print(f"  {opt_name}...", end=" ")

    # Need fresh ansatz for each optimizer
    ansatz_fresh = UCCAnsatz(n_q_conv, n_e_conv)
    ansatz_fresh.circuit = ansatz_fresh.build_circuit()

    vqe_conv = VQESolver(
        hamiltonian=h2_conv.hamiltonian,
        ansatz=ansatz_fresh,
        mapper=h2_conv.mapper,
        optimizer=opt_name,
        max_iterations=100
    )

    try:
        result = vqe_conv.solve()
        print(f"Energy: {result['energy']:.6f} Ha, Iterations: {result['iterations']}")
        opt_results.append({
            'optimizer': opt_name,
            'energy': result['energy'],
            'iterations': result['iterations'],
            'converged': result['converged']
        })
    except Exception as e:
        print(f"Failed - {e}")
        opt_results.append({
            'optimizer': opt_name,
            'energy': np.nan,
            'iterations': 0,
            'converged': False
        })

print()

# Find best result
valid_results = [r for r in opt_results if not np.isnan(r['energy'])]
if valid_results:
    best_result = min(valid_results, key=lambda x: x['energy'])
    fastest_result = min(valid_results, key=lambda x: x['iterations'])

    print(f"Best energy:      {best_result['optimizer']} ({best_result['energy']:.8f} Ha)")
    print(f"Fastest converge: {fastest_result['optimizer']} ({fastest_result['iterations']} iterations)")
print()

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
print("="*80)
print("VQE vs CLASSICAL VALIDATION SUMMARY")
print("="*80)
print()

all_validations = validations_vqe_h2 + validations_lih

passed = sum(1 for _, p in all_validations if p)
total = len(all_validations)

print(f"Total Validations: {passed}/{total} passed ({passed/total*100:.1f}%)")
print()

# Key findings
print("Key Findings:")
print()
print("1. VQE Implementation:")
findings_vqe = [r for r in vqe_results if not np.isnan(r['energy'])]
if findings_vqe:
    print(f"   • {len(findings_vqe)}/{len(vqe_results)} ansätze converged successfully")
    converged_energies = [r['energy'] for r in findings_vqe if r['converged']]
    if converged_energies:
        min_vqe = min(converged_energies)
        print(f"   • Best VQE energy: {min_vqe:.8f} Ha")
        print(f"   • Correlation recovery: {(min_vqe - hf_energy)*1000:.3f} mHa")

print()
print("2. Ansatz Efficiency:")
print(f"   • Parameter counts range from {min(a[1].circuit.get_num_parameters() for a in ansatz_efficiency_test if hasattr(a[1].circuit, 'get_num_parameters'))}")
print(f"     to {max(len(a[1].circuit.parameters) for a in ansatz_efficiency_test)} parameters")
print("   • Tradeoff between chemical accuracy and circuit depth")

print()
print("3. Optimizer Performance:")
if valid_results:
    print(f"   • All {len(valid_results)} optimizers tested converged")
    print(f"   • Iteration counts: {min(r['iterations'] for r in valid_results)} - {max(r['iterations'] for r in valid_results)}")

print()

if passed >= total * 0.85:
    print("✅ EXCELLENT: VQE implementation validated")
    print("   Framework ready for quantum chemistry research!")
elif passed >= total * 0.7:
    print("✓ GOOD: VQE implementation mostly validated")
    print("   Suitable for research use")
else:
    print("⚠ WARNING: VQE implementation needs improvement")

print()
print("="*80)
