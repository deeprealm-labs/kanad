"""Test configuration sampling for Hi-VQE"""
import numpy as np
from kanad.core.configuration import (
    Configuration,
    ConfigurationSubspace,
    sample_configurations_from_statevector
)

print("="*80)
print("CONFIGURATION SAMPLING TEST")
print("="*80)

# Test 1: Configuration basics
print(f"\n{'─'*80}")
print("Test 1: Configuration Creation and Validation")
print(f"{'─'*80}")

config1 = Configuration('1100', n_qubits=4)
print(f"\nConfig: {config1}")
print(f"  Bitstring: {config1.bitstring}")
print(f"  Electrons: {config1.n_electrons}")
print(f"  Integer: {config1.to_int()}")

config2 = Configuration.from_int(12, n_qubits=4)  # 12 = 1100
print(f"\nConfig from int(12): {config2}")
assert config1 == config2, "Conversion failed"
print("✅ Configuration creation PASSED")

# Test 2: Configuration subspace
print(f"\n{'─'*80}")
print("Test 2: Configuration Subspace Management")
print(f"{'─'*80}")

subspace = ConfigurationSubspace(n_qubits=4, n_electrons=2)
print(f"\nSubspace: {subspace.n_qubits} qubits, {subspace.n_electrons} electrons")

# Get HF configuration
hf = subspace.get_hf_configuration()
print(f"HF config: {hf}")
assert hf.bitstring == '1100', f"Expected '1100', got {hf.bitstring}"

# Add configurations
added = subspace.add_config(hf)
print(f"Added HF: {added}")
assert added == True, "Should add new config"

# Try adding again (should fail)
added = subspace.add_config(hf)
print(f"Added HF again: {added}")
assert added == False, "Should not add duplicate"

# Add invalid config (wrong electron count)
invalid = Configuration('1110', n_qubits=4)  # 3 electrons instead of 2
added = subspace.add_config(invalid)
print(f"Added invalid (3 electrons): {added}")
assert added == False, "Should reject invalid config"

print(f"Subspace size: {len(subspace)}")
assert len(subspace) == 1, "Should have 1 config"

print("✅ Subspace management PASSED")

# Test 3: Excitation generation
print(f"\n{'─'*80}")
print("Test 3: Single and Double Excitations")
print(f"{'─'*80}")

# Single excitations from HF (|1100⟩)
single_excs = subspace.generate_single_excitations(hf)
print(f"\nSingle excitations from {hf}:")
for exc in single_excs[:5]:  # Show first 5
    print(f"  {exc}")
print(f"Total single excitations: {len(single_excs)}")

# For 2 electrons in 4 orbitals: 2 occupied × 2 virtual = 4 single excitations
assert len(single_excs) == 4, f"Expected 4 single excitations, got {len(single_excs)}"

# Double excitations
double_excs = subspace.generate_double_excitations(hf)
print(f"\nDouble excitations from {hf}:")
for exc in double_excs[:3]:
    print(f"  {exc}")
print(f"Total double excitations: {len(double_excs)}")

# For 2 electrons in 4 orbitals: C(2,2) × C(2,2) = 1 × 1 = 1 double excitation
assert len(double_excs) == 1, f"Expected 1 double excitation, got {len(double_excs)}"

print("✅ Excitation generation PASSED")

# Test 4: Configuration sampling from statevector
print(f"\n{'─'*80}")
print("Test 4: Configuration Sampling from Statevector")
print(f"{'─'*80}")

# Create a simple state: equal superposition of HF and single excitation
# |ψ⟩ = 0.9|1100⟩ + 0.3|1010⟩ + 0.3|0110⟩
statevector = np.zeros(16, dtype=complex)
statevector[12] = 0.9  # |1100⟩ = 12
statevector[10] = 0.3  # |1010⟩ = 10
statevector[6] = 0.3   # |0110⟩ = 6

# Normalize
statevector = statevector / np.linalg.norm(statevector)

# Sample configurations
samples = sample_configurations_from_statevector(statevector, n_shots=100, n_electrons=2)

print(f"\nSampled configurations (n_shots=100, filter n_electrons=2):")
for bitstring, count in sorted(samples, key=lambda x: x[1], reverse=True):
    prob = count / 100
    print(f"  {bitstring}: {count:3d} shots ({prob:.1%})")

# Should have sampled mostly |1100⟩ (81% probability)
bitstrings = [bs for bs, _ in samples]
assert '1100' in bitstrings, "Should sample HF configuration"

# Total shots should approximately match (some may be filtered out if wrong electron count)
total_shots = sum(count for _, count in samples)
print(f"\nTotal valid shots (n_electrons=2): {total_shots}/100")

print("✅ Configuration sampling PASSED")

# Test 5: Subspace building and pruning
print(f"\n{'─'*80}")
print("Test 5: Subspace Growth and Pruning")
print(f"{'─'*80}")

# Start with HF
subspace2 = ConfigurationSubspace(n_qubits=6, n_electrons=3)
hf2 = subspace2.get_hf_configuration()
subspace2.add_config(hf2)

print(f"\nInitial subspace: {len(subspace2)} configs")

# Add single excitations
single_excs2 = subspace2.generate_single_excitations(hf2)
added = subspace2.add_configs(single_excs2)
print(f"After adding single excitations: {len(subspace2)} configs (+{added})")

# Add double excitations from HF
double_excs2 = subspace2.generate_double_excitations(hf2)
added = subspace2.add_configs(double_excs2)
print(f"After adding double excitations: {len(subspace2)} configs (+{added})")

# Simulate amplitudes (HF is large, excitations are small)
amplitudes = np.zeros(len(subspace2))
amplitudes[0] = 0.9  # HF
amplitudes[1:] = np.random.rand(len(subspace2) - 1) * 0.01  # Small excitations

# Prune low-amplitude configs
removed = subspace2.prune(amplitudes, threshold=0.005)
print(f"After pruning (threshold=0.005): {len(subspace2)} configs (-{removed})")

# HF should remain
assert len(subspace2) > 0, "Should keep at least HF"
assert subspace2[0] == hf2, "HF should be first config"

print("✅ Subspace growth and pruning PASSED")

print(f"\n{'='*80}")
print("ALL TESTS PASSED ✅")
print(f"{'='*80}")
print("\nKey Results:")
print(f"  ✅ Configuration creation and validation")
print(f"  ✅ Subspace management (add, filter, index)")
print(f"  ✅ Excitation generation (single and double)")
print(f"  ✅ Sampling from statevector (simulates Z measurement)")
print(f"  ✅ Subspace growth and pruning")
print(f"\n🎯 Configuration sampling ready for Hi-VQE integration!")
