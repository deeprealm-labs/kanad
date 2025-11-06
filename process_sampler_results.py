"""
Process IBM Sampler Results for Hi-VQE

Converts raw measurement counts into:
- Diagonal energies
- Ground state energy
- Rich visualization data
- Error analysis

Usage:
    export IBM_API='your_token'
    export IBM_CRN='your_crn'
    python process_sampler_results.py <job_id>
"""
import os
import sys
import json
from kanad.backends.ibm.sampler_backend import IBMSamplerBackend
from kanad.bonds import BondFactory
from kanad.core.configuration import ConfigurationSubspace
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol

if len(sys.argv) < 2:
    print("ERROR: Job ID required")
    print("Usage: python process_sampler_results.py <job_id>")
    exit(1)

job_id = sys.argv[1]

IBM_API = os.getenv('IBM_API')
IBM_CRN = os.getenv('IBM_CRN')

if not IBM_API:
    print("ERROR: IBM_API not set")
    exit(1)

print("="*80)
print("PROCESSING IBM SAMPLER RESULTS - HI-VQE")
print("="*80)
print(f"\nJob ID: {job_id}")

# Initialize backend
backend = IBMSamplerBackend(
    backend_name='ibm_torino',
    api_token=IBM_API,
    channel='ibm_cloud' if IBM_CRN else 'ibm_quantum_platform',
    crn=IBM_CRN
)

# Check job status
print(f"\n📋 Checking job status...")
status = backend.get_job_status(job_id)
print(f"   Status: {status}")

if status != 'DONE':
    print(f"\n⏳ Job not complete yet (status: {status})")
    print(f"   Check back later!")
    exit(0)

print(f"\n✅ Job completed! Processing results...")

# Rebuild subspace and Hamiltonian
print(f"\n🔧 Rebuilding Hi-VQE configuration...")
bond = BondFactory.create_bond('H', 'H', distance=0.74)
hamiltonian = bond.hamiltonian.to_sparse_hamiltonian(mapper='jordan_wigner')

protocol = CovalentGovernanceProtocol()
subspace = ConfigurationSubspace(4, 2, protocol=protocol)
hf_config = subspace.get_hf_configuration()
subspace.add_config(hf_config)
single_excs = subspace.generate_single_excitations(hf_config)
subspace.add_configs(single_excs)

print(f"   Subspace size: {len(subspace)}")

# Process results
print(f"\n📊 Processing measurement counts...")
results = backend.process_hivqe_results(
    job_id=job_id,
    observable=hamiltonian,
    subspace=subspace
)

# Display results
print(f"\n" + "="*80)
print("RESULTS")
print("="*80)

print(f"\n🎯 Diagonal Energies (⟨config|H|config⟩):")
for i, (config, energy) in enumerate(zip(subspace, results['diagonal_energies'])):
    print(f"   {config}: {energy:.8f} Ha")

print(f"\n🏆 Ground State Energy:")
print(f"   IBM Hardware:     {results['ground_energy']:.8f} Ha")
print(f"   Expected (FCI):   -1.13728383 Ha")
print(f"   Difference:       {abs(results['ground_energy'] - (-1.13728383)):.8f} Ha")

error = abs(results['ground_energy'] - (-1.13728383))
if error < 0.001:
    print(f"   ✅ EXCELLENT! Chemical accuracy achieved")
elif error < 0.01:
    print(f"   ✅ VERY GOOD! Within 0.01 Ha")
elif error < 0.1:
    print(f"   ✅ GOOD! Within 0.1 Ha")
else:
    print(f"   ⚠️  Large error (hardware noise)")

print(f"\n📈 Energy Spectrum:")
for i, eigenval in enumerate(results['eigenvalues'][:5]):  # Show first 5
    print(f"   State {i}: {eigenval:.8f} Ha")

print(f"\n🎨 Ground State Composition:")
viz_data = results['visualization_data']
for i, (config, amp, prob) in enumerate(zip(
    viz_data['configurations'],
    viz_data['ground_state_composition']['amplitudes'],
    viz_data['ground_state_composition']['probabilities']
)):
    # Convert to real if it's complex
    prob_real = prob.real if hasattr(prob, 'real') else prob
    amp_real = amp.real if hasattr(amp, 'real') else amp
    if prob_real > 0.01:  # Show significant contributions
        print(f"   {config}: amplitude={amp_real:+.4f}, probability={prob_real:.4f}")

print(f"\n📊 Measurement Statistics:")
for stat in viz_data['measurement_statistics']:
    print(f"   Config {stat['configuration']}:")
    print(f"     Most probable outcome: {stat['most_probable_outcome']} ({stat['max_probability']:.3f})")
    print(f"     Unique outcomes: {stat['num_outcomes']}")
    print(f"     Entropy: {stat['entropy']:.3f}")

print(f"\n🔍 Error Analysis:")
error_data = results['error_analysis']
print(f"   Shot noise:")
print(f"     Average: {error_data['shot_noise']['average']:.6f}")
print(f"     Max: {error_data['shot_noise']['max']:.6f}")
print(f"   Readout fidelity:")
print(f"     Average: {error_data['readout_fidelity']['average']:.4f}")
print(f"     Min: {error_data['readout_fidelity']['min']:.4f}")
print(f"   Mitigation applied:")
print(f"     Twirling: {error_data['mitigation_applied']['twirling']}")

# Save full results for web app
output_file = f'hivqe_results_{job_id}.json'
with open(output_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n✓ Full results saved to: {output_file}")

# Save summary for easy viewing
summary = {
    'job_id': job_id,
    'backend': results['backend'],
    'ground_energy': results['ground_energy'],
    'expected_energy': -1.13728383,
    'error': abs(results['ground_energy'] - (-1.13728383)),
    'relative_error_percent': abs(results['ground_energy'] - (-1.13728383)) / abs(-1.13728383) * 100,
    'diagonal_energies': results['diagonal_energies'],
    'dominant_configurations': viz_data['ground_state_composition']['dominant_configs']
}

summary_file = f'hivqe_summary_{job_id}.json'
with open(summary_file, 'w') as f:
    json.dump(summary, f, indent=2)

print(f"✓ Summary saved to: {summary_file}")

print(f"\n💡 For Web App:")
print(f"   Load {output_file} for full visualization data")
print(f"   Includes: counts, probabilities, energy spectrum, error analysis")

print("\n" + "="*80)
print("🎉 Hi-VQE results processed successfully!")
print("="*80)
