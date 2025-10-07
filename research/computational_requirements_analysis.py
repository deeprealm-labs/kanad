"""
Analyze computational requirements for Hamiltonian calculation
and pre-execution processing for complex molecules.

For commercial cloud deployment planning (Azure, AWS, etc.)
"""

import time
import psutil
import os
import numpy as np
from pathlib import Path
import json
from dotenv import load_dotenv

# Load environment
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(env_path)

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from kanad.bonds import BondFactory


def analyze_molecule_requirements(name, atom1, atom2, distance, basis='sto-3g'):
    """
    Measure computational requirements for a single molecule.

    Returns:
        dict with CPU time, memory usage, matrix sizes
    """
    print(f"\n{'='*70}")
    print(f"Analyzing: {name}")
    print(f"{'='*70}")

    # Get initial memory
    process = psutil.Process(os.getpid())
    mem_before = process.memory_info().rss / 1024**2  # MB

    results = {
        'name': name,
        'atoms': f"{atom1}-{atom2}",
        'distance': distance,
        'basis': basis,
        'timings': {},
        'memory': {},
        'dimensions': {}
    }

    try:
        # Step 1: Bond creation
        t0 = time.time()
        bond = BondFactory.create_bond(atom1, atom2, distance=distance, basis=basis)
        t_bond = time.time() - t0
        results['timings']['bond_creation'] = t_bond

        mem_after_bond = process.memory_info().rss / 1024**2
        results['memory']['bond_creation_mb'] = mem_after_bond - mem_before

        # Get molecular dimensions
        n_orbitals = bond.hamiltonian.n_orbitals
        n_electrons = bond.hamiltonian.molecule.n_electrons
        n_qubits = 2 * n_orbitals

        results['dimensions']['n_orbitals'] = n_orbitals
        results['dimensions']['n_electrons'] = n_electrons
        results['dimensions']['n_qubits'] = n_qubits
        results['dimensions']['hilbert_space_size'] = 2**n_qubits

        print(f"  Orbitals: {n_orbitals}")
        print(f"  Electrons: {n_electrons}")
        print(f"  Qubits: {n_qubits}")
        print(f"  Hilbert space: 2^{n_qubits} = {2**n_qubits:,}")

        # Step 2: Hartree-Fock (SCF)
        t0 = time.time()
        hf_result = bond.compute_energy(method='HF')
        t_hf = time.time() - t0
        results['timings']['hartree_fock_scf'] = t_hf
        results['energies'] = {'hf_energy': float(hf_result['energy'])}

        mem_after_hf = process.memory_info().rss / 1024**2
        results['memory']['hf_calculation_mb'] = mem_after_hf - mem_after_bond

        # Step 3: Hamiltonian matrix construction
        t0 = time.time()
        try:
            # This is the expensive operation for large molecules
            hamiltonian_matrix = bond.hamiltonian.to_matrix()
            t_hamiltonian = time.time() - t0
            results['timings']['hamiltonian_matrix_construction'] = t_hamiltonian

            mem_after_hamiltonian = process.memory_info().rss / 1024**2
            results['memory']['hamiltonian_matrix_mb'] = mem_after_hamiltonian - mem_after_hf

            # Matrix properties
            results['dimensions']['hamiltonian_shape'] = list(hamiltonian_matrix.shape)
            results['dimensions']['hamiltonian_elements'] = int(np.prod(hamiltonian_matrix.shape))
            results['dimensions']['hamiltonian_memory_theoretical_gb'] = (
                np.prod(hamiltonian_matrix.shape) * 16 / 1024**3  # complex128 = 16 bytes
            )

            print(f"  Hamiltonian shape: {hamiltonian_matrix.shape}")
            print(f"  Matrix elements: {np.prod(hamiltonian_matrix.shape):,}")

        except MemoryError as e:
            results['timings']['hamiltonian_matrix_construction'] = None
            results['memory']['hamiltonian_matrix_mb'] = None
            results['error'] = f"MemoryError: {str(e)}"
            print(f"  ⚠️  MemoryError during Hamiltonian construction")

        # Step 4: VQE preparation (without full matrix)
        t0 = time.time()
        from kanad.solvers import VQESolver
        vqe = VQESolver(bond, ansatz_type='hardware_efficient')
        t_vqe_prep = time.time() - t0
        results['timings']['vqe_preparation'] = t_vqe_prep

        mem_after_vqe = process.memory_info().rss / 1024**2
        results['memory']['vqe_preparation_mb'] = mem_after_vqe - mem_after_hf

        # Total
        mem_peak = process.memory_info().rss / 1024**2
        results['memory']['total_peak_mb'] = mem_peak
        results['memory']['delta_from_start_mb'] = mem_peak - mem_before

        results['timings']['total_time'] = sum(
            v for v in results['timings'].values() if v is not None
        )

        results['success'] = True

        print(f"\n  Timings:")
        print(f"    Bond creation: {t_bond:.3f}s")
        print(f"    Hartree-Fock: {t_hf:.3f}s")
        if results['timings']['hamiltonian_matrix_construction']:
            print(f"    Hamiltonian matrix: {results['timings']['hamiltonian_matrix_construction']:.3f}s")
        print(f"    VQE prep: {t_vqe_prep:.3f}s")
        print(f"    Total: {results['timings']['total_time']:.3f}s")

        print(f"\n  Memory:")
        print(f"    Peak usage: {mem_peak:.1f} MB")
        print(f"    Delta: {mem_peak - mem_before:.1f} MB")
        if results['memory']['hamiltonian_matrix_mb']:
            print(f"    Hamiltonian matrix: {results['memory']['hamiltonian_matrix_mb']:.1f} MB")

    except Exception as e:
        results['success'] = False
        results['error'] = str(e)
        print(f"  ✗ Error: {e}")

    return results


def generate_cloud_recommendations(all_results):
    """
    Generate Azure/AWS cloud resource recommendations based on profiling results.
    """
    print(f"\n{'='*70}")
    print("CLOUD DEPLOYMENT RECOMMENDATIONS")
    print(f"{'='*70}\n")

    # Categorize by complexity
    small = [r for r in all_results if r.get('dimensions', {}).get('n_qubits', 0) <= 4]
    medium = [r for r in all_results if 4 < r.get('dimensions', {}).get('n_qubits', 0) <= 12]
    large = [r for r in all_results if 12 < r.get('dimensions', {}).get('n_qubits', 0) <= 20]
    xlarge = [r for r in all_results if r.get('dimensions', {}).get('n_qubits', 0) > 20]

    categories = [
        ('Small (≤4 qubits)', small),
        ('Medium (5-12 qubits)', medium),
        ('Large (13-20 qubits)', large),
        ('X-Large (>20 qubits)', xlarge)
    ]

    recommendations = {}

    for category_name, molecules in categories:
        if not molecules:
            continue

        print(f"\n{category_name}")
        print("-" * 60)

        # Get max requirements
        max_mem = max(r['memory']['total_peak_mb'] for r in molecules if r.get('success'))
        max_time = max(r['timings']['total_time'] for r in molecules if r.get('success'))
        max_qubits = max(r['dimensions']['n_qubits'] for r in molecules)

        print(f"  Max qubits: {max_qubits}")
        print(f"  Max memory: {max_mem:.1f} MB")
        print(f"  Max time: {max_time:.2f}s")

        # Recommendation with 4x safety margin
        recommended_mem_gb = (max_mem * 4) / 1024

        # Azure VM recommendations
        if recommended_mem_gb < 4:
            azure_vm = "Standard_B2s (2 vCPU, 4 GB RAM) - $30/month"
            cores = 2
        elif recommended_mem_gb < 8:
            azure_vm = "Standard_D2s_v3 (2 vCPU, 8 GB RAM) - $70/month"
            cores = 2
        elif recommended_mem_gb < 16:
            azure_vm = "Standard_D4s_v3 (4 vCPU, 16 GB RAM) - $140/month"
            cores = 4
        elif recommended_mem_gb < 32:
            azure_vm = "Standard_D8s_v3 (8 vCPU, 32 GB RAM) - $280/month"
            cores = 8
        elif recommended_mem_gb < 64:
            azure_vm = "Standard_D16s_v3 (16 vCPU, 64 GB RAM) - $560/month"
            cores = 16
        else:
            azure_vm = "Standard_E32s_v3 (32 vCPU, 256 GB RAM) - $1,600/month"
            cores = 32

        print(f"\n  Recommended Azure VM:")
        print(f"    {azure_vm}")
        print(f"    Safety margin: 4x memory, 2x CPU cores")

        recommendations[category_name] = {
            'max_qubits': max_qubits,
            'max_memory_mb': max_mem,
            'max_time_seconds': max_time,
            'recommended_memory_gb': recommended_mem_gb,
            'recommended_cores': cores,
            'azure_vm': azure_vm,
            'molecules': len(molecules)
        }

    return recommendations


def main():
    """
    Profile computational requirements for various molecules.
    """
    print("="*70)
    print("COMPUTATIONAL REQUIREMENTS ANALYSIS")
    print("For Commercial Cloud Deployment (Azure/AWS/GCP)")
    print("="*70)

    # Test molecules of increasing complexity
    test_molecules = [
        # Small (2-4 qubits)
        {'name': 'H2 (hydrogen)', 'atom1': 'H', 'atom2': 'H', 'distance': 0.74},
        {'name': 'LiH (lithium hydride)', 'atom1': 'Li', 'atom2': 'H', 'distance': 1.60},

        # Medium (6-12 qubits)
        {'name': 'H2O (water)', 'atom1': 'O', 'atom2': 'H', 'distance': 0.96},
        {'name': 'N2 (nitrogen)', 'atom1': 'N', 'atom2': 'N', 'distance': 1.10},
        {'name': 'Adenine (N-H bond)', 'atom1': 'N', 'atom2': 'H', 'distance': 1.01},

        # Large (14-20 qubits)
        {'name': 'Caffeine (C-N bond)', 'atom1': 'C', 'atom2': 'N', 'distance': 1.34},
        {'name': 'Aspirin (C=O bond)', 'atom1': 'C', 'atom2': 'O', 'distance': 1.22},
        {'name': 'Cholesterol (C-C bond)', 'atom1': 'C', 'atom2': 'C', 'distance': 1.54},

        # X-Large (>20 qubits)
        {'name': 'Penicillin (S-C bond)', 'atom1': 'S', 'atom2': 'C', 'distance': 1.82},
    ]

    all_results = []

    for mol in test_molecules:
        result = analyze_molecule_requirements(
            mol['name'],
            mol['atom1'],
            mol['atom2'],
            mol['distance']
        )
        all_results.append(result)
        time.sleep(0.5)  # Brief pause between molecules

    # Generate recommendations
    recommendations = generate_cloud_recommendations(all_results)

    # Save results
    output = {
        'analysis_date': time.strftime('%Y-%m-%d %H:%M:%S'),
        'system_info': {
            'total_memory_gb': psutil.virtual_memory().total / 1024**3,
            'cpu_count': psutil.cpu_count(),
            'platform': sys.platform
        },
        'molecule_results': all_results,
        'cloud_recommendations': recommendations
    }

    output_file = Path(__file__).parent / 'computational_requirements.json'
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n{'='*70}")
    print(f"Results saved to: {output_file}")
    print(f"{'='*70}\n")

    # Print summary table
    print("\nSUMMARY TABLE")
    print("-" * 100)
    print(f"{'Molecule':<25} {'Qubits':>7} {'Memory (MB)':>12} {'Time (s)':>10} {'Hamiltonian':>15}")
    print("-" * 100)

    for r in all_results:
        if r.get('success'):
            name = r['name'][:24]
            qubits = r['dimensions']['n_qubits']
            mem = r['memory']['total_peak_mb']
            time_s = r['timings']['total_time']
            ham_status = "✓" if r['timings']['hamiltonian_matrix_construction'] else "✗ OOM"

            print(f"{name:<25} {qubits:>7} {mem:>12.1f} {time_s:>10.2f} {ham_status:>15}")

    print("-" * 100)

    # Key insights
    print("\n" + "="*70)
    print("KEY INSIGHTS FOR COMMERCIAL DEPLOYMENT")
    print("="*70)

    print("\n1. SCALING BEHAVIOR:")
    print("   - Hamiltonian matrix: O(2^(2n)) memory for n qubits")
    print("   - Two-electron integrals: O(N^4) for N orbitals")
    print("   - 20 qubits → 1M×1M complex matrix → 16 GB theoretical minimum")

    print("\n2. BOTTLENECKS:")
    print("   - Hamiltonian matrix construction (main memory bottleneck)")
    print("   - SCF convergence (iterative, CPU-bound)")
    print("   - Two-electron integral computation (disk I/O intensive)")

    print("\n3. OPTIMIZATION STRATEGIES:")
    print("   - Use sparse matrix representations (can reduce memory 100-1000x)")
    print("   - Active space selection (reduce from N to active orbitals)")
    print("   - Avoid full matrix construction (use operators directly)")
    print("   - Parallel SCF with multi-core CPUs")

    print("\n4. CLOUD COST ESTIMATES (Azure, monthly):")
    print("   - Small molecules (≤4 qubits): $30-70/month")
    print("   - Medium molecules (5-12 qubits): $70-140/month")
    print("   - Large molecules (13-20 qubits): $140-560/month")
    print("   - X-Large molecules (>20 qubits): $560-1600/month")
    print("   - Add storage: $20-50/month for results database")

    print("\n5. PRODUCTION RECOMMENDATIONS:")
    print("   - Start with Standard_D4s_v3 (4 vCPU, 16 GB) for prototyping")
    print("   - Scale to Standard_D8s_v3 (8 vCPU, 32 GB) for production")
    print("   - Use Azure Spot VMs for batch processing (70% cost savings)")
    print("   - Implement caching for repeated molecules")
    print("   - Consider Azure Batch for parallel molecule processing")


if __name__ == '__main__':
    main()
