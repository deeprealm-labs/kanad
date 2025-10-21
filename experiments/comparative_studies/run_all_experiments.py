"""
Comparative Studies Campaign - Master Runner
============================================
Runs all Option A comparative experiments sequentially

Experiments:
1. Bond Type Comparison
2. Ansatz Benchmarking
3. Mapper Efficiency Tests
"""

import sys
import os
import time
from datetime import datetime

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

def run_experiment(script_path, experiment_name):
    """Run a single experiment script"""
    print("\n\n" + "="*80)
    print("="*80)
    print(f"RUNNING: {experiment_name}")
    print("="*80)
    print("="*80)

    start_time = time.time()

    try:
        # Execute the experiment script
        with open(script_path, 'r') as f:
            code = compile(f.read(), script_path, 'exec')
            exec(code, {'__name__': '__main__'})

        elapsed = time.time() - start_time
        print(f"\n✓ {experiment_name} completed in {elapsed:.2f} seconds")
        return True, elapsed

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"\n✗ {experiment_name} failed after {elapsed:.2f} seconds")
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False, elapsed


def main():
    """Run all comparative studies experiments"""
    print("="*80)
    print("KANAD FRAMEWORK - COMPARATIVE STUDIES CAMPAIGN")
    print("="*80)
    print(f"\nStarted at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    experiments = [
        ("bond_types/01_bond_type_comparison.py", "Experiment 1: Bond Type Comparison"),
        ("ansatz_benchmark/02_ansatz_benchmark.py", "Experiment 2: Ansatz Benchmarking"),
        ("mapper_efficiency/03_mapper_efficiency.py", "Experiment 3: Mapper Efficiency"),
    ]

    results = []
    total_start = time.time()

    for script_path, experiment_name in experiments:
        full_path = os.path.join(os.path.dirname(__file__), script_path)
        success, elapsed = run_experiment(full_path, experiment_name)
        results.append({
            'name': experiment_name,
            'success': success,
            'time': elapsed
        })

    total_elapsed = time.time() - total_start

    # Print final summary
    print("\n\n" + "="*80)
    print("="*80)
    print("CAMPAIGN SUMMARY")
    print("="*80)
    print("="*80)

    successful = sum(1 for r in results if r['success'])
    failed = len(results) - successful

    print(f"\nTotal experiments: {len(results)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"\nTotal time: {total_elapsed:.2f} seconds ({total_elapsed/60:.2f} minutes)")

    print("\n" + "-"*80)
    print("Experiment Details:")
    print("-"*80)
    for r in results:
        status = "✓ PASS" if r['success'] else "✗ FAIL"
        print(f"{status} | {r['name']:<50} | {r['time']:>8.2f}s")

    print("\n" + "="*80)
    print("RESULTS LOCATIONS")
    print("="*80)
    print("\n1. Bond Type Comparison:")
    print("   experiments/comparative_studies/bond_types/results_bond_comparison.json")
    print("\n2. Ansatz Benchmarking:")
    print("   experiments/comparative_studies/ansatz_benchmark/results_ansatz_benchmark.json")
    print("\n3. Mapper Efficiency:")
    print("   experiments/comparative_studies/mapper_efficiency/results_mapper_efficiency.json")

    print("\n" + "="*80)
    print(f"Campaign completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)

    return results


if __name__ == "__main__":
    results = main()
