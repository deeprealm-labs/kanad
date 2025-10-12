#!/usr/bin/env python3
"""
MASTER BENCHMARK SCRIPT FOR SQD AND EXCITED STATE SOLVERS

Runs all SQD and excited state solver tests systematically and generates analysis.
"""

import subprocess
import json
import time
from datetime import datetime
from pathlib import Path

# Define all test scripts
SQD_TESTS = [
    "test_h2_sqd.py",
    "test_lih_sqd.py",
    "test_heh_sqd.py",
]

EXCITED_TESTS = [
    "test_h2_excited.py",
    "test_lih_excited.py",
    "test_heh_excited.py",
]

ALL_TESTS = SQD_TESTS + EXCITED_TESTS

def run_test(test_file: str, timeout: int = 180) -> dict:
    """
    Run a single test script and parse output.

    Args:
        test_file: Name of test file
        timeout: Maximum time in seconds

    Returns:
        Dictionary with test results
    """
    test_path = Path(__file__).parent / test_file

    result = {
        'test': test_file,
        'passed': False,
        'error': None,
        'output': '',
        'time': 0,
        'energy': None,
        'solver_type': 'sqd' if 'sqd' in test_file else 'excited',
        'molecule': test_file.split('_')[1].upper().replace('.PY', ''),
    }

    try:
        proc_result = subprocess.run(
            [f". env/bin/activate && python {test_path}"],
            shell=True,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=Path(__file__).parent.parent
        )

        result['output'] = proc_result.stdout + proc_result.stderr
        result['passed'] = proc_result.returncode == 0 and '✓ PASS' in result['output']

        # Extract energy
        for line in result['output'].split('\n'):
            if 'SQD Energy:' in line:
                result['energy'] = float(line.split(':')[1].split('Ha')[0].strip())
            elif 'Ground State:' in line and 'Ha' in line:
                result['energy'] = float(line.split(':')[1].split('Ha')[0].strip())

            if 'Time:' in line and 's' in line:
                result['time'] = float(line.split(':')[1].split('s')[0].strip())

        if not result['passed']:
            result['error'] = 'Test failed (assertion or other error)'

    except subprocess.TimeoutExpired:
        result['error'] = f'Timeout after {timeout}s'
    except Exception as e:
        result['error'] = str(e)

    return result

def main():
    """Run all tests and generate reports."""
    print("=" * 80)
    print("SQD AND EXCITED STATE SOLVER BENCHMARKS")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total tests: {len(ALL_TESTS)}")
    print("=" * 80)

    all_results = []
    passed = 0
    failed = 0

    # Run all tests
    for i, test in enumerate(ALL_TESTS, 1):
        print(f"\n[{i}/{len(ALL_TESTS)}] Running {test}...")
        result = run_test(test, timeout=180)
        all_results.append(result)

        if result['passed']:
            passed += 1
            status = "✓ PASS"
        else:
            failed += 1
            status = "✗ FAIL"

        print(f"  {status} ({result['time']:.1f}s)")
        if result['error']:
            print(f"  Error: {result['error']}")

    # Generate summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total tests:  {len(ALL_TESTS)}")
    print(f"Passed:       {passed} ({100*passed/len(ALL_TESTS):.0f}%)")
    print(f"Failed:       {failed}")

    # Summary by solver type
    sqd_results = [r for r in all_results if r['solver_type'] == 'sqd']
    excited_results = [r for r in all_results if r['solver_type'] == 'excited']

    sqd_passed = sum(1 for r in sqd_results if r['passed'])
    excited_passed = sum(1 for r in excited_results if r['passed'])

    print(f"\nSQD Solver:        {sqd_passed}/{len(sqd_results)} passed")
    print(f"Excited States:    {excited_passed}/{len(excited_results)} passed")

    # Detailed results
    print("\n" + "=" * 80)
    print("DETAILED RESULTS")
    print("=" * 80)

    print("\n--- SQD Solver Results ---")
    for r in sqd_results:
        status = "✓" if r['passed'] else "✗"
        energy_str = f"{r['energy']:.6f} Ha" if r['energy'] else "N/A"
        print(f"  {status} {r['molecule']:6s} - Energy: {energy_str:15s} Time: {r['time']:5.1f}s")

    print("\n--- Excited States Solver Results ---")
    for r in excited_results:
        status = "✓" if r['passed'] else "✗"
        energy_str = f"{r['energy']:.6f} Ha" if r['energy'] else "N/A"
        print(f"  {status} {r['molecule']:6s} - Ground: {energy_str:15s} Time: {r['time']:5.1f}s")

    # Save results to JSON
    json_path = Path(__file__).parent / 'sqd_excited_results.json'
    with open(json_path, 'w') as f:
        json.dump({
            'timestamp': datetime.now().isoformat(),
            'summary': {
                'total': len(ALL_TESTS),
                'passed': passed,
                'failed': failed,
                'success_rate': passed / len(ALL_TESTS)
            },
            'by_solver': {
                'sqd': {'passed': sqd_passed, 'total': len(sqd_results)},
                'excited': {'passed': excited_passed, 'total': len(excited_results)}
            },
            'results': all_results
        }, f, indent=2)

    print(f"\nResults saved to: {json_path}")

    # Save detailed log
    log_path = Path(__file__).parent / 'sqd_excited_log.txt'
    with open(log_path, 'w') as f:
        f.write("SQD AND EXCITED STATE SOLVER BENCHMARK LOG\n")
        f.write("=" * 80 + "\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Passed: {passed}/{len(ALL_TESTS)}\n\n")

        for r in all_results:
            f.write(f"\n{'=' * 80}\n")
            f.write(f"Test: {r['test']}\n")
            f.write(f"Molecule: {r['molecule']}\n")
            f.write(f"Solver: {r['solver_type'].upper()}\n")
            f.write(f"Passed: {r['passed']}\n")
            f.write(f"Time: {r['time']:.2f}s\n")
            if r['energy']:
                f.write(f"Energy: {r['energy']:.8f} Ha\n")
            if r['error']:
                f.write(f"Error: {r['error']}\n")
            f.write(f"\n--- Output ---\n{r['output']}\n")

    print(f"Detailed log saved to: {log_path}")

    print("\n" + "=" * 80)
    print("BENCHMARK COMPLETE")
    print("=" * 80)

    return 0 if failed == 0 else 1

if __name__ == "__main__":
    exit(main())
