"""
Master Validation Runner.

Runs all validation scripts and generates a comprehensive report.
"""

import subprocess
import sys
import time
from pathlib import Path

print("=" * 80)
print("KANAD FRAMEWORK - COMPREHENSIVE VALIDATION")
print("=" * 80)

validation_scripts = [
    "01_vqe_solver_validation.py",
    "02_sqd_solver_validation.py",
    "03_excited_states_validation.py",
    "04_mapper_comparison.py",
    "05_hamiltonian_comparison.py",
    "06_basis_set_validation.py",
    "07_analysis_modules_validation.py",
    "08_io_modules_validation.py",
    "09_bluequbit_cloud_validation.py",
    "10_complex_molecules_validation.py",
]

results = {}
total_start = time.time()

for script in validation_scripts:
    print(f"\n{'=' * 80}")
    print(f"Running {script}")
    print('=' * 80)

    script_path = Path(__file__).parent / script
    start_time = time.time()

    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout per script
        )

        elapsed = time.time() - start_time

        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)

        # Check if script passed
        passed = "PASSED" in result.stdout and result.returncode == 0

        results[script] = {
            'passed': passed,
            'returncode': result.returncode,
            'elapsed': elapsed,
            'stdout': result.stdout,
            'stderr': result.stderr
        }

        print(f"\n{'✓' if passed else '✗'} {script} completed in {elapsed:.1f}s")

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        print(f"✗ {script} TIMEOUT after {elapsed:.1f}s")
        results[script] = {
            'passed': False,
            'timeout': True,
            'elapsed': elapsed
        }

    except Exception as e:
        elapsed = time.time() - start_time
        print(f"✗ {script} EXCEPTION: {e}")
        results[script] = {
            'passed': False,
            'exception': str(e),
            'elapsed': elapsed
        }

total_elapsed = time.time() - total_start

# Generate summary report
print("\n" + "=" * 80)
print("VALIDATION SUMMARY")
print("=" * 80)

print(f"\nTotal time: {total_elapsed:.1f}s\n")

passed_count = sum(1 for r in results.values() if r.get('passed', False))
total_count = len(results)

print(f"{'Script':<50} {'Status':<10} {'Time':<10}")
print("-" * 80)

for script, result in results.items():
    status = "✓ PASS" if result.get('passed', False) else "✗ FAIL"
    elapsed = f"{result.get('elapsed', 0):.1f}s"
    print(f"{script:<50} {status:<10} {elapsed:<10}")

print("-" * 80)
print(f"{'TOTAL':<50} {passed_count}/{total_count} {total_elapsed:.1f}s")

# Overall status
print("\n" + "=" * 80)
if passed_count == total_count:
    print("✓✓✓ ALL VALIDATIONS PASSED ✓✓✓")
elif passed_count >= total_count * 0.8:
    print(f"⚠ MOSTLY PASSED ({passed_count}/{total_count})")
elif passed_count >= total_count * 0.5:
    print(f"⚠ PARTIALLY PASSED ({passed_count}/{total_count})")
else:
    print(f"✗✗✗ VALIDATION FAILED ({passed_count}/{total_count}) ✗✗✗")
print("=" * 80)

# Exit code
sys.exit(0 if passed_count == total_count else 1)
