#!/usr/bin/env python3
"""
MASTER BENCHMARK SCRIPT
Runs all individual test scripts and collects results
"""
import subprocess
import json
import time
from datetime import datetime
from pathlib import Path

print("="*80)
print("KANAD FRAMEWORK - MASTER BENCHMARK")
print("="*80)
print()

# All test scripts
tests = [
    "test_h2_hardwareefficient.py",
    "test_h2_twolocal.py",
    "test_h2_covalentgov.py",
    "test_h2_ionicgov.py",
    "test_lih_hardwareefficient.py",
    "test_lih_twolocal.py",
    "test_lih_covalentgov.py",
    "test_lih_ionicgov.py",
    "test_heh+_hardwareefficient.py",
    "test_heh+_twolocal.py",
    "test_heh+_covalentgov.py",
    "test_heh+_ionicgov.py",
]

results = []
log_lines = []

start_time = time.time()

for i, test in enumerate(tests, 1):
    test_path = Path("benchmarks") / test
    print(f"[{i}/{len(tests)}] Running {test}... ", end='', flush=True)
    log_lines.append(f"\n{'='*80}\n[{i}/{len(tests)}] {test}\n{'='*80}\n")

    try:
        result = subprocess.run(
            [". env/bin/activate && python " + str(test_path)],
            shell=True,
            capture_output=True,
            text=True,
            timeout=180
        )

        output = result.stdout + result.stderr
        log_lines.append(output)

        if result.returncode == 0 and "✓ PASS" in output:
            # Parse output
            for line in output.split('\n'):
                if '/SLSQP:' in line:
                    parts = line.split()
                    mol_ansatz = parts[0].split('/')
                    energy = float(parts[1].replace('Ha,', ''))
                    error = float(parts[2].replace('mHa,', ''))
                    t = float(parts[3].replace('s', ''))

                    results.append({
                        'test': test,
                        'molecule': mol_ansatz[0],
                        'ansatz': mol_ansatz[1],
                        'energy': energy,
                        'error_mha': error,
                        'time_sec': t,
                        'status': 'PASS'
                    })
                    print(f"✓ {error:+7.3f} mHa ({t:.1f}s)")
                    break
        else:
            results.append({
                'test': test,
                'status': 'FAIL',
                'error': output[-200:]
            })
            print(f"✗ FAILED")

    except subprocess.TimeoutExpired:
        results.append({'test': test, 'status': 'TIMEOUT'})
        log_lines.append("TIMEOUT\n")
        print("✗ TIMEOUT")
    except Exception as e:
        results.append({'test': test, 'status': 'ERROR', 'error': str(e)})
        log_lines.append(f"ERROR: {e}\n")
        print(f"✗ ERROR: {e}")

total_time = time.time() - start_time

# Analysis
print(f"\n{'='*80}\nRESULTS SUMMARY\n{'='*80}\n")

passed = [r for r in results if r.get('status') == 'PASS']
failed = [r for r in results if r.get('status') in ['FAIL', 'ERROR', 'TIMEOUT']]

print(f"Total tests: {len(tests)}")
print(f"Passed: {len(passed)}")
print(f"Failed: {len(failed)}")
print(f"Success rate: {len(passed)/len(tests)*100:.1f}%")
print(f"Total time: {total_time/60:.1f} minutes")

if passed:
    import numpy as np
    errors = [abs(r['error_mha']) for r in passed]
    times = [r['time_sec'] for r in passed]

    print(f"\nAccuracy:")
    print(f"  Mean error: {np.mean(errors):.3f} mHa")
    print(f"  Best error: {min(errors):.3f} mHa")
    print(f"  Worst error: {max(errors):.3f} mHa")

    chemical = len([e for e in errors if e < 1.6])
    nobel = len([e for e in errors if e < 1.0])
    print(f"  Chemical accuracy (< 1.6 mHa): {chemical}/{len(passed)} ({chemical/len(passed)*100:.1f}%)")
    print(f"  Nobel accuracy (< 1.0 mHa): {nobel}/{len(passed)} ({nobel/len(passed)*100:.1f}%)")

    print(f"\nPerformance:")
    print(f"  Mean time: {np.mean(times):.1f}s")
    print(f"  Fastest: {min(times):.1f}s")
    print(f"  Slowest: {max(times):.1f}s")

    print(f"\n{'-'*80}\nTOP 10 RESULTS\n{'-'*80}")
    for i, r in enumerate(sorted(passed, key=lambda x: abs(x['error_mha']))[:10], 1):
        print(f"{i:2d}. {r['molecule']:5s}/{r['ansatz']:18s}: {r['error_mha']:+7.3f} mHa ({r['time_sec']:5.1f}s)")

if failed:
    print(f"\n{'-'*80}\nFAILED TESTS\n{'-'*80}")
    for r in failed:
        print(f"  {r['test']}: {r['status']}")

# Save results
output = {
    'timestamp': datetime.now().isoformat(),
    'total_tests': len(tests),
    'passed': len(passed),
    'failed': len(failed),
    'total_time_minutes': total_time / 60,
    'results': results
}

with open('benchmarks/benchmark_results.json', 'w') as f:
    json.dump(output, f, indent=2)

with open('benchmarks/benchmark_log.txt', 'w') as f:
    f.write(''.join(log_lines))

print(f"\n✓ Results saved to:")
print(f"  - benchmarks/benchmark_results.json")
print(f"  - benchmarks/benchmark_log.txt")
print("="*80)
