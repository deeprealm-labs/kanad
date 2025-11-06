#!/usr/bin/env python3
"""Analyze partial results from Phase 1B test run."""

import re

# Parse the log file
with open('tests/phase_1b_output.log', 'r') as f:
    log_content = f.read()

# Extract test results
results = []
test_pattern = r'\[(\d+)/60\] Testing (\w+)\.\.\.\n.*?(?:✅|❌|⚠️)\s+(\w+):\s+([-\d.]+)\s+mHa'

matches = re.findall(test_pattern, log_content, re.DOTALL)

print("=" * 80)
print("PARTIAL RESULTS ANALYSIS - Phase 1B (26/60 tests)")
print("=" * 80)
print(f"\nTests Completed: {len(matches)}/60")
print()

# Group by molecule
molecule_pattern = r'Molecule:\s+(\w+)\s+\((\w+)\)\s+\|\s+Basis:\s+([\w-]+)'
molecules = re.findall(molecule_pattern, log_content)

current_mol = None
current_basis = None
mol_results = {}

for mol_name, mol_type, basis in molecules[:10]:  # First 10 molecule/basis combos
    key = f"{mol_name} ({basis})"
    if key not in mol_results:
        mol_results[key] = {'name': mol_name, 'type': mol_type, 'basis': basis, 'tests': []}
    current_mol = key

# Re-parse for structured results
for test_num, ansatz, status, correlation_str in matches:
    correlation = float(correlation_str)
    results.append({
        'test_num': int(test_num),
        'ansatz': ansatz,
        'status': status,
        'correlation_mha': correlation
    })

# Print results table
print(f"{'Test':<6s} {'Ansatz':<30s} {'Status':<12s} {'Correlation (mHa)':<18s}")
print("-" * 80)

for r in results:
    symbol = "✅" if r['status'] in ['EXCELLENT', 'GOOD'] else ("⚠️" if r['status'] == 'FAIR' else "❌")
    print(f"{symbol} {r['test_num']:<4d} {r['ansatz']:<30s} {r['status']:<12s} {r['correlation_mha']:>15.3f}")

# Summary by status
print("\n" + "=" * 80)
print("STATUS SUMMARY")
print("=" * 80)

status_counts = {}
for r in results:
    status = r['status']
    status_counts[status] = status_counts.get(status, 0) + 1

total = len(results)
for status in sorted(status_counts.keys()):
    count = status_counts[status]
    pct = (count / total) * 100 if total > 0 else 0
    print(f"{status:15s}: {count:3d} / {total} ({pct:5.1f}%)")

# Summary by ansatz
print("\n" + "=" * 80)
print("ANSATZ PERFORMANCE")
print("=" * 80)

ansatz_stats = {}
for r in results:
    ansatz = r['ansatz']
    if ansatz not in ansatz_stats:
        ansatz_stats[ansatz] = {'correlations': [], 'statuses': []}
    if r['status'] not in ['ABOVE_HF', 'STUCK_AT_HF']:
        ansatz_stats[ansatz]['correlations'].append(r['correlation_mha'])
    ansatz_stats[ansatz]['statuses'].append(r['status'])

for ansatz in sorted(ansatz_stats.keys()):
    stats = ansatz_stats[ansatz]
    corrs = stats['correlations']
    if corrs:
        avg_corr = sum(corrs) / len(corrs)
        min_corr = min(corrs)
        max_corr = max(corrs)
        success_rate = sum(1 for s in stats['statuses'] if s in ['EXCELLENT', 'GOOD']) / len(stats['statuses']) * 100

        print(f"\n{ansatz}:")
        print(f"  Avg correlation: {avg_corr:>8.3f} mHa")
        print(f"  Range:           {min_corr:>8.3f} to {max_corr:>8.3f} mHa")
        print(f"  Success rate:    {success_rate:>5.1f}%")

# Best results
print("\n" + "=" * 80)
print("TOP 10 BEST RESULTS (by correlation)")
print("=" * 80)

sorted_results = sorted([r for r in results if r['correlation_mha'] < 0],
                        key=lambda x: x['correlation_mha'])[:10]

for i, r in enumerate(sorted_results, 1):
    print(f"{i:2d}. {r['ansatz']:<30s}: {r['correlation_mha']:>10.3f} mHa")

print("\n" + "=" * 80)
print("KEY FINDINGS")
print("=" * 80)

excellent_count = status_counts.get('EXCELLENT', 0)
good_count = status_counts.get('GOOD', 0)
total_success = excellent_count + good_count

print(f"\n✅ Success Rate: {total_success}/{total} tests ({total_success/total*100:.1f}%)")
print(f"✅ Tests Completed: {len(results)}/60 (43.3%)")
print(f"\n⏱️  Test Duration: ~4.6 hours for 26 tests")
print(f"⏱️  Estimated Total: ~10-11 hours for all 60 tests")
print(f"\n💡 Recommendation: Focus on smaller molecules (H2, HeH+, LiH) with STO-3G")
print(f"💡 Larger molecules (H2O 10e, NH3 10e) require >12 hours per test")

print("\n" + "=" * 80)
