"""
Master Validation Suite Runner
===============================

Runs all validation scripts and generates a comprehensive report on
the Kanad framework's scientific credibility and readiness for research use.

This script:
1. Runs all validation modules
2. Collects results and metrics
3. Generates summary report
4. Provides pass/fail assessment for professional research use

Usage:
    python validation_suite/run_all_validations.py
"""

import sys
import os
import time
import subprocess
from pathlib import Path

print("="*80)
print("KANAD FRAMEWORK - COMPREHENSIVE VALIDATION SUITE")
print("="*80)
print()
print("Running complete scientific validation of the Kanad quantum chemistry framework")
print()

# Get validation suite directory
suite_dir = Path(__file__).parent
validation_scripts = [
    ("01_molecular_benchmarks.py", "Molecular Benchmarks - Literature Comparison"),
    ("02_vqe_vs_classical.py", "VQE vs Classical Methods"),
    ("03_physical_properties.py", "Physical Properties and Symmetries"),
    ("04_edge_cases_stress_test.py", "Edge Cases and Stress Testing"),
]

results = []
total_start_time = time.time()

# ==============================================================================
# Run Each Validation Script
# ==============================================================================

for script_name, description in validation_scripts:
    script_path = suite_dir / script_name

    if not script_path.exists():
        print(f"✗ Script not found: {script_name}")
        results.append({
            'name': description,
            'status': 'NOT FOUND',
            'time': 0,
            'output': '',
        })
        continue

    print(f"{'='*80}")
    print(f"Running: {description}")
    print(f"Script:  {script_name}")
    print(f"{'='*80}")
    print()

    # Run the validation script
    start_time = time.time()

    try:
        # Run script and capture output
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )

        elapsed_time = time.time() - start_time

        # Print output
        print(result.stdout)

        if result.stderr:
            print("STDERR:", file=sys.stderr)
            print(result.stderr, file=sys.stderr)

        # Determine status
        if result.returncode == 0:
            status = "COMPLETED"
        else:
            status = f"ERROR (exit code {result.returncode})"

        results.append({
            'name': description,
            'status': status,
            'time': elapsed_time,
            'output': result.stdout,
            'returncode': result.returncode,
        })

    except subprocess.TimeoutExpired:
        elapsed_time = time.time() - start_time
        print(f"✗ TIMEOUT after {elapsed_time:.1f} seconds")

        results.append({
            'name': description,
            'status': 'TIMEOUT',
            'time': elapsed_time,
            'output': '',
            'returncode': -1,
        })

    except Exception as e:
        elapsed_time = time.time() - start_time
        print(f"✗ EXCEPTION: {e}")

        results.append({
            'name': description,
            'status': f'EXCEPTION: {type(e).__name__}',
            'time': elapsed_time,
            'output': '',
            'returncode': -1,
        })

    print()

total_elapsed_time = time.time() - total_start_time

# ==============================================================================
# Generate Summary Report
# ==============================================================================

print("="*80)
print("COMPREHENSIVE VALIDATION SUMMARY REPORT")
print("="*80)
print()

print(f"Total Execution Time: {total_elapsed_time:.2f} seconds ({total_elapsed_time/60:.1f} minutes)")
print()

# Script-by-script summary
print("Validation Suite Results:")
print("-" * 80)
print(f"{'Test Suite':<50} {'Status':<15} {'Time (s)':<10}")
print("-" * 80)

for result in results:
    status_symbol = "✓" if result['status'] == "COMPLETED" else "✗"
    print(f"{status_symbol} {result['name']:<48} {result['status']:<15} {result['time']:>8.2f}")

print("-" * 80)
print()

# Overall statistics
completed = sum(1 for r in results if r['status'] == 'COMPLETED')
total_tests = len(results)

print(f"Overall Success Rate: {completed}/{total_tests} ({completed/total_tests*100:.1f}%)")
print()

# Parse outputs for validation counts
print("Detailed Validation Metrics:")
print("-" * 80)

total_validations = 0
total_passed = 0

for result in results:
    if result['status'] == 'COMPLETED' and result['output']:
        # Try to extract validation counts from output
        output = result['output']

        # Look for patterns like "X/Y passed" or "Validation score: X/Y"
        import re

        # Pattern 1: "X/Y passed"
        matches1 = re.findall(r'(\d+)/(\d+)\s+passed', output, re.IGNORECASE)

        # Pattern 2: "Validation score: X/Y"
        matches2 = re.findall(r'validation\s+score:\s*(\d+)/(\d+)', output, re.IGNORECASE)

        # Pattern 3: "Total Validations: X/Y passed"
        matches3 = re.findall(r'total\s+validations:\s*(\d+)/(\d+)', output, re.IGNORECASE)

        all_matches = matches1 + matches2 + matches3

        if all_matches:
            # Take the last match (usually the final summary)
            passed_str, total_str = all_matches[-1]
            passed_count = int(passed_str)
            total_count = int(total_str)

            total_validations += total_count
            total_passed += passed_count

            success_rate = (passed_count / total_count * 100) if total_count > 0 else 0

            print(f"{result['name']:<50} {passed_count:>3}/{total_count:<3} ({success_rate:>5.1f}%)")
        else:
            print(f"{result['name']:<50} {'N/A':>8}")

print("-" * 80)

if total_validations > 0:
    overall_pass_rate = (total_passed / total_validations * 100)
    print(f"{'TOTAL ACROSS ALL SUITES':<50} {total_passed:>3}/{total_validations:<3} ({overall_pass_rate:>5.1f}%)")
    print()

# ==============================================================================
# Scientific Credibility Assessment
# ==============================================================================

print()
print("="*80)
print("SCIENTIFIC CREDIBILITY ASSESSMENT")
print("="*80)
print()

# Criteria for research readiness
criteria = [
    ("All test suites completed", completed == total_tests),
    ("No timeouts or exceptions", all(r['status'] == 'COMPLETED' for r in results)),
]

if total_validations > 0:
    criteria.extend([
        ("Overall pass rate > 85%", overall_pass_rate > 85),
        ("Overall pass rate > 70%", overall_pass_rate > 70),
        ("Overall pass rate > 50%", overall_pass_rate > 50),
    ])

print("Readiness Criteria:")
for criterion, met in criteria:
    symbol = "✓" if met else "✗"
    print(f"  {symbol} {criterion}")

print()

# Final verdict
all_completed = completed == total_tests

if total_validations > 0:
    if all_completed and overall_pass_rate >= 85:
        verdict = "EXCELLENT - READY FOR PROFESSIONAL RESEARCH"
        color = "GREEN"
        details = [
            "The Kanad framework has passed comprehensive scientific validation.",
            "It demonstrates:",
            "  • Accurate energy calculations vs literature benchmarks",
            "  • Correct VQE implementation and quantum circuit handling",
            "  • Proper physical symmetries and conservation laws",
            "  • Robust error handling and numerical stability",
            "",
            "The framework is suitable for professional quantum chemistry research,",
            "publication-quality calculations, and educational purposes.",
        ]
    elif all_completed and overall_pass_rate >= 70:
        verdict = "GOOD - SUITABLE FOR RESEARCH WITH VALIDATION"
        color = "YELLOW"
        details = [
            "The Kanad framework performs well in most scientific tests.",
            "It is suitable for research use, with the following recommendations:",
            "  • Cross-validate critical results with established codes",
            "  • Document any limitations encountered",
            "  • Use for exploratory and educational purposes",
            "",
            "The framework shows solid scientific foundations but may have",
            "some edge cases or advanced features requiring further development.",
        ]
    elif overall_pass_rate >= 50:
        verdict = "ACCEPTABLE - RESEARCH USE WITH CAUTION"
        color = "ORANGE"
        details = [
            "The Kanad framework shows promise but has some limitations.",
            "Recommended usage:",
            "  • Educational and learning purposes",
            "  • Prototyping and exploration",
            "  • Development and testing",
            "",
            "Not recommended for production research without thorough validation",
            "of specific use cases and cross-checking with established methods.",
        ]
    else:
        verdict = "NEEDS IMPROVEMENT - NOT READY FOR RESEARCH"
        color = "RED"
        details = [
            "The Kanad framework requires further development before research use.",
            "Current limitations prevent confident scientific use.",
            "Recommended actions:",
            "  • Review failed validation tests",
            "  • Fix critical issues identified",
            "  • Enhance numerical stability and convergence",
            "  • Expand test coverage",
        ]
else:
    verdict = "INCOMPLETE - VALIDATION METRICS NOT AVAILABLE"
    color = "GRAY"
    details = [
        "Unable to extract validation metrics from test outputs.",
        "Manual review of test results required.",
    ]

print("="*80)
print(f"FINAL VERDICT: {verdict}")
print("="*80)
print()

for line in details:
    print(line)

print()

# ==============================================================================
# Recommendations
# ==============================================================================

print("="*80)
print("RECOMMENDATIONS FOR USERS")
print("="*80)
print()

if total_validations > 0 and overall_pass_rate >= 70:
    print("✓ The framework is validated for the following use cases:")
    print("  • Small molecule electronic structure calculations (< 10 atoms)")
    print("  • Hartree-Fock and variational quantum eigensolver methods")
    print("  • Bond analysis and molecular property prediction")
    print("  • Educational quantum chemistry demonstrations")
    print("  • Prototyping governance-based quantum algorithms")
    print()
    print("Best Practices:")
    print("  • Start with simple systems (H2, LiH, H2O)")
    print("  • Verify convergence of SCF calculations")
    print("  • Compare results with reference data when available")
    print("  • Use appropriate basis sets (STO-3G validated)")
    print("  • Monitor numerical stability for challenging systems")

else:
    print("⚠ Further development recommended before research use:")
    print("  • Review failed test cases in validation outputs")
    print("  • Focus on critical failures first")
    print("  • Enhance error handling and convergence algorithms")
    print("  • Expand basis set library if needed")
    print("  • Add more robust SCF convergence techniques")

print()
print("="*80)
print("VALIDATION SUITE COMPLETE")
print("="*80)
print()

# Save results to file
report_file = suite_dir / 'validation_report.txt'

with open(report_file, 'w') as f:
    f.write("KANAD FRAMEWORK VALIDATION REPORT\n")
    f.write("="*80 + "\n\n")
    f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Total Time: {total_elapsed_time:.2f} seconds\n\n")

    f.write("Test Results:\n")
    f.write("-" * 80 + "\n")
    for result in results:
        f.write(f"{result['name']}: {result['status']} ({result['time']:.2f}s)\n")

    f.write("\n")

    if total_validations > 0:
        f.write(f"Total Validations: {total_passed}/{total_validations} ({overall_pass_rate:.1f}%)\n\n")

    f.write(f"Verdict: {verdict}\n\n")

    for line in details:
        f.write(line + "\n")

print(f"Report saved to: {report_file}")
print()

# Exit with appropriate code
if all_completed and (total_validations == 0 or overall_pass_rate >= 70):
    sys.exit(0)  # Success
else:
    sys.exit(1)  # Failure
