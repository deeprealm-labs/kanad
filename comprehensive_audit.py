#!/usr/bin/env python3
"""
Comprehensive Framework Audit

Systematically finds ALL issues:
1. Hardcoded empirical factors
2. Result schema mismatches
3. Placeholders and TODOs
4. Missing governance integration
5. Inconsistent density matrix access
"""

import os
import re
from pathlib import Path
from collections import defaultdict

print("=" * 80)
print("COMPREHENSIVE KANAD FRAMEWORK AUDIT")
print("=" * 80)

issues = defaultdict(list)

# Search patterns
patterns = {
    'hardcoded_factors': [
        (r'\* 0\.[0-9]+', "Hardcoded multiplication factor"),
        (r'= [0-9]+\.[0-9]+ \*', "Hardcoded multiplication"),
        (r'n_electrons \* 0\.[0-9]+', "Electron count scaling"),
    ],
    'return_zero': [
        (r'return 0\.0(?!\))', "Return 0.0 placeholder"),
    ],
    'result_schema': [
        (r"result\['density_matrix'\]", "Old density_matrix key"),
        (r"result\['state_vector'\]", "Old state_vector key"),
        (r'results\[.*density.*\]', "Density access in results"),
    ],
    'placeholders': [
        (r'TODO|FIXME|HACK|XXX', "TODO/FIXME marker"),
        (r'placeholder|Placeholder|PLACEHOLDER', "Placeholder code"),
        (r'NotImplemented|NotImplementedError', "Not implemented"),
        (r'pass\s*#.*placeholder', "Placeholder pass statement"),
    ],
    'governance_missing': [
        (r'def _generate_subspace_basis', "Subspace generation (check governance)"),
        (r'excitations.*=.*\[\]', "Excitation generation"),
    ],
}

# Scan kanad directory
kanad_path = Path('kanad')
total_files = 0
total_lines = 0

for py_file in kanad_path.rglob('*.py'):
    total_files += 1
    try:
        with open(py_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            total_lines += len(lines)

            for line_num, line in enumerate(lines, 1):
                for category, pattern_list in patterns.items():
                    for pattern, description in pattern_list:
                        if re.search(pattern, line):
                            # Skip comments explaining the issue
                            if 'CRITICAL FIX' in line or 'Removed' in line:
                                continue

                            issues[category].append({
                                'file': str(py_file),
                                'line': line_num,
                                'description': description,
                                'code': line.strip()
                            })
    except Exception as e:
        print(f"Error reading {py_file}: {e}")

# Report findings
print(f"\n📊 Scanned: {total_files} files, {total_lines} lines\n")

print("=" * 80)
print("AUDIT RESULTS")
print("=" * 80)

# Critical issues
critical_categories = ['hardcoded_factors', 'return_zero', 'result_schema']
high_priority = ['governance_missing']
medium_priority = ['placeholders']

def print_issues(category, issue_list):
    if not issue_list:
        print(f"  ✅ No issues found")
        return

    print(f"  ❌ Found {len(issue_list)} issue(s):\n")

    # Group by file
    by_file = defaultdict(list)
    for issue in issue_list:
        by_file[issue['file']].append(issue)

    for file, file_issues in sorted(by_file.items()):
        print(f"  📁 {file}")
        for issue in file_issues[:5]:  # Show first 5 per file
            print(f"     Line {issue['line']}: {issue['description']}")
            print(f"     Code: {issue['code'][:80]}")
        if len(file_issues) > 5:
            print(f"     ... and {len(file_issues) - 5} more")
        print()

print("\n🔴 CRITICAL ISSUES\n")
for cat in critical_categories:
    if cat in issues:
        print(f"\n{cat.upper().replace('_', ' ')}:")
        print_issues(cat, issues[cat])

print("\n🟠 HIGH PRIORITY\n")
for cat in high_priority:
    if cat in issues:
        print(f"\n{cat.upper().replace('_', ' ')}:")
        print_issues(cat, issues[cat])

print("\n🟡 MEDIUM PRIORITY\n")
for cat in medium_priority:
    if cat in issues:
        print(f"\n{cat.upper().replace('_', ' ')}:")
        print_issues(cat, issues[cat])

# Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

total_critical = sum(len(issues[cat]) for cat in critical_categories if cat in issues)
total_high = sum(len(issues[cat]) for cat in high_priority if cat in issues)
total_medium = sum(len(issues[cat]) for cat in medium_priority if cat in issues)

print(f"\n🔴 Critical issues: {total_critical}")
print(f"🟠 High priority: {total_high}")
print(f"🟡 Medium priority: {total_medium}")
print(f"\n📊 Total issues: {total_critical + total_high + total_medium}")

# Specific critical file list
print("\n" + "=" * 80)
print("TOP PRIORITY FILES TO FIX")
print("=" * 80)

file_issue_count = defaultdict(int)
for category in critical_categories + high_priority:
    if category in issues:
        for issue in issues[category]:
            file_issue_count[issue['file']] += 1

top_files = sorted(file_issue_count.items(), key=lambda x: x[1], reverse=True)[:10]
for rank, (file, count) in enumerate(top_files, 1):
    print(f"{rank}. {file}: {count} issue(s)")

print("\n" + "=" * 80)
