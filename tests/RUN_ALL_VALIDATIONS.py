"""
Master Test Runner for Kanad Framework

Runs all validation and comparison tests in sequence:
1. Framework Validation Suite
2. Use Case Demonstrations
3. Quantum vs Traditional Comparison
4. MD Benchmarks

Generates comprehensive report at the end.

Run with: python tests/RUN_ALL_VALIDATIONS.py
"""

import sys
import os
import subprocess
import logging
import json
from datetime import datetime
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


class MasterValidator:
    """Run all validation tests."""

    def __init__(self):
        self.tests_dir = Path(__file__).parent
        self.results = {
            'timestamp': datetime.now().isoformat(),
            'test_suites': []
        }

    def run_all_tests(self):
        """Run all test suites."""
        logger.info("="*80)
        logger.info("KANAD FRAMEWORK - MASTER VALIDATION RUNNER")
        logger.info("="*80)
        logger.info(f"Timestamp: {datetime.now().isoformat()}")
        logger.info("="*80)

        # Test suites to run
        suites = [
            {
                'name': 'Framework Validation Suite',
                'script': 'FRAMEWORK_VALIDATION_SUITE.py',
                'timeout': 600
            },
            {
                'name': 'Use Case Demonstrations',
                'script': 'USE_CASE_DEMONSTRATIONS.py',
                'timeout': 600
            },
            {
                'name': 'Quantum vs Traditional Comparison',
                'script': 'QUANTUM_VS_TRADITIONAL_COMPARISON.py',
                'timeout': 600
            }
        ]

        for suite in suites:
            self.run_suite(suite)

        self.generate_report()

    def run_suite(self, suite):
        """Run a single test suite."""
        logger.info("\n" + "="*80)
        logger.info(f"RUNNING: {suite['name']}")
        logger.info("="*80)

        script_path = self.tests_dir / suite['script']

        if not script_path.exists():
            logger.error(f"Script not found: {script_path}")
            self.results['test_suites'].append({
                'name': suite['name'],
                'status': 'SKIPPED',
                'error': 'Script not found'
            })
            return

        try:
            # Run the test suite
            result = subprocess.run(
                [sys.executable, str(script_path)],
                capture_output=True,
                text=True,
                timeout=suite['timeout'],
                cwd=self.tests_dir.parent
            )

            # Log output
            if result.stdout:
                logger.info(result.stdout)

            if result.stderr:
                logger.error(result.stderr)

            # Record result
            self.results['test_suites'].append({
                'name': suite['name'],
                'status': 'PASS' if result.returncode == 0 else 'FAIL',
                'exit_code': result.returncode,
                'duration': 'N/A'  # Could add timing if needed
            })

            if result.returncode == 0:
                logger.info(f"✅ {suite['name']}: PASSED")
            else:
                logger.error(f"❌ {suite['name']}: FAILED (exit code {result.returncode})")

        except subprocess.TimeoutExpired:
            logger.error(f"⏱️  {suite['name']}: TIMEOUT")
            self.results['test_suites'].append({
                'name': suite['name'],
                'status': 'TIMEOUT',
                'timeout': suite['timeout']
            })

        except Exception as e:
            logger.error(f"❌ {suite['name']}: ERROR - {e}")
            self.results['test_suites'].append({
                'name': suite['name'],
                'status': 'ERROR',
                'error': str(e)
            })

    def generate_report(self):
        """Generate final validation report."""
        logger.info("\n" + "="*80)
        logger.info("MASTER VALIDATION REPORT")
        logger.info("="*80)

        total = len(self.results['test_suites'])
        passed = sum(1 for s in self.results['test_suites'] if s['status'] == 'PASS')
        failed = sum(1 for s in self.results['test_suites'] if s['status'] == 'FAIL')
        errors = sum(1 for s in self.results['test_suites'] if s['status'] in ['ERROR', 'TIMEOUT', 'SKIPPED'])

        logger.info(f"\nTotal Test Suites: {total}")
        logger.info(f"Passed: {passed} ✅")
        logger.info(f"Failed: {failed} ❌")
        logger.info(f"Errors: {errors} ⚠️")

        if total > 0:
            logger.info(f"Success Rate: {passed/total*100:.1f}%")

        logger.info("\nDetailed Results:")
        for suite in self.results['test_suites']:
            status_icon = {
                'PASS': '✅',
                'FAIL': '❌',
                'ERROR': '⚠️',
                'TIMEOUT': '⏱️',
                'SKIPPED': '⏭️'
            }.get(suite['status'], '❓')

            logger.info(f"  {status_icon} {suite['name']}: {suite['status']}")

        # Save report
        report_file = self.tests_dir / f"VALIDATION_REPORT_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(report_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        logger.info(f"\nFull report saved to: {report_file}")

        logger.info("\n" + "="*80)
        logger.info("VALIDATION COMPLETE")
        logger.info("="*80)

        # Exit code based on results
        if failed > 0 or errors > 0:
            logger.warning("\n⚠️  Some tests failed or encountered errors")
            sys.exit(1)
        else:
            logger.info("\n✅ All validations passed!")
            sys.exit(0)


def main():
    """Run master validation."""
    validator = MasterValidator()
    validator.run_all_tests()


if __name__ == '__main__':
    main()
