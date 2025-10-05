#!/usr/bin/env python3
"""
Scientific Test 10: New Solvers and Algorithms

This test rigorously validates the new solvers and algorithms
in the expanded framework.

Scientific Focus:
- Solver accuracy and convergence
- Algorithm performance
- Numerical stability
- Edge cases and error handling
- Cross-solver validation
"""

import numpy as np
import sys
import traceback
from typing import Dict, Any, List, Tuple

# Add the project root to the path
sys.path.insert(0, '/home/mk/deeprealm/kanad')

from kanad.core.atom import Atom
from kanad.bonds.metallic_bond import MetallicBond

class TestNewSolvers:
    """Comprehensive test of new solvers and algorithms."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_solver_availability(self):
        """Test 10.1: Solver availability and basic functionality."""
        print("="*70)
        print("TEST 10.1: SOLVER AVAILABILITY")
        print("="*70)
        
        issues = []
        
        # Test different solvers
        solvers = ['hartree_fock', 'dft', 'ccsd', 'mp2', 'cisd', 'fci', 'vqe', 'qpe']
        
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
        
        print(f"Testing solver availability:")
        
        for solver in solvers:
            print(f"\nSolver: {solver}")
            
            try:
                result = bond.compute_energy(method=solver)
                
                print(f"  Result: {result}")
                
                # Check if result is reasonable
                if not isinstance(result, dict):
                    issues.append(f"{solver}: Result not a dictionary: {type(result)}")
                
                if 'energy' in result:
                    energy = result['energy']
                    if np.isnan(energy) or np.isinf(energy):
                        issues.append(f"{solver}: Invalid energy (NaN or Inf): {energy}")
                
            except Exception as e:
                issues.append(f"{solver}: Exception during computation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ All solvers work")
            self.passed_tests += 1
        
        return issues
    
    def test_solver_accuracy(self):
        """Test 10.2: Solver accuracy comparison."""
        print("="*70)
        print("TEST 10.2: SOLVER ACCURACY COMPARISON")
        print("="*70)
        
        issues = []
        
        # Test solver accuracy
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
        
        print(f"Testing solver accuracy:")
        
        # Get reference energy (tight-binding)
        try:
            ref_result = bond.compute_energy(method='tight_binding')
            ref_energy = ref_result['energy']
            print(f"Reference energy (tight-binding): {ref_energy:.6f} eV")
        except Exception as e:
            issues.append(f"Reference energy calculation failed: {e}")
            return issues
        
        # Test other solvers
        solvers = ['hartree_fock', 'dft', 'ccsd', 'mp2', 'cisd', 'fci', 'vqe']
        
        for solver in solvers:
            print(f"\nSolver: {solver}")
            
            try:
                result = bond.compute_energy(method=solver)
                energy = result['energy']
                
                print(f"  Energy: {energy:.6f} eV")
                
                # Check accuracy
                error = abs(energy - ref_energy)
                relative_error = error / abs(ref_energy) if ref_energy != 0 else error
                
                print(f"  Absolute error: {error:.6f} eV")
                print(f"  Relative error: {relative_error:.6f}")
                
                # Check if error is reasonable
                if relative_error > 0.1:  # 10% tolerance
                    issues.append(f"{solver}: Large relative error: {relative_error:.4f}")
                
            except Exception as e:
                issues.append(f"{solver}: Exception during accuracy test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Solver accuracy is good")
            self.passed_tests += 1
        
        return issues
    
    def test_solver_convergence(self):
        """Test 10.3: Solver convergence behavior."""
        print("="*70)
        print("TEST 10.3: SOLVER CONVERGENCE BEHAVIOR")
        print("="*70)
        
        issues = []
        
        # Test solver convergence
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
        
        print(f"Testing solver convergence:")
        
        solvers = ['hartree_fock', 'dft', 'ccsd', 'mp2', 'cisd', 'fci', 'vqe']
        
        for solver in solvers:
            print(f"\nSolver: {solver}")
            
            try:
                result = bond.compute_energy(method=solver)
                
                print(f"  Result: {result}")
                
                # Check convergence
                if 'converged' in result:
                    converged = result['converged']
                    print(f"  Converged: {converged}")
                    
                    if not isinstance(converged, bool):
                        issues.append(f"{solver}: Converged not boolean: {converged}")
                
                # Check iterations
                if 'iterations' in result:
                    iterations = result['iterations']
                    print(f"  Iterations: {iterations}")
                    
                    if iterations < 0:
                        issues.append(f"{solver}: Negative iterations: {iterations}")
                    elif iterations > 1000:
                        issues.append(f"{solver}: Too many iterations: {iterations}")
                
            except Exception as e:
                issues.append(f"{solver}: Exception during convergence test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Solver convergence is good")
            self.passed_tests += 1
        
        return issues
    
    def test_solver_performance(self):
        """Test 10.4: Solver performance and timing."""
        print("="*70)
        print("TEST 10.4: SOLVER PERFORMANCE")
        print("="*70)
        
        issues = []
        
        # Test solver performance
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
        
        print(f"Testing solver performance:")
        
        solvers = ['hartree_fock', 'dft', 'ccsd', 'mp2', 'cisd', 'fci', 'vqe']
        
        for solver in solvers:
            print(f"\nSolver: {solver}")
            
            try:
                import time
                start_time = time.time()
                result = bond.compute_energy(method=solver)
                end_time = time.time()
                
                computation_time = end_time - start_time
                print(f"  Computation time: {computation_time:.4f} seconds")
                
                # Check if computation time is reasonable
                if computation_time > 60.0:  # More than 60 seconds
                    issues.append(f"{solver}: Computation too slow: {computation_time:.2f}s")
                
            except Exception as e:
                issues.append(f"{solver}: Exception during performance test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Solver performance is good")
            self.passed_tests += 1
        
        return issues
    
    def test_solver_edge_cases(self):
        """Test 10.5: Solver edge cases and error handling."""
        print("="*70)
        print("TEST 10.5: SOLVER EDGE CASES")
        print("="*70)
        
        issues = []
        
        # Test solver edge cases
        edge_cases = [
            # Single atom
            ([Atom('Na', position=np.array([0.0, 0, 0]))], "Single atom"),
            
            # Two atoms
            ([Atom('Na', position=np.array([0.0, 0, 0])),
              Atom('Na', position=np.array([3.0, 0, 0]))], "Two atoms"),
            
            # Large system
            ([Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(20)], "Large system"),
        ]
        
        print(f"Testing solver edge cases:")
        
        for atoms, name in edge_cases:
            print(f"\nCase: {name}")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Test different solvers
                solvers = ['hartree_fock', 'dft', 'ccsd', 'mp2', 'cisd', 'fci', 'vqe']
                
                for solver in solvers:
                    try:
                        result = bond.compute_energy(method=solver)
                        
                        if 'energy' in result:
                            energy = result['energy']
                            if np.isnan(energy) or np.isinf(energy):
                                issues.append(f"{name}, {solver}: Invalid energy (NaN or Inf): {energy}")
                        
                    except Exception as e:
                        issues.append(f"{name}, {solver}: Exception during computation: {e}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during edge case test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Solver edge cases are handled correctly")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all new solver tests."""
        print("="*70)
        print("SCIENTIFIC TEST 10: NEW SOLVERS AND ALGORITHMS")
        print("="*70)
        
        # Run all tests
        self.test_solver_availability()
        self.test_solver_accuracy()
        self.test_solver_convergence()
        self.test_solver_performance()
        self.test_solver_edge_cases()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 10 SUMMARY")
        print("="*70)
        
        print(f"Passed tests: {self.passed_tests}")
        print(f"Failed tests: {self.failed_tests}")
        print(f"Total issues found: {len(self.critical_issues)}")
        
        if self.critical_issues:
            print(f"\n🚨 CRITICAL ISSUES FOUND:")
            for i, issue in enumerate(self.critical_issues, 1):
                print(f"  {i}. {issue}")
        else:
            print(f"\n✅ NO CRITICAL ISSUES FOUND")
        
        return {
            'passed': self.passed_tests,
            'failed': self.failed_tests,
            'issues': self.critical_issues
        }

def main():
    """Run the new solver test suite."""
    tester = TestNewSolvers()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
