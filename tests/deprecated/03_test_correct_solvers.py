#!/usr/bin/env python3
"""
Scientific Test 3: Correct Solver Testing (Fixed API Usage)

This test uses the CORRECT API and focuses on REAL solver issues,
not test bugs.

Scientific Focus:
- Proper solver API usage
- Real solver capabilities
- Actual framework functionality
"""

import numpy as np
import sys
import traceback
from typing import Dict, Any, List, Tuple

# Add the project root to the path
sys.path.insert(0, '/home/mk/deeprealm/kanad')

from kanad.core.atom import Atom
from kanad.bonds.metallic_bond import MetallicBond
from kanad.solvers.vqe_solver import VQESolver
from kanad.solvers.excited_states_solver import ExcitedStatesSolver
from kanad.solvers.vibrational_solver import VibrationalSolver
from kanad.solvers.alloy_solver import AlloySolver

class TestCorrectSolvers:
    """Scientific test using CORRECT API calls for solvers."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_available_solvers(self):
        """Test 3.1: Test actually available solvers."""
        print("="*70)
        print("TEST 3.1: AVAILABLE SOLVERS (CORRECT API)")
        print("="*70)
        
        issues = []
        
        # Test actually available solvers (use 2 atoms for speed)
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(2)]
        
        print(f"Testing available solvers:")
        
        # Test VQE solver
        print(f"\nVQE Solver:")
        try:
            vqe_solver = VQESolver()
            result = vqe_solver.solve(atoms, method='vqe')
            print(f"  VQE result: {result}")
            
            if 'energy' in result:
                energy = result['energy']
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"VQE: Invalid energy (NaN or Inf): {energy}")
            
        except Exception as e:
            issues.append(f"VQE solver failed: {e}")
        
        # Test Excited States solver
        print(f"\nExcited States Solver:")
        try:
            excited_solver = ExcitedStatesSolver()
            result = excited_solver.solve(atoms, method='excited_states')
            print(f"  Excited states result: {result}")
            
        except Exception as e:
            issues.append(f"Excited states solver failed: {e}")
        
        # Test Vibrational solver
        print(f"\nVibrational Solver:")
        try:
            vib_solver = VibrationalSolver()
            result = vib_solver.solve(atoms, method='vibrational')
            print(f"  Vibrational result: {result}")
            
        except Exception as e:
            issues.append(f"Vibrational solver failed: {e}")
        
        # Test Alloy solver
        print(f"\nAlloy Solver:")
        try:
            alloy_solver = AlloySolver()
            result = alloy_solver.solve(atoms, method='alloy')
            print(f"  Alloy result: {result}")
            
        except Exception as e:
            issues.append(f"Alloy solver failed: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Available solvers work correctly")
            self.passed_tests += 1
        
        return issues
    
    def test_metallic_bond_solvers(self):
        """Test 3.2: Test metallic bond solvers."""
        print("="*70)
        print("TEST 3.2: METALLIC BOND SOLVERS")
        print("="*70)
        
        issues = []
        
        # Test metallic bond with different solvers
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing metallic bond solvers:")
        
        # Test tight-binding solver
        print(f"\nTight-Binding Solver:")
        try:
            bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
            result = bond.compute_energy(method='tight_binding')
            energy = result['energy']
            print(f"  Energy: {energy:.6f} eV")
            
            if np.isnan(energy) or np.isinf(energy):
                issues.append(f"Tight-binding: Invalid energy (NaN or Inf): {energy}")
            
        except Exception as e:
            issues.append(f"Tight-binding solver failed: {e}")
        
        # Test VQE solver
        print(f"\nVQE Solver:")
        try:
            bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
            result = bond.compute_energy(method='VQE')
            energy = result['energy']
            print(f"  Energy: {energy:.6f} eV")
            
            if np.isnan(energy) or np.isinf(energy):
                issues.append(f"VQE: Invalid energy (NaN or Inf): {energy}")
            
        except Exception as e:
            issues.append(f"VQE solver failed: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Metallic bond solvers work correctly")
            self.passed_tests += 1
        
        return issues
    
    def test_solver_accuracy_comparison(self):
        """Test 3.3: Compare solver accuracy."""
        print("="*70)
        print("TEST 3.3: SOLVER ACCURACY COMPARISON")
        print("="*70)
        
        issues = []
        
        # Test solver accuracy comparison
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing solver accuracy comparison:")
        
        try:
            bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
            
            # Get reference energy
            ref_result = bond.compute_energy(method='tight_binding')
            ref_energy = ref_result['energy']
            
            # Get VQE energy
            vqe_result = bond.compute_energy(method='VQE')
            vqe_energy = vqe_result['energy']
            
            print(f"  Reference energy (tight-binding): {ref_energy:.6f} eV")
            print(f"  VQE energy: {vqe_energy:.6f} eV")
            
            # Check accuracy
            error = abs(vqe_energy - ref_energy)
            relative_error = error / abs(ref_energy) if ref_energy != 0 else error
            
            print(f"  Absolute error: {error:.6f} eV")
            print(f"  Relative error: {relative_error:.6f}")
            
            # Check if error is reasonable
            if relative_error > 0.1:  # 10% tolerance
                issues.append(f"VQE accuracy issue: {relative_error:.4f} relative error")
            else:
                print(f"  ✅ VQE accuracy is good")
            
        except Exception as e:
            issues.append(f"Solver accuracy comparison failed: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Solver accuracy is acceptable")
            self.passed_tests += 1
        
        return issues
    
    def test_solver_performance(self):
        """Test 3.4: Test solver performance."""
        print("="*70)
        print("TEST 3.4: SOLVER PERFORMANCE")
        print("="*70)
        
        issues = []
        
        # Test solver performance
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing solver performance:")
        
        solvers = ['tight_binding', 'VQE']
        
        for solver in solvers:
            print(f"\nSolver: {solver}")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
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
        """Test 3.5: Test solver edge cases."""
        print("="*70)
        print("TEST 3.5: SOLVER EDGE CASES")
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
            ([Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(10)], "Large system"),
        ]
        
        print(f"Testing solver edge cases:")
        
        for atoms, name in edge_cases:
            print(f"\nCase: {name}")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Test different solvers
                solvers = ['tight_binding', 'VQE']
                
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
        """Run all corrected solver tests."""
        print("="*70)
        print("SCIENTIFIC TEST 3: CORRECT SOLVERS (FIXED API)")
        print("="*70)
        
        # Run all tests
        self.test_available_solvers()
        self.test_metallic_bond_solvers()
        self.test_solver_accuracy_comparison()
        self.test_solver_performance()
        self.test_solver_edge_cases()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 3 SUMMARY")
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
    """Run the corrected solver test suite."""
    tester = TestCorrectSolvers()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
