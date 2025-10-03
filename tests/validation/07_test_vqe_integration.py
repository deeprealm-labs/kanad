#!/usr/bin/env python3
"""
Scientific Test 7: VQE Integration and Quantum Algorithm Accuracy

This test rigorously validates the VQE integration and quantum
algorithm accuracy in the metallic bonding framework.

Scientific Focus:
- VQE algorithm implementation
- Quantum circuit construction
- Parameter optimization
- Ground state energy accuracy
- Convergence behavior
- Edge cases and numerical stability
"""

import numpy as np
import sys
import traceback
from typing import Dict, Any, List, Tuple
import matplotlib.pyplot as plt

# Add the project root to the path
sys.path.insert(0, '/home/mk/deeprealm/kanad')

from kanad.core.atom import Atom
from kanad.bonds.metallic_bond import MetallicBond
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian

class TestVQEIntegration:
    """Comprehensive test of VQE integration and quantum algorithm accuracy."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_vqe_basic_functionality(self):
        """Test 7.1: Basic VQE functionality."""
        print("="*70)
        print("TEST 7.1: BASIC VQE FUNCTIONALITY")
        print("="*70)
        
        issues = []
        
        # Test VQE on different system sizes
        system_sizes = [2, 4, 6, 8]
        
        print(f"Testing basic VQE functionality:")
        
        for n in system_sizes:
            print(f"\nSystem size N={n}:")
            
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Test VQE computation
                result = bond.compute_energy(method='VQE')
                
                print(f"  VQE result: {result}")
                
                # Check if result is a dictionary
                if not isinstance(result, dict):
                    issues.append(f"N={n}: VQE result not a dictionary: {type(result)}")
                
                # Check required keys
                required_keys = ['energy', 'method', 'converged', 'n_qubits', 'n_layers', 'entanglement', 'iterations']
                for key in required_keys:
                    if key not in result:
                        issues.append(f"N={n}: VQE result missing key '{key}'")
                
                # Check energy
                if 'energy' in result:
                    energy = result['energy']
                    if np.isnan(energy) or np.isinf(energy):
                        issues.append(f"N={n}: VQE energy invalid (NaN or Inf): {energy}")
                
                # Check convergence
                if 'converged' in result:
                    converged = result['converged']
                    if not isinstance(converged, bool):
                        issues.append(f"N={n}: VQE converged not boolean: {converged}")
                
            except Exception as e:
                issues.append(f"N={n}: Exception during VQE computation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Basic VQE functionality works")
            self.passed_tests += 1
        
        return issues
    
    def test_vqe_energy_accuracy(self):
        """Test 7.2: VQE energy accuracy against exact diagonalization."""
        print("="*70)
        print("TEST 7.2: VQE ENERGY ACCURACY")
        print("="*70)
        
        issues = []
        
        # Test VQE accuracy on different system sizes
        system_sizes = [2, 4, 6, 8]
        
        print(f"Testing VQE energy accuracy:")
        
        for n in system_sizes:
            print(f"\nSystem size N={n}:")
            
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get exact energy
                exact_result = bond.compute_energy(method='tight_binding')
                exact_energy = exact_result['energy']
                
                # Get VQE energy
                vqe_result = bond.compute_energy(method='VQE')
                vqe_energy = vqe_result['energy']
                
                print(f"  Exact energy: {exact_energy:.6f} eV")
                print(f"  VQE energy: {vqe_energy:.6f} eV")
                
                # Check energy accuracy
                energy_error = abs(vqe_energy - exact_energy)
                relative_error = energy_error / abs(exact_energy) if exact_energy != 0 else energy_error
                
                print(f"  Absolute error: {energy_error:.6f} eV")
                print(f"  Relative error: {relative_error:.6f}")
                
                # VQE should be accurate to within 1% for small systems
                if relative_error > 0.01:  # 1% tolerance
                    issues.append(f"N={n}: VQE energy inaccurate: {vqe_energy:.6f} vs {exact_energy:.6f} (error={relative_error:.4f})")
                
            except Exception as e:
                issues.append(f"N={n}: Exception during VQE accuracy test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ VQE energy accuracy is good")
            self.passed_tests += 1
        
        return issues
    
    def test_vqe_convergence_behavior(self):
        """Test 7.3: VQE convergence behavior."""
        print("="*70)
        print("TEST 7.3: VQE CONVERGENCE BEHAVIOR")
        print("="*70)
        
        issues = []
        
        # Test VQE convergence on different system sizes
        system_sizes = [2, 4, 6, 8]
        
        print(f"Testing VQE convergence behavior:")
        
        for n in system_sizes:
            print(f"\nSystem size N={n}:")
            
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get VQE result
                vqe_result = bond.compute_energy(method='VQE')
                
                print(f"  VQE result: {vqe_result}")
                
                # Check convergence
                if 'converged' in vqe_result:
                    converged = vqe_result['converged']
                    print(f"  Converged: {converged}")
                    
                    # For small systems, VQE should converge
                    if n <= 4 and not converged:
                        issues.append(f"N={n}: VQE should converge for small systems")
                
                # Check iterations
                if 'iterations' in vqe_result:
                    iterations = vqe_result['iterations']
                    print(f"  Iterations: {iterations}")
                    
                    # Check if iterations is reasonable
                    if iterations < 0:
                        issues.append(f"N={n}: VQE iterations negative: {iterations}")
                    elif iterations > 1000:  # Unreasonably many iterations
                        issues.append(f"N={n}: VQE iterations too many: {iterations}")
                
            except Exception as e:
                issues.append(f"N={n}: Exception during VQE convergence test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ VQE convergence behavior is good")
            self.passed_tests += 1
        
        return issues
    
    def test_vqe_parameter_sensitivity(self):
        """Test 7.4: VQE parameter sensitivity."""
        print("="*70)
        print("TEST 7.4: VQE PARAMETER SENSITIVITY")
        print("="*70)
        
        issues = []
        
        # Test VQE with different parameters
        test_params = [
            {'n_layers': 1, 'entanglement': 'linear'},
            {'n_layers': 2, 'entanglement': 'linear'},
            {'n_layers': 3, 'entanglement': 'linear'},
            {'n_layers': 2, 'entanglement': 'circular'},
            {'n_layers': 2, 'entanglement': 'full'},
        ]
        
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing VQE parameter sensitivity:")
        
        for params in test_params:
            print(f"\nParameters: {params}")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get VQE result with specific parameters
                vqe_result = bond.compute_energy(method='VQE', **params)
                
                print(f"  VQE result: {vqe_result}")
                
                # Check if parameters are respected
                if 'n_layers' in vqe_result and 'n_layers' in params:
                    if vqe_result['n_layers'] != params['n_layers']:
                        issues.append(f"VQE n_layers not respected: {vqe_result['n_layers']} vs {params['n_layers']}")
                
                if 'entanglement' in vqe_result and 'entanglement' in params:
                    if vqe_result['entanglement'] != params['entanglement']:
                        issues.append(f"VQE entanglement not respected: {vqe_result['entanglement']} vs {params['entanglement']}")
                
            except Exception as e:
                issues.append(f"Parameters {params}: Exception during VQE test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ VQE parameter sensitivity is good")
            self.passed_tests += 1
        
        return issues
    
    def test_vqe_quantum_circuit_construction(self):
        """Test 7.5: VQE quantum circuit construction."""
        print("="*70)
        print("TEST 7.5: VQE QUANTUM CIRCUIT CONSTRUCTION")
        print("="*70)
        
        issues = []
        
        # Test VQE quantum circuit construction
        system_sizes = [2, 4, 6, 8]
        
        print(f"Testing VQE quantum circuit construction:")
        
        for n in system_sizes:
            print(f"\nSystem size N={n}:")
            
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get VQE result
                vqe_result = bond.compute_energy(method='VQE')
                
                print(f"  VQE result: {vqe_result}")
                
                # Check quantum circuit properties
                if 'n_qubits' in vqe_result:
                    n_qubits = vqe_result['n_qubits']
                    print(f"  Number of qubits: {n_qubits}")
                    
                    # Check if number of qubits is reasonable
                    if n_qubits <= 0:
                        issues.append(f"N={n}: VQE n_qubits non-positive: {n_qubits}")
                    elif n_qubits > 2 * n:  # Should not need more than 2 qubits per atom
                        issues.append(f"N={n}: VQE n_qubits too large: {n_qubits} > 2*{n}")
                
                if 'n_layers' in vqe_result:
                    n_layers = vqe_result['n_layers']
                    print(f"  Number of layers: {n_layers}")
                    
                    # Check if number of layers is reasonable
                    if n_layers <= 0:
                        issues.append(f"N={n}: VQE n_layers non-positive: {n_layers}")
                    elif n_layers > 10:  # Unreasonably many layers
                        issues.append(f"N={n}: VQE n_layers too many: {n_layers}")
                
                if 'entanglement' in vqe_result:
                    entanglement = vqe_result['entanglement']
                    print(f"  Entanglement: {entanglement}")
                    
                    # Check if entanglement is valid
                    valid_entanglement = ['linear', 'circular', 'full']
                    if entanglement not in valid_entanglement:
                        issues.append(f"N={n}: VQE entanglement invalid: {entanglement}")
                
            except Exception as e:
                issues.append(f"N={n}: Exception during VQE circuit test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ VQE quantum circuit construction is good")
            self.passed_tests += 1
        
        return issues
    
    def test_vqe_edge_cases(self):
        """Test 7.6: VQE edge cases and numerical stability."""
        print("="*70)
        print("TEST 7.6: VQE EDGE CASES")
        print("="*70)
        
        issues = []
        
        # Test VQE edge cases
        edge_cases = [
            ([Atom('Na', position=np.array([0.0, 0, 0]))], "Single atom"),
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('Na', position=np.array([3.0, 0, 0]))], "Two atoms"),
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('K', position=np.array([3.0, 0, 0]))], "Different elements"),
        ]
        
        print(f"Testing VQE edge cases:")
        
        for atoms, name in edge_cases:
            print(f"\nCase: {name}")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get VQE result
                vqe_result = bond.compute_energy(method='VQE')
                
                print(f"  VQE result: {vqe_result}")
                
                # Check if result is reasonable
                if 'energy' in vqe_result:
                    energy = vqe_result['energy']
                    if np.isnan(energy) or np.isinf(energy):
                        issues.append(f"{name}: VQE energy invalid (NaN or Inf): {energy}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during VQE edge case test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ VQE edge cases are handled correctly")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all VQE integration tests."""
        print("="*70)
        print("SCIENTIFIC TEST 7: VQE INTEGRATION AND QUANTUM ALGORITHM ACCURACY")
        print("="*70)
        
        # Run all tests
        self.test_vqe_basic_functionality()
        self.test_vqe_energy_accuracy()
        self.test_vqe_convergence_behavior()
        self.test_vqe_parameter_sensitivity()
        self.test_vqe_quantum_circuit_construction()
        self.test_vqe_edge_cases()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 7 SUMMARY")
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
    """Run the VQE integration test suite."""
    tester = TestVQEIntegration()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
