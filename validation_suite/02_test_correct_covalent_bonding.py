#!/usr/bin/env python3
"""
Scientific Test 2: Correct Covalent Bonding (Fixed API Usage)

This test uses the CORRECT API for covalent bonding and focuses on
REAL physics issues, not test bugs.

Scientific Focus:
- Proper CovalentBond API usage (atom_1, atom_2)
- Molecule class for multi-atom systems
- Real physics validation
- Actual framework capabilities
"""

import numpy as np
import sys
import traceback
from typing import Dict, Any, List, Tuple

# Add the project root to the path
sys.path.insert(0, '/home/mk/deeprealm/kanad')

from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.core.molecule import Molecule

class TestCorrectCovalentBonding:
    """Scientific test using CORRECT API calls for covalent bonding."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_covalent_bond_basic_functionality(self):
        """Test 2.1: Basic covalent bonding with correct API."""
        print("="*70)
        print("TEST 2.1: BASIC COVALENT BONDING (CORRECT API)")
        print("="*70)
        
        issues = []
        
        # Test basic covalent bonding with CORRECT API
        test_bonds = [
            # H2 molecule
            (Atom('H', position=np.array([0.0, 0.0, 0.0])), 
             Atom('H', position=np.array([0.74, 0.0, 0.0])), "H2"),
            
            # HCl molecule
            (Atom('H', position=np.array([0.0, 0.0, 0.0])), 
             Atom('Cl', position=np.array([1.27, 0.0, 0.0])), "HCl"),
            
            # CO molecule
            (Atom('C', position=np.array([0.0, 0.0, 0.0])), 
             Atom('O', position=np.array([1.13, 0.0, 0.0])), "CO"),
        ]
        
        print(f"Testing basic covalent bonding with correct API:")
        
        for atom1, atom2, name in test_bonds:
            print(f"\nBond: {name}")
            
            try:
                # Use CORRECT API - CovalentBond takes two atoms
                bond = CovalentBond(atom1, atom2, distance=1.5)
                
                # Test basic functionality
                result = bond.compute_energy(method='hartree_fock')
                energy = result['energy']
                
                print(f"  Atoms: {atom1.symbol}-{atom2.symbol}")
                print(f"  Energy: {energy:.6f} eV")
                print(f"  Method: {result.get('method', 'N/A')}")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"{name}: Invalid energy (NaN or Inf): {energy}")
                
                # Check if energy is negative (bound system)
                if energy > 0:
                    issues.append(f"{name}: Positive energy (unbound system): {energy}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during covalent bonding: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Basic covalent bonding works correctly")
            self.passed_tests += 1
        
        return issues
    
    def test_molecule_class_for_multi_atom_systems(self):
        """Test 2.2: Molecule class for multi-atom systems."""
        print("="*70)
        print("TEST 2.2: MOLECULE CLASS FOR MULTI-ATOM SYSTEMS")
        print("="*70)
        
        issues = []
        
        # Test multi-atom systems with Molecule class
        test_molecules = [
            # Methane (CH4)
            ([Atom('C', position=np.array([0.0, 0.0, 0.0])),
              Atom('H', position=np.array([1.09, 1.09, 1.09])),
              Atom('H', position=np.array([-1.09, -1.09, 1.09])),
              Atom('H', position=np.array([-1.09, 1.09, -1.09])),
              Atom('H', position=np.array([1.09, -1.09, -1.09]))], "CH4"),
            
            # Water (H2O)
            ([Atom('O', position=np.array([0.0, 0.0, 0.0])),
              Atom('H', position=np.array([0.96, 0.0, 0.0])),
              Atom('H', position=np.array([-0.24, 0.93, 0.0]))], "H2O"),
        ]
        
        print(f"Testing multi-atom systems with Molecule class:")
        
        for atoms, name in test_molecules:
            print(f"\nMolecule: {name}")
            
            try:
                # Use CORRECT API - Molecule class for multi-atom systems
                molecule = Molecule(atoms)
                
                # Test basic functionality
                result = molecule.compute_energy(method='hartree_fock')
                energy = result['energy']
                
                print(f"  Atoms: {[atom.symbol for atom in atoms]}")
                print(f"  Energy: {energy:.6f} eV")
                print(f"  Method: {result.get('method', 'N/A')}")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"{name}: Invalid energy (NaN or Inf): {energy}")
                
                # Check if energy is negative (bound system)
                if energy > 0:
                    issues.append(f"{name}: Positive energy (unbound system): {energy}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during molecule computation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Multi-atom systems work correctly")
            self.passed_tests += 1
        
        return issues
    
    def test_covalent_bond_properties(self):
        """Test 2.3: Covalent bond properties."""
        print("="*70)
        print("TEST 2.3: COVALENT BOND PROPERTIES")
        print("="*70)
        
        issues = []
        
        # Test covalent bond properties
        atom1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        atom2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
        
        print(f"Testing covalent bond properties:")
        
        try:
            bond = CovalentBond(atom1, atom2, distance=0.74)
            
            # Test bond properties
            bond_length = bond.get_bond_length()
            bond_energy = bond.get_bond_energy()
            bond_order = bond.get_bond_order()
            
            print(f"  Bond length: {bond_length:.6f} Å")
            print(f"  Bond energy: {bond_energy:.6f} eV")
            print(f"  Bond order: {bond_order:.6f}")
            
            # Check if properties are reasonable
            if bond_length <= 0 or bond_length > 5.0:
                issues.append(f"Unreasonable bond length: {bond_length}")
            
            if bond_energy > 0:
                issues.append(f"Positive bond energy (should be negative): {bond_energy}")
            
            if bond_order <= 0 or bond_order > 3:
                issues.append(f"Unreasonable bond order: {bond_order}")
            
        except Exception as e:
            issues.append(f"Exception during bond properties test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Covalent bond properties work correctly")
            self.passed_tests += 1
        
        return issues
    
    def test_covalent_bond_optimization(self):
        """Test 2.4: Covalent bond optimization."""
        print("="*70)
        print("TEST 2.4: COVALENT BOND OPTIMIZATION")
        print("="*70)
        
        issues = []
        
        # Test bond optimization
        atom1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        atom2 = Atom('H', position=np.array([1.5, 0.0, 0.0]))  # Start with wrong distance
        
        print(f"Testing covalent bond optimization:")
        
        try:
            bond = CovalentBond(atom1, atom2, distance=1.5)
            
            # Get initial energy
            initial_result = bond.compute_energy(method='hartree_fock')
            initial_energy = initial_result['energy']
            
            print(f"  Initial energy: {initial_energy:.6f} eV")
            print(f"  Initial distance: {bond.get_bond_length():.6f} Å")
            
            # Optimize bond
            optimized_bond = bond.optimize_geometry()
            
            # Get optimized energy
            optimized_result = optimized_bond.compute_energy(method='hartree_fock')
            optimized_energy = optimized_result['energy']
            
            print(f"  Optimized energy: {optimized_energy:.6f} eV")
            print(f"  Optimized distance: {optimized_bond.get_bond_length():.6f} Å")
            
            # Check if optimization improved energy
            if optimized_energy >= initial_energy:
                issues.append(f"Bond optimization did not improve energy: {initial_energy:.6f} -> {optimized_energy:.6f}")
            
        except Exception as e:
            issues.append(f"Exception during bond optimization test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Covalent bond optimization works correctly")
            self.passed_tests += 1
        
        return issues
    
    def test_covalent_bond_edge_cases(self):
        """Test 2.5: Covalent bond edge cases."""
        print("="*70)
        print("TEST 2.5: COVALENT BOND EDGE CASES")
        print("="*70)
        
        issues = []
        
        # Test edge cases
        edge_cases = [
            # Very short distance
            (Atom('H', position=np.array([0.0, 0.0, 0.0])), 
             Atom('H', position=np.array([0.1, 0.0, 0.0])), "Very short distance"),
            
            # Very long distance
            (Atom('H', position=np.array([0.0, 0.0, 0.0])), 
             Atom('H', position=np.array([5.0, 0.0, 0.0])), "Very long distance"),
            
            # Same atom (should fail)
            (Atom('H', position=np.array([0.0, 0.0, 0.0])), 
             Atom('H', position=np.array([0.0, 0.0, 0.0])), "Same atom"),
        ]
        
        print(f"Testing covalent bond edge cases:")
        
        for atom1, atom2, name in edge_cases:
            print(f"\nCase: {name}")
            
            try:
                bond = CovalentBond(atom1, atom2, distance=1.0)
                
                # Test basic functionality
                result = bond.compute_energy(method='hartree_fock')
                energy = result['energy']
                
                print(f"  Energy: {energy:.6f} eV")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"{name}: Invalid energy (NaN or Inf): {energy}")
                
            except Exception as e:
                if "same atom" in name.lower():
                    print(f"  ✅ Expected error for same atom: {e}")
                else:
                    issues.append(f"{name}: Exception during edge case test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Covalent bond edge cases are handled correctly")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all corrected covalent bonding tests."""
        print("="*70)
        print("SCIENTIFIC TEST 2: CORRECT COVALENT BONDING (FIXED API)")
        print("="*70)
        
        # Run all tests
        self.test_covalent_bond_basic_functionality()
        self.test_molecule_class_for_multi_atom_systems()
        self.test_covalent_bond_properties()
        self.test_covalent_bond_optimization()
        self.test_covalent_bond_edge_cases()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 2 SUMMARY")
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
    """Run the corrected covalent bonding test suite."""
    tester = TestCorrectCovalentBonding()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
