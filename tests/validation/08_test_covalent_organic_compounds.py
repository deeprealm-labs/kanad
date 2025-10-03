#!/usr/bin/env python3
"""
Scientific Test 8: Covalent Bonding with Organic Compounds

This test rigorously validates covalent bonding implementations
with organic compounds and molecular systems.

Scientific Focus:
- Covalent bond formation
- Organic molecule structures
- Molecular geometry
- Electronic structure
- Bond angles and lengths
- Aromaticity and conjugation
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
from kanad.bonds.covalent_bond import CovalentBond
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

class TestCovalentOrganicCompounds:
    """Comprehensive test of covalent bonding with organic compounds."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_simple_organic_molecules(self):
        """Test 8.1: Simple organic molecules."""
        print("="*70)
        print("TEST 8.1: SIMPLE ORGANIC MOLECULES")
        print("="*70)
        
        issues = []
        
        # Test simple organic molecules
        test_molecules = [
            # Methane (CH4)
            ([Atom('C', position=np.array([0.0, 0.0, 0.0])),
              Atom('H', position=np.array([1.09, 1.09, 1.09])),
              Atom('H', position=np.array([-1.09, -1.09, 1.09])),
              Atom('H', position=np.array([-1.09, 1.09, -1.09])),
              Atom('H', position=np.array([1.09, -1.09, -1.09]))], "CH4"),
            
            # Ethane (C2H6)
            ([Atom('C', position=np.array([0.0, 0.0, 0.0])),
              Atom('C', position=np.array([1.54, 0.0, 0.0])),
              Atom('H', position=np.array([-0.51, 0.89, 0.89])),
              Atom('H', position=np.array([-0.51, -0.89, 0.89])),
              Atom('H', position=np.array([-0.51, 0.89, -0.89])),
              Atom('H', position=np.array([-0.51, -0.89, -0.89])),
              Atom('H', position=np.array([2.05, 0.89, 0.89])),
              Atom('H', position=np.array([2.05, -0.89, 0.89])),
              Atom('H', position=np.array([2.05, 0.89, -0.89])),
              Atom('H', position=np.array([2.05, -0.89, -0.89]))], "C2H6"),
            
            # Benzene (C6H6)
            ([Atom('C', position=np.array([1.39, 0.0, 0.0])),
              Atom('C', position=np.array([0.70, 1.21, 0.0])),
              Atom('C', position=np.array([-0.70, 1.21, 0.0])),
              Atom('C', position=np.array([-1.39, 0.0, 0.0])),
              Atom('C', position=np.array([-0.70, -1.21, 0.0])),
              Atom('C', position=np.array([0.70, -1.21, 0.0])),
              Atom('H', position=np.array([2.48, 0.0, 0.0])),
              Atom('H', position=np.array([1.25, 2.16, 0.0])),
              Atom('H', position=np.array([-1.25, 2.16, 0.0])),
              Atom('H', position=np.array([-2.48, 0.0, 0.0])),
              Atom('H', position=np.array([-1.25, -2.16, 0.0])),
              Atom('H', position=np.array([1.25, -2.16, 0.0]))], "C6H6"),
        ]
        
        print(f"Testing simple organic molecules:")
        
        for atoms, name in test_molecules:
            print(f"\nMolecule: {name}")
            
            try:
                # Create molecule (not single bond - molecules have multiple bonds)
                from kanad.core.representations.base_representation import Molecule
                from kanad.core.representations.lcao_representation import LCAORepresentation
                from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

                molecule = Molecule(atoms)
                representation = LCAORepresentation(molecule)
                hamiltonian = CovalentHamiltonian(molecule, representation)

                # Get energy from Hamiltonian
                energy = hamiltonian.total_energy()['total']
                
                print(f"  Atoms: {[atom.symbol for atom in atoms]}")
                print(f"  Energy: {energy:.6f} eV")
                
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
            print(f"\n✅ Simple organic molecules work")
            self.passed_tests += 1
        
        return issues
    
    def test_bond_lengths_and_angles(self):
        """Test 8.2: Bond lengths and angles."""
        print("="*70)
        print("TEST 8.2: BOND LENGTHS AND ANGLES")
        print("="*70)
        
        issues = []
        
        # Test bond lengths and angles
        test_molecules = [
            # Methane (CH4) - tetrahedral geometry
            ([Atom('C', position=np.array([0.0, 0.0, 0.0])),
              Atom('H', position=np.array([1.09, 1.09, 1.09])),
              Atom('H', position=np.array([-1.09, -1.09, 1.09])),
              Atom('H', position=np.array([-1.09, 1.09, -1.09])),
              Atom('H', position=np.array([1.09, -1.09, -1.09]))], "CH4"),
            
            # Water (H2O) - bent geometry
            ([Atom('O', position=np.array([0.0, 0.0, 0.0])),
              Atom('H', position=np.array([0.96, 0.0, 0.0])),
              Atom('H', position=np.array([-0.24, 0.93, 0.0]))], "H2O"),
            
            # Ammonia (NH3) - trigonal pyramidal
            ([Atom('N', position=np.array([0.0, 0.0, 0.0])),
              Atom('H', position=np.array([1.01, 0.0, 0.0])),
              Atom('H', position=np.array([-0.34, 0.94, 0.0])),
              Atom('H', position=np.array([-0.34, -0.47, 0.82]))], "NH3"),
        ]
        
        print(f"Testing bond lengths and angles:")
        
        for atoms, name in test_molecules:
            print(f"\nMolecule: {name}")
            
            try:
                # Create molecule
                from kanad.core.representations.base_representation import Molecule
                from kanad.core.representations.lcao_representation import LCAORepresentation
                from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

                molecule = Molecule(atoms)
                representation = LCAORepresentation(molecule)
                hamiltonian = CovalentHamiltonian(molecule, representation)
                
                # Get bond analysis
                bond_analysis = bond.analyze_bonds()
                
                print(f"  Bond analysis: {bond_analysis}")
                
                # Check if bond analysis is reasonable
                if not isinstance(bond_analysis, dict):
                    issues.append(f"{name}: Bond analysis not a dictionary: {type(bond_analysis)}")
                
                # Check for specific bond types
                if 'bond_lengths' in bond_analysis:
                    bond_lengths = bond_analysis['bond_lengths']
                    print(f"  Bond lengths: {bond_lengths}")
                    
                    # Check if bond lengths are reasonable (in Angstroms)
                    for bond_type, length in bond_lengths.items():
                        if length <= 0 or length > 5.0:  # Unreasonable bond length
                            issues.append(f"{name}: Unreasonable bond length {bond_type}: {length}")
                
                if 'bond_angles' in bond_analysis:
                    bond_angles = bond_analysis['bond_angles']
                    print(f"  Bond angles: {bond_angles}")
                    
                    # Check if bond angles are reasonable (in degrees)
                    for angle_type, angle in bond_angles.items():
                        if angle < 0 or angle > 180:  # Unreasonable bond angle
                            issues.append(f"{name}: Unreasonable bond angle {angle_type}: {angle}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during bond analysis: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Bond lengths and angles are correct")
            self.passed_tests += 1
        
        return issues
    
    def test_aromaticity_and_conjugation(self):
        """Test 8.3: Aromaticity and conjugation."""
        print("="*70)
        print("TEST 8.3: AROMATICITY AND CONJUGATION")
        print("="*70)
        
        issues = []
        
        # Test aromatic and conjugated systems
        test_molecules = [
            # Benzene (C6H6) - aromatic
            ([Atom('C', position=np.array([1.39, 0.0, 0.0])),
              Atom('C', position=np.array([0.70, 1.21, 0.0])),
              Atom('C', position=np.array([-0.70, 1.21, 0.0])),
              Atom('C', position=np.array([-1.39, 0.0, 0.0])),
              Atom('C', position=np.array([-0.70, -1.21, 0.0])),
              Atom('C', position=np.array([0.70, -1.21, 0.0])),
              Atom('H', position=np.array([2.48, 0.0, 0.0])),
              Atom('H', position=np.array([1.25, 2.16, 0.0])),
              Atom('H', position=np.array([-1.25, 2.16, 0.0])),
              Atom('H', position=np.array([-2.48, 0.0, 0.0])),
              Atom('H', position=np.array([-1.25, -2.16, 0.0])),
              Atom('H', position=np.array([1.25, -2.16, 0.0]))], "C6H6"),
            
            # Butadiene (C4H6) - conjugated
            ([Atom('C', position=np.array([0.0, 0.0, 0.0])),
              Atom('C', position=np.array([1.34, 0.0, 0.0])),
              Atom('C', position=np.array([2.68, 0.0, 0.0])),
              Atom('C', position=np.array([4.02, 0.0, 0.0])),
              Atom('H', position=np.array([-0.51, 0.89, 0.0])),
              Atom('H', position=np.array([-0.51, -0.89, 0.0])),
              Atom('H', position=np.array([1.85, 0.89, 0.0])),
              Atom('H', position=np.array([1.85, -0.89, 0.0])),
              Atom('H', position=np.array([2.68, 0.89, 0.0])),
              Atom('H', position=np.array([2.68, -0.89, 0.0]))], "C4H6"),
        ]
        
        print(f"Testing aromaticity and conjugation:")
        
        for atoms, name in test_molecules:
            print(f"\nMolecule: {name}")
            
            try:
                # Create molecule
                from kanad.core.representations.base_representation import Molecule
                from kanad.core.representations.lcao_representation import LCAORepresentation
                from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

                molecule = Molecule(atoms)
                representation = LCAORepresentation(molecule)
                hamiltonian = CovalentHamiltonian(molecule, representation)
                
                # Get aromaticity analysis
                aromaticity = bond.is_aromatic()
                conjugation = bond.get_conjugation_energy()
                
                print(f"  Is aromatic: {aromaticity}")
                print(f"  Conjugation energy: {conjugation:.6f} eV")
                
                # Check if aromaticity makes sense
                if name == "C6H6" and not aromaticity:
                    issues.append(f"{name}: Benzene should be aromatic")
                
                # Check if conjugation energy is reasonable
                if np.isnan(conjugation) or np.isinf(conjugation):
                    issues.append(f"{name}: Invalid conjugation energy (NaN or Inf): {conjugation}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during aromaticity test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Aromaticity and conjugation work")
            self.passed_tests += 1
        
        return issues
    
    def test_electronic_structure(self):
        """Test 8.4: Electronic structure calculations."""
        print("="*70)
        print("TEST 8.4: ELECTRONIC STRUCTURE CALCULATIONS")
        print("="*70)
        
        issues = []
        
        # Test electronic structure
        test_molecules = [
            # Methane (CH4)
            ([Atom('C', position=np.array([0.0, 0.0, 0.0])),
              Atom('H', position=np.array([1.09, 1.09, 1.09])),
              Atom('H', position=np.array([-1.09, -1.09, 1.09])),
              Atom('H', position=np.array([-1.09, 1.09, -1.09])),
              Atom('H', position=np.array([1.09, -1.09, -1.09]))], "CH4"),
            
            # Ethylene (C2H4)
            ([Atom('C', position=np.array([0.0, 0.0, 0.0])),
              Atom('C', position=np.array([1.34, 0.0, 0.0])),
              Atom('H', position=np.array([-0.51, 0.89, 0.0])),
              Atom('H', position=np.array([-0.51, -0.89, 0.0])),
              Atom('H', position=np.array([1.85, 0.89, 0.0])),
              Atom('H', position=np.array([1.85, -0.89, 0.0]))], "C2H4"),
        ]
        
        print(f"Testing electronic structure calculations:")
        
        for atoms, name in test_molecules:
            print(f"\nMolecule: {name}")
            
            try:
                # Create molecule
                from kanad.core.representations.base_representation import Molecule
                from kanad.core.representations.lcao_representation import LCAORepresentation
                from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

                molecule = Molecule(atoms)
                representation = LCAORepresentation(molecule)
                hamiltonian = CovalentHamiltonian(molecule, representation)
                
                # Get electronic structure
                electronic_structure = bond.get_electronic_structure()
                
                print(f"  Electronic structure: {electronic_structure}")
                
                # Check if electronic structure is reasonable
                if not isinstance(electronic_structure, dict):
                    issues.append(f"{name}: Electronic structure not a dictionary: {type(electronic_structure)}")
                
                # Check for specific properties
                if 'homo_energy' in electronic_structure:
                    homo_energy = electronic_structure['homo_energy']
                    print(f"  HOMO energy: {homo_energy:.6f} eV")
                    
                    if np.isnan(homo_energy) or np.isinf(homo_energy):
                        issues.append(f"{name}: Invalid HOMO energy (NaN or Inf): {homo_energy}")
                
                if 'lumo_energy' in electronic_structure:
                    lumo_energy = electronic_structure['lumo_energy']
                    print(f"  LUMO energy: {lumo_energy:.6f} eV")
                    
                    if np.isnan(lumo_energy) or np.isinf(lumo_energy):
                        issues.append(f"{name}: Invalid LUMO energy (NaN or Inf): {lumo_energy}")
                
                if 'band_gap' in electronic_structure:
                    band_gap = electronic_structure['band_gap']
                    print(f"  Band gap: {band_gap:.6f} eV")
                    
                    if band_gap < 0:
                        issues.append(f"{name}: Negative band gap: {band_gap}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during electronic structure test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Electronic structure calculations work")
            self.passed_tests += 1
        
        return issues
    
    def test_molecular_geometry_optimization(self):
        """Test 8.5: Molecular geometry optimization."""
        print("="*70)
        print("TEST 8.5: MOLECULAR GEOMETRY OPTIMIZATION")
        print("="*70)
        
        issues = []
        
        # Test geometry optimization
        test_molecules = [
            # Water (H2O) - start with wrong geometry
            ([Atom('O', position=np.array([0.0, 0.0, 0.0])),
              Atom('H', position=np.array([1.0, 0.0, 0.0])),  # Wrong distance
              Atom('H', position=np.array([1.0, 1.0, 0.0]))], "H2O"),  # Wrong angle
        ]
        
        print(f"Testing molecular geometry optimization:")
        
        for atoms, name in test_molecules:
            print(f"\nMolecule: {name}")
            
            try:
                # Create molecule
                from kanad.core.representations.base_representation import Molecule
                from kanad.core.representations.lcao_representation import LCAORepresentation
                from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

                molecule = Molecule(atoms)
                representation = LCAORepresentation(molecule)
                hamiltonian = CovalentHamiltonian(molecule, representation)
                
                # Get initial geometry
                initial_positions = [atom.position for atom in atoms]
                print(f"  Initial positions: {initial_positions}")
                
                # Optimize geometry
                optimized_bond = bond.optimize_geometry()
                
                # Get optimized geometry
                optimized_positions = [atom.position for atom in optimized_bond.atoms]
                print(f"  Optimized positions: {optimized_positions}")
                
                # Check if geometry changed
                geometry_changed = not np.allclose(initial_positions, optimized_positions, atol=1e-6)
                
                if not geometry_changed:
                    issues.append(f"{name}: Geometry optimization did not change positions")
                
                # Check if optimized geometry is reasonable
                for i, pos in enumerate(optimized_positions):
                    if np.any(np.isnan(pos)) or np.any(np.isinf(pos)):
                        issues.append(f"{name}: Invalid optimized position {i}: {pos}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during geometry optimization: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Molecular geometry optimization works")
            self.passed_tests += 1
        
        return issues
    
    def test_covalent_edge_cases(self):
        """Test 8.6: Covalent bonding edge cases."""
        print("="*70)
        print("TEST 8.6: COVALENT BONDING EDGE CASES")
        print("="*70)
        
        issues = []
        
        # Test edge cases
        edge_cases = [
            # Single atom
            ([Atom('C', position=np.array([0.0, 0.0, 0.0]))], "Single C"),
            
            # Two atoms (no bonds)
            ([Atom('C', position=np.array([0.0, 0.0, 0.0])),
              Atom('C', position=np.array([10.0, 0.0, 0.0]))], "Two C far apart"),
            
            # Very large molecule
            ([Atom('C', position=np.array([i*1.5, 0.0, 0.0])) for i in range(20)], "C20 chain"),
        ]
        
        print(f"Testing covalent bonding edge cases:")
        
        for atoms, name in edge_cases:
            print(f"\nCase: {name}")
            
            try:
                # Create molecule (not single bond - molecules have multiple bonds)
                from kanad.core.representations.base_representation import Molecule
                from kanad.core.representations.lcao_representation import LCAORepresentation
                from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

                molecule = Molecule(atoms)
                representation = LCAORepresentation(molecule)
                hamiltonian = CovalentHamiltonian(molecule, representation)

                # Get energy from Hamiltonian
                energy = hamiltonian.total_energy()['total']
                
                print(f"  Atoms: {[atom.symbol for atom in atoms]}")
                print(f"  Energy: {energy:.6f} eV")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"{name}: Invalid energy (NaN or Inf): {energy}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during covalent edge case test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Covalent bonding edge cases are handled correctly")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all covalent organic compound tests."""
        print("="*70)
        print("SCIENTIFIC TEST 8: COVALENT BONDING WITH ORGANIC COMPOUNDS")
        print("="*70)
        
        # Run all tests
        self.test_simple_organic_molecules()
        self.test_bond_lengths_and_angles()
        self.test_aromaticity_and_conjugation()
        self.test_electronic_structure()
        self.test_molecular_geometry_optimization()
        self.test_covalent_edge_cases()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 8 SUMMARY")
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
    """Run the covalent organic compound test suite."""
    tester = TestCovalentOrganicCompounds()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
