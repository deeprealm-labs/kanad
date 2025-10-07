"""
PySCF Reference Benchmarks (Classical)

Provides classical chemistry baseline for comparison:
- Hartree-Fock (HF)
- MP2 (second-order perturbation theory)
- CCSD (Coupled Cluster Singles Doubles)
- FCI (Full Configuration Interaction - exact within basis)

Same molecules as Kanad benchmarks: H2, HeH+, LiH
"""

import sys
import time
import json
from pathlib import Path
from datetime import datetime
import numpy as np

try:
    from pyscf import gto, scf, mp, cc, fci
    PYSCF_AVAILABLE = True
except ImportError:
    print("⚠️  PySCF not installed. Run: pip install pyscf")
    PYSCF_AVAILABLE = False
    sys.exit(1)


class PySCFBenchmark:
    """Classical chemistry benchmarks using PySCF."""

    def __init__(self, output_dir='benchmarks/results'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {
            'framework': 'PySCF',
            'version': 'Classical Reference',
            'timestamp': datetime.now().isoformat(),
            'molecules': {}
        }

    def benchmark_molecule(self, name, atom1, atom2, distance, charge=0):
        """
        Benchmark molecule with classical methods.

        Args:
            name: Molecule name
            atom1, atom2: Atom symbols
            distance: Bond distance in Angstroms
            charge: Molecular charge
        """
        print(f"\n{'='*80}")
        print(f"PySCF Benchmark: {name}")
        print(f"{'='*80}")

        mol_results = {
            'name': name,
            'atoms': [atom1, atom2],
            'distance': distance,
            'charge': charge,
            'methods': {}
        }

        # Build PySCF molecule
        mol = gto.M(
            atom=f'{atom1} 0 0 0; {atom2} 0 0 {distance}',
            basis='sto-3g',
            charge=charge,
            spin=0 if (sum([self._get_atomic_number(atom1), self._get_atomic_number(atom2)]) - charge) % 2 == 0 else 1
        )

        print(f"  Atoms: {atom1}-{atom2}")
        print(f"  Distance: {distance} Å")
        print(f"  Charge: {charge}")
        print(f"  Basis: sto-3g")
        print(f"  Orbitals: {mol.nao}")
        print(f"  Electrons: {mol.nelectron}")

        # Method 1: Hartree-Fock
        print(f"\n--- Hartree-Fock (HF) ---")
        t_start = time.time()
        mf = scf.RHF(mol) if mol.spin == 0 else scf.UHF(mol)
        hf_energy = mf.kernel()
        t_hf = time.time() - t_start

        print(f"  Energy: {hf_energy:.6f} Ha")
        print(f"  Time: {t_hf:.3f}s")
        print(f"  Converged: {mf.converged}")

        mol_results['methods']['HF'] = {
            'energy': float(hf_energy),
            'time_seconds': float(t_hf),
            'converged': bool(mf.converged)
        }

        # Method 2: MP2
        print(f"\n--- MP2 (Møller-Plesset 2nd order) ---")
        try:
            t_start = time.time()
            mp2 = mp.MP2(mf)
            mp2_energy, _ = mp2.kernel()
            t_mp2 = time.time() - t_start

            mp2_total = hf_energy + mp2_energy
            mp2_corr = mp2_energy * 1000  # mHa

            print(f"  Total Energy: {mp2_total:.6f} Ha")
            print(f"  Correlation: {mp2_corr:.3f} mHa")
            print(f"  Time: {t_mp2:.3f}s")

            mol_results['methods']['MP2'] = {
                'total_energy': float(mp2_total),
                'correlation_energy': float(mp2_energy),
                'correlation_mha': float(mp2_corr),
                'time_seconds': float(t_mp2)
            }
        except Exception as e:
            print(f"  ✗ MP2 failed: {e}")
            mol_results['methods']['MP2'] = {'error': str(e)}

        # Method 3: CCSD
        print(f"\n--- CCSD (Coupled Cluster Singles Doubles) ---")
        try:
            t_start = time.time()
            ccsd_solver = cc.CCSD(mf)
            ccsd_energy, _ = ccsd_solver.kernel()
            t_ccsd = time.time() - t_start

            ccsd_total = hf_energy + ccsd_energy
            ccsd_corr = ccsd_energy * 1000  # mHa

            print(f"  Total Energy: {ccsd_total:.6f} Ha")
            print(f"  Correlation: {ccsd_corr:.3f} mHa")
            print(f"  Time: {t_ccsd:.3f}s")

            mol_results['methods']['CCSD'] = {
                'total_energy': float(ccsd_total),
                'correlation_energy': float(ccsd_energy),
                'correlation_mha': float(ccsd_corr),
                'time_seconds': float(t_ccsd)
            }
        except Exception as e:
            print(f"  ✗ CCSD failed: {e}")
            mol_results['methods']['CCSD'] = {'error': str(e)}

        # Method 4: FCI (exact within basis)
        print(f"\n--- FCI (Full Configuration Interaction - Exact) ---")
        try:
            t_start = time.time()
            fci_solver = fci.FCI(mf)
            fci_energy = fci_solver.kernel()[0]
            t_fci = time.time() - t_start

            fci_corr = (fci_energy - hf_energy) * 1000  # mHa

            print(f"  Energy: {fci_energy:.6f} Ha")
            print(f"  Correlation: {fci_corr:.3f} mHa")
            print(f"  Time: {t_fci:.3f}s")

            mol_results['methods']['FCI'] = {
                'energy': float(fci_energy),
                'correlation_energy': float(fci_energy - hf_energy),
                'correlation_mha': float(fci_corr),
                'time_seconds': float(t_fci)
            }
            mol_results['exact_energy'] = float(fci_energy)

        except Exception as e:
            print(f"  ✗ FCI failed: {e}")
            mol_results['methods']['FCI'] = {'error': str(e)}

        self.results['molecules'][name] = mol_results
        return mol_results

    def _get_atomic_number(self, symbol):
        """Get atomic number from symbol."""
        atomic_numbers = {
            'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5,
            'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10
        }
        return atomic_numbers.get(symbol, 0)

    def run_all_benchmarks(self):
        """Run benchmarks for all molecules."""
        print("="*80)
        print("PYSCF CLASSICAL CHEMISTRY BENCHMARKS")
        print("="*80)
        print(f"Timestamp: {self.results['timestamp']}")

        molecules = [
            {'name': 'H2', 'atom1': 'H', 'atom2': 'H', 'distance': 0.74, 'charge': 0},
            {'name': 'HeH+', 'atom1': 'He', 'atom2': 'H', 'distance': 0.772, 'charge': 1},
            {'name': 'LiH', 'atom1': 'Li', 'atom2': 'H', 'distance': 1.60, 'charge': 0}
        ]

        for mol in molecules:
            try:
                self.benchmark_molecule(
                    mol['name'],
                    mol['atom1'],
                    mol['atom2'],
                    mol['distance'],
                    mol.get('charge', 0)
                )
            except Exception as e:
                print(f"\n✗ Failed: {e}")
                import traceback
                traceback.print_exc()

        self.save_results()
        self.print_summary()

    def save_results(self):
        """Save results to JSON."""
        output_file = self.output_dir / 'pyscf_results.json'
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\n✓ Results saved to: {output_file}")

    def print_summary(self):
        """Print summary table."""
        print("\n" + "="*80)
        print("PYSCF BENCHMARK SUMMARY")
        print("="*80)

        print(f"\n{'Molecule':<10} {'Method':<8} {'Energy (Ha)':<15} {'Correlation (mHa)':<18} {'Time (s)':<10}")
        print("-" * 80)

        for mol_name, mol_data in self.results['molecules'].items():
            for method, data in mol_data['methods'].items():
                if 'energy' in data or 'total_energy' in data:
                    energy = data.get('energy') or data.get('total_energy')
                    corr = data.get('correlation_mha', 0.0)
                    time_s = data.get('time_seconds', 0.0)
                    print(f"{mol_name:<10} {method:<8} {energy:<15.6f} {corr:<18.3f} {time_s:<10.3f}")

        print("-" * 80)


def main():
    """Run PySCF benchmarks."""
    if not PYSCF_AVAILABLE:
        return

    benchmark = PySCFBenchmark()
    benchmark.run_all_benchmarks()

    print("\n" + "="*80)
    print("✓ PySCF benchmarking complete!")
    print("="*80)


if __name__ == '__main__':
    main()
