"""
DFT (Density Functional Theory) Benchmarks

Classical DFT using PySCF with multiple functionals.
Provides another classical baseline for comparison.

DFT Functionals tested:
- LDA (Local Density Approximation)
- PBE (Perdew-Burke-Ernzerhof GGA)
- B3LYP (Hybrid functional)
"""

import sys
import time
import json
from pathlib import Path
from datetime import datetime
import numpy as np

try:
    from pyscf import gto, scf, dft
    PYSCF_AVAILABLE = True
except ImportError:
    print("⚠️  PySCF not installed. Run: pip install pyscf")
    PYSCF_AVAILABLE = False
    sys.exit(1)


class DFTBenchmark:
    """DFT benchmarks using PySCF."""

    def __init__(self, output_dir='benchmarks/results'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {
            'framework': 'PySCF DFT',
            'version': 'Classical DFT',
            'timestamp': datetime.now().isoformat(),
            'molecules': {}
        }

    def benchmark_molecule(self, name, atom1, atom2, distance, charge=0, exact_energy=None):
        """
        Benchmark molecule with DFT methods.

        Args:
            name: Molecule name
            atom1, atom2: Atom symbols
            distance: Bond distance in Angstroms
            charge: Molecular charge
            exact_energy: FCI energy for comparison
        """
        print(f"\n{'='*80}")
        print(f"DFT Benchmark: {name}")
        print(f"{'='*80}")

        mol_results = {
            'name': name,
            'atoms': [atom1, atom2],
            'distance': distance,
            'charge': charge,
            'exact_energy': exact_energy,
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
        print(f"  Basis: sto-3g")
        print(f"  Orbitals: {mol.nao}")
        print(f"  Electrons: {mol.nelectron}")

        # Reference HF energy
        print(f"\n--- Hartree-Fock Reference ---")
        t_start = time.time()
        mf_hf = scf.RHF(mol) if mol.spin == 0 else scf.UHF(mol)
        hf_energy = mf_hf.kernel()
        t_hf = time.time() - t_start

        print(f"  Energy: {hf_energy:.6f} Ha")
        print(f"  Time: {t_hf:.3f}s")

        mol_results['hf_energy'] = float(hf_energy)
        mol_results['hf_time'] = float(t_hf)

        # DFT Functionals
        functionals = [
            ('LDA', 'lda,vwn'),  # Local Density Approximation
            ('PBE', 'pbe,pbe'),  # Perdew-Burke-Ernzerhof
            ('B3LYP', 'b3lyp')   # Hybrid functional
        ]

        for functional_name, xc_code in functionals:
            print(f"\n--- DFT: {functional_name} ---")
            try:
                t_start = time.time()

                # Create DFT object
                mf_dft = dft.RKS(mol) if mol.spin == 0 else dft.UKS(mol)
                mf_dft.xc = xc_code
                dft_energy = mf_dft.kernel()

                t_dft = time.time() - t_start

                # Calculate correlation relative to HF
                corr_vs_hf = (dft_energy - hf_energy) * 1000  # mHa

                # Error vs exact (if available)
                if exact_energy:
                    error_vs_exact = abs(dft_energy - exact_energy) * 1000
                else:
                    error_vs_exact = None

                print(f"  Energy: {dft_energy:.6f} Ha")
                print(f"  Correlation vs HF: {corr_vs_hf:.3f} mHa")
                print(f"  Time: {t_dft:.3f}s")
                if error_vs_exact:
                    print(f"  Error vs Exact: {error_vs_exact:.3f} mHa")

                mol_results['methods'][functional_name] = {
                    'energy': float(dft_energy),
                    'correlation_vs_hf_mha': float(corr_vs_hf),
                    'error_vs_exact_mha': error_vs_exact,
                    'time_seconds': float(t_dft),
                    'converged': bool(mf_dft.converged),
                    'functional': xc_code
                }

            except Exception as e:
                print(f"  ✗ {functional_name} failed: {e}")
                mol_results['methods'][functional_name] = {'error': str(e)}

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
        print("DFT BENCHMARKS (Classical)")
        print("="*80)
        print(f"Timestamp: {self.results['timestamp']}")

        molecules = [
            {'name': 'H2', 'atom1': 'H', 'atom2': 'H', 'distance': 0.74, 'charge': 0, 'exact_energy': -1.1373},
            {'name': 'HeH+', 'atom1': 'He', 'atom2': 'H', 'distance': 0.772, 'charge': 1, 'exact_energy': -2.8510},
            {'name': 'LiH', 'atom1': 'Li', 'atom2': 'H', 'distance': 1.60, 'charge': 0, 'exact_energy': -7.8823}
        ]

        for mol in molecules:
            try:
                self.benchmark_molecule(
                    mol['name'],
                    mol['atom1'],
                    mol['atom2'],
                    mol['distance'],
                    mol.get('charge', 0),
                    mol.get('exact_energy')
                )
            except Exception as e:
                print(f"\n✗ Failed: {e}")
                import traceback
                traceback.print_exc()

        self.save_results()
        self.print_summary()

    def save_results(self):
        """Save results to JSON."""
        output_file = self.output_dir / 'dft_results.json'
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\n✓ Results saved to: {output_file}")

    def print_summary(self):
        """Print summary table."""
        print("\n" + "="*80)
        print("DFT BENCHMARK SUMMARY")
        print("="*80)

        print(f"\n{'Molecule':<10} {'Method':<10} {'Energy (Ha)':<15} {'Corr vs HF (mHa)':<18} {'Time (s)':<10}")
        print("-" * 80)

        for mol_name, mol_data in self.results['molecules'].items():
            # HF reference
            if 'hf_energy' in mol_data:
                print(f"{mol_name:<10} {'HF':<10} {mol_data['hf_energy']:<15.6f} {0.0:<18.3f} {mol_data['hf_time']:<10.3f}")

            # DFT methods
            for method, data in mol_data.get('methods', {}).items():
                if 'energy' in data:
                    energy = data['energy']
                    corr = data.get('correlation_vs_hf_mha', 0.0)
                    time_s = data.get('time_seconds', 0.0)
                    print(f"{'':<10} {method:<10} {energy:<15.6f} {corr:<18.3f} {time_s:<10.3f}")

        print("-" * 80)

        # Accuracy summary
        print("\nACCURACY SUMMARY (vs Exact/FCI):")
        for mol_name, mol_data in self.results['molecules'].items():
            print(f"\n{mol_name}:")
            if 'exact_energy' in mol_data and mol_data['exact_energy']:
                exact = mol_data['exact_energy']
                for method, data in mol_data.get('methods', {}).items():
                    if 'error_vs_exact_mha' in data and data['error_vs_exact_mha']:
                        error = data['error_vs_exact_mha']
                        print(f"  {method}: {error:.3f} mHa")


def main():
    """Run DFT benchmarks."""
    if not PYSCF_AVAILABLE:
        return

    benchmark = DFTBenchmark()
    benchmark.run_all_benchmarks()

    print("\n" + "="*80)
    print("✓ DFT benchmarking complete!")
    print("="*80)


if __name__ == '__main__':
    main()
