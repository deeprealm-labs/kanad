"""
Verify CIS Results Against PySCF's Built-in TDA/CIS

This confirms that our CIS implementation gives same results as PySCF.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from pyscf import gto, scf, tdscf


def test_pyscf_tda_h2():
    """Test PySCF's built-in TDA for H2 to get reference value."""
    print("=" * 80)
    print("PySCF Reference: TDA for H2 (STO-3G basis)")
    print("=" * 80)

    # H2 molecule
    mol = gto.Mole()
    mol.atom = [['H', [0, 0, 0]], ['H', [0, 0, 0.74]]]
    mol.basis = 'sto-3g'
    mol.build()

    # HF calculation
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()

    print(f"\nGround State Energy (HF): {mf.e_tot:.8f} Ha")

    # TDA (Tamm-Dancoff Approximation = CIS for closed-shell HF)
    td = tdscf.TDA(mf)
    td.nstates = 3
    td.verbose = 0
    td.kernel()

    print(f"\nExcited States (TDA/CIS):")
    for i, E_ha in enumerate(td.e, 1):
        E_ev = E_ha * 27.2114
        print(f"  S{i}: {E_ev:8.4f} eV  ({E_ha:.8f} Ha)")

    return td.e


def test_kanad_cis_h2():
    """Test Kanad's CIS for H2."""
    print("\n" * 2)
    print("=" * 80)
    print("Kanad CIS for H2 (STO-3G basis)")
    print("=" * 80)

    from kanad.bonds import BondFactory
    from kanad.solvers import ExcitedStatesSolver

    bond = BondFactory.create_bond('H', 'H', distance=0.74)
    solver = ExcitedStatesSolver(bond, method='cis', n_states=3, enable_analysis=False)
    result = solver.solve()

    print(f"\nGround State Energy (HF): {result['ground_state_energy']:.8f} Ha")
    print(f"\nExcited States (CIS):")
    for i, (E_ha, E_ev) in enumerate(zip(result['excitation_energies_ha'], result['excitation_energies_ev']), 1):
        print(f"  S{i}: {E_ev:8.4f} eV  ({E_ha:.8f} Ha)")

    return result['excitation_energies_ha']


if __name__ == '__main__':
    print("\n🔬 VERIFICATION: Kanad CIS vs PySCF TDA\n")

    pyscf_exc = test_pyscf_tda_h2()
    kanad_exc = test_kanad_cis_h2()

    print("\n" * 2)
    print("=" * 80)
    print("COMPARISON")
    print("=" * 80)

    print(f"\n{'State':<8} {'PySCF TDA':<15} {'Kanad CIS':<15} {'Difference':<15} {'Match':<10}")
    print("-" * 80)

    all_match = True
    for i, (pyscf_e, kanad_e) in enumerate(zip(pyscf_exc, kanad_exc), 1):
        pyscf_ev = pyscf_e * 27.2114
        kanad_ev = kanad_e * 27.2114
        diff_ev = abs(pyscf_ev - kanad_ev)

        # Consider match if within 0.1 eV
        match = "✅" if diff_ev < 0.1 else "❌"
        if diff_ev >= 0.1:
            all_match = False

        print(f"S{i:<7} {pyscf_ev:>12.4f} eV  {kanad_ev:>12.4f} eV  {diff_ev:>12.4f} eV  {match}")

    print("=" * 80)

    if all_match:
        print("\n✅ ALL STATES MATCH!")
        print("   Kanad CIS implementation is correct!")
        print("   Now using real ERIs (not placeholders)")
        print("\n📝 Note: H2 CIS excitation energy (~26 eV) is higher than experimental")
        print("         This is EXPECTED with minimal STO-3G basis")
        print("         Experimental value (~11.4 eV) requires larger basis sets")
        sys.exit(0)
    else:
        print("\n❌ SOME STATES DON'T MATCH")
        print("   Check CIS implementation")
        sys.exit(1)