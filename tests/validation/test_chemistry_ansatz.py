"""
Validation Test: Chemistry-Efficient Ansatz
============================================
Compare new chemistry-aware ansatz with old hardware-efficient
"""

from kanad import BondFactory
from kanad.solvers import VQESolver

def test_chemistry_ansatz():
    """Test chemistry-efficient ansatz vs hardware-efficient."""
    print("\n" + "="*60)
    print("TEST: Chemistry-Efficient Ansatz")
    print("="*60)

    # Create H2 bond
    bond = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
    print(f"✓ Created H2 bond")

    # Reference: UCC (gold standard)
    print("\n1. UCC Ansatz (reference):")
    solver_ucc = VQESolver(bond, ansatz_type='ucc', mapper_type='jordan_wigner', backend='statevector')
    result_ucc = solver_ucc.solve()
    print(f"   Energy: {result_ucc['energy']:.8f} Ha")
    print(f"   HF Energy: {result_ucc['hf_energy']:.8f} Ha")
    print(f"   Correlation: {result_ucc['correlation_energy']:.8f} Ha")

    # Test: Hardware-Efficient (problematic)
    print("\n2. Hardware-Efficient Ansatz (old):")
    try:
        solver_hw = VQESolver(bond, ansatz_type='hardware_efficient', mapper_type='jordan_wigner', backend='statevector')
        result_hw = solver_hw.solve()
        print(f"   Energy: {result_hw['energy']:.8f} Ha")
        if result_hw['energy'] > result_hw.get('hf_energy', -999):
            print(f"   ⚠️  WARNING: Energy above HF! (unphysical)")
        error_hw = abs(result_hw['energy'] - result_ucc['energy']) * 1000
        print(f"   Error vs UCC: {error_hw:.3f} mHa")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        result_hw = None

    # Test: Chemistry-Efficient (NEW!)
    print("\n3. Chemistry-Efficient Ansatz (NEW):")
    try:
        # Import the new ansatz
        from kanad.ansatze.chemistry_efficient_ansatz import ChemistryEfficientAnsatz

        # Create solver with chemistry-efficient ansatz
        n_qubits = 2 * bond.hamiltonian.n_orbitals
        n_electrons = bond.molecule.n_electrons

        ansatz_chem = ChemistryEfficientAnsatz(
            n_qubits=n_qubits,
            n_electrons=n_electrons,
            n_layers=2,
            include_singles=True,
            include_doubles=True
        )

        from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
        mapper = JordanWignerMapper()

        solver_chem = VQESolver(
            bond=None,  # Low-level API
            hamiltonian=bond.hamiltonian,
            ansatz=ansatz_chem,
            mapper=mapper,
            molecule=bond.molecule,
            backend='statevector'
        )

        result_chem = solver_chem.solve()
        print(f"   Energy: {result_chem['energy']:.8f} Ha")
        print(f"   HF Energy: {result_chem.get('hf_energy', 'N/A')}")

        if 'hf_energy' in result_chem:
            if result_chem['energy'] < result_chem['hf_energy']:
                print(f"   ✓ Energy below HF (physical)")
            else:
                print(f"   ⚠️  Energy above HF")

        error_chem = abs(result_chem['energy'] - result_ucc['energy']) * 1000
        print(f"   Error vs UCC: {error_chem:.3f} mHa")

    except Exception as e:
        print(f"   ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        result_chem = None

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"UCC (reference):              {result_ucc['energy']:.8f} Ha")
    if result_hw:
        print(f"Hardware-Efficient (old):     {result_hw['energy']:.8f} Ha")
    if result_chem:
        print(f"Chemistry-Efficient (NEW):    {result_chem['energy']:.8f} Ha")

    print("\nConclusion:")
    if result_chem and abs(result_chem['energy'] - result_ucc['energy']) < 0.1:
        print("✅ Chemistry-efficient ansatz works reasonably well!")
    else:
        print("⚠️  Chemistry-efficient ansatz needs more tuning")


if __name__ == "__main__":
    test_chemistry_ansatz()
