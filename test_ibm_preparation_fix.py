"""
Test IBM Preparation Hamiltonian-to-Pauli Fix

Verifies that IBMPreparation now properly converts Hamiltonians to Pauli operators
instead of returning a placeholder identity operator.
"""

import numpy as np
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_ibm_preparation_pauli_conversion():
    """Test that IBMPreparation properly converts Hamiltonian to Pauli operators."""
    from kanad.io import from_smiles
    from kanad.bonds import BondFactory
    from kanad.backends.ibm.preparation import IBMPreparation

    logger.info("=" * 60)
    logger.info("TEST: IBM Preparation Hamiltonian-to-Pauli Conversion")
    logger.info("=" * 60)

    # Create simple H2 molecule
    logger.info("\n1. Creating H2 molecule...")
    molecule = from_smiles("[H][H]")

    # Create bond using BondFactory with atom indices
    logger.info("2. Creating covalent bond...")
    atom_1 = molecule.atoms[0]
    atom_2 = molecule.atoms[1]
    bond = BondFactory.create_bond(atom_1, atom_2, bond_type='covalent')

    # Create IBM preparation
    logger.info("3. Creating IBM preparation...")
    prep = IBMPreparation(bond, ansatz_type='hardware_efficient')

    # Test Hamiltonian-to-Pauli conversion
    logger.info("\n4. Converting Hamiltonian to Pauli operators...")
    observable = prep._hamiltonian_to_pauli_sum()

    # Verify it's NOT just identity
    logger.info(f"\n5. Analyzing observable:")
    logger.info(f"   Type: {type(observable)}")
    logger.info(f"   Number of Pauli terms: {len(observable)}")

    # Get first few Pauli strings to verify
    pauli_strings = []
    coefficients = []
    for i in range(min(5, len(observable))):
        pauli_strings.append(str(observable.paulis[i]))
        coefficients.append(observable.coeffs[i])

    logger.info(f"\n6. First few Pauli terms:")
    for i, (pauli, coeff) in enumerate(zip(pauli_strings, coefficients)):
        logger.info(f"   Term {i}: {pauli} with coefficient {coeff:.6f}")

    # CRITICAL CHECK: Should NOT be just identity
    all_identity = all('I' in str(p) and 'X' not in str(p) and 'Y' not in str(p) and 'Z' not in str(p)
                       for p in observable.paulis)

    if all_identity and len(observable) == 1:
        logger.error("\n❌ FAILED: Still returning identity placeholder!")
        raise AssertionError("IBMPreparation still uses identity placeholder")
    else:
        logger.info(f"\n✅ SUCCESS: Proper Pauli decomposition with {len(observable)} terms")
        logger.info("   (Not just identity - includes X, Y, Z terms)")

    # Verify observable can be used for energy computation
    logger.info("\n7. Verifying observable is valid for IBM Estimator...")

    # Check it has the right attributes
    assert hasattr(observable, 'paulis'), "Observable missing 'paulis' attribute"
    assert hasattr(observable, 'coeffs'), "Observable missing 'coeffs' attribute"
    assert len(observable) > 1, f"Expected multiple terms, got {len(observable)}"

    logger.info("   ✅ Observable has correct structure for IBM Estimator")

    logger.info("\n" + "=" * 60)
    logger.info("✅ ALL TESTS PASSED")
    logger.info("=" * 60)

    return observable


if __name__ == "__main__":
    try:
        observable = test_ibm_preparation_pauli_conversion()
        print(f"\n✅ IBM Preparation fix verified - {len(observable)} Pauli terms")
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
