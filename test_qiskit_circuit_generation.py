"""
Test script to verify real Qiskit circuit generation for all solvers.
"""

from kanad.visualization import create_preview_circuit

# Test molecules
test_configs = [
    {
        "name": "H2 with VQE/UCC",
        "molecule": {
            "atoms": [
                {"symbol": "H", "x": 0.0, "y": 0.0, "z": 0.0},
                {"symbol": "H", "x": 0.0, "y": 0.0, "z": 0.74}
            ],
            "basis": "sto-3g",
            "charge": 0,
            "multiplicity": 1
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "ucc",
            "mapper": "jordan_wigner"
        }
    },
    {
        "name": "H2 with VQE/Hardware-Efficient",
        "molecule": {
            "atoms": [
                {"symbol": "H", "x": 0.0, "y": 0.0, "z": 0.0},
                {"symbol": "H", "x": 0.0, "y": 0.0, "z": 0.74}
            ],
            "basis": "sto-3g",
            "charge": 0,
            "multiplicity": 1
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "hardware_efficient",
            "mapper": "jordan_wigner"
        }
    },
    {
        "name": "H2 with SQD",
        "molecule": {
            "atoms": [
                {"symbol": "H", "x": 0.0, "y": 0.0, "z": 0.0},
                {"symbol": "H", "x": 0.0, "y": 0.0, "z": 0.74}
            ],
            "basis": "sto-3g",
            "charge": 0,
            "multiplicity": 1
        },
        "configuration": {
            "method": "SQD",
            "mapper": "jordan_wigner"
        }
    },
    {
        "name": "LiH with VQE/UCC",
        "molecule": {
            "atoms": [
                {"symbol": "Li", "x": 0.0, "y": 0.0, "z": 0.0},
                {"symbol": "H", "x": 0.0, "y": 0.0, "z": 1.5}
            ],
            "basis": "sto-3g",
            "charge": 0,
            "multiplicity": 1
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "ucc",
            "mapper": "jordan_wigner"
        }
    }
]

def test_circuit_generation():
    """Test circuit generation for all configurations."""
    print("=" * 80)
    print("TESTING REAL QISKIT CIRCUIT GENERATION")
    print("=" * 80)

    all_passed = True

    for test_config in test_configs:
        print(f"\n{'=' * 80}")
        print(f"Test: {test_config['name']}")
        print(f"{'=' * 80}")

        try:
            result = create_preview_circuit(
                test_config["molecule"],
                test_config["configuration"]
            )

            # Check if circuit image was generated
            has_image = result.get('circuit_image') is not None
            has_diagram = result.get('diagram') is not None
            has_stats = result.get('stats') is not None

            print(f"\n✅ SUCCESS:")
            print(f"   - Circuit image generated: {has_image}")
            print(f"   - ASCII diagram generated: {has_diagram}")
            print(f"   - Statistics generated: {has_stats}")
            print(f"   - Qubits: {result.get('n_qubits', 'N/A')}")
            print(f"   - Electrons: {result.get('n_electrons', 'N/A')}")
            print(f"   - Parameters: {result.get('n_parameters', 'N/A')}")

            if has_image:
                image_size = len(result['circuit_image'])
                print(f"   - Image size: {image_size:,} bytes (base64)")
                if image_size < 1000:
                    print("   ⚠️  WARNING: Image size is very small, may be incorrect")
                    all_passed = False
            else:
                print("   ❌ ERROR: No circuit image generated!")
                all_passed = False

            # Show ASCII diagram (first 10 lines)
            if has_diagram:
                lines = result['diagram'].split('\n')[:5]
                print(f"\n   ASCII Preview (first 5 lines):")
                for line in lines:
                    print(f"   {line}")
                if len(result['diagram'].split('\n')) > 5:
                    print(f"   ... ({len(result['diagram'].split('\n')) - 5} more lines)")

        except Exception as e:
            print(f"\n❌ FAILED: {e}")
            import traceback
            traceback.print_exc()
            all_passed = False

    print(f"\n{'=' * 80}")
    if all_passed:
        print("✅ ALL TESTS PASSED")
    else:
        print("❌ SOME TESTS FAILED")
    print(f"{'=' * 80}")

    return all_passed


if __name__ == "__main__":
    test_circuit_generation()
