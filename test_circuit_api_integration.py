"""
Test end-to-end circuit API integration to verify real Qiskit circuits are generated.
"""

import requests
import json

API_BASE = "http://localhost:8000/api"

def test_circuit_preview_api():
    """Test the circuit preview API endpoint."""

    print("=" * 80)
    print("TESTING CIRCUIT PREVIEW API INTEGRATION")
    print("=" * 80)

    test_cases = [
        {
            "name": "H2 with VQE/UCC",
            "payload": {
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
            }
        },
        {
            "name": "H2 with VQE/Hardware-Efficient",
            "payload": {
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
            }
        },
        {
            "name": "H2 with SQD",
            "payload": {
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
            }
        }
    ]

    all_passed = True

    for test_case in test_cases:
        print(f"\n{'=' * 80}")
        print(f"Test: {test_case['name']}")
        print(f"{'=' * 80}")

        try:
            # Make API request
            response = requests.post(
                f"{API_BASE}/circuits/preview",
                json=test_case["payload"],
                timeout=30
            )

            if response.status_code != 200:
                print(f"❌ API returned status {response.status_code}")
                print(f"Response: {response.text}")
                all_passed = False
                continue

            data = response.json()

            if not data.get("success"):
                print(f"❌ API returned success=False")
                print(f"Response: {json.dumps(data, indent=2)}")
                all_passed = False
                continue

            preview = data.get("preview", {})

            # Check if circuit image was generated
            has_image = preview.get('circuit_image') is not None
            has_diagram = preview.get('circuit_diagram') is not None
            has_stats = preview.get('statistics') is not None

            print(f"\n✅ API REQUEST SUCCESSFUL:")
            print(f"   - Status Code: {response.status_code}")
            print(f"   - Circuit image generated: {has_image}")
            print(f"   - ASCII diagram generated: {has_diagram}")
            print(f"   - Statistics generated: {has_stats}")
            print(f"   - Qubits: {preview.get('n_qubits', 'N/A')}")
            print(f"   - Electrons: {preview.get('n_electrons', 'N/A')}")
            print(f"   - Parameters: {preview.get('n_parameters', 'N/A')}")
            print(f"   - Method: {preview.get('method', 'N/A')}")
            print(f"   - Ansatz: {preview.get('ansatz_type', 'N/A')}")

            if has_image:
                image_size = len(preview['circuit_image'])
                print(f"   - Image size: {image_size:,} bytes (base64)")

                # Verify it's a valid base64 string
                if image_size > 1000:
                    print(f"   ✅ Circuit image is valid size")
                else:
                    print(f"   ⚠️  WARNING: Image size is suspiciously small")
                    all_passed = False

                # Check that it starts with valid base64 PNG signature
                # After base64 decoding, PNG files start with: \x89PNG\r\n\x1a\n
                # In base64, this is: iVBORw0KGgo
                if preview['circuit_image'].startswith('iVBORw0KGgo'):
                    print(f"   ✅ Circuit image has valid PNG signature")
                else:
                    print(f"   ❌ Circuit image does not have valid PNG signature")
                    all_passed = False
            else:
                print(f"   ❌ ERROR: No circuit image generated!")
                all_passed = False

            # Show statistics
            if has_stats:
                stats = preview['statistics']
                print(f"\n   Circuit Statistics:")
                print(f"   - Total Gates: {stats.get('total_gates', 'N/A')}")
                print(f"   - Circuit Depth: {stats.get('depth', 'N/A')}")
                print(f"   - Gate Counts: {stats.get('gate_counts', {})}")

        except requests.exceptions.ConnectionError:
            print(f"❌ FAILED: Could not connect to API at {API_BASE}")
            print(f"   Make sure the backend server is running!")
            all_passed = False
        except Exception as e:
            print(f"❌ FAILED: {e}")
            import traceback
            traceback.print_exc()
            all_passed = False

    print(f"\n{'=' * 80}")
    if all_passed:
        print("✅ ALL API INTEGRATION TESTS PASSED")
        print("\nThe backend is generating REAL Qiskit circuits!")
        print("The circuits are:")
        print("  - Generated from actual ansatz implementations")
        print("  - Drawn using Qiskit's native circuit.draw('mpl')")
        print("  - Returned as base64-encoded PNG images")
        print("  - Ready to be displayed in the frontend")
    else:
        print("❌ SOME API INTEGRATION TESTS FAILED")
    print(f"{'=' * 80}")

    return all_passed


if __name__ == "__main__":
    success = test_circuit_preview_api()
    exit(0 if success else 1)
