import requests
import json

response = requests.post(
    'http://localhost:8000/api/circuits/preview',
    json={
        'molecule': {
            'smiles': 'O',
            'basis': 'sto-3g',
            'charge': 0,
            'multiplicity': 1
        },
        'configuration': {
            'method': 'VQE',
            'ansatz': 'hardware_efficient',
            'mapper': 'jordan_wigner'
        }
    }
)

data = response.json()
preview = data.get('preview', {})

print(f"Success: {data.get('success')}")
print(f"Has circuit_image: {'circuit_image' in preview}")
print(f"circuit_image is not None: {preview.get('circuit_image') is not None}")
if preview.get('circuit_image'):
    print(f"circuit_image length: {len(preview['circuit_image'])} chars")
    print(f"Starts with PNG signature: {preview['circuit_image'][:20]}")
else:
    print("❌ circuit_image is None or missing!")
    print(f"Preview keys: {list(preview.keys())}")
