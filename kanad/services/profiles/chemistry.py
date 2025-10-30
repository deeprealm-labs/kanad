"""
Chemistry Profile - Standard Quantum Chemistry Analysis

For general molecular chemistry research including:
- Reaction energetics
- Molecular structure and stability
- Thermochemistry
- Ground state properties

Recommended for: Organic chemistry, inorganic chemistry, computational chemistry
"""

CHEMISTRY_PROFILE = {
    'name': 'Standard Chemistry',
    'description': 'Comprehensive quantum chemistry analysis for ground state molecules',
    'recommended_for': [
        'Reaction energetics and barriers',
        'Molecular structure optimization',
        'Thermochemical calculations',
        'Ground state properties',
        'Organic and inorganic molecules'
    ],
    'use_cases': [
        'Computing reaction enthalpies and free energies',
        'Predicting molecular stability',
        'Analyzing bonding patterns',
        'Calculating molecular properties (dipole, polarizability)'
    ],
    'analyses': [
        'energy',           # Energy decomposition and convergence
        'bonding',          # Bond orders, charges, HOMO-LUMO
        'correlation',      # Correlation energy analysis
        'properties',       # Dipole, polarizability
        'thermochemistry',  # H, S, G at specified temperature
    ],
    'default_parameters': {
        'thermochemistry': {
            'temperature': 298.15,  # K
            'pressure': 101325.0,   # Pa
        },
        'properties': {
            'compute_polarizability': True,
            'field_strength': 0.001  # a.u. for finite field
        }
    },
    'required_data': [
        'final_energy',
        'geometry',
        'nuclear_repulsion'
    ],
    'optional_data': [
        'rdm1',           # For bonding analysis
        'hf_energy',      # For correlation analysis
        'frequencies'     # For thermochemistry
    ]
}
