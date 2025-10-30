"""
Materials Science Profile - Solid State and Materials Properties

For materials science and solid-state physics including:
- Electronic structure (band structure, DOS)
- Mechanical properties
- Surface properties
- Alloys and metallurgy
- Semiconductor properties

Recommended for: Materials science, metallurgy, solid-state physics, semiconductor research
"""

MATERIALS_PROFILE = {
    'name': 'Materials Science',
    'description': 'Electronic structure and materials properties for solids and surfaces',
    'recommended_for': [
        'Electronic band structure',
        'Density of states analysis',
        'Mechanical properties (elastic constants)',
        'Alloy composition analysis',
        'Semiconductor materials',
        'Surface chemistry'
    ],
    'use_cases': [
        'Computing band gaps for semiconductors',
        'Analyzing density of states',
        'Predicting mechanical properties',
        'Alloy thermodynamics',
        'Surface adsorption studies',
        'Material stability analysis'
    ],
    'analyses': [
        'dos',              # Density of states - UNIQUE
        'frequencies',      # Vibrational frequencies (required for thermochemistry)
        'thermochemistry',  # Temperature-dependent properties for alloys
    ],
    'default_parameters': {
        'dos': {
            'broadening': 0.1,  # eV
            'n_points': 500,
            'compute_pdos': True,  # Projected DOS
            'fermi_level': 'auto'
        },
        'thermochemistry': {
            'temperature': 298.15,
            'pressure': 101325.0,
            'compute_phase_diagram': False  # Future feature
        },
        'mechanical': {  # Future feature
            'compute_elastic_tensor': True,
            'compute_bulk_modulus': True,
            'compute_shear_modulus': True
        }
    },
    'required_data': [
        'final_energy',
        'geometry'
    ],
    'optional_data': [
        'orbital_energies',
        'kpoints',
        'rdm1',
        'band_structure',
        'cell_parameters'
    ],
    'future_features': [
        'Elastic tensor computation',
        'Bulk and shear moduli',
        'Phase diagram prediction',
        'Surface energy and work function',
        'Alloy formation energies'
    ]
}
