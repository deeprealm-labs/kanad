"""
Spectroscopy Profile - Optical and Vibrational Spectroscopy

For spectroscopic analysis including:
- UV-Vis absorption spectroscopy
- Excited state calculations
- Vibrational spectroscopy (IR, Raman)
- Fluorescence and phosphorescence
- NMR spectroscopy (future: ECHOES algorithm)

Recommended for: Photochemistry, spectroscopy, photophysics
"""

SPECTROSCOPY_PROFILE = {
    'name': 'Spectroscopy',
    'description': 'NMR, UV-Vis, IR spectroscopy and excited state analysis',
    'recommended_for': [
        'NMR spectroscopy with Google ECHOES',
        'UV-Vis absorption spectroscopy',
        'IR and Raman spectroscopy',
        'Excited state calculations',
        'Photochemistry and photophysics'
    ],
    'use_cases': [
        'Predicting absorption spectra',
        'Computing excited state energies',
        'Analyzing vibrational frequencies',
        'Understanding photochemical reactions',
        'Chromophore design'
    ],
    'analyses': [
        'nmr',             # NMR spectroscopy with Google ECHOES - UNIQUE
        'frequencies',      # Vibrational modes (IR) - UNIQUE
        'excited_states',   # Electronic excited states - UNIQUE
        'uv_vis',          # UV-Vis absorption spectrum - UNIQUE
        'vibronic',        # Vibronic coupling - UNIQUE (future)
    ],
    'default_parameters': {
        'excited_states': {
            'n_states': 5,
            'method': 'tda',  # Tamm-Dancoff approximation
            'triplets': False
        },
        'uv_vis': {
            'n_states': 5,
            'broadening': 0.3,  # eV
            'wavelength_range': [200, 800]  # nm
        },
        'frequencies': {
            'compute_intensities': True,
            'temperature': 298.15
        }
    },
    'required_data': [
        'final_energy',
        'geometry'
    ],
    'optional_data': [
        'rdm1',
        'excited_states',
        'transition_dipoles',
        'frequencies',
        'hessian'
    ],
    'future_features': [
        'NMR spectroscopy with Google ECHOES algorithm',
        'Raman intensities',
        'Time-resolved spectroscopy',
        'Fluorescence quantum yields'
    ]
}
