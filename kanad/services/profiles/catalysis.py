"""
Catalysis Profile - Catalytic Reactions and Transition States

For catalysis research including:
- Reaction pathways
- Transition state analysis
- Activation energies
- Reaction mechanisms
- Surface catalysis

Recommended for: Catalysis research, reaction mechanism studies, enzyme catalysis
"""

CATALYSIS_PROFILE = {
    'name': 'Catalysis',
    'description': 'Transition state analysis and reaction pathway studies for catalytic systems',
    'recommended_for': [
        'Heterogeneous catalysis',
        'Homogeneous catalysis',
        'Enzyme catalysis',
        'Reaction mechanism elucidation',
        'Activation barrier calculations'
    ],
    'use_cases': [
        'Computing activation energies',
        'Locating transition states',
        'Analyzing reaction pathways',
        'Catalyst design and optimization',
        'Understanding reaction mechanisms'
    ],
    'analyses': [
        'frequencies',      # Transition state verification (imaginary frequency) - UNIQUE
        # 'transition_state', # TS analysis (to be implemented) - UNIQUE
    ],
    'default_parameters': {
        'thermochemistry': {
            'temperature': 298.15,
            'pressure': 101325.0,
            'compute_rate_constant': True  # Future: TST rate constants
        },
        'frequencies': {
            'verify_ts': True,  # Check for single imaginary frequency
            'compute_irc': False  # Future: Intrinsic reaction coordinate
        },
        'transition_state': {  # Future module
            'compute_reaction_pathway': True,
            'compute_activation_energy': True,
            'compute_rate_constant': True,
            'temperature_range': [273, 373]  # K
        },
        'bonding': {
            'track_bond_changes': True,
            'compute_reaction_coordinate': True
        }
    },
    'required_data': [
        'final_energy',
        'geometry',
        'frequencies'  # Essential for TS verification
    ],
    'optional_data': [
        'rdm1',
        'hessian',
        'reaction_coordinate',
        'reactant_energy',
        'product_energy',
        'irc_path'
    ],
    'future_features': [
        'Transition state location (TS optimization)',
        'Intrinsic reaction coordinate (IRC) following',
        'Rate constant calculation (Eyring equation)',
        'Reaction pathway analysis',
        'Nudged elastic band (NEB) method',
        'Kinetic isotope effects',
        'Tunneling corrections',
        'Microkinetic modeling'
    ],
    'reaction_analysis': [
        'Activation energy (Ea)',
        'Activation enthalpy (ΔH‡)',
        'Activation entropy (ΔS‡)',
        'Activation free energy (ΔG‡)',
        'Rate constant (k)',
        'Reaction energy (ΔE)',
        'Reaction enthalpy (ΔH)',
        'Reaction free energy (ΔG)',
        'Imaginary frequency (ν‡)'
    ]
}
