"""
Drug Discovery Profile - Pharmaceutical and ADME Properties

For drug discovery and pharmaceutical research including:
- ADME properties (Absorption, Distribution, Metabolism, Excretion)
- Drug-likeness (Lipinski's Rule of Five)
- Molecular descriptors
- Pharmacokinetic properties
- Bioavailability prediction

Recommended for: Medicinal chemistry, pharmaceutical research, drug design
"""

DRUG_DISCOVERY_PROFILE = {
    'name': 'Drug Discovery',
    'description': 'ADME properties and drug-likeness analysis for pharmaceutical molecules',
    'recommended_for': [
        'Drug design and optimization',
        'ADME property prediction',
        'Lead compound screening',
        'Pharmacokinetic profiling',
        'Drug-likeness assessment'
    ],
    'use_cases': [
        'Predicting oral bioavailability',
        'Computing lipophilicity (logP)',
        'Assessing blood-brain barrier permeability',
        'Evaluating Lipinski Rule of Five compliance',
        'Optimizing molecular properties for drug candidates'
    ],
    'analyses': [
        'adme',             # ADME properties (UNIQUE - not in inline analysis)
    ],
    'default_parameters': {
        'properties': {
            'compute_polarizability': True,
            'solvation_model': 'pcm',  # Future: Implicit solvation
            'solvent': 'water'
        },
        'adme': {  # Future module
            'compute_logp': True,
            'compute_pka': True,
            'compute_solubility': True,
            'compute_permeability': True,
            'lipinski_rules': True,
            'ghose_filter': True,
            'veber_rules': True
        },
        'thermochemistry': {
            'temperature': 310.15,  # Body temperature (37°C)
            'pressure': 101325.0,
            'solvation': True  # Future: Solvation free energy
        }
    },
    'required_data': [
        'final_energy',
        'geometry',
        'molecular_formula'
    ],
    'optional_data': [
        'rdm1',
        'molecular_weight',
        'smiles',
        'solvation_energy'
    ],
    'future_features': [
        'ADME calculator module',
        'Lipophilicity (logP/logD) prediction',
        'pKa calculation',
        'Aqueous solubility prediction',
        'Membrane permeability (Caco-2)',
        'Blood-brain barrier penetration',
        'Plasma protein binding',
        'Metabolic stability',
        'Cytochrome P450 interaction prediction',
        'Toxicity prediction (hERG, Ames)'
    ],
    'molecular_descriptors': [
        'Molecular weight',
        'Heavy atom count',
        'Hydrogen bond donors/acceptors',
        'Rotatable bonds',
        'Topological polar surface area (TPSA)',
        'Molar refractivity',
        'Number of rings',
        'Aromatic rings'
    ]
}
