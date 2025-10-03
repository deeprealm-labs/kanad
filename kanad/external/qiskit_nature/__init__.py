"""
Qiskit Nature Components for Kanad Framework.

================================================================================
ATTRIBUTION
================================================================================

This module contains code adapted from Qiskit Nature:
- Repository: https://github.com/qiskit-community/qiskit-nature
- Copyright: (C) Copyright IBM 2021-2025
- License: Apache License 2.0
- Original Authors: Qiskit Nature Development Team

Qiskit Nature is licensed under the Apache License, Version 2.0.
You may obtain a copy of the license at:
http://www.apache.org/licenses/LICENSE-2.0

================================================================================
MODIFICATIONS
================================================================================

Modified by Kanad Framework developers for integration with the Kanad
governance-based quantum chemistry framework. Key adaptations:

1. Standalone operation (removed Qiskit Nature internal dependencies)
2. Compatibility with Qiskit 2.2+
3. Integration with Kanad's bond-type specific Hamiltonians
4. cuQuantum GPU backend support

All modifications maintain compatibility with Apache License 2.0.

================================================================================
USAGE IN KANAD
================================================================================

These components provide the foundational quantum chemistry operators:
- FermionicOp: Fermionic creation/annihilation operators
- ElectronicIntegrals: One- and two-body integral containers
- JordanWignerMapper: Correct fermionic-to-qubit transformation
- ElectronicEnergy: Electronic Hamiltonian builder

Kanad's innovation layers on top:
- Governance protocols (covalent/ionic/metallic)
- Bond-type specific Hamiltonians
- cuQuantum GPU acceleration
- Material science applications

================================================================================
"""

from .operators import FermionicOp, ElectronicIntegrals, PolynomialTensor
from .mappers import JordanWignerMapper, BravyiKitaevMapper
from .hamiltonians import ElectronicEnergy

__all__ = [
    'FermionicOp',
    'ElectronicIntegrals',
    'PolynomialTensor',
    'JordanWignerMapper',
    'BravyiKitaevMapper',
    'ElectronicEnergy',
]
