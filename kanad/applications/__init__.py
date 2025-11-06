"""
Domain-Specific Applications Layer

Production-ready tools for real-world users (NOT quantum developers):
- Drug Discovery: Compete with SwissADME, Schrödinger
- Alloy Designer: Compete with CALPHAD, Thermo-Calc
- Catalyst Optimizer: Compete with Materials Project
- Materials Scout: Compete with Schrödinger Materials

Each application leverages:
- Quantum accuracy (SQD, Hi-VQE)
- Kanad governance (fast, physically-aware)
- Environmental effects (T, P, pH, solvent)
- Real-time configuration exploration

COMPETITIVE POSITIONING:
========================
We don't compete with Qiskit/Quanto. We compete with SwissADME/Schrödinger.
Target: Domain experts (chemists, materials scientists), NOT quantum developers.

TOTAL ADDRESSABLE MARKET:
==========================
- Drug Discovery: $50-100M/year
- Alloy Design: $50-75M/year
- Catalysis: $30-50M/year
- Materials: $40-60M/year
Total: $170-285M/year
"""

from kanad.applications.drug_discovery import DrugDiscoveryPlatform, DrugCandidate, BindingResult
from kanad.applications.alloy_designer import AlloyDesigner, AlloyCandidate, PhaseDiagram
from kanad.applications.catalyst_optimizer import CatalystOptimizer, CatalystCandidate, ActivityResult
from kanad.applications.materials_scout import MaterialsScout, MaterialCandidate, BandStructure, OpticalSpectrum

__all__ = [
    # Drug Discovery
    'DrugDiscoveryPlatform',
    'DrugCandidate',
    'BindingResult',
    # Alloy Design
    'AlloyDesigner',
    'AlloyCandidate',
    'PhaseDiagram',
    # Catalysis
    'CatalystOptimizer',
    'CatalystCandidate',
    'ActivityResult',
    # Materials
    'MaterialsScout',
    'MaterialCandidate',
    'BandStructure',
    'OpticalSpectrum',
]
