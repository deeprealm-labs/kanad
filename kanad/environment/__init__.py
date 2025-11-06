"""
Environmental Effects Module

Modulates molecular Hamiltonians based on external conditions:
- Temperature: Thermal effects on bonding
- Pressure: Compression and phase transitions
- pH: Protonation states
- Solvent: Solvation and dielectric effects
- Electric fields: External perturbations

These effects integrate with Kanad's governance system to provide
realistic molecular simulations under various conditions.
"""

from kanad.environment.temperature import TemperatureModulator
from kanad.environment.solvent import SolventModulator
from kanad.environment.ph_effects import pHModulator, ProtonationSite
from kanad.environment.pressure import PressureModulator

__all__ = [
    'TemperatureModulator',
    'SolventModulator',
    'pHModulator',
    'ProtonationSite',
    'PressureModulator',
]
