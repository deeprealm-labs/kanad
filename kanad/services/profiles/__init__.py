"""
Domain-specific analysis profiles for different research fields.

Each profile defines:
- A set of analysis modules to run
- Default parameters for each analysis
- Recommended use cases
- Required and optional data
"""

from kanad.services.profiles.chemistry import CHEMISTRY_PROFILE
from kanad.services.profiles.spectroscopy import SPECTROSCOPY_PROFILE
from kanad.services.profiles.materials import MATERIALS_PROFILE
from kanad.services.profiles.drug_discovery import DRUG_DISCOVERY_PROFILE
from kanad.services.profiles.catalysis import CATALYSIS_PROFILE

__all__ = [
    'CHEMISTRY_PROFILE',
    'SPECTROSCOPY_PROFILE',
    'MATERIALS_PROFILE',
    'DRUG_DISCOVERY_PROFILE',
    'CATALYSIS_PROFILE',
]
