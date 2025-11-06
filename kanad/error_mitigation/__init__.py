"""
Kanad Error Mitigation Module

Cost-effective error mitigation for achieving chemical accuracy on NISQ devices.

Implements:
- Lightweight M3: Matrix-free measurement mitigation (4-qubit readout errors)
- REM: Reference-state Error Mitigation (Hi-VQE specific)
- Calibration caching (amortized cost ~0%)

Expected improvement: 15-20x error reduction with only 15% overhead
Target: <1 mHa error (chemical accuracy) on IBM quantum hardware
"""

from .lite_mitigation import LiteMitigator, create_lite_mitigator

__all__ = [
    'LiteMitigator',
    'create_lite_mitigator',
]
