"""
Services module for Kanad - Decoupled business logic.

This module provides service classes that operate independently of
the solver execution pipeline, enabling on-demand analysis and
other post-processing tasks.
"""

from kanad.services.analysis_service import AnalysisService, create_analysis_service

__all__ = [
    'AnalysisService',
    'create_analysis_service',
]
