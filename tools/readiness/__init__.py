"""Readiness label extraction package for Cystic Fibrosis machine learning papers.

Extracts cohort size, split type, external validation presence, and headline metric
with verbatim source sentence grounding to prevent overclaiming of AI readiness.
"""

from tools.readiness.extractor import (
    ReadinessExtractor,
    ReadinessRecord,
    extract_readiness,
    run_benchmark,
)

__all__ = [
    "ReadinessExtractor",
    "ReadinessRecord",
    "extract_readiness",
    "run_benchmark",
]
