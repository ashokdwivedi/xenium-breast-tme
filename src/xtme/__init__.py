"""xtme - Xenium tumour microenvironment analysis toolkit.

A small, tested library behind a reproducible three-stage pipeline:
QC -> cell typing against an scRNA-seq reference -> spatial organisation.
"""
__version__ = "0.1.0"

from . import annotate, io, plotting, qc, spatial  # noqa: F401

__all__ = ["io", "qc", "annotate", "spatial", "plotting", "__version__"]
