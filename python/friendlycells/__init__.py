"""
friendlycells computes similarity between spatial transcriptomics samples based
on which labels (e.g. cell types) border one another.
"""

from .similarity import spatial_similarity

__all__ = ["spatial_similarity"]
