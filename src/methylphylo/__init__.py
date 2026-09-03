"""Public API for MethylPhylo."""

from .core import (
    AnalysisResult,
    Matrix,
    TransitionProbabilities,
    TreeNode,
    analyze,
    binary_state_probabilities,
    distance_matrix,
    jukes_cantor_distance,
    jukes_cantor_probabilities,
    preprocess,
    read_matrix,
    to_newick,
    upgma,
)

__all__ = [
    "AnalysisResult",
    "Matrix",
    "TransitionProbabilities",
    "TreeNode",
    "analyze",
    "binary_state_probabilities",
    "distance_matrix",
    "jukes_cantor_distance",
    "jukes_cantor_probabilities",
    "preprocess",
    "read_matrix",
    "to_newick",
    "upgma",
]

__version__ = "0.2.0"
