"""Scientific core for methylation-based exploratory sample phylogenies."""

from __future__ import annotations

import csv
import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray

DistanceMetric = Literal["euclidean", "manhattan", "correlation"]


@dataclass(frozen=True)
class Matrix:
    """A CpG-by-sample beta-value matrix."""

    loci: tuple[str, ...]
    samples: tuple[str, ...]
    values: NDArray[np.float64]

    def __post_init__(self) -> None:
        if self.values.ndim != 2:
            raise ValueError("values must be a two-dimensional array")
        if self.values.shape != (len(self.loci), len(self.samples)):
            raise ValueError("matrix dimensions do not match locus and sample labels")
        if len(set(self.samples)) != len(self.samples):
            raise ValueError("sample names must be unique")


@dataclass(frozen=True)
class PreprocessedMatrix:
    matrix: Matrix
    input_loci: int
    missing_values: int
    missing_rate: float
    dropped_loci: int


@dataclass
class TreeNode:
    """A rooted UPGMA node whose height is measured from its descendant leaves."""

    name: str | None
    height: float
    members: tuple[int, ...]
    children: tuple[TreeNode, TreeNode] | None = None


@dataclass(frozen=True)
class AnalysisResult:
    processed: PreprocessedMatrix
    distances: NDArray[np.float64]
    tree: TreeNode
    newick: str
    metric: DistanceMetric


@dataclass(frozen=True)
class TransitionProbabilities:
    """Forward transition probabilities at a non-negative model distance."""

    distance: float
    same: float
    change: float
    specific_change: float | None = None


def read_matrix(path: str | Path) -> Matrix:
    """Read a CSV/TSV matrix with CpG IDs in column one and samples in columns 2..n."""

    path = Path(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        preview = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(preview, delimiters=",\t")
        except csv.Error as exc:
            raise ValueError("could not determine whether the input is CSV or TSV") from exc
        rows = list(csv.reader(handle, dialect))

    if len(rows) < 2 or len(rows[0]) < 4:
        raise ValueError("matrix requires a header, CpG IDs, and at least three samples")

    header = [cell.strip() for cell in rows[0]]
    samples = tuple(header[1:])
    if any(not sample for sample in samples):
        raise ValueError("sample names must not be empty")

    loci: list[str] = []
    values: list[list[float]] = []
    missing_tokens = {"", "na", "nan", "null", "."}
    for line_number, row in enumerate(rows[1:], start=2):
        if not row or not any(cell.strip() for cell in row):
            continue
        if len(row) != len(header):
            raise ValueError(f"row {line_number} has {len(row)} columns; expected {len(header)}")
        locus = row[0].strip()
        if not locus:
            raise ValueError(f"row {line_number} has no CpG identifier")
        parsed: list[float] = []
        for cell in row[1:]:
            token = cell.strip()
            if token.lower() in missing_tokens:
                parsed.append(float("nan"))
                continue
            try:
                parsed.append(float(token))
            except ValueError as exc:
                raise ValueError(f"non-numeric beta value {token!r} at {locus}") from exc
        loci.append(locus)
        values.append(parsed)

    if len(loci) < 5:
        raise ValueError("matrix requires at least five CpG loci")
    return Matrix(tuple(loci), samples, np.asarray(values, dtype=np.float64))


def preprocess(matrix: Matrix, top_loci: int = 500) -> PreprocessedMatrix:
    """Validate beta values, median-impute each locus, and retain high-variance loci."""

    if len(matrix.samples) < 3:
        raise ValueError("at least three samples are required")
    observed = matrix.values[~np.isnan(matrix.values)]
    if observed.size == 0:
        raise ValueError("matrix contains no observed beta values")
    if np.any((observed < 0) | (observed > 1)):
        raise ValueError("observed beta values must lie within 0 and 1")

    missing_values = int(np.isnan(matrix.values).sum())
    keep_mask = ~np.all(np.isnan(matrix.values), axis=1)
    kept = matrix.values[keep_mask].copy()
    kept_loci = np.asarray(matrix.loci, dtype=object)[keep_mask]
    dropped_loci = int((~keep_mask).sum())
    if kept.shape[0] < 5:
        raise ValueError("at least five loci must contain observed beta values")

    row_medians = np.nanmedian(kept, axis=1)
    missing_rows, missing_columns = np.where(np.isnan(kept))
    kept[missing_rows, missing_columns] = row_medians[missing_rows]

    variances = np.var(kept, axis=1)
    order = np.argsort(-variances, kind="stable")
    count = min(max(int(top_loci), 5), len(order))
    selected = order[:count]
    processed = Matrix(
        tuple(str(locus) for locus in kept_loci[selected]),
        matrix.samples,
        kept[selected],
    )
    return PreprocessedMatrix(
        matrix=processed,
        input_loci=len(matrix.loci),
        missing_values=missing_values,
        missing_rate=missing_values / matrix.values.size,
        dropped_loci=dropped_loci,
    )


def distance_matrix(matrix: Matrix, metric: DistanceMetric = "euclidean") -> NDArray[np.float64]:
    """Calculate a symmetric sample-by-sample distance matrix."""

    samples = matrix.values.T
    if np.isnan(samples).any():
        raise ValueError("distance calculation requires a matrix without missing values")
    count = samples.shape[0]
    result = np.zeros((count, count), dtype=np.float64)
    for left in range(count):
        for right in range(left + 1, count):
            a = samples[left]
            b = samples[right]
            if metric == "euclidean":
                value = float(np.sqrt(np.mean((a - b) ** 2)))
            elif metric == "manhattan":
                value = float(np.mean(np.abs(a - b)))
            elif metric == "correlation":
                std_a = float(np.std(a))
                std_b = float(np.std(b))
                value = 1.0 - float(np.corrcoef(a, b)[0, 1]) if std_a and std_b else 1.0
            else:
                raise ValueError(f"unsupported distance metric: {metric}")
            result[left, right] = result[right, left] = value
    return result


def upgma(samples: Sequence[str], distances: NDArray[np.float64]) -> TreeNode:
    """Build an average-linkage UPGMA tree from a precomputed distance matrix."""

    sample_names = tuple(samples)
    if distances.shape != (len(sample_names), len(sample_names)):
        raise ValueError("distance matrix dimensions must match sample names")
    if not np.allclose(distances, distances.T) or not np.allclose(np.diag(distances), 0):
        raise ValueError("distance matrix must be symmetric with a zero diagonal")

    clusters = [TreeNode(name, 0.0, (index,)) for index, name in enumerate(sample_names)]

    def average_distance(left: TreeNode, right: TreeNode) -> float:
        block = distances[np.ix_(left.members, right.members)]
        return float(np.mean(block))

    while len(clusters) > 1:
        best: tuple[float, int, int] | None = None
        for left in range(len(clusters)):
            for right in range(left + 1, len(clusters)):
                candidate = (average_distance(clusters[left], clusters[right]), left, right)
                if best is None or candidate < best:
                    best = candidate
        assert best is not None
        value, left_index, right_index = best
        left = clusters[left_index]
        right = clusters[right_index]
        merged = TreeNode(
            name=None,
            height=value / 2.0,
            members=left.members + right.members,
            children=(left, right),
        )
        clusters = [
            node for index, node in enumerate(clusters) if index not in (left_index, right_index)
        ]
        clusters.append(merged)
    return clusters[0]


def to_newick(node: TreeNode, parent_height: float | None = None) -> str:
    """Serialize a TreeNode as Newick with non-negative branch lengths."""

    branch = "" if parent_height is None else f":{max(0.0, parent_height - node.height):.6f}"
    if node.children is None:
        if node.name is None:
            raise ValueError("leaf nodes require names")
        name = (
            node.name
            if _is_safe_newick_label(node.name)
            else f"'{node.name.replace(chr(39), chr(39) * 2)}'"
        )
        return f"{name}{branch}"
    children = ",".join(to_newick(child, node.height) for child in node.children)
    return f"({children}){branch}"


def _is_safe_newick_label(value: str) -> bool:
    return bool(value) and all(character.isalnum() or character in "_.-" for character in value)


def jukes_cantor_probabilities(distance: float) -> TransitionProbabilities:
    """Return JC69 probabilities for a nucleotide after evolutionary distance ``d``.

    ``distance`` is the expected number of substitutions per site. ``change`` is
    the probability of any different nucleotide; ``specific_change`` is the
    probability of one named substitution, such as A→G.
    """

    _validate_distance(distance)
    decay = math.exp(-4.0 * distance / 3.0)
    same = 0.25 + 0.75 * decay
    specific_change = 0.25 - 0.25 * decay
    return TransitionProbabilities(
        distance=distance,
        same=same,
        change=3.0 * specific_change,
        specific_change=specific_change,
    )


def jukes_cantor_distance(p_distance: float) -> float:
    """Correct an observed nucleotide mismatch proportion under JC69.

    The correction is defined for ``0 <= p < 0.75``. At and above 0.75 the
    finite JC69 distance is not identifiable because of substitution saturation.
    """

    if not math.isfinite(p_distance) or not 0.0 <= p_distance < 0.75:
        raise ValueError("p-distance must be finite and satisfy 0 <= p < 0.75")
    return -0.75 * math.log(1.0 - 4.0 * p_distance / 3.0)


def binary_state_probabilities(distance: float) -> TransitionProbabilities:
    """Return symmetric methylated↔unmethylated transition probabilities.

    ``distance`` is the product of the per-state transition rate and time.
    """

    _validate_distance(distance)
    decay = math.exp(-2.0 * distance)
    return TransitionProbabilities(
        distance=distance,
        same=0.5 + 0.5 * decay,
        change=0.5 - 0.5 * decay,
    )


def _validate_distance(distance: float) -> None:
    if not math.isfinite(distance) or distance < 0:
        raise ValueError("model distance must be a finite, non-negative number")


def analyze(
    matrix: Matrix, *, top_loci: int = 500, metric: DistanceMetric = "euclidean"
) -> AnalysisResult:
    """Run the complete validation, feature-selection, distance, and tree workflow."""

    processed = preprocess(matrix, top_loci=top_loci)
    distances = distance_matrix(processed.matrix, metric=metric)
    tree = upgma(processed.matrix.samples, distances)
    return AnalysisResult(processed, distances, tree, f"{to_newick(tree)};", metric)
