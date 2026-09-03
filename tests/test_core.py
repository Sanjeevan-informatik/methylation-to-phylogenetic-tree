from __future__ import annotations

import numpy as np
import pytest

from methylphylo import (
    Matrix,
    analyze,
    binary_state_probabilities,
    distance_matrix,
    jukes_cantor_distance,
    jukes_cantor_probabilities,
    preprocess,
    to_newick,
    upgma,
)


def example_matrix() -> Matrix:
    values = np.array(
        [
            [0.10, 0.12, 0.11, 0.80, 0.78, 0.82],
            [0.20, 0.21, 0.19, 0.70, 0.72, 0.69],
            [0.80, 0.78, 0.81, 0.22, 0.20, 0.24],
            [0.65, 0.66, 0.64, 0.40, 0.42, 0.39],
            [0.30, np.nan, 0.31, 0.60, 0.58, 0.62],
            [0.48, 0.49, 0.47, 0.51, 0.52, 0.50],
        ]
    )
    return Matrix(
        loci=tuple(f"cg{index:08d}" for index in range(values.shape[0])),
        samples=("Control_1", "Control_2", "Control_3", "Case_1", "Case_2", "Case_3"),
        values=values,
    )


def test_preprocess_imputes_and_selects_variable_loci() -> None:
    processed = preprocess(example_matrix(), top_loci=5)
    assert processed.matrix.values.shape == (5, 6)
    assert not np.isnan(processed.matrix.values).any()
    assert processed.missing_values == 1
    assert processed.missing_rate == pytest.approx(1 / 36)


@pytest.mark.parametrize("metric", ["euclidean", "manhattan", "correlation"])
def test_distance_matrix_is_symmetric(metric: str) -> None:
    matrix = preprocess(example_matrix(), top_loci=6).matrix
    distances = distance_matrix(matrix, metric=metric)  # type: ignore[arg-type]
    assert np.allclose(distances, distances.T)
    assert np.allclose(np.diag(distances), 0)
    assert distances[0, 1] < distances[0, 3]


def test_upgma_keeps_every_sample_in_newick() -> None:
    matrix = preprocess(example_matrix(), top_loci=6).matrix
    distances = distance_matrix(matrix)
    tree = upgma(matrix.samples, distances)
    serialized = f"{to_newick(tree)};"
    assert serialized.endswith(";")
    assert all(sample in serialized for sample in matrix.samples)
    assert tree.height > 0


def test_complete_analysis_is_reproducible() -> None:
    first = analyze(example_matrix(), top_loci=5, metric="euclidean")
    second = analyze(example_matrix(), top_loci=5, metric="euclidean")
    assert first.newick == second.newick
    assert np.array_equal(first.distances, second.distances)


def test_out_of_range_beta_values_are_rejected() -> None:
    matrix = example_matrix()
    invalid = Matrix(matrix.loci, matrix.samples, matrix.values.copy())
    invalid.values[0, 0] = 1.2
    with pytest.raises(ValueError, match="within 0 and 1"):
        preprocess(invalid)


def test_jc69_probabilities_sum_to_one() -> None:
    probabilities = jukes_cantor_probabilities(0.1)
    assert probabilities.same == pytest.approx(0.9063799893)
    assert probabilities.change == pytest.approx(0.0936200107)
    assert probabilities.specific_change == pytest.approx(probabilities.change / 3)
    assert probabilities.same + probabilities.change == pytest.approx(1.0)


def test_jc69_distance_inverts_observed_difference() -> None:
    corrected = jukes_cantor_distance(0.1)
    probabilities = jukes_cantor_probabilities(corrected)
    assert corrected == pytest.approx(0.1073256327)
    assert probabilities.change == pytest.approx(0.1)


def test_binary_state_probabilities_have_expected_limits() -> None:
    at_zero = binary_state_probabilities(0.0)
    at_large_distance = binary_state_probabilities(20.0)
    assert (at_zero.same, at_zero.change) == pytest.approx((1.0, 0.0))
    assert (at_large_distance.same, at_large_distance.change) == pytest.approx((0.5, 0.5))
    assert at_large_distance.specific_change is None


@pytest.mark.parametrize("invalid", [-0.1, float("nan"), float("inf")])
def test_transition_models_reject_invalid_distance(invalid: float) -> None:
    with pytest.raises(ValueError, match="non-negative"):
        jukes_cantor_probabilities(invalid)


@pytest.mark.parametrize("invalid", [-0.1, 0.75, 1.0, float("nan")])
def test_jc69_rejects_unidentifiable_p_distance(invalid: float) -> None:
    with pytest.raises(ValueError, match="0 <= p < 0.75"):
        jukes_cantor_distance(invalid)
