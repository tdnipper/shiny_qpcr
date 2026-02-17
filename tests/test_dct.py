import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from calculations.dct import calculate_dct_between_groups


def test_identical_groups_give_fc_one():
    """If experimental == control CTs, fold change should be 1.0."""
    result = calculate_dct_between_groups(
        ct_experimental=np.array([20.0, 20.0, 20.0]),
        ct_control=np.array([20.0, 20.0, 20.0]),
    )
    assert result["fold_change"] == pytest.approx(1.0)
    assert result["dCT"] == pytest.approx(0.0)
    assert result["pct_control"] == pytest.approx(100.0)


def test_one_cycle_higher_gives_half_fc():
    """1 cycle higher CT in experimental = 2^(-1) = 0.5 fold change."""
    result = calculate_dct_between_groups(
        ct_experimental=np.array([21.0, 21.0, 21.0]),
        ct_control=np.array([20.0, 20.0, 20.0]),
    )
    assert result["fold_change"] == pytest.approx(0.5)
    assert result["dCT"] == pytest.approx(1.0)
    assert result["pct_control"] == pytest.approx(50.0)


def test_one_cycle_lower_gives_double_fc():
    """1 cycle lower CT in experimental = 2^(1) = 2.0 fold change."""
    result = calculate_dct_between_groups(
        ct_experimental=np.array([19.0, 19.0, 19.0]),
        ct_control=np.array([20.0, 20.0, 20.0]),
    )
    assert result["fold_change"] == pytest.approx(2.0)
    assert result["dCT"] == pytest.approx(-1.0)
    assert result["pct_control"] == pytest.approx(200.0)


def test_error_propagation_zero_variance():
    """With zero variance in both groups, SEM should be 0."""
    result = calculate_dct_between_groups(
        ct_experimental=np.array([20.0, 20.0, 20.0]),
        ct_control=np.array([20.0, 20.0, 20.0]),
    )
    assert result["dCT_SEM"] == pytest.approx(0.0)
    assert result["fold_change_SE"] == pytest.approx(0.0)


def test_error_propagation_with_variance():
    """SEM should be nonzero when there is variance in the data."""
    result = calculate_dct_between_groups(
        ct_experimental=np.array([19.0, 20.0, 21.0]),
        ct_control=np.array([20.0, 20.0, 20.0]),
    )
    assert result["dCT_SEM"] > 0
    assert result["fold_change_SE"] > 0


def test_single_replicate():
    """Single replicate should work (SEM = 0)."""
    result = calculate_dct_between_groups(
        ct_experimental=np.array([21.0]),
        ct_control=np.array([20.0]),
    )
    assert result["fold_change"] == pytest.approx(0.5)
    assert result["dCT_SEM"] == pytest.approx(0.0)
