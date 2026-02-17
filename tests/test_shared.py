import numpy as np
import pandas as pd
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from shared import (
    filter_targets,
    group_ct_data,
    welch_ttest,
    one_sample_ttest,
    correct_pvalues,
    p_to_star,
    export_to_csv,
    export_to_excel,
)


@pytest.fixture
def sample_df():
    return pd.DataFrame({
        "Sample Name": ["A", "A", "B", "B"],
        "Target Name": ["Gene1", "Gene1", "Gene1", "Gene1"],
        "CT": [20.0, 21.0, 22.0, 23.0],
    })


class TestFilterTargets:
    def test_exact_match(self, sample_df):
        result = filter_targets(sample_df, "Gene1")
        assert len(result) == 4

    def test_no_match(self, sample_df):
        result = filter_targets(sample_df, "GeneX")
        assert len(result) == 0


class TestGroupCtData:
    def test_groups_correctly(self, sample_df):
        result = group_ct_data(sample_df)
        assert "Gene1_A" in result
        assert "Gene1_B" in result
        assert len(result["Gene1_A"]) == 2
        assert len(result["Gene1_B"]) == 2


class TestWelchTtest:
    def test_identical_groups(self):
        g1 = np.array([20.0, 20.0, 20.0])
        g2 = np.array([20.0, 20.0, 20.0])
        t, p = welch_ttest(g1, g2)
        assert np.isnan(t)  # 0/0 with no variance

    def test_different_groups(self):
        g1 = np.array([10.0, 11.0, 12.0])
        g2 = np.array([20.0, 21.0, 22.0])
        t, p = welch_ttest(g1, g2)
        assert p < 0.05

    def test_single_replicate_returns_nan(self):
        t, p = welch_ttest(np.array([10.0]), np.array([20.0]))
        assert np.isnan(p)


class TestOneSampleTtest:
    def test_mean_zero(self):
        vals = np.array([0.0, 0.0, 0.0])
        t, p = one_sample_ttest(vals, popmean=0)
        assert np.isnan(t)  # 0/0

    def test_mean_nonzero(self):
        vals = np.array([5.0, 5.1, 4.9])
        t, p = one_sample_ttest(vals, popmean=0)
        assert p < 0.05


class TestCorrectPvalues:
    def test_none_correction(self):
        p = np.array([0.01, 0.04, 0.1])
        result = correct_pvalues(p, method="none")
        np.testing.assert_array_equal(result, p)

    def test_bonferroni(self):
        p = np.array([0.01, 0.04, 0.1])
        result = correct_pvalues(p, method="bonferroni")
        np.testing.assert_allclose(result, [0.03, 0.12, 0.3])

    def test_bonferroni_caps_at_one(self):
        p = np.array([0.5, 0.8])
        result = correct_pvalues(p, method="bonferroni")
        assert all(r <= 1.0 for r in result)

    def test_bh_correction(self):
        p = np.array([0.01, 0.04, 0.1])
        result = correct_pvalues(p, method="bh")
        # BH should be less aggressive than Bonferroni
        bonf = correct_pvalues(p, method="bonferroni")
        assert all(result <= bonf + 1e-10)


class TestPToStar:
    def test_significant(self):
        assert p_to_star(0.0001) == "***"
        assert p_to_star(0.005) == "**"
        assert p_to_star(0.03) == "*"
        assert p_to_star(0.1) == "ns"

    def test_nan(self):
        assert p_to_star(np.nan) == ""


class TestExport:
    def test_csv_export(self, sample_df):
        result = export_to_csv(sample_df)
        assert isinstance(result, bytes)
        assert b"Sample Name" in result

    def test_excel_export(self, sample_df):
        result = export_to_excel(sample_df)
        assert isinstance(result, bytes)
        assert len(result) > 0
