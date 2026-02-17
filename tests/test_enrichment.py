import numpy as np
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from calculations.enrichment import (
    calculate_percent_input,
    calculate_fold_enrichment,
)


class TestPercentInput:
    def test_basic_calculation(self):
        """Known: 10% input, CT_input=20, CT_IP=18.
        Adjusted_CT_input = 20 - log2(10) = 20 - 3.3219 = 16.678
        % Input = 100 * 2^(16.678 - 18) = 100 * 2^(-1.322) = ~40%
        """
        result = calculate_percent_input(
            ct_ip=np.array([18.0, 18.0, 18.0]),
            ct_input=np.array([20.0, 20.0, 20.0]),
            dilution_factor=10.0,
        )
        expected = 100 * 2 ** (20.0 - np.log2(10) - 18.0)
        assert result["pct_input"] == pytest.approx(expected, rel=1e-4)

    def test_100_percent_dilution_factor_1(self):
        """With dilution_factor=1 and equal CTs, should get 100%."""
        result = calculate_percent_input(
            ct_ip=np.array([20.0, 20.0]),
            ct_input=np.array([20.0, 20.0]),
            dilution_factor=1.0,
        )
        assert result["pct_input"] == pytest.approx(100.0, rel=1e-4)

    def test_higher_ip_ct_gives_lower_pct(self):
        """Higher IP CT (less DNA) should give lower % input."""
        result_low = calculate_percent_input(
            ct_ip=np.array([18.0]), ct_input=np.array([20.0]),
            dilution_factor=10.0,
        )
        result_high = calculate_percent_input(
            ct_ip=np.array([22.0]), ct_input=np.array([20.0]),
            dilution_factor=10.0,
        )
        assert result_low["pct_input"] > result_high["pct_input"]

    def test_error_zero_variance(self):
        """With no variance, SE should be 0."""
        result = calculate_percent_input(
            ct_ip=np.array([18.0, 18.0, 18.0]),
            ct_input=np.array([20.0, 20.0, 20.0]),
            dilution_factor=10.0,
        )
        assert result["pct_input_SE"] == pytest.approx(0.0)


class TestFoldEnrichment:
    def test_equal_ip_and_igg_give_fc_one(self):
        """When IP and IgG have same CTs, fold enrichment = 1."""
        result = calculate_fold_enrichment(
            ct_ip=np.array([20.0, 20.0]),
            ct_input_for_ip=np.array([15.0, 15.0]),
            ct_igg=np.array([20.0, 20.0]),
            ct_input_for_igg=np.array([15.0, 15.0]),
            dilution_factor=10.0,
        )
        assert result["fold_enrichment"] == pytest.approx(1.0, rel=1e-4)
        assert result["ddCT"] == pytest.approx(0.0, abs=1e-10)

    def test_ip_lower_ct_than_igg_gives_enrichment(self):
        """IP with lower CT (more DNA) than IgG should give enrichment > 1."""
        result = calculate_fold_enrichment(
            ct_ip=np.array([18.0, 18.0]),
            ct_input_for_ip=np.array([15.0, 15.0]),
            ct_igg=np.array([22.0, 22.0]),
            ct_input_for_igg=np.array([15.0, 15.0]),
            dilution_factor=10.0,
        )
        assert result["fold_enrichment"] > 1.0

    def test_ip_higher_ct_than_igg_gives_depletion(self):
        """IP with higher CT (less DNA) than IgG should give enrichment < 1."""
        result = calculate_fold_enrichment(
            ct_ip=np.array([25.0, 25.0]),
            ct_input_for_ip=np.array([15.0, 15.0]),
            ct_igg=np.array([20.0, 20.0]),
            ct_input_for_igg=np.array([15.0, 15.0]),
            dilution_factor=10.0,
        )
        assert result["fold_enrichment"] < 1.0

    def test_known_fold_enrichment_value(self):
        """dCT_IP = 18 - (15-3.322) = 18 - 11.678 = 6.322
        dCT_IgG = 22 - (15-3.322) = 22 - 11.678 = 10.322
        ddCT = 6.322 - 10.322 = -4
        fold = 2^4 = 16
        """
        result = calculate_fold_enrichment(
            ct_ip=np.array([18.0]),
            ct_input_for_ip=np.array([15.0]),
            ct_igg=np.array([22.0]),
            ct_input_for_igg=np.array([15.0]),
            dilution_factor=10.0,
        )
        assert result["fold_enrichment"] == pytest.approx(16.0, rel=1e-3)
