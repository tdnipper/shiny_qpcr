import numpy as np
import pandas as pd
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from calculations.ddct import calculate_ddct, summarize_ddct


def _make_results(gene_group_ct, ctrl_gene_group_ct,
                  gene_ctrl_cond_ct, ctrl_gene_ctrl_cond_ct):
    """Helper to build the results dict expected by calculate_ddct."""
    return {
        "GeneA_Grp_exp": pd.DataFrame({"CT": gene_group_ct}),
        "HK_Grp_exp": pd.DataFrame({"CT": ctrl_gene_group_ct}),
        "GeneA_Grp_ctrl": pd.DataFrame({"CT": gene_ctrl_cond_ct}),
        "HK_Grp_ctrl": pd.DataFrame({"CT": ctrl_gene_ctrl_cond_ct}),
    }


def test_self_referential_returns_zero():
    """When gene == control_gene and condition == control_condition, ddCT should be 0."""
    results = _make_results(
        [20.0, 20.0], [10.0, 10.0], [20.0, 20.0], [10.0, 10.0]
    )
    ddct = calculate_ddct(
        results, gene="HK", group="Grp", condition="ctrl",
        control_gene="HK", control_condition="ctrl",
    )
    assert ddct == 0


def test_identical_cts_give_zero_ddct():
    """When all CTs are the same, ddCT should be 0."""
    ct = [20.0, 20.0, 20.0]
    results = _make_results(ct, ct, ct, ct)
    ddct = calculate_ddct(
        results, gene="GeneA", group="Grp", condition="exp",
        control_gene="HK", control_condition="ctrl",
    )
    np.testing.assert_allclose(ddct, 0.0, atol=1e-10)


def test_known_ddct_values():
    """Manually verify ddCT computation.

    gene_group = [25], hk_group = [15]  -> dCT_sample = 10
    gene_control = [22], hk_control = [12] -> dCT_control = 10
    ddCT = 10 - 10 = 0
    """
    results = _make_results([25.0], [15.0], [22.0], [12.0])
    ddct = calculate_ddct(
        results, gene="GeneA", group="Grp", condition="exp",
        control_gene="HK", control_condition="ctrl",
    )
    np.testing.assert_allclose(ddct, 0.0, atol=1e-10)


def test_nonzero_ddct():
    """gene_group=25, hk_group=15 -> dCT=10
    gene_control=20, hk_control=12 -> dCT_ctrl=8
    ddCT = 10 - 8 = 2
    """
    results = _make_results([25.0], [15.0], [20.0], [12.0])
    ddct = calculate_ddct(
        results, gene="GeneA", group="Grp", condition="exp",
        control_gene="HK", control_condition="ctrl",
    )
    np.testing.assert_allclose(ddct, 2.0, atol=1e-10)


def test_summarize_ddct_mean_and_foldchange():
    """Mean ddCT and geometric mean fold change are correct."""
    ddct = np.array([1.0, 2.0, 3.0])
    result = summarize_ddct(ddct)
    np.testing.assert_allclose(result["ddct_mean"], 2.0)
    np.testing.assert_allclose(result["foldchange_mean"], 2.0 ** -2.0)


def test_summarize_ddct_sem():
    """SEM is std(ddof=1)/sqrt(n) of per-rep ddCT values."""
    ddct = np.array([1.0, 2.0, 3.0])
    result = summarize_ddct(ddct)
    expected_sem = np.std(ddct, ddof=1) / np.sqrt(3)
    np.testing.assert_allclose(result["ddct_SEM"], expected_sem)


def test_summarize_ddct_foldchange_se():
    """Fold change SE uses delta method: ln(2) * fc_mean * ddct_SEM."""
    ddct = np.array([1.0, 2.0, 3.0])
    result = summarize_ddct(ddct)
    expected_se = np.log(2) * result["foldchange_mean"] * result["ddct_SEM"]
    np.testing.assert_allclose(result["foldchange_se"], expected_se)


def test_summarize_ddct_single_replicate_sem_zero():
    """With a single replicate SE should be 0, not NaN."""
    result = summarize_ddct(np.array([2.0]))
    assert result["ddct_SEM"] == 0.0
    assert result["foldchange_se"] == 0.0


def test_missing_key_raises():
    """If a required key is missing, KeyError should be raised."""
    results = {"GeneA_Grp_exp": pd.DataFrame({"CT": [20.0]})}
    with pytest.raises(KeyError):
        calculate_ddct(
            results, gene="GeneA", group="Grp", condition="exp",
            control_gene="HK", control_condition="ctrl",
        )
