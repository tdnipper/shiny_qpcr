import numpy as np
import pandas as pd
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from calculations.ddct import calculate_ddct


def _make_results(gene_group_ct, ctrl_gene_group_ct,
                  gene_ctrl_group_ct, ctrl_gene_ctrl_group_ct):
    """Helper to build the results dict expected by calculate_ddct."""
    return {
        "GeneA_Exp": pd.DataFrame({"CT": gene_group_ct}),
        "HK_Exp": pd.DataFrame({"CT": ctrl_gene_group_ct}),
        "GeneA_Control": pd.DataFrame({"CT": gene_ctrl_group_ct}),
        "HK_Control": pd.DataFrame({"CT": ctrl_gene_ctrl_group_ct}),
    }


def test_self_referential_returns_zero():
    """When gene == control_gene and group == control_group, ddCT should be 0."""
    results = _make_results(
        [20.0, 20.0], [10.0, 10.0], [20.0, 20.0], [10.0, 10.0]
    )
    ddct, error = calculate_ddct(
        results, gene="HK", group="Control",
        control_gene="HK", control_group="Control",
    )
    assert ddct == 0
    assert error == 0


def test_identical_cts_give_zero_ddct():
    """When all CTs are the same, ddCT should be 0."""
    ct = [20.0, 20.0, 20.0]
    results = _make_results(ct, ct, ct, ct)
    ddct, error = calculate_ddct(
        results, gene="GeneA", group="Exp",
        control_gene="HK", control_group="Control",
    )
    np.testing.assert_allclose(ddct, 0.0, atol=1e-10)


def test_known_ddct_values():
    """Manually verify ddCT computation.

    gene_group = [25], hk_group = [15]  -> dCT_sample = 10
    gene_control = [22], hk_control = [12] -> dCT_control = 10
    ddCT = 10 - 10 = 0
    """
    results = _make_results([25.0], [15.0], [22.0], [12.0])
    ddct, error = calculate_ddct(
        results, gene="GeneA", group="Exp",
        control_gene="HK", control_group="Control",
    )
    np.testing.assert_allclose(ddct, 0.0, atol=1e-10)


def test_nonzero_ddct():
    """gene_group=25, hk_group=15 -> dCT=10
    gene_control=20, hk_control=12 -> dCT_ctrl=8
    ddCT = 10 - 8 = 2
    """
    results = _make_results([25.0], [15.0], [20.0], [12.0])
    ddct, error = calculate_ddct(
        results, gene="GeneA", group="Exp",
        control_gene="HK", control_group="Control",
    )
    np.testing.assert_allclose(ddct, 2.0, atol=1e-10)


def test_missing_key_raises():
    """If a required key is missing, KeyError should be raised."""
    results = {"GeneA_Exp": pd.DataFrame({"CT": [20.0]})}
    with pytest.raises(KeyError):
        calculate_ddct(
            results, gene="GeneA", group="Exp",
            control_gene="HK", control_group="Control",
        )
