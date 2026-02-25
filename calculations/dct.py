import numpy as np


def calculate_dct_between_groups(
    ct_experimental: np.ndarray,
    ct_control: np.ndarray,
) -> dict:
    """Simple fold-change between two groups for the same target.

    Per experimental bio rep:
        dCT_rep = CT_experimental_rep - mean(CT_control)
        fold_change_rep = 2^(-dCT_rep)

    Summary values are the mean and SE of the per-rep fold change values.
    Groups are treated as unpaired; each experimental bio rep is independently
    compared to the grand mean of the control group.

    Parameters
    ----------
    ct_experimental : array-like
        CT values for the experimental group, one per biological replicate.
    ct_control : array-like
        CT values for the control group, one per biological replicate.

    Returns
    -------
    dict with keys:
        dCT, dCT_SEM,
        fold_change, fold_change_SE,
        pct_control, pct_control_SE
    """
    ct_experimental = np.asarray(ct_experimental, dtype=float)
    ct_control = np.asarray(ct_control, dtype=float)

    mean_ctrl = np.mean(ct_control)

    # Per-rep dCT and fold change for experimental group
    dct_per_rep = ct_experimental - mean_ctrl
    fc_per_rep = 2.0 ** (-dct_per_rep)

    n = len(fc_per_rep)
    dct = np.mean(dct_per_rep)
    dct_sem = np.std(dct_per_rep, ddof=1) / np.sqrt(n) if n > 1 else 0.0
    fold_change = np.mean(fc_per_rep)
    fold_change_se = np.std(fc_per_rep, ddof=1) / np.sqrt(n) if n > 1 else 0.0

    pct_control = fold_change * 100.0
    pct_control_se = fold_change_se * 100.0

    return {
        "dCT": dct,
        "dCT_SEM": dct_sem,
        "fold_change": fold_change,
        "fold_change_SE": fold_change_se,
        "pct_control": pct_control,
        "pct_control_SE": pct_control_se,
    }
