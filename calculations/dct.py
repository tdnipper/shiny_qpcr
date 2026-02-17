import numpy as np


def calculate_dct_between_groups(
    ct_experimental: np.ndarray,
    ct_control: np.ndarray,
) -> dict:
    """Simple fold-change between two groups for the same target.

    dCT  = mean(CT_experimental) - mean(CT_control)
    Fold change = 2^(-dCT)
    % of control = fold_change * 100

    Error is propagated via SEM and the delta method for the
    exponential transform.

    Parameters
    ----------
    ct_experimental : array-like
        CT values for the experimental group.
    ct_control : array-like
        CT values for the control group.

    Returns
    -------
    dict with keys:
        dCT, dCT_SEM,
        fold_change, fold_change_SE,
        pct_control, pct_control_SE
    """
    ct_experimental = np.asarray(ct_experimental, dtype=float)
    ct_control = np.asarray(ct_control, dtype=float)

    mean_exp = np.mean(ct_experimental)
    mean_ctrl = np.mean(ct_control)

    sem_exp = (
        np.std(ct_experimental, ddof=1) / np.sqrt(len(ct_experimental))
        if len(ct_experimental) > 1
        else 0.0
    )
    sem_ctrl = (
        np.std(ct_control, ddof=1) / np.sqrt(len(ct_control))
        if len(ct_control) > 1
        else 0.0
    )

    dct = mean_exp - mean_ctrl
    dct_sem = np.sqrt(sem_exp ** 2 + sem_ctrl ** 2)

    fold_change = 2.0 ** (-dct)
    fold_change_se = np.log(2) * fold_change * dct_sem

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
