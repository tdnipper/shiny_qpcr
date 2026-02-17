import numpy as np


def calculate_percent_input(
    ct_ip: np.ndarray,
    ct_input: np.ndarray,
    dilution_factor: float = 10.0,
) -> dict:
    """Percent-Input method for RIP-qPCR / ChIP-qPCR.

    Adjusted_CT_input = CT_input - log2(dilution_factor)
    % Input = 100 * 2^(Adjusted_CT_input - CT_IP)

    Parameters
    ----------
    ct_ip : array-like
        CT values for the IP (immunoprecipitation) sample.
    ct_input : array-like
        CT values for the Input sample.
    dilution_factor : float
        Dilution factor of the input fraction (e.g. 10 for 10 % input).

    Returns
    -------
    dict with keys:
        pct_input, pct_input_SE, delta, delta_SEM
    """
    ct_ip = np.asarray(ct_ip, dtype=float)
    ct_input = np.asarray(ct_input, dtype=float)

    adjustment = np.log2(dilution_factor)
    adjusted_ct_input_mean = np.mean(ct_input) - adjustment

    sem_input = (
        np.std(ct_input, ddof=1) / np.sqrt(len(ct_input))
        if len(ct_input) > 1
        else 0.0
    )
    sem_ip = (
        np.std(ct_ip, ddof=1) / np.sqrt(len(ct_ip))
        if len(ct_ip) > 1
        else 0.0
    )

    delta = adjusted_ct_input_mean - np.mean(ct_ip)
    delta_sem = np.sqrt(sem_input ** 2 + sem_ip ** 2)

    pct_input = 100.0 * (2.0 ** delta)
    pct_input_se = np.log(2) * pct_input * delta_sem

    return {
        "pct_input": pct_input,
        "pct_input_SE": pct_input_se,
        "delta": delta,
        "delta_SEM": delta_sem,
    }


def calculate_fold_enrichment(
    ct_ip: np.ndarray,
    ct_input_for_ip: np.ndarray,
    ct_igg: np.ndarray,
    ct_input_for_igg: np.ndarray,
    dilution_factor: float = 10.0,
) -> dict:
    """Fold enrichment normalized to IgG for RIP-qPCR / ChIP-qPCR.

    dCT_IP  = CT_IP  - Adjusted_CT_input(IP)
    dCT_IgG = CT_IgG - Adjusted_CT_input(IgG)
    ddCT    = dCT_IP - dCT_IgG
    Fold enrichment = 2^(-ddCT)

    Parameters
    ----------
    ct_ip : array-like
        CT values for the antibody-of-interest IP.
    ct_input_for_ip : array-like
        CT values for the Input sample paired with the IP.
    ct_igg : array-like
        CT values for the IgG (negative control) IP.
    ct_input_for_igg : array-like
        CT values for the Input sample paired with IgG.
    dilution_factor : float
        Dilution factor of the input fraction.

    Returns
    -------
    dict with keys:
        ddCT, ddCT_SEM, fold_enrichment, fold_enrichment_SE
    """
    arrays = [ct_ip, ct_input_for_ip, ct_igg, ct_input_for_igg]
    arrays = [np.asarray(a, dtype=float) for a in arrays]
    ct_ip, ct_input_for_ip, ct_igg, ct_input_for_igg = arrays

    adjustment = np.log2(dilution_factor)

    dct_ip = np.mean(ct_ip) - (np.mean(ct_input_for_ip) - adjustment)
    dct_igg = np.mean(ct_igg) - (np.mean(ct_input_for_igg) - adjustment)

    ddct = dct_ip - dct_igg

    sems = []
    for arr in arrays:
        sem = np.std(arr, ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0.0
        sems.append(sem)
    ddct_sem = np.sqrt(sum(s ** 2 for s in sems))

    fold_enrichment = 2.0 ** (-ddct)
    fold_enrichment_se = np.log(2) * fold_enrichment * ddct_sem

    return {
        "ddCT": ddct,
        "ddCT_SEM": ddct_sem,
        "fold_enrichment": fold_enrichment,
        "fold_enrichment_SE": fold_enrichment_se,
    }
