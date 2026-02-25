import numpy as np


def calculate_percent_input(
    ct_ip: np.ndarray,
    ct_input: np.ndarray,
    dilution_factor: float = 10.0,
) -> dict:
    """Percent-Input method for RIP-qPCR / ChIP-qPCR.

    ct_ip and ct_input must be paired arrays of equal length — one element per
    biological replicate, already averaged over technical replicates.

    Per bio rep:
        dCT = CT_IP - (CT_input - log2(dilution_factor))
        % Input = 100 * 2^(-dCT)

    Summary values are the mean and SE of the per-bio-rep % Input values.

    Parameters
    ----------
    ct_ip : array-like
        CT values for the IP sample, one per biological replicate.
    ct_input : array-like
        CT values for the Input sample, paired to ct_ip by biological replicate.
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

    # Per-bio-rep dCT and % Input
    delta_per_rep = (ct_input - adjustment) - ct_ip
    pct_input_per_rep = 100.0 * (2.0 ** delta_per_rep)

    n = len(pct_input_per_rep)
    pct_input_mean = np.mean(pct_input_per_rep)
    pct_input_se = np.std(pct_input_per_rep, ddof=1) / np.sqrt(n) if n > 1 else 0.0
    delta_mean = np.mean(delta_per_rep)
    delta_sem = np.std(delta_per_rep, ddof=1) / np.sqrt(n) if n > 1 else 0.0

    return {
        "pct_input": pct_input_mean,
        "pct_input_SE": pct_input_se,
        "delta": delta_mean,
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

    All four arrays must be paired by biological replicate (equal length),
    already averaged over technical replicates.

    Per bio rep:
        dCT_IP  = CT_IP  - (CT_input - log2(dilution_factor))
        dCT_IgG = CT_IgG - (CT_input - log2(dilution_factor))
        ddCT    = dCT_IP - dCT_IgG
        Fold enrichment = 2^(-ddCT)

    Summary values are the mean and SE of the per-bio-rep fold enrichment values.

    Parameters
    ----------
    ct_ip : array-like
        CT values for the antibody-of-interest IP, one per biological replicate.
    ct_input_for_ip : array-like
        CT values for the Input paired with ct_ip by biological replicate.
    ct_igg : array-like
        CT values for the IgG (negative control) IP, one per biological replicate.
    ct_input_for_igg : array-like
        CT values for the Input paired with ct_igg by biological replicate.
    dilution_factor : float
        Dilution factor of the input fraction.

    Returns
    -------
    dict with keys:
        ddCT, ddCT_SEM, fold_enrichment, fold_enrichment_SE
    """
    ct_ip = np.asarray(ct_ip, dtype=float)
    ct_input_for_ip = np.asarray(ct_input_for_ip, dtype=float)
    ct_igg = np.asarray(ct_igg, dtype=float)
    ct_input_for_igg = np.asarray(ct_input_for_igg, dtype=float)

    adjustment = np.log2(dilution_factor)

    # Per-bio-rep dCT for IP and IgG
    dct_ip_per_rep = ct_ip - (ct_input_for_ip - adjustment)
    dct_igg_per_rep = ct_igg - (ct_input_for_igg - adjustment)

    # Per-bio-rep ddCT and fold enrichment
    ddct_per_rep = dct_ip_per_rep - dct_igg_per_rep
    fe_per_rep = 2.0 ** (-ddct_per_rep)

    n = len(fe_per_rep)
    fold_enrichment = np.mean(fe_per_rep)
    fold_enrichment_se = np.std(fe_per_rep, ddof=1) / np.sqrt(n) if n > 1 else 0.0
    ddct = np.mean(ddct_per_rep)
    ddct_sem = np.std(ddct_per_rep, ddof=1) / np.sqrt(n) if n > 1 else 0.0

    return {
        "ddCT": ddct,
        "ddCT_SEM": ddct_sem,
        "fold_enrichment": fold_enrichment,
        "fold_enrichment_SE": fold_enrichment_se,
    }
