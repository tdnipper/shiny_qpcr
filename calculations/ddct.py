import numpy as np


def calculate_ddct(results: dict, gene: str, group: str, condition: str,
                   control_gene: str, control_condition: str):
    """Compute per-biological-replicate ddCT values for a gene/group/condition combination.

    Parameters
    ----------
    results : dict
        Mapping of ``"TargetName_GroupName_Condition"`` to DataFrames that each
        contain a ``"CT"`` column.
    gene : str
        Target gene name.
    group : str
        Experiment arm / group name (stays constant; only condition varies).
    condition : str
        Sample condition (fraction/treatment) being measured.
    control_gene : str
        Housekeeping (reference) gene name.
    control_condition : str
        Control (baseline) condition name.

    Returns
    -------
    ddct : ndarray
        Per-biological-replicate delta-delta-Ct values.
    """
    if gene == control_gene and condition == control_condition:
        return 0

    try:
        gene_group_ct = results[f"{gene}_{group}_{condition}"]["CT"].values
        control_gene_group_ct = results[f"{control_gene}_{group}_{condition}"]["CT"].values
        gene_control_condition_ct = results[f"{gene}_{group}_{control_condition}"]["CT"].values
        control_gene_control_condition_ct = results[f"{control_gene}_{group}_{control_condition}"]["CT"].values
    except KeyError as e:
        raise KeyError(
            f"Missing key in results: {e}. "
            "Check if all input values exist in the data."
        )

    dct_control_average = np.mean(
        gene_control_condition_ct - control_gene_control_condition_ct
    )

    ddct = (gene_group_ct - control_gene_group_ct) - dct_control_average

    return ddct
