import numpy as np


def calculate_ddct(results: dict, gene: str, group: str, condition: str,
                   control_gene: str, control_condition: str):
    """Compute ddCT and its propagated error for a gene/group/condition combination.

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
    ddct : float or ndarray
        The delta-delta-Ct value(s).
    ddct_error : float
        Propagated SEM of the control dCT.
    """
    if gene == control_gene and condition == control_condition:
        return 0, 0

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

    gene_control_condition_ct_sem = (
        np.std(gene_control_condition_ct) / np.sqrt(len(gene_control_condition_ct))
    )
    control_gene_control_condition_ct_sem = (
        np.std(control_gene_control_condition_ct)
        / np.sqrt(len(control_gene_control_condition_ct))
    )

    dct_control_average = np.mean(
        gene_control_condition_ct - control_gene_control_condition_ct
    )
    dct_control_average_sem = np.sqrt(
        gene_control_condition_ct_sem ** 2
        + control_gene_control_condition_ct_sem ** 2
    )

    ddct = (gene_group_ct - control_gene_group_ct) - dct_control_average
    ddct_error = dct_control_average_sem

    return ddct, ddct_error
