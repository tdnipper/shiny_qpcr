import pandas as pd
import numpy as np
import io
from scipy import stats


# ---------------------------------------------------------------------------
# File import
# ---------------------------------------------------------------------------

def _is_quantstudio_export(file_path: str) -> bool:
    if not file_path.endswith((".xls", ".xlsx")):
        return False
    import openpyxl
    wb = openpyxl.load_workbook(file_path, read_only=True, data_only=True)
    result = "Results" in wb.sheetnames
    wb.close()
    return result


def _split_sample_name(sample_name: str) -> tuple[str, str]:
    if "_" not in sample_name:
        raise ValueError(
            f"Sample Name '{sample_name}' has no underscore; "
            "cannot split into Group and Condition."
        )
    group, _, condition = sample_name.rpartition("_")
    return group, condition


def _transform_quantstudio_df(df: pd.DataFrame) -> pd.DataFrame:
    split = df["Sample Name"].apply(_split_sample_name)
    df = df.copy()
    df["Group"] = [pair[0] for pair in split]
    df["Condition"] = [pair[1] for pair in split]
    df = df.rename(columns={"Task": "Biological Replicate"})
    df = df.drop(columns=["Sample Name"])
    cols = ["Group", "Condition", "Target Name", "CT", "Biological Replicate"]
    extra = [c for c in df.columns if c not in cols]
    return df[cols + extra]


def import_file(file_path: str) -> pd.DataFrame:
    """Read CSV/XLSX, remove Undetermined/NTC, validate required columns."""
    if file_path.endswith(".csv"):
        df = pd.read_csv(file_path, na_values=["Undetermined", "NTC"])
    elif file_path.endswith((".xls", ".xlsx")):
        if _is_quantstudio_export(file_path):
            from qpcr_importer import FileImporter
            raw = FileImporter(file_path).import_file()
            df = _transform_quantstudio_df(raw)
        else:
            df = pd.read_excel(file_path, na_values=["Undetermined", "NTC"])
    else:
        raise ValueError(f"Unsupported file format: {file_path}")

    required = {"Group", "Condition", "Target Name", "CT"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df["CT"] = pd.to_numeric(df["CT"], errors="coerce")
    df = df.dropna(subset=["CT"]).reset_index(drop=True)
    return df


def average_technical_replicates(df: pd.DataFrame) -> pd.DataFrame:
    """Average technical replicates within each biological replicate.

    Groups by Sample Name, Target Name, Biological Replicate and returns
    one row per biological replicate with mean CT. Other columns (e.g.
    Dilution Factor) are preserved by taking the first value per group.
    """
    group_cols = ["Group", "Condition", "Target Name", "Biological Replicate"]
    other_cols = [c for c in df.columns if c not in group_cols + ["CT"]]
    agg = {"CT": "mean"}
    for col in other_cols:
        agg[col] = "first"
    return df.groupby(group_cols).agg(agg).reset_index()


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

def filter_targets(df: pd.DataFrame, target: str) -> pd.DataFrame:
    """Filter DataFrame rows where Target Name contains *target* (literal)."""
    return df[df["Target Name"].str.contains(target, regex=False)].copy()


def group_ct_data(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Split a DataFrame into a dict keyed by 'TargetName_GroupName_Condition'."""
    results = {}
    for (target, group, condition), group_df in df.groupby(["Target Name", "Group", "Condition"]):
        results[f"{target}_{group}_{condition}"] = group_df.copy()
    return results


# ---------------------------------------------------------------------------
# Statistical tests
# ---------------------------------------------------------------------------

def welch_ttest(group1: np.ndarray, group2: np.ndarray):
    """Welch's two-sample t-test (unequal variance). Returns (t_stat, p_value)."""
    if len(group1) < 2 or len(group2) < 2:
        return np.nan, np.nan
    t_stat, p_val = stats.ttest_ind(group1, group2, equal_var=False)
    return t_stat, p_val


def one_sample_ttest(values: np.ndarray, popmean: float = 0.0):
    """One-sample t-test. Returns (t_stat, p_value)."""
    if len(values) < 2:
        return np.nan, np.nan
    t_stat, p_val = stats.ttest_1samp(values, popmean)
    return t_stat, p_val


def correct_pvalues(p_values: np.ndarray, method: str = "none") -> np.ndarray:
    """Apply multiple testing correction.

    Methods: 'none', 'bonferroni', 'bh' (Benjamini-Hochberg).
    """
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0 or method == "none":
        return p

    if method == "bonferroni":
        return np.minimum(p * n, 1.0)

    if method == "bh":
        order = np.argsort(p)
        sorted_p = p[order]
        adjusted = np.empty(n)
        cummin = 1.0
        for i in range(n - 1, -1, -1):
            val = min(sorted_p[i] * n / (i + 1), 1.0)
            cummin = min(cummin, val)
            adjusted[order[i]] = cummin
        return adjusted

    return p


def p_to_star(p: float) -> str:
    """Convert p-value to significance star annotation."""
    if np.isnan(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# ---------------------------------------------------------------------------
# Advanced statistical tests (ANOVA, Kruskal-Wallis, post-hoc)
# ---------------------------------------------------------------------------

def one_way_anova(*groups) -> tuple[float, float]:
    """One-way ANOVA. Returns (F_stat, p_value) or (nan, nan) if < 2 valid groups."""
    valid = [g for g in groups if len(g) >= 2]
    if len(valid) < 2:
        return np.nan, np.nan
    try:
        F, p = stats.f_oneway(*valid)
        return float(F), float(p)
    except Exception:
        return np.nan, np.nan


def kruskal_wallis(*groups) -> tuple[float, float]:
    """Kruskal-Wallis test. Returns (H_stat, p_value) or (nan, nan)."""
    valid = [g for g in groups if len(g) >= 1]
    if len(valid) < 2:
        return np.nan, np.nan
    try:
        H, p = stats.kruskal(*valid)
        return float(H), float(p)
    except Exception:
        return np.nan, np.nan


def tukey_hsd_posthoc(groups_dict: dict) -> pd.DataFrame:
    """Tukey HSD post-hoc test.

    Parameters
    ----------
    groups_dict : dict[str, np.ndarray]
        Mapping of group label to array of values.

    Returns
    -------
    pd.DataFrame with columns: group1, group2, p_adj, significance
    """
    labels = list(groups_dict.keys())
    arrays = [np.asarray(groups_dict[k]) for k in labels]
    valid_mask = [len(a) >= 2 for a in arrays]
    if sum(valid_mask) < 2:
        return pd.DataFrame(columns=["group1", "group2", "p_adj", "significance"])

    valid_labels = [l for l, v in zip(labels, valid_mask) if v]
    valid_arrays = [a for a, v in zip(arrays, valid_mask) if v]

    try:
        result = stats.tukey_hsd(*valid_arrays)
    except Exception:
        return pd.DataFrame(columns=["group1", "group2", "p_adj", "significance"])

    rows = []
    n = len(valid_labels)
    for i in range(n):
        for j in range(i + 1, n):
            p = float(result.pvalue[i, j])
            rows.append({
                "group1": valid_labels[i],
                "group2": valid_labels[j],
                "p_adj": p,
                "significance": p_to_star(p),
            })
    return pd.DataFrame(rows)


def pairwise_welch_posthoc(groups_dict: dict) -> pd.DataFrame:
    """Pairwise Welch t-tests with Bonferroni correction.

    Parameters
    ----------
    groups_dict : dict[str, np.ndarray]
        Mapping of group label to array of values.

    Returns
    -------
    pd.DataFrame with columns: group1, group2, p_adj, significance
    """
    labels = list(groups_dict.keys())
    raw_p = []
    pair_labels = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            g1 = np.asarray(groups_dict[labels[i]])
            g2 = np.asarray(groups_dict[labels[j]])
            _, p = welch_ttest(g1, g2)
            raw_p.append(p)
            pair_labels.append((labels[i], labels[j]))

    adjusted = correct_pvalues(np.array(raw_p, dtype=float), method="bonferroni")
    rows = []
    for (g1, g2), p in zip(pair_labels, adjusted):
        rows.append({
            "group1": g1,
            "group2": g2,
            "p_adj": float(p),
            "significance": p_to_star(float(p)),
        })
    return pd.DataFrame(rows)


def dunns_posthoc(groups_dict: dict, correction: str = "bonferroni") -> pd.DataFrame:
    """Dunn's post-hoc test (pairwise Mann-Whitney U with correction).

    Parameters
    ----------
    groups_dict : dict[str, np.ndarray]
        Mapping of group label to array of values.
    correction : str
        Correction method: 'bonferroni', 'bh', 'dunn_bonf', or 'dunn_bh'.

    Returns
    -------
    pd.DataFrame with columns: group1, group2, p_adj, significance
    """
    corr_method = "bh" if "bh" in correction else "bonferroni"
    labels = list(groups_dict.keys())
    raw_p = []
    pair_labels = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            g1 = np.asarray(groups_dict[labels[i]])
            g2 = np.asarray(groups_dict[labels[j]])
            if len(g1) >= 1 and len(g2) >= 1:
                try:
                    _, p = stats.mannwhitneyu(g1, g2, alternative="two-sided")
                    raw_p.append(float(p))
                except Exception:
                    raw_p.append(np.nan)
            else:
                raw_p.append(np.nan)
            pair_labels.append((labels[i], labels[j]))

    adjusted = correct_pvalues(np.array(raw_p, dtype=float), method=corr_method)
    rows = []
    for (g1, g2), p in zip(pair_labels, adjusted):
        rows.append({
            "group1": g1,
            "group2": g2,
            "p_adj": float(p),
            "significance": p_to_star(float(p)),
        })
    return pd.DataFrame(rows)


def two_way_anova(df: pd.DataFrame, value_col: str, factor1: str, factor2: str) -> pd.DataFrame:
    """Two-way ANOVA using statsmodels OLS formula API.

    Parameters
    ----------
    df : pd.DataFrame
        Data containing value_col, factor1, and factor2 columns.
    value_col : str
        Name of the response variable column.
    factor1 : str
        Name of the first factor column.
    factor2 : str
        Name of the second factor column.

    Returns
    -------
    pd.DataFrame with columns: Factor, F, p_value, significance
    """
    try:
        import statsmodels.formula.api as smf
        from statsmodels.stats.anova import anova_lm

        tmp = df[[value_col, factor1, factor2]].copy()
        tmp.columns = ["value", "f1", "f2"]
        tmp = tmp.dropna()

        if len(tmp) < 4:
            return pd.DataFrame(columns=["Factor", "F", "p_value", "significance"])

        model = smf.ols("value ~ C(f1) + C(f2) + C(f1):C(f2)", data=tmp).fit()
        anova_table = anova_lm(model, typ=2)

        name_map = {
            "C(f1)": factor1,
            "C(f2)": factor2,
            "C(f1):C(f2)": f"{factor1}:{factor2}",
        }

        rows = []
        for idx, row in anova_table.iterrows():
            if idx == "Residual":
                continue
            factor_name = name_map.get(idx, idx)
            p = float(row["PR(>F)"])
            rows.append({
                "Factor": factor_name,
                "F": round(float(row["F"]), 4),
                "p_value": round(p, 4),
                "significance": p_to_star(p),
            })
        return pd.DataFrame(rows)
    except Exception:
        return pd.DataFrame(columns=["Factor", "F", "p_value", "significance"])


def _lookup_posthoc_pval(ph_df: pd.DataFrame, group_a: str, group_b: str) -> float:
    """Find p_adj from post-hoc DataFrame for a given pair (order-independent)."""
    if ph_df is None or ph_df.empty:
        return np.nan
    mask = (
        ((ph_df["group1"] == group_a) & (ph_df["group2"] == group_b))
        | ((ph_df["group1"] == group_b) & (ph_df["group2"] == group_a))
    )
    matches = ph_df[mask]
    if matches.empty:
        return np.nan
    return float(matches.iloc[0]["p_adj"])


def _run_posthoc(groups_dict: dict, stat_test: str, posthoc_test: str) -> pd.DataFrame:
    """Dispatch to the appropriate post-hoc function.

    Parameters
    ----------
    groups_dict : dict[str, np.ndarray]
    stat_test : str
        'anova1', 'anova2', or 'kruskal'
    posthoc_test : str
        'tukey', 'welch_bonf', 'dunn_bonf', or 'dunn_bh'
    """
    if stat_test == "kruskal":
        return dunns_posthoc(groups_dict, correction=posthoc_test)
    if posthoc_test == "tukey":
        return tukey_hsd_posthoc(groups_dict)
    return pairwise_welch_posthoc(groups_dict)


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def apply_classic_theme(fig):
    """Apply a classic theme (white background, axis lines, no gridlines)."""
    fig.update_layout(
        template="simple_white",
        plot_bgcolor="white",
        xaxis=dict(
            showline=True, linewidth=1, linecolor="black",
            mirror=False, showgrid=False,
        ),
        yaxis=dict(
            showline=True, linewidth=1, linecolor="black",
            mirror=False, showgrid=False,
        ),
        toImageButtonOptions=dict(
            format="png",
            scale=3,
        ),
    )
    return fig


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

def export_to_excel(df: pd.DataFrame) -> bytes:
    buf = io.BytesIO()
    df.to_excel(buf, index=False)
    return buf.getvalue()


def export_to_csv(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")
