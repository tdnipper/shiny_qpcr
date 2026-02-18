import pandas as pd
import numpy as np
import io
from scipy import stats


# ---------------------------------------------------------------------------
# File import
# ---------------------------------------------------------------------------

def import_file(file_path: str) -> pd.DataFrame:
    """Read CSV/XLSX, remove Undetermined/NTC, validate required columns."""
    if file_path.endswith(".csv"):
        df = pd.read_csv(file_path, na_values=["Undetermined", "NTC"])
    elif file_path.endswith((".xls", ".xlsx")):
        df = pd.read_excel(file_path, na_values=["Undetermined", "NTC"])
    else:
        raise ValueError(f"Unsupported file format: {file_path}")

    required = {"Sample Name", "Target Name", "CT"}
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
    group_cols = ["Sample Name", "Target Name", "Biological Replicate"]
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
    """Split a DataFrame into a dict keyed by 'TargetName_SampleName'."""
    results = {}
    for (target, sample), group_df in df.groupby(["Target Name", "Sample Name"]):
        results[f"{target}_{sample}"] = group_df.copy()
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
