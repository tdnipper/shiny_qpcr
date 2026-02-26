import pandas as pd
import numpy as np
from shiny import module, ui, render, reactive
from shinywidgets import output_widget, render_widget
import plotly.express as px

from calculations.enrichment import (
    calculate_percent_input,
    calculate_fold_enrichment,
)
from shared import (
    welch_ttest,
    correct_pvalues,
    p_to_star,
    export_to_excel,
    export_to_csv,
    apply_classic_theme,
    make_plot_widget,
    _run_posthoc,
    _lookup_posthoc_pval,
    two_way_anova,
)


@module.ui
def enrichment_ui():
    return ui.TagList(
        ui.panel_conditional(
            "input.stat_test === 'anova2'",
            ui.card(
                ui.card_header("Two-Way ANOVA Summary (Group × Condition)"),
                ui.output_data_frame("anova2_table"),
            ),
        ),
        ui.card(
            ui.card_header("Enrichment Results"),
            ui.output_data_frame("results_table"),
        ),
        ui.layout_columns(
            ui.download_button("download_xlsx", "Download Excel (.xlsx)"),
            ui.download_button("download_csv", "Download CSV (.csv)"),
            col_widths=[3, 3],
        ),
        ui.card(
            ui.card_header("Enrichment Plot"),
            output_widget("enrichment_plot"),
        ),
    )


@module.server
def enrichment_server(
    input,
    output,
    session,
    df_reactive,
    input_group_reactive,
    dilution_factor_reactive,
    normalize_igg_reactive,
    igg_group_reactive,
    plot_targets_reactive,
    correction_method_reactive,
    stat_test_reactive,
    posthoc_reactive,
):
    """Server logic for RIP-qPCR / ChIP-qPCR enrichment analysis.

    Parameters
    ----------
    df_reactive : reactive.Calc
        Returns the uploaded DataFrame.
    input_group_reactive : reactive.Calc
        Returns the name of the Input / reference group.
    dilution_factor_reactive : reactive.Calc
        Returns the dilution factor (numeric).
    normalize_igg_reactive : reactive.Calc
        Returns True if IgG normalization is enabled.
    igg_group_reactive : reactive.Calc
        Returns the name of the IgG / negative-control group.
    plot_targets_reactive : reactive.Calc
        Returns the list of target names selected for plotting.
    correction_method_reactive : reactive.Calc
        Returns the p-value correction method string.
    stat_test_reactive : reactive.Calc
        Returns the selected statistical test string.
    posthoc_reactive : reactive.Calc
        Returns the selected post-hoc test string.
    """

    @reactive.calc
    def get_dilution(self=None):
        """Get the dilution factor from the data column or sidebar input."""
        df = df_reactive()
        sidebar_val = dilution_factor_reactive()
        if df is not None and "Dilution Factor" in df.columns:
            # Use per-row dilution factor from data (take first non-null)
            vals = df["Dilution Factor"].dropna().unique()
            if len(vals) > 0:
                return float(vals[0])
        return float(sidebar_val) if sidebar_val else 10.0

    @reactive.calc
    def pct_input_results():
        """Calculate % Input for each IP group and target."""
        df = df_reactive()
        input_group = input_group_reactive()
        if df is None or not input_group:
            return None

        stat_test = stat_test_reactive()
        posthoc_test = posthoc_reactive()
        dilution = get_dilution()
        targets = df["Target Name"].unique().tolist()
        groups = df["Group"].unique().tolist()
        ip_conditions = [
            c for c in df["Condition"].unique() if c != input_group
        ]

        # If IgG normalization is on, also exclude IgG from IP conditions
        normalize_igg = normalize_igg_reactive()
        igg_group = igg_group_reactive()
        if normalize_igg and igg_group:
            ip_conditions = [c for c in ip_conditions if c != igg_group]

        rows = []
        for target in targets:
            for grp in groups:
                input_rows = df[
                    (df["Target Name"] == target)
                    & (df["Group"] == grp)
                    & (df["Condition"] == input_group)
                ]
                ct_input = input_rows["CT"].values

                if len(ct_input) == 0:
                    continue

                # For ANOVA/KW: compute post-hoc once per target+group
                if stat_test not in ("ttest", "none"):
                    groups_data = {input_group: ct_input}
                    for ip_cond in ip_conditions:
                        ct_ip = df[
                            (df["Target Name"] == target)
                            & (df["Group"] == grp)
                            & (df["Condition"] == ip_cond)
                        ]["CT"].values
                        if len(ct_ip) > 0:
                            groups_data[ip_cond] = ct_ip
                    ph_df = _run_posthoc(groups_data, stat_test, posthoc_test)
                else:
                    ph_df = None

                for ip_cond in ip_conditions:
                    ip_rows = df[
                        (df["Target Name"] == target)
                        & (df["Group"] == grp)
                        & (df["Condition"] == ip_cond)
                    ]

                    if ip_rows.empty:
                        continue

                    # Pair input and IP by biological replicate when possible
                    if "Biological Replicate" in df.columns:
                        merged = (
                            input_rows[["Biological Replicate", "CT"]]
                            .rename(columns={"CT": "ct_input"})
                            .merge(
                                ip_rows[["Biological Replicate", "CT"]].rename(
                                    columns={"CT": "ct_ip"}
                                ),
                                on="Biological Replicate",
                            )
                        )
                        if merged.empty:
                            continue
                        ct_input_paired = merged["ct_input"].values
                        ct_ip = merged["ct_ip"].values
                    else:
                        ct_input_paired = ct_input
                        ct_ip = ip_rows["CT"].values

                    result = calculate_percent_input(ct_ip, ct_input_paired, dilution)

                    if stat_test == "none":
                        p_val = np.nan
                    elif stat_test == "ttest":
                        _, p_val = welch_ttest(ct_ip, ct_input)
                    else:
                        p_val = _lookup_posthoc_pval(ph_df, input_group, ip_cond)

                    rows.append(
                        {
                            "Target Name": target,
                            "Group": grp,
                            "IP Condition": ip_cond,
                            "Input Condition": input_group,
                            "id": f"{target}_{grp}_{ip_cond}",
                            "pct_input": result["pct_input"],
                            "pct_input_SE": result["pct_input_SE"],
                            "delta": result["delta"],
                            "delta_SEM": result["delta_SEM"],
                            "dilution_factor": dilution,
                            "p_value": p_val,
                        }
                    )

        if not rows:
            return None

        result_df = pd.DataFrame(rows)
        method = correction_method_reactive()
        result_df["p_value"] = correct_pvalues(
            result_df["p_value"].values, method=method
        )
        result_df["significance"] = result_df["p_value"].apply(p_to_star)
        return result_df

    @reactive.calc
    def fold_enrichment_results():
        """Calculate fold enrichment over IgG for each IP group and target."""
        df = df_reactive()
        input_group = input_group_reactive()
        igg_group = igg_group_reactive()
        if df is None or not input_group or not igg_group:
            return None

        stat_test = stat_test_reactive()
        posthoc_test = posthoc_reactive()
        dilution = get_dilution()
        targets = df["Target Name"].unique().tolist()
        groups = df["Group"].unique().tolist()
        ip_conditions = [
            c
            for c in df["Condition"].unique()
            if c not in (input_group, igg_group)
        ]

        rows = []
        for target in targets:
            for grp in groups:
                input_rows_fe = df[
                    (df["Target Name"] == target)
                    & (df["Group"] == grp)
                    & (df["Condition"] == input_group)
                ]
                ct_input = input_rows_fe["CT"].values

                igg_rows = df[
                    (df["Target Name"] == target)
                    & (df["Group"] == grp)
                    & (df["Condition"] == igg_group)
                ]
                ct_igg = igg_rows["CT"].values

                if len(ct_input) == 0 or len(ct_igg) == 0:
                    continue

                # For ANOVA/KW: compute post-hoc once per target+group
                # Reference is IgG; compare each IP condition to IgG
                if stat_test not in ("ttest", "none"):
                    groups_data = {igg_group: ct_igg}
                    for ip_cond in ip_conditions:
                        ct_ip = df[
                            (df["Target Name"] == target)
                            & (df["Group"] == grp)
                            & (df["Condition"] == ip_cond)
                        ]["CT"].values
                        if len(ct_ip) > 0:
                            groups_data[ip_cond] = ct_ip
                    ph_df = _run_posthoc(groups_data, stat_test, posthoc_test)
                else:
                    ph_df = None

                for ip_cond in ip_conditions:
                    ip_rows_fe = df[
                        (df["Target Name"] == target)
                        & (df["Group"] == grp)
                        & (df["Condition"] == ip_cond)
                    ]

                    if ip_rows_fe.empty:
                        continue

                    # Three-way pair: input, IP, IgG by biological replicate
                    if "Biological Replicate" in df.columns:
                        merged_fe = (
                            input_rows_fe[["Biological Replicate", "CT"]]
                            .rename(columns={"CT": "ct_input"})
                            .merge(
                                ip_rows_fe[["Biological Replicate", "CT"]].rename(
                                    columns={"CT": "ct_ip"}
                                ),
                                on="Biological Replicate",
                            )
                            .merge(
                                igg_rows[["Biological Replicate", "CT"]].rename(
                                    columns={"CT": "ct_igg"}
                                ),
                                on="Biological Replicate",
                            )
                        )
                        if merged_fe.empty:
                            continue
                        ct_ip_paired = merged_fe["ct_ip"].values
                        ct_input_paired_fe = merged_fe["ct_input"].values
                        ct_igg_paired = merged_fe["ct_igg"].values
                    else:
                        ct_ip_paired = ip_rows_fe["CT"].values
                        ct_input_paired_fe = ct_input
                        ct_igg_paired = ct_igg

                    result = calculate_fold_enrichment(
                        ct_ip_paired, ct_input_paired_fe,
                        ct_igg_paired, ct_input_paired_fe,
                        dilution,
                    )

                    if stat_test == "none":
                        p_val = np.nan
                    elif stat_test == "ttest":
                        _, p_val = welch_ttest(ct_ip_paired, ct_igg_paired)
                    else:
                        p_val = _lookup_posthoc_pval(ph_df, igg_group, ip_cond)

                    rows.append(
                        {
                            "Target Name": target,
                            "Group": grp,
                            "IP Condition": ip_cond,
                            "IgG Condition": igg_group,
                            "Input Condition": input_group,
                            "id": f"{target}_{grp}_{ip_cond}",
                            "ddCT": result["ddCT"],
                            "ddCT_SEM": result["ddCT_SEM"],
                            "fold_enrichment": result["fold_enrichment"],
                            "fold_enrichment_SE": result["fold_enrichment_SE"],
                            "dilution_factor": dilution,
                            "p_value": p_val,
                        }
                    )

        if not rows:
            return None

        result_df = pd.DataFrame(rows)
        method = correction_method_reactive()
        result_df["p_value"] = correct_pvalues(
            result_df["p_value"].values, method=method
        )
        result_df["significance"] = result_df["p_value"].apply(p_to_star)
        return result_df

    @reactive.calc
    def anova2_results():
        if stat_test_reactive() != "anova2":
            return None
        df = df_reactive()
        if df is None:
            return None

        targets = df["Target Name"].unique().tolist()
        all_rows = []
        for target in targets:
            target_df = df[df["Target Name"] == target].copy()
            anova_df = two_way_anova(target_df, "CT", "Group", "Condition")
            if not anova_df.empty:
                anova_df.insert(0, "Target Name", target)
                all_rows.append(anova_df)

        if not all_rows:
            return None
        return pd.concat(all_rows, ignore_index=True)

    @reactive.calc
    def active_results():
        """Return the appropriate results based on IgG toggle."""
        if normalize_igg_reactive():
            return fold_enrichment_results()
        return pct_input_results()

    @reactive.calc
    def enrichment_individual_reps():
        df = df_reactive()
        summary = active_results()
        input_group = input_group_reactive()
        if df is None or summary is None or not input_group:
            return None

        normalize_igg = normalize_igg_reactive()
        igg_group = igg_group_reactive()
        dilution = get_dilution()
        adjustment = np.log2(dilution)

        targets = df["Target Name"].unique().tolist()
        groups = df["Group"].unique().tolist()
        ip_conditions = [c for c in df["Condition"].unique() if c != input_group]
        if normalize_igg and igg_group:
            ip_conditions = [c for c in ip_conditions if c != igg_group]

        rows = []
        for target in targets:
            for grp in groups:
                ct_input_rows = df[
                    (df["Target Name"] == target)
                    & (df["Group"] == grp)
                    & (df["Condition"] == input_group)
                ]
                if ct_input_rows.empty:
                    continue

                if normalize_igg and igg_group:
                    ct_igg_rows = df[
                        (df["Target Name"] == target)
                        & (df["Group"] == grp)
                        & (df["Condition"] == igg_group)
                    ]
                    if ct_igg_rows.empty:
                        continue

                for ip_cond in ip_conditions:
                    ip_rows = df[
                        (df["Target Name"] == target)
                        & (df["Group"] == grp)
                        & (df["Condition"] == ip_cond)
                    ]
                    if ip_rows.empty:
                        continue
                    row_id = f"{target}_{grp}_{ip_cond}"
                    for _, r in ip_rows.iterrows():
                        # Look up this bio rep's own input CT
                        bio_rep = r.get("Biological Replicate", None)
                        if (
                            bio_rep is not None
                            and "Biological Replicate" in ct_input_rows.columns
                        ):
                            rep_input = ct_input_rows[
                                ct_input_rows["Biological Replicate"] == bio_rep
                            ]["CT"]
                            input_ct = (
                                rep_input.iloc[0]
                                if len(rep_input) > 0
                                else ct_input_rows["CT"].mean()
                            )
                        else:
                            input_ct = ct_input_rows["CT"].mean()

                        if normalize_igg and igg_group:
                            if (
                                bio_rep is not None
                                and "Biological Replicate" in ct_igg_rows.columns
                            ):
                                rep_igg = ct_igg_rows[
                                    ct_igg_rows["Biological Replicate"] == bio_rep
                                ]["CT"]
                                igg_ct = (
                                    rep_igg.iloc[0]
                                    if len(rep_igg) > 0
                                    else ct_igg_rows["CT"].mean()
                                )
                            else:
                                igg_ct = ct_igg_rows["CT"].mean()
                            dct_ip = r["CT"] - (input_ct - adjustment)
                            dct_igg_rep = igg_ct - (input_ct - adjustment)
                            ddct = dct_ip - dct_igg_rep
                            val = 2.0 ** -ddct
                            val_col = "fold_enrichment_rep"
                        else:
                            val = 100.0 * (2.0 ** ((input_ct - adjustment) - r["CT"]))
                            val_col = "pct_input_rep"
                        rows.append({
                            "id": row_id,
                            "Target Name": target,
                            "Group": grp,
                            "IP Condition": ip_cond,
                            "Biological Replicate": r.get("Biological Replicate", ""),
                            val_col: val,
                        })

        if not rows:
            return None

        reps_df = pd.DataFrame(rows)
        sig = summary[["id", "significance"]].drop_duplicates()
        return reps_df.merge(sig, on="id", how="left")

    # ---- Downloads ----

    @render.download(filename="enrichment_results.xlsx")
    def download_xlsx():
        data = active_results()
        if data is not None:
            yield export_to_excel(data)

    @render.download(filename="enrichment_results.csv")
    def download_csv():
        data = active_results()
        if data is not None:
            yield export_to_csv(data)

    # ---- Rendering ----

    @reactive.calc
    def plot_data():
        data = active_results()
        targets = plot_targets_reactive()
        if data is None or not targets:
            return None
        mask = data["Target Name"].isin(list(targets))
        return data[mask]

    @reactive.calc
    def plot_data_reps():
        data = enrichment_individual_reps()
        targets = plot_targets_reactive()
        if data is None or not targets:
            return None
        return data[data["Target Name"].isin(list(targets))]

    @render.data_frame
    def results_table():
        return active_results()

    @render.data_frame
    def anova2_table():
        return anova2_results()

    @render_widget
    def enrichment_plot():
        reps = plot_data_reps()
        summary = plot_data()
        if reps is None or reps.empty:
            return make_plot_widget(apply_classic_theme(px.scatter(title="No data to display")))

        normalize_igg = normalize_igg_reactive()

        if normalize_igg:
            value_col = "fold_enrichment_rep"
            y_label = "Fold Enrichment over IgG"
            title = "Fold Enrichment over IgG"
            ref_line = 1.0
        else:
            value_col = "pct_input_rep"
            y_label = "% Input"
            title = "Percent Input"
            ref_line = None

        fig = px.strip(
            reps,
            x="id",
            y=value_col,
            color="Target Name",
            stripmode="overlay",
            labels={value_col: y_label, "id": ""},
            title=title,
        )
        fig.update_traces(jitter=0.3, marker=dict(size=8, opacity=0.8))
        if ref_line is not None:
            fig.add_hline(y=ref_line, line_dash="dash", line_color="gray")
        fig.update_layout(showlegend=False)
        apply_classic_theme(fig)

        if summary is not None and not summary.empty:
            summary_y_col = "fold_enrichment" if normalize_igg else "pct_input"
            max_per_id = reps.groupby("id")[value_col].max()
            for _, row in summary.iterrows():
                star = row.get("significance", "")
                if star and star not in ("", "ns"):
                    y_pos = max_per_id.get(row["id"], row[summary_y_col]) * 1.15
                    fig.add_annotation(
                        x=row["id"],
                        y=y_pos,
                        text=star,
                        showarrow=False,
                        font=dict(size=14),
                    )
        return make_plot_widget(fig)

    return {"results": active_results}
