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
)


@module.ui
def enrichment_ui():
    return ui.TagList(
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

        dilution = get_dilution()
        targets = df["Target Name"].unique().tolist()
        ip_groups = [
            g for g in df["Sample Name"].unique() if g != input_group
        ]

        # If IgG normalization is on, also exclude IgG from IP groups
        normalize_igg = normalize_igg_reactive()
        igg_group = igg_group_reactive()
        if normalize_igg and igg_group:
            ip_groups = [g for g in ip_groups if g != igg_group]

        rows = []
        for target in targets:
            ct_input = df[
                (df["Target Name"] == target)
                & (df["Sample Name"] == input_group)
            ]["CT"].values

            if len(ct_input) == 0:
                continue

            for ip in ip_groups:
                ct_ip = df[
                    (df["Target Name"] == target)
                    & (df["Sample Name"] == ip)
                ]["CT"].values

                if len(ct_ip) == 0:
                    continue

                result = calculate_percent_input(ct_ip, ct_input, dilution)

                # Stat test: Welch's t-test on CT values (IP vs Input)
                _, p_val = welch_ttest(ct_ip, ct_input)

                rows.append(
                    {
                        "Target Name": target,
                        "IP Group": ip,
                        "Input Group": input_group,
                        "id": f"{target}_{ip}",
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

        dilution = get_dilution()
        targets = df["Target Name"].unique().tolist()
        ip_groups = [
            g
            for g in df["Sample Name"].unique()
            if g not in (input_group, igg_group)
        ]

        rows = []
        for target in targets:
            ct_input_ip = df[
                (df["Target Name"] == target)
                & (df["Sample Name"] == input_group)
            ]["CT"].values

            ct_igg = df[
                (df["Target Name"] == target)
                & (df["Sample Name"] == igg_group)
            ]["CT"].values

            # Input for IgG is the same Input sample
            ct_input_igg = ct_input_ip

            if len(ct_input_ip) == 0 or len(ct_igg) == 0:
                continue

            for ip in ip_groups:
                ct_ip = df[
                    (df["Target Name"] == target)
                    & (df["Sample Name"] == ip)
                ]["CT"].values

                if len(ct_ip) == 0:
                    continue

                result = calculate_fold_enrichment(
                    ct_ip, ct_input_ip, ct_igg, ct_input_igg, dilution
                )

                # Stat test: Welch's t-test comparing IP to IgG CTs
                _, p_val = welch_ttest(ct_ip, ct_igg)

                rows.append(
                    {
                        "Target Name": target,
                        "IP Group": ip,
                        "IgG Group": igg_group,
                        "Input Group": input_group,
                        "id": f"{target}_{ip}",
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
    def active_results():
        """Return the appropriate results based on IgG toggle."""
        if normalize_igg_reactive():
            return fold_enrichment_results()
        return pct_input_results()

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

    @render.data_frame
    def results_table():
        return active_results()

    @render_widget
    def enrichment_plot():
        data = plot_data()
        if data is None or data.empty:
            return apply_classic_theme(px.scatter(title="No data to display"))

        normalize_igg = normalize_igg_reactive()

        if normalize_igg:
            y_col = "fold_enrichment"
            y_err = "fold_enrichment_SE"
            y_label = "Fold Enrichment over IgG (mean +/- SE)"
            title = "Fold Enrichment over IgG"
            ref_line = 1.0
        else:
            y_col = "pct_input"
            y_err = "pct_input_SE"
            y_label = "% Input (mean +/- SE)"
            title = "Percent Input"
            ref_line = None

        fig = px.scatter(
            data,
            x="id",
            y=y_col,
            color="Target Name",
            error_y=y_err,
            labels={y_col: y_label, "id": ""},
            title=title,
        )
        fig.update_traces(marker=dict(size=10))
        if ref_line is not None:
            fig.add_hline(y=ref_line, line_dash="dash", line_color="gray")

        fig.update_layout(showlegend=False)
        apply_classic_theme(fig)

        # Add significance annotations
        for _, row in data.iterrows():
            star = row.get("significance", "")
            if star and star not in ("", "ns"):
                fig.add_annotation(
                    x=row["id"],
                    y=row[y_col] + row[y_err] * 1.3,
                    text=star,
                    showarrow=False,
                    font=dict(size=14),
                )
        return fig

    return {"results": active_results}
