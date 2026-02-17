import pandas as pd
import numpy as np
from shiny import module, ui, render, reactive
from shinywidgets import output_widget, render_widget
import plotly.express as px
import plotly.graph_objects as go

from calculations.dct import calculate_dct_between_groups
from shared import (
    welch_ttest,
    correct_pvalues,
    p_to_star,
    export_to_excel,
    export_to_csv,
)


@module.ui
def dct_ui():
    return ui.TagList(
        ui.card(
            ui.card_header("dCT Results (Fold Change Between Groups)"),
            ui.output_data_frame("results_table"),
        ),
        ui.layout_columns(
            ui.download_button("download_xlsx", "Download Excel (.xlsx)"),
            ui.download_button("download_csv", "Download CSV (.csv)"),
            col_widths=[3, 3],
        ),
        ui.card(
            ui.card_header("Fold Change Bar Chart"),
            output_widget("foldchange_bar"),
        ),
        ui.card(
            ui.card_header("Percent of Control"),
            output_widget("pct_control_bar"),
        ),
    )


@module.server
def dct_server(
    input,
    output,
    session,
    df_reactive,
    control_reactive,
    plot_targets_reactive,
    correction_method_reactive,
):
    """Server logic for the dCT (fold-change between groups) analysis mode.

    Parameters
    ----------
    df_reactive : reactive.Calc
        Returns the uploaded DataFrame.
    control_reactive : reactive.Calc
        Returns the selected control group name.
    plot_targets_reactive : reactive.Calc
        Returns the list of target names selected for plotting.
    correction_method_reactive : reactive.Calc
        Returns the p-value correction method string.
    """

    @reactive.calc
    def dct_results():
        df = df_reactive()
        control = control_reactive()
        if df is None or not control:
            return None

        targets = df["Target Name"].unique().tolist()
        groups = df["Sample Name"].unique().tolist()

        rows = []
        for target in targets:
            ct_control = df[
                (df["Target Name"] == target) & (df["Sample Name"] == control)
            ]["CT"].values

            if len(ct_control) == 0:
                continue

            for group in groups:
                if group == control:
                    continue

                ct_exp = df[
                    (df["Target Name"] == target)
                    & (df["Sample Name"] == group)
                ]["CT"].values

                if len(ct_exp) == 0:
                    continue

                result = calculate_dct_between_groups(ct_exp, ct_control)

                # Statistical test: Welch's t-test on raw CT values
                _, p_val = welch_ttest(ct_exp, ct_control)

                rows.append(
                    {
                        "Target Name": target,
                        "Sample Name": group,
                        "Control": control,
                        "id": f"{target}_{group}",
                        "dCT": result["dCT"],
                        "dCT_SEM": result["dCT_SEM"],
                        "fold_change": result["fold_change"],
                        "fold_change_SE": result["fold_change_SE"],
                        "pct_control": result["pct_control"],
                        "pct_control_SE": result["pct_control_SE"],
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

    # ---- Downloads ----

    @render.download(filename="dct_foldchange.xlsx")
    def download_xlsx():
        data = dct_results()
        if data is not None:
            yield export_to_excel(data)

    @render.download(filename="dct_foldchange.csv")
    def download_csv():
        data = dct_results()
        if data is not None:
            yield export_to_csv(data)

    # ---- Plots ----

    @reactive.calc
    def plot_data():
        data = dct_results()
        targets = plot_targets_reactive()
        if data is None or not targets:
            return None
        mask = data["Target Name"].isin(list(targets))
        return data[mask]

    @render.data_frame
    def results_table():
        return dct_results()

    @render_widget
    def foldchange_bar():
        data = plot_data()
        if data is None or data.empty:
            return px.bar(title="No data to display")

        fig = px.bar(
            data,
            x="id",
            y="fold_change",
            color="Target Name",
            error_y="fold_change_SE",
            labels={
                "fold_change": "Fold Change (mean +/- SE)",
                "id": "",
            },
            title="Fold Change vs Control (2^-dCT)",
            template="simple_white",
        )
        fig.add_hline(y=1.0, line_dash="dash", line_color="gray")
        fig.update_layout(showlegend=False)

        # Add significance annotations
        for _, row in data.iterrows():
            star = row.get("significance", "")
            if star and star not in ("", "ns"):
                fig.add_annotation(
                    x=row["id"],
                    y=row["fold_change"] + row["fold_change_SE"] * 1.3,
                    text=star,
                    showarrow=False,
                    font=dict(size=14),
                )
        return fig

    @render_widget
    def pct_control_bar():
        data = plot_data()
        if data is None or data.empty:
            return px.bar(title="No data to display")

        fig = px.bar(
            data,
            x="id",
            y="pct_control",
            color="Target Name",
            error_y="pct_control_SE",
            labels={
                "pct_control": "% of Control (mean +/- SE)",
                "id": "",
            },
            title="Expression as Percent of Control",
            template="simple_white",
        )
        fig.add_hline(y=100.0, line_dash="dash", line_color="gray")
        fig.update_layout(showlegend=False)
        return fig

    return {"results": dct_results}
