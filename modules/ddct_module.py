import pandas as pd
import numpy as np
from shiny import module, ui, render, reactive
from shinywidgets import output_widget, render_widget
import plotly.express as px

from calculations.ddct import calculate_ddct
from shared import (
    filter_targets,
    group_ct_data,
    one_sample_ttest,
    correct_pvalues,
    p_to_star,
    export_to_excel,
    export_to_csv,
    apply_classic_theme,
    _run_posthoc,
    _lookup_posthoc_pval,
    two_way_anova,
)


@module.ui
def ddct_ui():
    return ui.TagList(
        ui.panel_conditional(
            "input.stat_test === 'anova2'",
            ui.card(
                ui.card_header("Two-Way ANOVA Summary (Group × Condition)"),
                ui.output_data_frame("anova2_table"),
            ),
        ),
        ui.card(
            ui.card_header("ddCT Results (Fold Change Summary)"),
            ui.output_data_frame("results_table"),
        ),
        ui.layout_columns(
            ui.download_button("download_xlsx", "Download Excel (.xlsx)"),
            ui.download_button("download_csv", "Download CSV (.csv)"),
            col_widths=[3, 3],
        ),
        ui.card(
            ui.card_header("Fold Change Plot"),
            output_widget("foldchange_plot"),
        ),
        ui.card(
            ui.card_header("ddCT Plot"),
            output_widget("ddct_plot"),
        ),
    )


@module.server
def ddct_server(
    input,
    output,
    session,
    df_reactive,
    housekeeping_reactive,
    control_reactive,
    plot_targets_reactive,
    correction_method_reactive,
    stat_test_reactive,
    posthoc_reactive,
):
    """Server logic for the ddCT (relative expression) analysis mode.

    Parameters
    ----------
    df_reactive : reactive.Calc
        Returns the uploaded DataFrame.
    housekeeping_reactive : reactive.Calc
        Returns the selected housekeeping gene name.
    control_reactive : reactive.Calc
        Returns the selected control condition name.
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
    def processed():
        dataframe = df_reactive()
        if dataframe is None:
            return None
        df = dataframe.copy()
        df["mean"] = (
            df.groupby(["Group", "Condition", "Target Name"])["CT"]
            .transform("mean")
            .reset_index(drop=True)
        )
        df["id"] = df["Group"] + "_" + df["Condition"]
        return df

    @reactive.calc
    def results_dict():
        df = processed()
        if df is None:
            return None
        return group_ct_data(df)

    @reactive.calc
    def ddct_results():
        results = results_dict()
        df = processed()
        if results is None or df is None:
            return None

        housekeeping = housekeeping_reactive()
        control = control_reactive()
        if not housekeeping or not control:
            return None

        genes = df["Target Name"].unique().tolist()
        groups = df["Group"].unique().tolist()
        conditions = df["Condition"].unique().tolist()

        for gene in genes:
            for grp in groups:
                for cond in conditions:
                    if cond == control:
                        continue
                    key = f"{gene}_{grp}_{cond}"
                    if key in results:
                        results[key]["ddct"], results[key]["ddct_error"] = (
                            calculate_ddct(
                                results=results,
                                gene=gene,
                                group=grp,
                                condition=cond,
                                control_gene=housekeeping,
                                control_condition=control,
                            )
                        )

        ddct_df = pd.concat(results.values()).reset_index(drop=True)
        ddct_df = ddct_df.dropna(subset=["ddct"]).reset_index(drop=True)
        return ddct_df.sort_values(
            by=["Group", "Condition", "Target Name"]
        ).reset_index(drop=True)

    @reactive.calc
    def foldchange():
        df_dd = ddct_results()
        if df_dd is None:
            return None

        stats = (
            df_dd.groupby(["Group", "Condition", "Target Name", "id"])
            .agg(
                ddct_mean=("ddct", "mean"),
                ddct_error=(
                    "ddct_error",
                    lambda x: np.sqrt(np.sum(x**2)),
                ),
            )
            .reset_index()
        )

        stats["foldchange_mean"] = 2 ** -stats["ddct_mean"]
        stats["foldchange_se"] = (
            np.log(2) * stats["foldchange_mean"] * stats["ddct_error"]
        )
        return stats

    @reactive.calc
    def foldchange_with_stats():
        fc = foldchange()
        df_dd = ddct_results()
        if fc is None or df_dd is None:
            return fc

        control = control_reactive()
        housekeeping = housekeeping_reactive()
        if not control or not housekeeping:
            return fc

        stat_test = stat_test_reactive()
        posthoc_test = posthoc_reactive()

        p_values = []
        for _, row in fc.iterrows():
            target = row["Target Name"]
            grp = row["Group"]
            cond = row["Condition"]

            if target == housekeeping or cond == control:
                p_values.append(np.nan)
                continue

            exp_data = df_dd[
                (df_dd["Target Name"] == target)
                & (df_dd["Group"] == grp)
                & (df_dd["Condition"] == cond)
            ]

            if exp_data.empty:
                p_values.append(np.nan)
                continue

            if stat_test == "none":
                p_values.append(np.nan)
            elif stat_test == "ttest":
                _, p = one_sample_ttest(exp_data["ddct"].values, popmean=0.0)
                p_values.append(p)
            else:
                # Build groups_dict: control = zeros (n = # control CT replicates)
                # Each experimental condition = its ddct values
                df_raw = df_reactive()
                if df_raw is None:
                    p_values.append(np.nan)
                    continue

                conditions = df_dd["Condition"].unique().tolist()
                groups_ddct = {}

                # Control group: zeros with n = # control CT replicates for this target/group
                ct_control = df_raw[
                    (df_raw["Target Name"] == target)
                    & (df_raw["Group"] == grp)
                    & (df_raw["Condition"] == control)
                ]["CT"].values
                n_control = len(ct_control)
                if n_control > 0:
                    groups_ddct[control] = np.zeros(n_control)

                for c in conditions:
                    if c == control:
                        continue
                    exp_c = df_dd[
                        (df_dd["Target Name"] == target)
                        & (df_dd["Group"] == grp)
                        & (df_dd["Condition"] == c)
                    ]
                    if not exp_c.empty:
                        groups_ddct[c] = exp_c["ddct"].values

                ph_df = _run_posthoc(groups_ddct, stat_test, posthoc_test)
                p = _lookup_posthoc_pval(ph_df, control, cond)
                p_values.append(p)

        fc = fc.copy()
        method = correction_method_reactive()
        fc["p_value"] = correct_pvalues(np.array(p_values), method=method)
        fc["significance"] = fc["p_value"].apply(p_to_star)
        return fc

    @reactive.calc
    def anova2_results():
        if stat_test_reactive() != "anova2":
            return None
        df_dd = ddct_results()
        if df_dd is None:
            return None

        housekeeping = housekeeping_reactive()
        control = control_reactive()
        targets = [t for t in df_dd["Target Name"].unique() if t != housekeeping]

        all_rows = []
        for target in targets:
            target_df = df_dd[df_dd["Target Name"] == target][
                ["Group", "Condition", "ddct"]
            ].copy()
            # Include control rows as ddct=0
            df_raw = df_reactive()
            if df_raw is not None and control:
                groups_in_data = df_raw["Group"].unique().tolist()
                control_rows = []
                for grp in groups_in_data:
                    ct_control = df_raw[
                        (df_raw["Target Name"] == target)
                        & (df_raw["Group"] == grp)
                        & (df_raw["Condition"] == control)
                    ]["CT"].values
                    for _ in range(len(ct_control)):
                        control_rows.append({
                            "Group": grp,
                            "Condition": control,
                            "ddct": 0.0,
                        })
                if control_rows:
                    ctrl_df = pd.DataFrame(control_rows)
                    target_df = pd.concat([target_df, ctrl_df], ignore_index=True)

            anova_df = two_way_anova(target_df, "ddct", "Group", "Condition")
            if not anova_df.empty:
                anova_df.insert(0, "Target Name", target)
                all_rows.append(anova_df)

        if not all_rows:
            return None
        return pd.concat(all_rows, ignore_index=True)

    # ---- Downloads ----

    @render.download(filename="ddct_foldchange.xlsx")
    def download_xlsx():
        data = foldchange_with_stats()
        if data is not None:
            yield export_to_excel(data)

    @render.download(filename="ddct_foldchange.csv")
    def download_csv():
        data = foldchange_with_stats()
        if data is not None:
            yield export_to_csv(data)

    # ---- Plots ----

    @reactive.calc
    def foldchange_individual_reps():
        dd = ddct_results()
        summary = foldchange_with_stats()
        if dd is None or summary is None:
            return None
        reps = dd[["id", "Target Name", "Group", "Condition", "ddct"]].copy()
        if "Biological Replicate" in dd.columns:
            reps["Biological Replicate"] = dd["Biological Replicate"]
        reps["fold_change_rep"] = 2.0 ** -reps["ddct"]
        sig = summary[["id", "significance"]].drop_duplicates()
        return reps.merge(sig, on="id", how="left")

    @reactive.calc
    def plot_data_foldchange():
        data = foldchange_with_stats()
        targets = plot_targets_reactive()
        if data is None or not targets:
            return None
        mask = data["Target Name"].isin(list(targets))
        return data[mask]

    @reactive.calc
    def plot_data_foldchange_reps():
        data = foldchange_individual_reps()
        targets = plot_targets_reactive()
        if data is None or not targets:
            return None
        return data[data["Target Name"].isin(list(targets))]

    @reactive.calc
    def plot_data_ddct():
        data = ddct_results()
        targets = plot_targets_reactive()
        if data is None or not targets:
            return None
        mask = data["Target Name"].isin(list(targets))
        return data[mask]

    @render.data_frame
    def results_table():
        return foldchange_with_stats()

    @render.data_frame
    def anova2_table():
        return anova2_results()

    @render_widget
    def foldchange_plot():
        reps = plot_data_foldchange_reps()
        summary = plot_data_foldchange()
        if reps is None or reps.empty:
            return apply_classic_theme(px.scatter(title="No data to display"))

        fig = px.strip(
            reps,
            x="id",
            y="fold_change_rep",
            color="Target Name",
            stripmode="overlay",
            labels={"fold_change_rep": "Fold Change", "id": ""},
            title="Relative Expression (2^-ddCt)",
        )
        fig.update_traces(jitter=0.3, marker=dict(size=8, opacity=0.8))
        fig.update_layout(showlegend=False)
        apply_classic_theme(fig)

        if summary is not None and not summary.empty:
            max_per_id = reps.groupby("id")["fold_change_rep"].max()
            for _, row in summary.iterrows():
                star = row.get("significance", "")
                if star and star not in ("", "ns"):
                    y_pos = max_per_id.get(row["id"], row["foldchange_mean"]) * 1.15
                    fig.add_annotation(
                        x=row["id"],
                        y=y_pos,
                        text=star,
                        showarrow=False,
                        font=dict(size=14),
                    )
        return fig

    @render_widget
    def ddct_plot():
        data = plot_data_ddct()
        if data is None or data.empty:
            return apply_classic_theme(px.scatter(title="No data to display"))

        fig = px.strip(
            data,
            x="id",
            y="ddct",
            color="Target Name",
            labels={"ddct": "ddCT", "id": ""},
            title="Relative Expression (ddCt)",
            stripmode="overlay",
        )
        fig.update_traces(
            jitter=0.3, marker=dict(size=8, opacity=0.7),
            showlegend=False,
        )
        apply_classic_theme(fig)
        return fig

    # Expose results for parent access
    return {"foldchange": foldchange_with_stats, "ddct_full": ddct_results}
