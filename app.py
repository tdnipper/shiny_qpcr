import pandas as pd
from shiny import App, ui, render, reactive
from shinywidgets import render_widget, output_widget
import plotly.graph_objects as go
import plotly.express as px

# UI
app_ui = ui.page_fluid(
    ui.input_file("file", "Upload Data File", accept=[".csv", ".xlsx"], multiple=False),
    ui.card(ui.output_data_frame("table")),
    ui.input_selectize(
        "housekeeping", "Select housekeeping gene:", choices=[], selected="RNA18S1"
    ),
    ui.input_selectize("control", "Select control group:", choices=[], selected="mock"),
    ui.card(ui.output_data_frame("ddct")),
    output_widget("log2fc_plot"),
)


# Server
def calculate_ddct(results, gene, group, control_gene, control_group):
    return (
        results[f"{gene}_{group}"]["CT"].values
        - results[f"{control_gene}_{group}"]["CT"].values
    ) - (
        results[f"{gene}_{control_group}"]["CT"].values
        - results[f"{control_gene}_{control_group}"]["CT"].values
    )


def filter_targets(df, target):
    filtered_df = df[df["Target Name"].str.contains(target)].copy()
    return filtered_df


def server(input, output, session):

    @reactive.calc
    def df():
        """Read in a csv, xlsx, or xls file and return the dataframe for display and analysis."""
        file_info = input.file()
        if file_info is not None:
            file_path = file_info[0]["datapath"]
            if file_path.endswith(".csv"):
                df = (
                    pd.read_csv(file_path, na_values=["Undetermined", "NTC"])
                    .reset_index()
                    .dropna()
                )
            elif file_path.endswith(".xlsx"):
                df = (
                    pd.read_excel(file_path, na_values=["Undetermined", "NTC"])
                    .reset_index()
                    .dropna()
                )
            elif file_path.endswith(".xls"):
                df = (
                    pd.read_excel(file_path, na_values=["Undetermined", "NTC"])
                    .reset_index()
                    .dropna()
                )
            return df

    @reactive.calc
    def mean():
        dataframe = df()
        if dataframe is not None:
            dataframe["mean"] = dataframe.groupby(["Sample Name", "Target Name"])[
                "CT"
            ].transform("mean")
            dataframe["id"] = dataframe[["Target Name", "Sample Name"]].agg(
                "_".join, axis=1
            )
            return dataframe

    @reactive.calc
    def store_filtered_dfs():
        old_df = mean()
        if old_df is not None:
            targets = old_df["Target Name"].unique()
            groups = old_df["Sample Name"].unique()
            target_dict = {}
            for target in targets:
                target_dict[target] = filter_targets(old_df, target)
            results = {}
            for key, data in target_dict.items():
                for group in groups:
                    results[f"{key}_{group}"] = data[data["Sample Name"] == group]
            return results

    @reactive.calc
    def ddCT():
        """Calculate ddCT given table of CTs"""
        results = store_filtered_dfs()
        df = mean()
        housekeeping = input.housekeeping()
        control = input.control()
        if results is not None:
            groups = df["Sample Name"].unique()
            genes = df["Target Name"].unique()
            # control_val = sum()
            ddct_results = {}
            for gene in genes:
                for group in groups:
                    results[f"{gene}_{group}"]["ddct"] = calculate_ddct(
                        results,
                        gene,
                        group=group,
                        control_gene=housekeeping,
                        control_group=control,
                    )
            return results

    @reactive.calc
    def ddct_dfs():
        results = ddCT()
        ddct_df_list = []
        if results is not None:
            for key, result in results.items():
                ddct_df_list.append(result)
            ddct_df = pd.concat(ddct_df_list)
            return ddct_df

    @reactive.calc
    def log2fc():
        if ddct_dfs() is not None:
            data = ddct_dfs()
            data["log2fc"] = 2 ** -data["ddct"]
            data["log2fc_mean"] = data.groupby(["Sample Name", "Target Name"])["log2fc"].transform("mean")
            data["log2fc_std"] = data.groupby(["Sample Name", "Target Name"])["log2fc"].transform("std")
            data["ddct_mean"] = data.groupby(["Sample Name", "Target Name"])["ddct"].transform("mean")
            data["ddct_std"] = data.groupby(["Sample Name", "Target Name"])["ddct"].transform("std")

            new_data = data.drop(["CT", "index", "ddct", "log2fc"], axis=1).reset_index(drop=True).drop_duplicates()
            return new_data

    @reactive.Effect
    # Update when df changes
    @reactive.event(df)
    # Update housekeeping selectize options based on Target Name from df
    def update_selectize_housekeeping():
        dataframe = df()  # Call the reactive function to get the DataFrame
        if not dataframe.empty and "Target Name" in dataframe.columns:
            # Get unique target names
            target_names = sorted(dataframe["Target Name"].dropna().unique())
            # Update the choices in the selectize input
            ui.update_selectize(
                "housekeeping", choices=target_names, selected="RNA18S1"
            )

    @reactive.Effect
    @reactive.event(df)
    def update_selectize_control():
        dataframe = df()
        if not dataframe.empty and "Sample Name" in dataframe.columns:
            group_names = sorted(dataframe["Sample Name"].dropna().unique())
            ui.update_selectize("control", choices=group_names, selected="mock")

    @render.data_frame
    def table():
        return mean()

    @render.data_frame
    def ddct():
        return log2fc()

    @render_widget
    def log2fc_plot():
        data = log2fc()
        if log2fc() is not None:
            plot = px.bar(data, x="id", y="log2fc_mean", color="Target Name", error_y='log2fc_std')
            return plot


app = App(app_ui, server)
