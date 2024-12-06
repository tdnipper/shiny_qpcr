import pandas as pd
from shiny import App, ui, render, reactive
from shinywidgets import render_widget, output_widget

# import plotly.graph_objects as go
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
    ui.input_selectize("plot_groups", "Select groups to plot:", choices=[], multiple=True),
    ui.card(output_widget("log2fc_plot"))
)


# Functions to be used in server calcs

def calculate_ddct(results, gene, group, control_gene, control_group):
    try:
        if gene == control_gene and group == control_group:
            # Simplified calculation since the comparison is self-referential
            return 0

        gene_group_ct = results[f"{gene}_{group}"]["CT"].values
        control_gene_group_ct = results[f"{control_gene}_{group}"]["CT"].values
        gene_control_group_ct = results[f"{gene}_{control_group}"]["CT"].values
        control_gene_control_group_ct = results[f"{control_gene}_{control_group}"]["CT"].values
    except KeyError as e:
        raise KeyError(f"Missing key in results: {e}. Check if all input values exist in the data.")

    return (gene_group_ct - control_gene_group_ct) - (gene_control_group_ct - control_gene_control_group_ct)


def filter_targets(df, target):
    filtered_df = df[df["Target Name"].str.contains(target)].copy()
    return filtered_df


# Server
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
                    .reset_index(drop=True)
                    .dropna()
                )
            elif file_path.endswith(".xlsx"):
                df = (
                    pd.read_excel(file_path, na_values=["Undetermined", "NTC"])
                    .reset_index(drop=True)
                    .dropna()
                )
            elif file_path.endswith(".xls"):
                df = (
                    pd.read_excel(file_path, na_values=["Undetermined", "NTC"])
                    .reset_index(drop=True)
                    .dropna()
                )
            return df

    @reactive.calc
    def mean():
        dataframe = df()
        if dataframe is not None:
            dataframe["mean"] = dataframe.groupby(["Sample Name", "Target Name"])[
                "CT"
            ].transform("mean").reset_index(drop=True)
            dataframe["id"] = dataframe["Target Name"] + "_" + dataframe["Sample Name"]
            return dataframe

    @reactive.calc
    def store_filtered_dfs():
        """Filter a df by gene and group into a dict of dfs with key as 'gene_group'."""
        old_df = mean()
        target_dict = {}
        results = {}
        if old_df is not None:
            genes = old_df["Target Name"].unique().tolist()
            groups = old_df["Sample Name"].unique().tolist()
            # Sort results by gene in target_dict
            for gene in genes:
                target_dict[gene] = filter_targets(old_df, gene).reset_index(drop=True)
            # Sort gene dfs by group in results dict
            for target, data in target_dict.items():
                for group in groups:
                    results[f"{target}_{group}"] = data[data["Sample Name"] == group]
            return results

    @reactive.calc
    def ddCT():
        """Calculate ddCT given table of CTs"""
        results = store_filtered_dfs()
        df = mean()

        if results is not None and df is not None:
            groups = df["Sample Name"].unique().tolist()
            genes = df["Target Name"].unique().tolist()
            housekeeping = input.housekeeping()
            control = input.control()
            for gene in genes:
                for group in groups:
                    results[f"{gene}_{group}"]["ddct"] = calculate_ddct(
                        results=results,
                        gene=gene,
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
            data["log2fc_mean"] = data.groupby(["Sample Name", "Target Name"])[
                "log2fc"
            ].transform("mean")
            data["log2fc_std"] = data.groupby(["Sample Name", "Target Name"])[
                "log2fc"
            ].transform("std")
            data["ddct_mean"] = data.groupby(["Sample Name", "Target Name"])[
                "ddct"
            ].transform("mean")
            data["ddct_std"] = data.groupby(["Sample Name", "Target Name"])[
                "ddct"
            ].transform("std")

            new_data = (
                data.drop(["CT", "ddct", "log2fc"], axis=1)
                .reset_index(drop=True)
                .drop_duplicates()
                .reset_index(drop=True)
            )
            return new_data
    
    @reactive.calc
    def select_for_plot():
        data = log2fc()
        groups = list(input.plot_groups())
        if data is not None:
            new_data_list = []
            for group in groups:
                data_sub = filter_targets(data, group)
                new_data_list.append(data_sub)
            new_data = pd.concat(new_data_list, axis=0)
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
    # Update control group selectize options based on Sample Name from df
    def update_selectize_control():
        dataframe = df()
        if not dataframe.empty and "Sample Name" in dataframe.columns:
            group_names = sorted(dataframe["Sample Name"].dropna().unique())
            ui.update_selectize("control", choices=group_names, selected="mock")

    @reactive.Effect
    @reactive.event(mean)
    def update_selectize_plot_groups():
        dataframe = mean()
        if dataframe is not None and not dataframe.empty:
            group_names = dataframe["Target Name"].unique().tolist()
            ui.update_selectize("plot_groups", choices=group_names, selected=group_names[1] if group_names else None)

    # Rendering
    @render.data_frame
    def table():
        return mean()

    @render.data_frame
    def ddct():
        return log2fc()

    @render_widget
    def log2fc_plot():
        data = select_for_plot()
        if select_for_plot() is not None:
            plot = px.bar(
                data, x="id", y="log2fc_mean", color="Target Name", error_y="log2fc_std"
            )
            return plot    


app = App(app_ui, server)
