import pandas as pd
from shiny import App, ui, render, reactive
from shinywidgets import render_widget, output_widget
import io
import numpy as np

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
    ui.download_button("foldchange_data_download", "Download foldchange data"),
    ui.download_button("full_data_download", "Download full ddCT data"),
    ui.input_selectize("plot_groups", "Select groups to plot:", choices=[], multiple=True),
    ui.card(output_widget("foldchange_plot")),
    ui.card(output_widget("ddct_plot"))
)


# Functions to be used in server calcs

def import_file(file):
    if file.endswith(".csv"):
        df = pd.read_csv(file, na_values=["Undetermined", "NTC"]).dropna().reset_index(drop=True)
    elif (file.endswith(".xls") or file.endswith(".xlsx")):
        df = pd.read_excel(file, na_values=["Undetermined", "NTC"]).dropna().reset_index(drop=True)
    return df

def calculate_ddct(results, gene, group, control_gene, control_group):
    try:
        if gene == control_gene and group == control_group:
            # Simplified calculation since the comparison is self-referential
            return 0, 0  # ddCT and its propagated error

        # Retrieve CT values and SEMs
        # Experimental gene and group
        gene_group_ct = results[f"{gene}_{group}"]["CT"].values
        # gene_group_ct_sem = np.std(gene_group_ct) / np.sqrt(len(gene_group_ct))
        # Control gene and experimental group
        control_gene_group_ct = results[f"{control_gene}_{group}"]["CT"].values
        # control_gene_group_ct_sem = np.std(control_gene_group_ct) / np.sqrt(len(control_gene_group_ct))
        # Experimental gene and control group
        gene_control_group_ct = results[f"{gene}_{control_group}"]["CT"].values
        gene_control_group_ct_sem = np.std(gene_control_group_ct) / np.sqrt(len(gene_control_group_ct))
        # Control gene and control group
        control_gene_control_group_ct = results[f"{control_gene}_{control_group}"]["CT"].values
        # control_gene_control_group_ct_mean = np.mean(control_gene_control_group_ct)
        control_gene_control_group_ct_sem = np.std(control_gene_control_group_ct) / np.sqrt(len(control_gene_control_group_ct))
    except KeyError as e:
        raise KeyError(f"Missing key in results: {e}. Check if all input values exist in the data.")

    # Calculate ddCT
    dct_control_average = np.mean(gene_control_group_ct - control_gene_control_group_ct)
    dct_control_average_sem = np.sqrt(gene_control_group_ct_sem**2 + control_gene_control_group_ct_sem**2)

    # ddct = (gene_group_ct - control_gene_group_ct) - (gene_control_group_ct - control_gene_control_group_ct_mean)
    ddct = (gene_group_ct - control_gene_group_ct) - dct_control_average

    # Propagate the error
    # sem_dct_sample = np.sqrt(
    #     gene_group_ct_sem**2 +
    #     control_gene_group_ct_sem**2
    # )

    # sem_dct_control = np.sqrt(
    #     gene_control_group_ct_sem**2 +
    #     control_gene_control_group_ct_sem**2
    # )
    # Calculate the propagated error
    # ddct_error = np.sqrt(
    #     sem_dct_sample**2
    #     sem_dct_control**2
    # )
    ddct_error = dct_control_average_sem

    return ddct, ddct_error


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
            df = import_file(file_path)
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
                    results[f"{gene}_{group}"]["ddct"], results[f"{gene}_{group}"]["ddct_error"] = calculate_ddct(
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
            ddct_df = pd.concat(ddct_df_list).reset_index(drop=True)
            ddct_df_sorted = ddct_df.sort_values(by=["Sample Name", "Target Name"]).reset_index(drop=True)
            return ddct_df_sorted

    @reactive.calc
    def foldchange():
        """Calculate mean fold‑change and its propagated error (standard error)."""
        df_dd = ddct_dfs()
        if df_dd is None:
            return None

        #1) summarize ΔΔCt per Sample+Target
        stats = (
            df_dd
            .groupby(["Sample Name", "Target Name", "id"])
            .agg(
                ddct_mean=("ddct", "mean"),
                # Propagated error of ddCT 
                ddct_error=("ddct_error", lambda x: np.sqrt(np.sum(x**2)))
                )
            .reset_index()
        )

        # 2) fold‑change and its error via delta‑method:
        stats['foldchange_mean'] = 2 ** -stats['ddct_mean']
        stats['foldchange_se']   = np.log(2) * stats['foldchange_mean'] * stats['ddct_error']

        return stats
    
    @render.download(filename = 'foldchange_data.xlsx')
    def foldchange_data_download():
        if foldchange() is not None:
            data = foldchange()
            with io.BytesIO() as buf:
                data.to_excel(buf)
                yield buf.getvalue()
    
    @render.download(filename = 'ddCT_data.xlsx')
    def full_data_download():
        if ddct_dfs() is not None:
            data = ddct_dfs()
            with io.BytesIO() as buf:
                data.to_excel(buf)
                yield buf.getvalue()
    
    @reactive.calc
    def select_for_ddct_plot():
        data = ddct_dfs()
        groups = list(input.plot_groups())
        if data is not None:
            new_data_list = []
            for group in groups:
                data_sub = filter_targets(data, group)
                new_data_list.append(data_sub)
            new_data = pd.concat(new_data_list, axis=0)
            return new_data

    @reactive.calc
    def select_for_plot():
        data = foldchange()
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
        return foldchange()
    
    @render_widget
    def ddct_plot():
        data = select_for_ddct_plot()
        if data is not None and not data.empty:
            # use propagated SE instead of raw SD
            plot = px.box(
                data,
                x="id",
                y="ddct",
                # box=True,
                # error_y="ddct_error",
                # stripmode="overlay",
                color="Target Name",
                labels={
                    "ddct": "ddCT",
                    "id": ""
                },
                title="Relative Expression (ΔΔCt)",
                template="simple_white",
            )
            plot.update_traces(showlegend=False)
            return plot

    @render_widget
    def foldchange_plot():
        data = select_for_plot()
        if data is not None and not data.empty:
            # use propagated SE instead of raw SD
            plot = px.scatter(
                data,
                x="id",
                y="foldchange_mean",
                symbol="Target Name",
                color="Target Name",
                error_y="foldchange_se",
                labels={
                    "foldchange_mean": "Fold Change (mean ± SE)",
                    "id": ""
                },
                title="Relative Expression (2^⁻ΔΔCt) with Propagated Error",
                template="simple_white",
            )
            plot.update_layout(showlegend=False)
            return plot


app = App(app_ui, server)
