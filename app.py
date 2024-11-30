import pandas as pd
from shiny import App, ui, render, reactive

# UI
app_ui = ui.page_fluid(
    ui.input_file("file", "Upload Data File", accept=[".csv", ".xlsx"], multiple=False),
    ui.card(ui.output_data_frame("table")),
    ui.input_selectize("housekeeping", "Select housekeeping gene:", choices=[]),
    ui.input_selectize("control", "Select control group:", choices=[]),
    ui.card(ui.output_data_frame("ddct"))
)

# Server
def calculate_ddct(results, gene, group, control_gene, control_group):
    return(
        results[f'{gene}_{group}']['mean'].values - results[f'{control_gene}_{group}']['mean'].values) - (results[f'{gene}_{control_group}']['mean'].values - results[f'{control_gene}_{control_group}']['mean'].values) 

def filter_targets(df, target):
    filtered_df = df[df['Target Name'].str.contains(target)].copy()
    return filtered_df

def server(input, output, session):
    
    @reactive.calc
    def df():
        """ Read in a csv, xlsx, or xls file and return the dataframe for display and analysis."""
        file_info = input.file()
        if file_info is not None:
            file_path = file_info[0]["datapath"]
            if file_path.endswith(".csv"):
                df = pd.read_csv(file_path, na_values=["Undetermined", "NTC"]).reset_index().dropna()
            elif file_path.endswith(".xlsx"):
                df = pd.read_excel(file_path, na_values=["Undetermined", "NTC"]).reset_index().dropna()
            elif file_path.endswith(".xls"):
                df = pd.read_excel(file_path, na_values=["Undetermined", "NTC"]).reset_index().dropna()
            return df
   
    @reactive.calc
    def mean():
        dataframe = df()
        if dataframe is not None:   
            dataframe["mean"] = dataframe.groupby(["Sample Name", "Target Name"])["CT"].transform("mean")
            dataframe['id'] = dataframe[['Target Name', 'Sample Name']].agg('_'.join, axis=1)
            return dataframe
        
    @reactive.calc
    def ddCT():
        """Calculate ddCT given table of CTs"""
        results = store_filtered_dfs()
        df = mean()
        housekeeping = input.housekeeping()
        control = input.control()
        if results is not None:
            groups = df['Sample Name'].unique()
            genes = df['Target Name'].unique()
            ddct_results = {}
            for gene in genes:
                for group in groups:
                    ddct_results[f'{gene}_{group}'] = calculate_ddct(results, gene, group=group, control_gene=housekeeping, control_group=control)
            return ddct_results

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
            ui.update_selectize("housekeeping", choices=target_names)
    @reactive.Effect
    @reactive.event(df)
    def update_selectize_control():
        dataframe = df()
        if not dataframe.empty and "Sample Name" in dataframe.columns:
            group_names = sorted(dataframe['Sample Name'].dropna().unique())
            ui.update_selectize("control", choices=group_names)

    @reactive.calc
    def store_filtered_dfs():
        old_df = mean()
        if old_df is not None:
            targets = [old_df['Target Name'].unique()]
            groups = [old_df['Sample Name'].unique()]        
            target_dict = {}
            for target in targets:
                target_dict[target] = filter_targets(old_df, target)
            results = {}
            for key, data in target_dict.items():
                for group in groups:
                    results[f'{key}_{group}'] = data[data['Sample Name'] == group]
            return results

    @output
    
    @render.data_frame
    def table():
        return mean()

    # @render.text
    # def ddct():
    #     return ddCT()

app = App(app_ui, server)