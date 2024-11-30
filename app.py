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
            return dataframe
        
    @reactive.calc
    def ddCT():
        """Calculate ddCT given table of CTs"""
        dataframe = df()
        dct_results = pd.DataFrame()
        housekeeping = input.housekeeping()
        control = input.control()
        if dataframe is not None:
            groups = dataframe['Sample Name'].unique()
            dataframe['id'] = dataframe[['Sample Name', 'Target Name']].agg('_'.join, axis=1)
            return dct_results

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

    @output
    
    @render.data_frame
    def table():
        return df()

    @render.data_frame
    def ddct():
        return ddCT()

app = App(app_ui, server)