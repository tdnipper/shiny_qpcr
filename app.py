import pandas as pd
from shiny import App, ui, render, reactive

# UI
app_ui = ui.page_fluid(
    ui.input_file("file", "Upload Data File", accept=[".csv", ".xlsx"], multiple=False),
    ui.card(ui.output_data_frame("table")),
    ui.input_selectize("housekeeping", "Select housekeeping gene:", choices=[])
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
                df = pd.read_csv(file_path)
            elif file_path.endswith(".xlsx"):
                df = pd.read_excel(file_path)
            elif file_path.endswith(".xls"):
                df = pd.read_excel(file_path)
            return df

    @reactive.Effect  
    # Update when df changes
    @reactive.event(df)   
    # Update housekeeping selectize options based on Target Name from df   
    def update_selectize():
        dataframe = df()  # Call the reactive function to get the DataFrame
        if not dataframe.empty and "Target Name" in dataframe.columns:
            # Get unique target names
            target_names = sorted(dataframe["Target Name"].dropna().unique())
            # Update the choices in the selectize input
            ui.update_selectize("housekeeping", choices=target_names)
    
    @output
    @render.data_frame
    def table():
        return df()


app = App(app_ui, server)