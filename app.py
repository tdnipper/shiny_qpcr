import pandas as pd
from shiny import App, ui, render, reactive

# UI
app_ui = ui.page_fluid(
    ui.input_file("file", "Upload Data File", accept=[".csv", ".xlsx"], multiple=False),
    ui.output_data_frame("table")
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
    @output
    @render.data_frame
    def table():
        return df()


app = App(app_ui, server)