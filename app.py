import pandas as pd
from shiny import App, ui, render

# UI
app_ui = ui.page_fluid(
    ui.input_file("file", "Upload CSV File", accept=[".csv", ".xlsx"], multiple=False),
    ui.output_table("table")
)

# Server
def server(input, output, session):
    @render.data_frame
    def summary():
        if ".csv" in input.file:
            df = pd.read_csv(input.file)
        elif ".xlsx" in input.file:
            df = pd.read_excel(input.file)
        return render.data_frame


app = App(app_ui, server)