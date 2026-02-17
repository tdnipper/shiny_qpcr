from shiny import App, ui, render, reactive
from shared import import_file
from modules.ddct_module import ddct_ui, ddct_server
from modules.dct_module import dct_ui, dct_server
from modules.enrichment_module import enrichment_ui, enrichment_server


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_file(
            "file", "Upload Data File",
            accept=[".csv", ".xlsx"], multiple=False,
        ),
        ui.input_select(
            "analysis_mode",
            "Analysis Mode",
            choices={
                "ddct": "ddCT (Relative Expression)",
                "dct": "dCT (Fold Change Between Groups)",
                "enrichment": "RIP-qPCR Enrichment",
            },
            selected="ddct",
        ),
        # --- ddCT-specific inputs ---
        ui.panel_conditional(
            "input.analysis_mode === 'ddct'",
            ui.input_selectize(
                "housekeeping", "Housekeeping gene:", choices=[],
            ),
            ui.input_selectize(
                "control_ddct", "Control group:", choices=[],
            ),
        ),
        # --- dCT-specific inputs ---
        ui.panel_conditional(
            "input.analysis_mode === 'dct'",
            ui.input_selectize(
                "control_dct", "Control group:", choices=[],
            ),
        ),
        # --- Enrichment-specific inputs ---
        ui.panel_conditional(
            "input.analysis_mode === 'enrichment'",
            ui.input_selectize(
                "input_group", "Input/Reference group:", choices=[],
            ),
            ui.input_numeric(
                "dilution_factor", "Input dilution factor:",
                value=10, min=1,
            ),
            ui.input_checkbox(
                "normalize_igg", "Normalize to IgG control", value=False,
            ),
            ui.panel_conditional(
                "input.normalize_igg",
                ui.input_selectize(
                    "igg_group", "IgG/Negative control group:", choices=[],
                ),
            ),
        ),
        # --- Shared controls ---
        ui.input_selectize(
            "plot_targets", "Select targets to plot:",
            choices=[], multiple=True,
        ),
        ui.accordion(
            ui.accordion_panel(
                "Statistical Options",
                ui.input_select(
                    "correction_method",
                    "Multiple testing correction:",
                    choices={
                        "none": "None",
                        "bonferroni": "Bonferroni",
                        "bh": "Benjamini-Hochberg",
                    },
                    selected="none",
                ),
            ),
            open=False,
        ),
        width=350,
    ),
    # --- Main content: raw data preview ---
    ui.card(
        ui.card_header("Uploaded Data"),
        ui.output_data_frame("raw_table"),
    ),
    # --- Mode-specific outputs (conditionally shown) ---
    ui.panel_conditional(
        "input.analysis_mode === 'ddct'",
        ddct_ui("ddct"),
    ),
    ui.panel_conditional(
        "input.analysis_mode === 'dct'",
        dct_ui("dct"),
    ),
    ui.panel_conditional(
        "input.analysis_mode === 'enrichment'",
        enrichment_ui("enrichment"),
    ),
    title="qPCR Analysis",
)


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

def server(input, output, session):

    # --- Shared reactive: uploaded data ---

    @reactive.calc
    def df():
        file_info = input.file()
        if file_info is not None:
            file_path = file_info[0]["datapath"]
            return import_file(file_path)
        return None

    @render.data_frame
    def raw_table():
        return df()

    # --- Wrapper reactives for sidebar inputs ---
    # These let us pass sidebar values into modules as reactive callables.

    @reactive.calc
    def housekeeping():
        return input.housekeeping()

    @reactive.calc
    def control_ddct():
        return input.control_ddct()

    @reactive.calc
    def control_dct():
        return input.control_dct()

    @reactive.calc
    def input_group():
        return input.input_group()

    @reactive.calc
    def dilution_factor():
        return input.dilution_factor()

    @reactive.calc
    def normalize_igg():
        return input.normalize_igg()

    @reactive.calc
    def igg_group():
        return input.igg_group()

    @reactive.calc
    def plot_targets():
        return input.plot_targets()

    @reactive.calc
    def correction_method():
        return input.correction_method()

    # --- Populate selectize dropdowns when data changes ---

    @reactive.effect
    @reactive.event(df)
    def update_selectize_inputs():
        dataframe = df()
        if dataframe is None or dataframe.empty:
            return

        target_names = sorted(
            dataframe["Target Name"].dropna().unique().tolist()
        )
        group_names = sorted(
            dataframe["Sample Name"].dropna().unique().tolist()
        )

        # ddCT mode
        ui.update_selectize(
            "housekeeping", choices=target_names,
            selected=target_names[0] if target_names else None,
        )
        ui.update_selectize(
            "control_ddct", choices=group_names,
            selected=group_names[0] if group_names else None,
        )

        # dCT mode
        ui.update_selectize(
            "control_dct", choices=group_names,
            selected=group_names[0] if group_names else None,
        )

        # Enrichment mode
        ui.update_selectize(
            "input_group", choices=group_names,
            selected=group_names[0] if group_names else None,
        )
        ui.update_selectize(
            "igg_group", choices=group_names,
            selected=group_names[1] if len(group_names) > 1 else None,
        )

        # Shared: targets for plotting
        ui.update_selectize(
            "plot_targets", choices=target_names,
            selected=target_names[1] if len(target_names) > 1 else (
                target_names[0] if target_names else None
            ),
        )

    # --- Wire up modules ---

    ddct_server(
        "ddct",
        df_reactive=df,
        housekeeping_reactive=housekeeping,
        control_reactive=control_ddct,
        plot_targets_reactive=plot_targets,
        correction_method_reactive=correction_method,
    )

    dct_server(
        "dct",
        df_reactive=df,
        control_reactive=control_dct,
        plot_targets_reactive=plot_targets,
        correction_method_reactive=correction_method,
    )

    enrichment_server(
        "enrichment",
        df_reactive=df,
        input_group_reactive=input_group,
        dilution_factor_reactive=dilution_factor,
        normalize_igg_reactive=normalize_igg,
        igg_group_reactive=igg_group,
        plot_targets_reactive=plot_targets,
        correction_method_reactive=correction_method,
    )


app = App(app_ui, server)
