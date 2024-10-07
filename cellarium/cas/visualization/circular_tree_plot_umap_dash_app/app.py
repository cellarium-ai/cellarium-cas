import os
import tempfile
import typing as t
from collections import OrderedDict

import dash_bootstrap_components as dbc
import numpy as np
import plotly.graph_objects as go
from anndata import AnnData
from Bio import Phylo
from dash import Dash, State, dcc, html
from dash.dependencies import Input, Output
from dash.development.base_component import Component
from plotly.express.colors import sample_colorscale

from cellarium.cas.logging import logger
from cellarium.cas.postprocessing import (
    CAS_CL_SCORES_ANNDATA_OBSM_KEY,
    CAS_METADATA_ANNDATA_UNS_KEY,
    CellOntologyScoresAggregationDomain,
    CellOntologyScoresAggregationOp,
    convert_aggregated_cell_ontology_scores_to_rooted_tree,
    generate_phyloxml_from_scored_cell_ontology_tree,
    get_aggregated_cas_ontology_aware_scores,
    get_obs_indices_for_cluster,
)
from cellarium.cas.postprocessing.cell_ontology import CL_CELL_ROOT_NODE, CellOntologyCache
from cellarium.cas.visualization._components.circular_tree_plot import CircularTreePlot
from cellarium.cas.visualization.ui_utils import ConfigValue, find_and_kill_process

# cell type ontology terms (and all descendents) to hide from the visualization
DEFAULT_HIDDEN_CL_NAMES_SET = {}

# header component ID -> default title mapping for the panels
DEFAULT_PANEL_TITLES = {
    "cell-selection-title-tree": "Cell Type Ontology View",
    "cell-selection-title-umap": "UMAP View",
}


class DomainSelectionConstants:
    NONE = 0
    USER_SELECTION = 1
    SEPARATOR = 2


# cell type ontology terms to always show as text labels in the visualization
DEFAULT_SHOWN_CL_NAMES_SET = {
    "CL_0000236",
    "CL_0000084",
    "CL_0000789",
    "CL_0000798",
    "CL_0002420",
    "CL_0002419",
    "CL_0000786",
    "CL_0000576",
    "CL_0001065",
    "CL_0000451",
    "CL_0000094",
    "CL_0000235",
    "CL_0000097",
    "CL_0000814",
    "CL_0000827",
    "CL_0000066",
    "CL_0000163",
    "CL_0000151",
    "CL_0000064",
    "CL_0000322",
    "CL_0000076",
    "CL_0005006",
    "CL_0000148",
    "CL_0000646",
    "CL_0009004",
    "CL_0000115",
    "CL_0000125",
    "CL_0002319",
    "CL_0000187",
    "CL_0000057",
    "CL_0008034",
    "CL_0000092",
    "CL_0000058",
    "CL_0000060",
    "CL_0000136",
    "CL_0000499",
    "CL_0000222",
    "CL_0007005",
    "CL_0000039",
    "CL_0000019",
    "CL_0000223",
    "CL_0008019",
    "CL_0005026",
    "CL_0000182",
    "CL_0000023",
    "CL_0000679",
    "CL_0000126",
    "CL_0000540",
    "CL_0000127",
    "CL_0011005",
}


class CASCircularTreePlotUMAPDashApp:
    """
    A Dash app for visualizing the results of a Cellarium CAS cell type ontology-aware analysis.

    :param adata: The AnnData object containing the cell type ontology-aware analysis results.
    :param cas_ontology_aware_response: The response from the Cellarium CAS cell type ontology-aware analysis. |br|
        `Default:` ``None``
    :param cluster_label_obs_column: The name of the observation column containing the cluster labels. |br|
        `Default:` ``None``
    :param aggregation_op: The aggregation operation to apply to the cell type ontology-aware scores. |br|
        `Default:` ``CellOntologyScoresAggregationOp.MEAN``
    :param aggregation_domain: The domain over which to aggregate the cell type ontology-aware scores. |br|
        `Default:` ``CellOntologyScoresAggregationDomain.OVER_THRESHOLD``
    :param score_threshold: The threshold for the cell type ontology-aware scores. |br|
        `Default:` ``0.05``
    :param min_cell_fraction: The minimum fraction of cells that must have a cell type ontology-aware score above the threshold. |br|
        `Default:` ``0.01``
    :param umap_marker_size: The size of the markers in the UMAP scatter plot. |br|
        `Default:` ``3.0``
    :param umap_padding: The padding to apply to the UMAP scatter plot bounds. |br|
        `Default:` ``0.15``
    :param umap_min_opacity: The minimum opacity for the UMAP scatter plot markers. |br|
        `Default:` ``0.1``
    :param umap_max_opacity: The maximum opacity for the UMAP scatter plot markers. |br|
        `Default:` ``1.0``
    :param umap_inactive_cell_color: The color for inactive cells in the UMAP scatter plot. |br|
        `Default:` ``"rgb(180,180,180)"``
    :param umap_inactive_cell_opacity: The opacity for inactive cells in the UMAP scatter plot. |br|
        `Default:` ``0.5``
    :param umap_active_cell_color: The color for active cells in the UMAP scatter plot. |br|
        `Default:` ``"rgb(250,50,50)"``
    :param umap_default_cell_color: The default color for cells in the UMAP scatter plot. |br|
        `Default:` ``"rgb(180,180,180)"``
    :param umap_default_opacity: The default opacity for cells in the UMAP scatter plot. |br|
        `Default:` ``0.9``
    :param circular_tree_plot_linecolor: The line color for the circular tree plot. |br|
        `Default:` ``"rgb(200,200,200)"``
    :param circular_tree_start_angle: The start angle for the circular tree plot. |br|
        `Default:` ``180``
    :param circular_tree_end_angle: The end angle for the circular tree plot. |br|
        `Default:` ``360``
    :param figure_height: The height of the figures in the Dash app. |br|
        `Default:` ``400``
    :param hidden_cl_names_set: The set of cell type ontology terms to hide from the visualization. |br|
        `Default:` ``DEFAULT_HIDDEN_CL_NAMES_SET``
    :param shown_cl_names_set: The set of cell type ontology terms to always show as text labels in the
        visualization. |br|
        `Default:` ``DEFAULT_SHOWN_CL_NAMES_SET``
    :param score_colorscale: The colorscale to use for the cell type ontology-aware scores. |br|
        `Default:` ``"Viridis"``

    Example:
    ________
        >>> from cellarium.cas._io import suppress_stderr
        >>> from cellarium.cas.visualization import CASCircularTreePlotUMAPDashApp
        >>> DASH_SERVER_PORT = 8050
        >>> adata = ... # get your matrix
        >>> cas_ontology_aware_response = cas.annotate_matrix_cell_type_ontology_aware_strategy(
        >>>     matrix=adata,
        >>>     chunk_size=500
        >>> )
        >>> with suppress_stderr():
        >>>     CASCircularTreePlotUMAPDashApp(
        >>>         adata,
        >>>         cas_ontology_aware_response,
        >>>         cluster_label_obs_column="cluster_label",
        >>>     ).run(port=DASH_SERVER_PORT, debug=False, jupyter_width="100%")
    """

    ALL_CELLS_DOMAIN_KEY = "all cells"
    CLUSTER_PREFIX_DOMAIN_KEY = "cluster "

    def __init__(
        self,
        adata: AnnData,
        cluster_label_obs_column: t.Optional[str] = None,
        aggregation_op: CellOntologyScoresAggregationOp = CellOntologyScoresAggregationOp.MEAN,
        aggregation_domain: CellOntologyScoresAggregationDomain = CellOntologyScoresAggregationDomain.OVER_THRESHOLD,
        score_threshold: float = 0.05,
        min_cell_fraction: float = 0.01,
        umap_marker_size: float = 3.0,
        umap_padding: float = 0.15,
        umap_min_opacity: float = 0.1,
        umap_max_opacity: float = 1.0,
        umap_inactive_cell_color: str = "rgb(180,180,180)",
        umap_inactive_cell_opacity: float = 0.5,
        umap_active_cell_color: str = "rgb(250,50,50)",
        umap_default_cell_color: str = "rgb(180,180,180)",
        umap_default_opacity: float = 0.9,
        circular_tree_plot_linecolor: str = "rgb(200,200,200)",
        circular_tree_start_angle: int = 180,
        circular_tree_end_angle: int = 360,
        figure_height: int = 400,
        root_node: str = CL_CELL_ROOT_NODE,
        hidden_cl_names_set: set[str] = DEFAULT_HIDDEN_CL_NAMES_SET,
        shown_cl_names_set: set[str] = DEFAULT_SHOWN_CL_NAMES_SET,
        score_colorscale: t.Union[str, list] = "Viridis",
    ):
        self.adata = adata
        self.aggregation_op = aggregation_op
        self.aggregation_domain = aggregation_domain
        self.score_threshold = ConfigValue(score_threshold)
        self.min_cell_fraction = ConfigValue(min_cell_fraction)
        self.umap_min_opacity = umap_min_opacity
        self.umap_max_opacity = umap_max_opacity
        self.umap_marker_size = umap_marker_size
        self.umap_padding = umap_padding
        self.umap_inactive_cell_color = umap_inactive_cell_color
        self.umap_inactive_cell_opacity = umap_inactive_cell_opacity
        self.umap_active_cell_color = umap_active_cell_color
        self.umap_default_cell_color = umap_default_cell_color
        self.umap_default_opacity = umap_default_opacity
        self.circular_tree_plot_linecolor = circular_tree_plot_linecolor
        self.circular_tree_start_angle = circular_tree_start_angle
        self.circular_tree_end_angle = circular_tree_end_angle
        self.height = figure_height
        self.root_node = root_node
        self.hidden_cl_names_set = hidden_cl_names_set
        self.shown_cl_names_set = shown_cl_names_set
        self.score_colorscale = score_colorscale

        assert "X_umap" in adata.obsm, (
            "UMAP coordinates not found in adata.obsm['X_umap']. "
            "This visualisation requires precomputed UMAP coordinates."
        )
        assert (CAS_CL_SCORES_ANNDATA_OBSM_KEY in adata.obsm) and (CAS_METADATA_ANNDATA_UNS_KEY in adata.uns), (
            "Cell type ontology scores not found in the provided AnnData file. Please please run "
            "`cellarium.cas.insert_cas_ontology_aware_response_into_adata` prior to running this visualisation."
        )

        # setup cell domains
        self.cell_domain_map = OrderedDict()
        self.cell_domain_map[self.ALL_CELLS_DOMAIN_KEY] = np.arange(adata.n_obs)
        if cluster_label_obs_column is not None:
            assert cluster_label_obs_column in adata.obs
            for cluster_label in adata.obs[cluster_label_obs_column].cat.categories:
                self.cell_domain_map[self.CLUSTER_PREFIX_DOMAIN_KEY + cluster_label] = get_obs_indices_for_cluster(
                    adata, cluster_label_obs_column, cluster_label
                )

        # default cell domain
        self.selected_cell_domain_key = ConfigValue(DomainSelectionConstants.NONE)
        # Selected cells (from UMAP chart)
        self.selected_cells = []
        # Selected cell class (from tree diagram)
        self.selected_cl_name = None

        # instantiate the cell type ontology cache
        self.cl = CellOntologyCache()

        # instantiate the Dash app
        self.app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP])
        self.server = self.app.server
        self.app.layout = self.__create_layout()
        self.__setup_initialization()
        self.__setup_callbacks()

    def run(self, port: int = 8050, jupyter_mode: str = "inline", **kwargs):
        """
        Run the Dash application on the specified port.

        :param port: The port on which to run the Dash application. |br|
            `Default:` ``8050``
        """
        logger.info(f"Starting Dash application on port {port}...")
        try:
            self.app.run_server(port=port, jupyter_mode=jupyter_mode, jupyter_height=self.height + 100, **kwargs)
        except OSError:  # Dash raises OSError if the port is already in use
            find_and_kill_process(port)
            self.app.run_server(port=port, jupyter_mode=jupyter_mode, jupyter_height=self.height + 100, **kwargs)

    def __instantiate_circular_tree_plot(self) -> CircularTreePlot:
        # reduce scores over the provided cells
        selected_cells = self.__get_effective_selected_cells()
        aggregated_scores = get_aggregated_cas_ontology_aware_scores(
            self.adata,
            obs_indices=(
                self.cell_domain_map[self.ALL_CELLS_DOMAIN_KEY] if len(selected_cells) == 0 else selected_cells
            ),
            aggregation_op=self.aggregation_op,
            aggregation_domain=self.aggregation_domain,
            threshold=self.score_threshold.get(),
        )

        # generate a Phylo tree
        rooted_tree = convert_aggregated_cell_ontology_scores_to_rooted_tree(
            aggregated_scores=aggregated_scores,
            cl=self.cl,
            root_cl_name=self.root_node,
            min_fraction=self.min_cell_fraction.get(),
            hidden_cl_names_set=self.hidden_cl_names_set,
        )
        phyloxml_string = generate_phyloxml_from_scored_cell_ontology_tree(
            rooted_tree, "Scored cell type ontology tree", self.cl, indent=3
        )

        with tempfile.NamedTemporaryFile(delete=False, mode="w+t") as temp_file:
            temp_file_name = temp_file.name
            temp_file.write(phyloxml_string)
            temp_file.flush()

        try:
            phyloxml_tree = Phylo.read(temp_file_name, "phyloxml")
        finally:
            os.remove(temp_file_name)

        return CircularTreePlot(
            tree=phyloxml_tree,
            score_colorscale=self.score_colorscale,
            linecolor=self.circular_tree_plot_linecolor,
            start_angle=self.circular_tree_start_angle,
            end_angle=self.circular_tree_end_angle,
            shown_cl_names_set=self.shown_cl_names_set,
        )

    def __get_padded_umap_bounds(self, umap_padding: float) -> t.Tuple[float, float, float, float]:
        actual_min_x = np.min(self.adata.obsm["X_umap"][:, 0])
        actual_max_x = np.max(self.adata.obsm["X_umap"][:, 0])
        actual_min_y = np.min(self.adata.obsm["X_umap"][:, 1])
        actual_max_y = np.max(self.adata.obsm["X_umap"][:, 1])
        padded_min_x = actual_min_x - umap_padding * (actual_max_x - actual_min_x)
        padded_max_x = actual_max_x + umap_padding * (actual_max_x - actual_min_x)
        padded_min_y = actual_min_y - umap_padding * (actual_max_y - actual_min_y)
        padded_max_y = actual_max_y + umap_padding * (actual_max_y - actual_min_y)

        return padded_min_x, padded_max_x, padded_min_y, padded_max_y

    def __get_scores_for_cl_name(self, cl_name: str) -> np.ndarray:
        cl_index = self.cl.cl_names_to_idx_map[cl_name]
        return self.adata.obsm[CAS_CL_SCORES_ANNDATA_OBSM_KEY][:, cl_index].toarray().flatten()

    def __get_scatter_plot_opacity_from_scores(self, scores: np.ndarray) -> np.ndarray:
        min_score = np.min(scores)
        max_score = np.max(scores)
        normalized_scores = (scores - min_score) / (1e-6 + max_score - min_score)
        return np.maximum(
            scores, self.umap_min_opacity + (self.umap_max_opacity - self.umap_min_opacity) * normalized_scores
        )

    def __create_layout(self):
        layout = html.Div(
            [
                dbc.Row(dbc.Col(className="gr-spacer", width=12)),
                dbc.Row(
                    dbc.Col(
                        [
                            html.H3(self.__render_breadcrumb(), id="selected-domain-label", className="gr-breadcrumb"),
                            html.Div(
                                [
                                    dbc.ButtonGroup(
                                        [
                                            dbc.Button(
                                                html.I(className="bi bi-gear-fill"),
                                                id="settings-button",
                                                n_clicks=0,
                                                size="sm",
                                            ),
                                        ]
                                    )
                                ],
                                className="gr-settings-buttons",
                            ),
                        ],
                        className="gr-title",
                        width=12,
                    )
                ),
                dbc.Row(
                    self.__render_cell_selection_panes(),
                    id="panel-titles",
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            html.Div(
                                [
                                    dcc.Graph(
                                        id="circular-tree-plot",
                                        style={
                                            "width": "100%",
                                            "display": "inline-block",
                                            "height": f"{self.height}px",
                                        },
                                        config={"scrollZoom": True},
                                    ),
                                ]
                            ),
                            width=6,
                        ),
                        dbc.Col(
                            html.Div(
                                [
                                    dcc.Graph(
                                        id="umap-scatter-plot",
                                        style={
                                            "width": "100%",
                                            "display": "inline-block",
                                            "height": f"{self.height-10}px",
                                        },
                                        # Zoom is very choppy with this enabled.  Users should use selection zoom
                                        config={"scrollZoom": False},
                                    ),
                                ]
                            ),
                            width=6,
                        ),
                    ]
                ),
                dbc.Offcanvas(
                    id="settings-pane", title="Settings", is_open=False, children=self.__render_closed_settings_pane()
                ),
                html.Div(id="init", style={"display": "none"}),
                html.Div(id="no-action", style={"display": "none"}),
            ],
        )

        return layout

    def __initialize_umap_scatter_plot(self) -> go.Figure:
        # calculate static bounds for the UMAP scatter plot
        self.umap_min_x, self.umap_max_x, self.umap_min_y, self.umap_max_y = self.__get_padded_umap_bounds(
            self.umap_padding
        )

        fig = go.Figure()

        color = self.umap_default_cell_color

        selected_cells = self.__get_effective_selected_cells()

        if len(selected_cells) > 0:
            color = [self.umap_inactive_cell_color] * self.adata.n_obs
            for i_obs in selected_cells:
                color[i_obs] = self.umap_active_cell_color

        fig.add_trace(
            go.Scatter(
                x=self.adata.obsm["X_umap"][:, 0],
                y=self.adata.obsm["X_umap"][:, 1],
                mode="markers",
                marker=dict(
                    color=color,
                    size=self.umap_marker_size,
                    opacity=self.umap_default_opacity,
                ),
                unselected=dict(marker=dict(opacity=self.umap_default_opacity)),
            )
        )
        self.__update_umap_scatter_plot_layout(fig)
        self._umap_scatter_plot_figure = fig
        return fig

    def __initialize_circular_tree_plot(self) -> go.Figure:
        self.circular_tree_plot = self.__instantiate_circular_tree_plot()
        fig = self.circular_tree_plot.plotly_figure
        return fig

    def __setup_initialization(self):
        @self.app.callback(Output("umap-scatter-plot", "figure"), Input("init", "children"))
        def __initialize_umap_scatter_plot(init):
            return self.__initialize_umap_scatter_plot()

        @self.app.callback(Output("circular-tree-plot", "figure"), Input("init", "children"))
        def __initialize_circular_tree_plot(init):
            return self.__initialize_circular_tree_plot()

    def __update_umap_scatter_plot_layout(self, umap_scatter_plot_fig):
        umap_scatter_plot_fig.update_layout(
            # uirevision is needed to maintain pan/zoom state.  It must be updated to trigger a refresh
            uirevision="true",
            plot_bgcolor="white",
            margin=dict(l=0, r=25, t=50, b=0),
            xaxis=dict(
                title="UMAP 1",
                showgrid=False,
                zeroline=False,  # Keep zero line enabled
                zerolinecolor="black",
                range=[self.umap_min_x, self.umap_max_x],  # Set x-axis limits
                showline=True,  # Show axis line
                linecolor="black",
                linewidth=1,
                tickmode="linear",
                tick0=-10,
                dtick=5,
            ),
            yaxis=dict(
                title="UMAP 2",
                showgrid=False,
                zeroline=False,  # Keep zero line enabled
                zerolinecolor="black",
                range=[self.umap_min_y, self.umap_max_y],  # Set y-axis limits
                showline=True,  # Show axis line
                linecolor="black",
                linewidth=1,
                tickmode="linear",
                tick0=-10,
                dtick=5,
            ),
            dragmode="pan",
        )

    def __render_breadcrumb(self) -> Component:
        selected_cells = self.__get_effective_selected_cells()
        if len(selected_cells) == 0 and self.selected_cell_domain_key.get() == DomainSelectionConstants.NONE:
            label = "Viewing results for all cells"
            show_clear = False
        elif (
            len(selected_cells) == 1 and self.selected_cell_domain_key.get() == DomainSelectionConstants.USER_SELECTION
        ):
            label = f"Selected cell index {selected_cells[0]}"
            show_clear = True
        elif len(selected_cells) > 1 and self.selected_cell_domain_key.get() == DomainSelectionConstants.USER_SELECTION:
            label = f"Selected {len(selected_cells)} cells"
            show_clear = True
        else:
            modifier = "cell" if len(selected_cells) == 1 else "cells"
            label = f"Selected cell domain {self.selected_cell_domain_key.get()} ({len(selected_cells)} {modifier})"
            show_clear = True
        children = [html.B(label, className="gr-breadcrumb-label")]

        if show_clear:
            children.append(
                html.Div(
                    [html.I(className="bi bi-x-circle")],
                    id="reset-selection-button",
                    n_clicks=0,
                    className="btn btn-link",
                    title="Clear selection",
                )
            )
        return html.Div(children)

    def __render_cell_selection_panes(self) -> Component:
        return [
            dbc.Col(self.__render_cell_selection_title(panel_id="cell-selection-title-tree"), width=6),
            dbc.Col(self.__render_cell_selection_title(panel_id="cell-selection-title-umap"), width=6),
        ]

    def __render_cell_selection_title(self, panel_id: str) -> Component:
        if self.selected_cl_name is not None:
            title = f"Selected cell class: {self.cl.cl_names_to_labels_map[self.selected_cl_name]}"
        else:
            title = DEFAULT_PANEL_TITLES[panel_id]

        return [html.Div(title, className="gr-header", id=panel_id)]

    def __render_closed_settings_pane(self) -> Component:
        return [
            html.Div(
                [
                    html.Label("Cell selection:", style={"margin-bottom": "5px"}),
                    self.__render_domain_dropdown(),
                ],
                className="gr-form-item",
            ),
            html.Div(
                [
                    dbc.Label("Evidence threshold:", html_for="evidence-threshold"),
                    dcc.Slider(
                        id="evidence-threshold",
                        min=0,
                        max=1,
                        value=self.score_threshold.get(dirty_read=True),
                        marks={
                            0: "0",
                            0.25: "0.25",
                            0.5: "0.5",
                            0.75: "0.75",
                            1: "1",
                        },
                        tooltip={"placement": "bottom", "always_visible": True, "style": {"margin": "0 5px"}},
                    ),
                ],
                className="gr-form-item",
            ),
            html.Div(
                [
                    dbc.Label("Minimum cell fraction:", html_for="cell-fraction"),
                    dcc.Slider(
                        id="cell-fraction",
                        min=0,
                        max=1,
                        value=self.min_cell_fraction.get(dirty_read=True),
                        marks={
                            0: "0",
                            0.25: "0.25",
                            0.5: "0.5",
                            0.75: "0.75",
                            1: "1",
                        },
                        tooltip={"placement": "bottom", "always_visible": True, "style": {"margin": "0 5px"}},
                    ),
                ],
                className="gr-form-item",
            ),
            html.Div(
                [
                    dbc.Button(
                        "Cancel",
                        id="cancel-button",
                        title="Cancel the changes and close the settings pane",
                        n_clicks=0,
                    ),
                    dbc.Button(
                        "Update",
                        id="update-button",
                        title="Update the graphs based on the specified configuration",
                        n_clicks=0,
                    ),
                ],
                className="gr-settings-button-bar",
            ),
            html.A(
                html.Img(
                    src="assets/cellarium-powered-400px.png",
                ),
                href="https://cellarium.ai",
                className="gr-powered-by",
                target="_blank",
            ),
        ]

    def __render_domain_dropdown(self) -> Component:
        labels = [{"label": "None selected", "value": DomainSelectionConstants.NONE}]
        if len(self.selected_cells) > 0:
            labels.append({"label": "User selection", "value": DomainSelectionConstants.USER_SELECTION})

        if len(self.cell_domain_map.keys()) > 1:
            labels.append({"label": "________________", "value": None, "disabled": True})
            labels.append({"label": html.Span("Provided domains"), "value": None, "disabled": True})

            for k in list(self.cell_domain_map.keys())[1:]:
                labels.append({"label": k, "value": k})

        return dcc.Dropdown(
            id="domain-dropdown",
            options=labels,
            value=self.selected_cell_domain_key.get(),  # default to no selection
            className="gr-custom-dropdown",
            clearable=False,
        )

    def __get_effective_selected_cells(self) -> list:
        # User has chosen not to show any highlighted cells
        if self.selected_cell_domain_key.get() == DomainSelectionConstants.NONE:
            return []

        # User has chose to highlight explicitly selected cells
        if self.selected_cell_domain_key.get() == DomainSelectionConstants.USER_SELECTION:
            return self.selected_cells

        # User has chose to highlight pre-calculated domain cells
        if self.selected_cell_domain_key.get() is not None:
            return self.cell_domain_map[self.selected_cell_domain_key.get()]

    def __clear_cell_selection(self):
        self.selected_cells = []
        self.selected_cell_domain_key.reset()

    def __clear_cell_class_selection(self):
        self.selected_cl_name = None
        self.circular_tree_plot.update_selected_nodes(selected_cl_path=[])

    def __setup_callbacks(self) -> None:
        # Cell selection callbacks
        @self.app.callback(
            Output("circular-tree-plot", "figure", allow_duplicate=True),
            Output("umap-scatter-plot", "figure", allow_duplicate=True),
            Output("panel-titles", "children", allow_duplicate=True),
            Input("circular-tree-plot", "clickData"),
            prevent_initial_call=True,
        )
        def __update_umap_scatter_plot_based_on_circular_tree_plot(clickData):
            if clickData is None or "points" not in clickData:
                return (
                    self.circular_tree_plot.plotly_figure,
                    self._umap_scatter_plot_figure,
                    self.__render_cell_selection_panes(),
                )

            point = clickData["points"][0]
            if "pointIndex" not in point:
                return (
                    self._umap_scatter_plot_figure,
                    self._umap_scatter_plot_figure,
                    self.__render_cell_selection_panes(),
                )

            node_index = point["pointIndex"]
            selected_cl_name = self.circular_tree_plot.clade_index_to_cl_name_map.get(node_index)
            if selected_cl_name is not None:
                # Toggle selection if the selected cell class was clicked
                if selected_cl_name == self.selected_cl_name:
                    self.__clear_cell_class_selection()
                else:
                    selected_cl_path = self.circular_tree_plot.get_clade_path_from_index(selected_cl_idx=node_index)
                    self.circular_tree_plot.update_selected_nodes(selected_cl_path)
                    self.selected_cl_name = selected_cl_name

            selected_cells_set = set(self.__get_effective_selected_cells())

            if self.selected_cl_name is None:
                # Optimizing for lower code complexity over performance by always treating color and opacity as arrays in this case
                color = [self.umap_default_cell_color] * self.adata.n_obs
                opacity = [self.umap_default_opacity] * self.adata.n_obs
            else:
                scores = self.__get_scores_for_cl_name(self.selected_cl_name)
                opacity = self.__get_scatter_plot_opacity_from_scores(scores)
                color = sample_colorscale(self.circular_tree_plot.score_colorscale, scores)

            selected_cells_set = set(self.__get_effective_selected_cells())
            # if no cells are selected but a cell class is, highlight all cells
            if len(selected_cells_set) == 0 and self.selected_cl_name is not None:
                selected_cells_set = set(self.cell_domain_map[self.ALL_CELLS_DOMAIN_KEY])

            for i_obs in range(self.adata.n_obs):
                if i_obs not in selected_cells_set:
                    color[i_obs] = self.umap_inactive_cell_color
                    opacity[i_obs] = self.umap_inactive_cell_opacity

            self._umap_scatter_plot_figure.update_traces(
                marker=dict(
                    color=color,
                    colorscale=self.circular_tree_plot.score_colorscale,
                    opacity=opacity,
                    cmin=0.0,
                    cmax=1.0,
                ),
                text=[f"{score:.5f}" for score in scores] if self.selected_cl_name is not None else None,
                hovertemplate=(
                    "<b>Evidence score: %{text}</b><extra></extra>" if self.selected_cl_name is not None else None
                ),
            )

            self._umap_scatter_plot_figure.update_layout(
                plot_bgcolor="white",
                margin=dict(l=0, r=25, t=50, b=0),
                # uirevision is needed to maintain pan/zoom state.  It must be updated to trigger a refresh
                uirevision=self._umap_scatter_plot_figure["layout"]["uirevision"],
            )
            return (
                self.circular_tree_plot.plotly_figure,
                self._umap_scatter_plot_figure,
                self.__render_cell_selection_panes(),
            )

        @self.app.callback(
            Output("circular-tree-plot", "figure", allow_duplicate=True),
            Output("umap-scatter-plot", "figure", allow_duplicate=True),
            Output("selected-domain-label", "children", allow_duplicate=True),
            Output("panel-titles", "children", allow_duplicate=True),
            Input("umap-scatter-plot", "clickData"),
            State("umap-scatter-plot", "selectedData"),
            prevent_initial_call=True,
        )
        def __update_circular_tree_plot_based_on_umap_scatter_plot(clickData, selectedData):
            self._umap_scatter_plot_figure.update_selections()
            if clickData is None or "points" not in clickData:
                return (
                    self.circular_tree_plot.plotly_figure,
                    self._umap_scatter_plot_figure,
                    self.__render_breadcrumb(),
                    self.__render_cell_selection_panes(),
                )

            point = clickData["points"][0]

            if "pointIndex" not in point:
                return (
                    self.circular_tree_plot.plotly_figure,
                    self._umap_scatter_plot_figure,
                    self.__render_breadcrumb(),
                    self.__render_cell_selection_panes(),
                )

            node_index = point["pointIndex"]
            self.selected_cells = [node_index]
            self.selected_cell_domain_key.set(DomainSelectionConstants.USER_SELECTION).commit()
            self.__clear_cell_class_selection()
            self.__initialize_circular_tree_plot()
            self.__initialize_umap_scatter_plot()

            return (
                self.circular_tree_plot.plotly_figure,
                self._umap_scatter_plot_figure,
                self.__render_breadcrumb(),
                self.__render_cell_selection_panes(),
            )

        @self.app.callback(
            Output("circular-tree-plot", "figure", allow_duplicate=True),
            Output("umap-scatter-plot", "figure", allow_duplicate=True),
            Output("umap-scatter-plot", "selectedData", allow_duplicate=True),
            Output("selected-domain-label", "children", allow_duplicate=True),
            Output("panel-titles", "children", allow_duplicate=True),
            Input("umap-scatter-plot", "selectedData"),
            prevent_initial_call=True,
        )
        def __update_circular_tree_plot_based_on_umap_scatter_plot_select(selectedData):
            # A selection event is firing on initialization. Ignore it by only accepting selectedData with a range field or lasso field
            if (
                selectedData is None
                or "points" not in selectedData
                or ("range" not in selectedData and "lassoPoints" not in selectedData)
            ):
                return (
                    self.circular_tree_plot.plotly_figure,
                    self._umap_scatter_plot_figure,
                    None,
                    self.__render_breadcrumb(),
                    self.__render_cell_selection_panes(),
                )

            points = selectedData["points"]

            node_indexes = [point["pointIndex"] for point in points]
            self.selected_cells = node_indexes
            self.selected_cell_domain_key.set(DomainSelectionConstants.USER_SELECTION).commit()
            self.__clear_cell_class_selection()
            self.__initialize_circular_tree_plot()
            self.__initialize_umap_scatter_plot()
            return (
                self.circular_tree_plot.plotly_figure,
                self._umap_scatter_plot_figure,
                selectedData,
                self.__render_breadcrumb(),
                self.__render_cell_selection_panes(),
            )

        @self.app.callback(
            Output("circular-tree-plot", "figure", allow_duplicate=True),
            Output("umap-scatter-plot", "figure", allow_duplicate=True),
            Output("selected-domain-label", "children", allow_duplicate=True),
            Output("settings-pane", "children", allow_duplicate=True),
            Output("panel-titles", "children", allow_duplicate=True),
            Input("reset-selection-button", "n_clicks"),
            prevent_initial_call=True,
        )
        def __reset_selection(n_clicks):
            if n_clicks != 0:
                self.__clear_cell_selection()
                self.__clear_cell_class_selection()

                # update the figures
                self.__initialize_circular_tree_plot()
                self.__initialize_umap_scatter_plot()

            return (
                self.circular_tree_plot.plotly_figure,
                self._umap_scatter_plot_figure,
                self.__render_breadcrumb(),
                self.__render_closed_settings_pane(),
                self.__render_cell_selection_panes(),
            )

        # Settings callbacks
        @self.app.callback(
            Output("circular-tree-plot", "figure", allow_duplicate=True),
            Output("umap-scatter-plot", "figure", allow_duplicate=True),
            Output("selected-domain-label", "children", allow_duplicate=True),
            Output("settings-pane", "children", allow_duplicate=True),
            Output("panel-titles", "children", allow_duplicate=True),
            Output("settings-pane", "is_open", allow_duplicate=True),
            Input("update-button", "n_clicks"),
            prevent_initial_call=True,
        )
        def __save_settings(n_clicks):
            if n_clicks > 0:
                # If a domain selection was changed and set to None, clear all selections
                if (
                    self.selected_cell_domain_key.is_dirty()
                    and self.selected_cell_domain_key.get(dirty_read=True) is DomainSelectionConstants.NONE
                ):
                    self.__clear_cell_selection()

                self.selected_cell_domain_key.commit()
                self.score_threshold.commit()
                self.min_cell_fraction.commit()

                self.__clear_cell_class_selection()

                # update the figures
                self.__initialize_circular_tree_plot()
                self.__initialize_umap_scatter_plot()

            return (
                self.circular_tree_plot.plotly_figure,
                self._umap_scatter_plot_figure,
                self.__render_breadcrumb(),
                self.__render_closed_settings_pane(),
                self.__render_cell_selection_panes(),
                False,
            )

        @self.app.callback(
            Output("settings-pane", "children", allow_duplicate=True),
            Output("settings-pane", "is_open", allow_duplicate=True),
            Input("cancel-button", "n_clicks"),
            prevent_initial_call=True,
        )
        def __cancel_settings(n_clicks):
            self.selected_cell_domain_key.rollback()
            self.score_threshold.rollback()
            self.min_cell_fraction.rollback()

            return self.__render_closed_settings_pane(), False

        @self.app.callback(
            Output("settings-pane", "is_open", allow_duplicate=True),
            Input("settings-button", "n_clicks"),
            [State("settings-pane", "is_open")],
            prevent_initial_call=True,
        )
        def __toggle_settings(n_clicks, is_open):
            if n_clicks:
                return not is_open
            return is_open

        @self.app.callback(
            Output("no-action", "children", allow_duplicate=True),
            Input("domain-dropdown", "value"),
            prevent_initial_call=True,
        )
        def __update_domain(domain):
            # set the domain
            self.selected_cell_domain_key.set(domain)

        @self.app.callback(
            Output("no-action", "children", allow_duplicate=True),
            Input("evidence-threshold", "value"),
            prevent_initial_call=True,
        )
        def __update_evidence_threshold(input_value):
            try:
                self.score_threshold.set(float(input_value))
            except ValueError:
                pass
            return input_value

        @self.app.callback(
            Output("no-action", "children", allow_duplicate=True),
            Input("cell-fraction", "value"),
            prevent_initial_call=True,
        )
        def __update_cell_fraction(input_value):
            try:
                self.min_cell_fraction.set(float(input_value))
            except ValueError:
                pass
            return input_value
