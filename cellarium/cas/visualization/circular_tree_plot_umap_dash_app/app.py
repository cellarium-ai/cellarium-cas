import os
import tempfile
import numpy as np
from collections import OrderedDict
from typing import Sequence, Tuple

from Bio import Phylo
from anndata import AnnData
from dash import Dash, dcc, html
from dash.dependencies import Input, Output
from plotly.express.colors import sample_colorscale
import plotly.graph_objects as go

from logging import log, INFO

from cellarium.cas.postprocessing import (
    get_obs_indices_for_cluster,
    get_aggregated_cas_ontology_aware_scores,
    generate_phyloxml_from_scored_cell_ontology_tree,
    convert_aggregated_cell_ontology_scores_to_rooted_tree,
    insert_cas_ontology_aware_response_into_adata,
    CellOntologyScoresAggregationDomain,
    CellOntologyScoresAggregationOp,
)

from cellarium.cas.postprocessing import CAS_CL_SCORES_ANNDATA_OBSM_KEY
from cellarium.cas.postprocessing.cell_ontology import CellOntologyCache, CL_CELL_ROOT_NODE
from cellarium.cas.visualization._components.circular_tree_plot import CircularTreePlot

# cell type ontology terms (and all descendents) to hide from the visualization
DEFAULT_HIDDEN_CL_NAMES_SET = {}


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
}


class CASCircularTreePlotUMAPDashApp:
    ALL_CELLS_DOMAIN_KEY = "all cells"
    CLUSTER_PREFIX_DOMAIN_KEY = "cluster "

    def __init__(
        self,
        adata: AnnData,
        cas_ontology_aware_response: list,
        cluster_label_obs_column: str | None = None,
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
        hidden_cl_names_set: set[str] = DEFAULT_HIDDEN_CL_NAMES_SET,
        shown_cl_names_set: set[str] = DEFAULT_SHOWN_CL_NAMES_SET,
        score_colorscale: str | list = "Viridis",
    ):
        self.adata = adata
        self.aggregation_op = aggregation_op
        self.aggregation_domain = aggregation_domain
        self.score_threshold = score_threshold
        self.min_cell_fraction = min_cell_fraction
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
        self.hidden_cl_names_set = hidden_cl_names_set
        self.shown_cl_names_set = shown_cl_names_set
        self.score_colorscale = score_colorscale

        assert "X_umap" in adata.obsm, "UMAP coordinates not found in adata.obsm['X_umap']"

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
        self.selected_cell_domain_key = self.ALL_CELLS_DOMAIN_KEY

        # instantiate the cell type ontology cache
        self.cl = CellOntologyCache()

        # insert CA ontology-aware response into adata
        insert_cas_ontology_aware_response_into_adata(cas_ontology_aware_response, adata, self.cl)

        # instantiate the Dash app
        self.app = Dash(__name__)
        self.server = self.app.server
        self.app.layout = self._create_layout()
        self._setup_initialization()
        self._setup_callbacks()

    def run(self, **kwargs):
        log(INFO, "Starting Dash application...")
        self.app.run_server(jupyter_mode="inline", jupyter_width="100%", jupyter_height=self.height + 50, **kwargs)

    def _instantiate_circular_tree_plot(
        self, obs_indices_override: Sequence | None = None, title_override: str | None = None
    ) -> CircularTreePlot:
        # reduce scores over the provided cells
        aggregated_scores = get_aggregated_cas_ontology_aware_scores(
            self.adata,
            obs_indices=self.cell_domain_map[self.selected_cell_domain_key]
            if obs_indices_override is None
            else obs_indices_override,
            aggregation_op=self.aggregation_op,
            aggregation_domain=self.aggregation_domain,
            threshold=self.score_threshold,
        )

        # generate a Phylo tree
        rooted_tree = convert_aggregated_cell_ontology_scores_to_rooted_tree(
            aggregated_scores=aggregated_scores,
            cl=self.cl,
            root_cl_name=CL_CELL_ROOT_NODE,
            min_fraction=self.min_cell_fraction,
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

        domain_label = self.selected_cell_domain_key if title_override is None else title_override
        return CircularTreePlot(
            tree=phyloxml_tree,
            title=f"Score summary of {domain_label}",
            score_colorscale=self.score_colorscale,
            linecolor=self.circular_tree_plot_linecolor,
            start_angle=self.circular_tree_start_angle,
            end_angle=self.circular_tree_end_angle,
            shown_cl_names_set=self.shown_cl_names_set,
        )

    def _get_padded_umap_bounds(self, umap_padding: float) -> Tuple[float, float, float, float]:
        actual_min_x = np.min(self.adata.obsm["X_umap"][:, 0])
        actual_max_x = np.max(self.adata.obsm["X_umap"][:, 0])
        actual_min_y = np.min(self.adata.obsm["X_umap"][:, 1])
        actual_max_y = np.max(self.adata.obsm["X_umap"][:, 1])
        padded_min_x = actual_min_x - umap_padding * (actual_max_x - actual_min_x)
        padded_max_x = actual_max_x + umap_padding * (actual_max_x - actual_min_x)
        padded_min_y = actual_min_y - umap_padding * (actual_max_y - actual_min_y)
        padded_max_y = actual_max_y + umap_padding * (actual_max_y - actual_min_y)

        return padded_min_x, padded_max_x, padded_min_y, padded_max_y

    def _get_scores_for_cl_name(self, cl_name: str) -> np.ndarray:
        cl_index = self.cl.cl_names_to_idx_map[cl_name]
        return self.adata.obsm[CAS_CL_SCORES_ANNDATA_OBSM_KEY][:, cl_index].toarray().flatten()

    def _get_scatter_plot_opacity_from_scores(self, scores: np.ndarray) -> np.ndarray:
        min_score = np.min(scores)
        max_score = np.max(scores)
        normalized_scores = (scores - min_score) / (1e-6 + max_score - min_score)
        return np.maximum(
            scores, self.umap_min_opacity + (self.umap_max_opacity - self.umap_min_opacity) * normalized_scores
        )

    def _create_layout(self):
        # Custom JavaScript for increasing scroll zoom sensitivity (doesn't seem to work)
        scroll_zoom_js = """
        function(graph) {
            var plot = document.getElementById(graph.id);
            plot.on('wheel', function(event) {
                event.deltaY *= 0.2;  // Adjust this value to change sensitivity
            });
            return graph;
        }
        """

        layout = html.Div(
            id="main-container",
            children=[
                html.Div(
                    [
                        html.H3("Control Panel", style={"margin-top": "0px", "margin-bottom": "50px"}),
                        html.Div(
                            [
                                html.Label("Cell selection:", style={"margin-bottom": "5px"}),
                                dcc.Dropdown(
                                    id="domain-dropdown",
                                    options=list(dict(label=k, value=k) for k in self.cell_domain_map.keys()),
                                    value=self.ALL_CELLS_DOMAIN_KEY,  # default to all cells
                                    className="custom-dropdown",
                                ),
                            ],
                            style={"margin-bottom": "20px"},
                        ),
                        html.Div(
                            [
                                html.Label("Evidence threshold:", style={"margin-bottom": "5px"}),
                                dcc.Input(id="evidence-threshold", type="text", value=self.score_threshold),
                            ],
                            style={"margin-bottom": "20px"},
                        ),
                        html.Div(
                            [
                                html.Label("Minimum cell fraction:"),
                                dcc.Input(id="cell-fraction", type="text", value=self.min_cell_fraction),
                            ],
                            style={"margin-bottom": "20px"},
                        ),
                        html.Div(
                            [
                                html.Button(
                                    "Update",
                                    id="update-button",
                                    n_clicks=0,
                                    style={"margin-top": "20px", "margin-left": "0px"},
                                )
                            ],
                            style={"margin-top": "20px"},
                        ),
                        html.Img(
                            src="assets/cellarium-powered-400px.png",
                            style={"width": "150px", "position": "absolute", "bottom": "20px"},
                        ),
                    ],
                    style={
                        "width": "20%",
                        "backgroundColor": "white",
                        "padding": "20px",
                        "boxSizing": "border-box",
                        "height": "100vh",
                    },
                ),
                html.Div(
                    [
                        dcc.Graph(
                            id="circular-tree-plot",
                            style={"width": "100%", "display": "inline-block", "height": f"{self.height}x"},
                            config={"scrollZoom": True},
                        )
                    ],
                    style={
                        "width": "40%",
                        "display": "inline-block",
                        "backgroundColor": "white",
                        "padding": "0px",
                        "boxSizing": "border-box",
                        "height": "100vh",
                        "float": "left",
                    },
                ),
                html.Div(
                    [
                        dcc.Graph(
                            id="umap-scatter-plot",
                            style={"width": "100%", "display": "inline-block", "height": f"{self.height}px"},
                            config={"scrollZoom": True},
                        )
                    ],
                    style={
                        "width": "40%",
                        "backgroundColor": "white",
                        "padding": "0px",
                        "boxSizing": "border-box",
                        "height": "100vh",
                    },
                ),
                html.Div(id="init", style={"display": "none"}),
                html.Div(id="no-action", style={"display": "none"}),
                html.Script(scroll_zoom_js),
            ],
            style={
                "display": "flex",
                "flexDirection": "row",
                "backgroundColor": "white",
                "margin": "0",
                "padding": "0",
                "height": "100vh",
            },
        )

        return layout

    def _initialize_umap_scatter_plot(self, highlight_active: bool = False) -> go.Figure:
        # calculate static bounds for the UMAP scatter plot
        self.umap_min_x, self.umap_max_x, self.umap_min_y, self.umap_max_y = self._get_padded_umap_bounds(
            self.umap_padding
        )

        fig = go.Figure()

        color = self.umap_default_cell_color
        if highlight_active:
            if self.selected_cell_domain_key != self.ALL_CELLS_DOMAIN_KEY:
                color = [self.umap_inactive_cell_color] * self.adata.n_obs
                for i_obs in self.cell_domain_map[self.selected_cell_domain_key]:
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
            )
        )
        self._update_umap_scatter_plot_layout(fig)
        self._umap_scatter_plot_figure = fig
        return fig

    def _initialize_circular_tree_plot(
        self, obs_indices_override: Sequence | None = None, title_override: str | None = None
    ) -> go.Figure:
        self.circular_tree_plot = self._instantiate_circular_tree_plot(obs_indices_override, title_override)
        fig = self.circular_tree_plot.plotly_figure
        self._circular_tree_plot_figure = fig
        return fig

    def _setup_initialization(self):
        @self.app.callback(Output("umap-scatter-plot", "figure"), Input("init", "children"))
        def _initialize_umap_scatter_plot(init):
            return self._initialize_umap_scatter_plot()

        @self.app.callback(Output("circular-tree-plot", "figure"), Input("init", "children"))
        def _initialize_circular_tree_plot(init):
            return self._initialize_circular_tree_plot()

    def _update_umap_scatter_plot_layout(self, umap_scatter_plot_fig):
        umap_scatter_plot_fig.update_layout(
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

    def _setup_callbacks(self):
        @self.app.callback(
            Output("umap-scatter-plot", "figure", allow_duplicate=True),
            Input("circular-tree-plot", "clickData"),
            prevent_initial_call=True,
        )
        def _update_umap_scatter_plot_based_on_circular_tree_plot(clickData):
            if clickData is None or "points" not in clickData:
                return self._umap_scatter_plot_figure

            point = clickData["points"][0]
            if "pointIndex" not in point:
                return self._umap_scatter_plot_figure

            node_index = point["pointIndex"]
            cl_name = self.circular_tree_plot.clade_index_to_cl_name_map.get(node_index)
            if cl_name is None:
                return self._umap_scatter_plot_figure

            scores = self._get_scores_for_cl_name(cl_name)
            opacity = self._get_scatter_plot_opacity_from_scores(scores)
            color = sample_colorscale(self.circular_tree_plot.score_colorscale, scores)
            selected_cells_set = set(self.cell_domain_map[self.selected_cell_domain_key])
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
                text=[f"{score:.5f}" for score in scores],
                hovertemplate="<b>Evidence score: %{text}</b><extra></extra>",
            )
            self._umap_scatter_plot_figure.update_layout(title=self.cl.cl_names_to_labels_map[cl_name])

            return self._umap_scatter_plot_figure

        @self.app.callback(
            Output("circular-tree-plot", "figure", allow_duplicate=True),
            Output("umap-scatter-plot", "figure", allow_duplicate=True),
            Input("domain-dropdown", "value"),
            prevent_initial_call=True,
        )
        def _update_figures_based_on_dropdown(domain):
            # set the domain
            self.selected_cell_domain_key = domain

            # update the figures
            self._initialize_umap_scatter_plot(highlight_active=True)
            self._initialize_circular_tree_plot()
            return self._circular_tree_plot_figure, self._umap_scatter_plot_figure

        @self.app.callback(
            Output("circular-tree-plot", "figure", allow_duplicate=True),
            Input("umap-scatter-plot", "clickData"),
            prevent_initial_call=True,
        )
        def _update_circular_tree_plot_based_on_umap_scatter_plot(clickData):
            if clickData is None or "points" not in clickData:
                return self._circular_tree_plot_figure

            point = clickData["points"][0]

            if "pointIndex" not in point:
                return self._circular_tree_plot_figure

            node_index = point["pointIndex"]
            self._initialize_circular_tree_plot([node_index], f"cell index {node_index}")
            return self._circular_tree_plot_figure

        @self.app.callback(
            Output("circular-tree-plot", "figure", allow_duplicate=True),
            Input("update-button", "n_clicks"),
            prevent_initial_call=True,
        )
        def _update_circular_tree_plot_based_on_update_button(n_clicks):
            self._initialize_circular_tree_plot()
            return self._circular_tree_plot_figure

        @self.app.callback(
            Output("no-action", "children", allow_duplicate=True),
            Input("evidence-threshold", "value"),
            prevent_initial_call=True,
        )
        def _update_evidence_threshold(input_value):
            try:
                self.score_threshold = float(input_value)
            except ValueError:
                pass
            return input_value

        @self.app.callback(
            Output("no-action", "children", allow_duplicate=True),
            Input("cell-fraction", "value"),
            prevent_initial_call=True,
        )
        def _update_cell_fraction(input_value):
            try:
                self.min_cell_fraction = float(input_value)
            except ValueError:
                pass
            return input_value
