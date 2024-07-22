import typing as t
from functools import lru_cache

import numpy as np
import plotly.graph_objs as go
from Bio import Phylo

from cellarium.cas.postprocessing.ontology_aware import (
    CAS_CL_LABEL_PROPERTY_REF,
    CAS_FRACTION_PROPERTY_REF,
    CAS_SCORE_PROPERTY_REF,
)

DEFAULT_SCORE_COLORSCALE = [
    [0.0, "rgb(214, 47, 38)"],
    [0.1, "rgb(214, 47, 38)"],
    [0.2, "rgb(244, 109, 67)"],
    [0.3, "rgb(252, 172, 96)"],
    [0.4, "rgb(254, 224, 139)"],
    [0.5, "rgb(254, 254, 189)"],
    [0.6, "rgb(217, 239, 139)"],
    [0.7, "rgb(164, 216, 105)"],
    [0.8, "rgb(102, 189, 99)"],
    [0.9, "rgb(63, 170, 89)"],
    [1.0, "rgb(25, 151, 79)"],
]


class CircularTreePlot:
    def __init__(
        self,
        tree: Phylo.BaseTree.Tree,
        title: t.Optional[str] = None,
        order: str = "preorder",
        dist: float = 1.0,
        start_angle: float = 0.0,
        end_angle: float = 360.0,
        start_leaf: str = "first",
        min_node_radius: float = 6.0,
        max_node_radius: float = 16.0,
        linecolor: str = "rgb(20,20,20)",
        score_colorscale: list = DEFAULT_SCORE_COLORSCALE,
        shown_cl_names_set: set[str] = set(),
    ):
        """
        Initialize the CircularTreePlot class.

        :param tree: An instance of Bio.Phylo.Newick.Tree or Bio.Phylo.PhyloXML.Phylogeny.
        :param order: Tree traversal method to associate polar coordinates to its nodes.
        :param dist: Vertical distance between two consecutive leaves in the rectangular tree layout.
        :param start_angle: Angle in degrees representing the angle of the first leaf mapped to a circle.
        :param end_angle: Angle in degrees representing the angle of the last leaf.
        :param start_leaf: Keyword with two possible values:
            'first' - to map the leaves in the list tree.get_terminals() onto a circle in the counter-clockwise direction.
            'last' - to map the leaves in the reversed list tree.get_terminals()[::-1].
        :param min_node_radius: Minimum radius of a node (the limit of 0. fraction of cells having that node in CAS response)
        :param max_node_radius: Minimum radius of a node (the limit of 1. fraction of cells having that node in CAS response)
        :param score_colorscale: Colorscale for scores
        """
        self.tree = tree
        self.title = title
        self.order = order
        self.dist = dist
        self.start_angle = np.deg2rad(start_angle)
        self.end_angle = np.deg2rad(end_angle)
        self.start_leaf = start_leaf
        self.min_node_radius = min_node_radius
        self.max_node_radius = max_node_radius
        self.linecolor = linecolor
        self.score_colorscale = score_colorscale
        self.shown_cl_names_set = shown_cl_names_set

        # generate circular tree coordinate data
        (
            self.x_nodes,
            self.y_nodes,
            self.x_lines,
            self.y_lines,
            self.x_arcs,
            self.y_arcs,
        ) = self._get_circular_tree_data()

        # generate assets for visualization
        self.tooltip_string_list = []
        self.score_list = []
        self.cl_name_list = []
        self.cl_label_list = []
        self.size_list = []
        self.clade_index_to_cl_name_map = dict()
        self.clade_to_index = dict()
        self.cl_name_to_cl_label_map = dict()
        for idx, clade in enumerate(self.tree.find_clades(order="preorder")):
            cl_name = clade.name
            self.clade_index_to_cl_name_map[idx] = cl_name
            self.clade_to_index[clade] = idx
            parsed_property_dict = dict()
            for property in clade.properties:
                if property.ref == CAS_SCORE_PROPERTY_REF:
                    parsed_property_dict["score"] = float(property.value)
                elif property.ref == CAS_FRACTION_PROPERTY_REF:
                    parsed_property_dict["fraction"] = float(property.value)
                elif property.ref == CAS_CL_LABEL_PROPERTY_REF:
                    parsed_property_dict["cl_label"] = property.value
                else:
                    raise ValueError
            assert "score" in parsed_property_dict
            assert "fraction" in parsed_property_dict
            assert "cl_label" in parsed_property_dict
            cl_label = parsed_property_dict["cl_label"]
            score = parsed_property_dict["score"]
            fraction = parsed_property_dict["fraction"]
            tooltip_string = f"ID: {cl_name}<br>Label: {cl_label}<br>Evidence score: {score:.5f}<br>Cell fraction: {fraction * 100:.3f}%"
            size = self.min_node_radius + (self.max_node_radius - self.min_node_radius) * fraction
            self.tooltip_string_list.append(tooltip_string)
            self.score_list.append(score)
            self.size_list.append(size)
            self.cl_name_list.append(cl_name)
            self.cl_label_list.append(cl_label)
            self.cl_name_to_cl_label_map[cl_name] = cl_label

        self.num_nodes = len(self.clade_to_index.keys())
        # text labels
        self.text_labels = [
            self.cl_name_to_cl_label_map[cl_name] if cl_name in shown_cl_names_set else ""
            for cl_name in self.cl_name_list
        ]

    def _get_radius(self) -> t.Dict:
        """
        Associates each clade root with its radius, equal to the distance from that clade to the tree root.

        :return: Dictionary {clade: node_radius}
        """
        node_radius = self.tree.depths()

        if not np.count_nonzero(list(node_radius.values())):
            node_radius = self.tree.depths(unit_branch_lengths=True)
        return node_radius

    def _get_vertical_position(self) -> t.Dict:
        """
        Returns a dictionary {clade: y_coord}, where y_coord is the Cartesian y-coordinate of a clade root in a rectangular phylogram.

        :return: Dictionary {clade: y_coord}
        """
        if self.start_leaf == "first":
            node_ycoord = {leaf: k for k, leaf in enumerate(self.tree.get_terminals())}
        elif self.start_leaf == "last":
            node_ycoord = {leaf: k for k, leaf in enumerate(reversed(self.tree.get_terminals()))}
        else:
            raise ValueError("start_leaf can be only 'first' or 'last'")

        def assign_ycoord(clade):
            for subclade in clade:
                if subclade not in node_ycoord:
                    assign_ycoord(subclade)
            node_ycoord[clade] = 0.5 * (node_ycoord[clade.clades[0]] + node_ycoord[clade.clades[-1]])

        if self.tree.root.clades:
            assign_ycoord(self.tree.root)
        return node_ycoord

    def _ycoord_to_theta(self, y: float, ymin: float, ymax: float) -> float:
        """
        Maps a y-coordinate to an angle in radians.

        :param y: y-coordinate.
        :param ymin: Minimum y-coordinate.
        :param ymax: Maximum y-coordinate.
        :return: Corresponding angle in radians.
        """
        return self.start_angle + (self.end_angle - self.start_angle) * (y - ymin) / float(ymax - ymin)

    def _get_points_on_lines(
        self,
        line_type: str,
        x_left: float = 0,
        x_right: float = 0,
        y_right: float = 0,
        y_bottom: float = 0,
        y_top: float = 0,
    ) -> t.Tuple[t.List[float], t.List[float]]:
        """
        Define the points that generate a radial branch and the circular arcs, perpendicular to that branch.

        :param line_type: Type of line ('radial' or 'angular').
        :param x_left: Left x-coordinate.
        :param x_right: Right x-coordinate.
        :param y_right: Right y-coordinate.
        :param y_bottom: Bottom y-coordinate.
        :param y_top: Top y-coordinate.
        :return: Tuple containing lists of x and y coordinates of the line representative points.
        """
        if line_type == "radial":
            theta = self._ycoord_to_theta(y_right, self.ymin, self.ymax)
            X = [x_left * np.cos(theta), x_right * np.cos(theta), None]
            Y = [x_left * np.sin(theta), x_right * np.sin(theta), None]
        elif line_type == "angular":
            theta_bottom = self._ycoord_to_theta(y_bottom, self.ymin, self.ymax)
            theta_top = self._ycoord_to_theta(y_top, self.ymin, self.ymax)
            t = np.linspace(0, 1, 10)
            theta = (1 - t) * theta_bottom + t * theta_top
            X = list(x_right * np.cos(theta)) + [None]
            Y = list(x_right * np.sin(theta)) + [None]
        else:
            raise ValueError("line_type can be only 'radial' or 'angular'")

        return X, Y

    def _get_line_lists(
        self,
        clade: Phylo.BaseTree.Clade,
        x_left: float,
        x_lines: t.List[float],
        y_lines: t.List[float],
        x_arcs: t.List[float],
        y_arcs: t.List[float],
    ):
        """
        Recursively compute the lists of points that span the tree branches.

        :param clade: Clade of the tree.
        :param x_left: Left x-coordinate.
        :param x_lines: List of x-coordinates of radial edge ends.
        :param y_lines: List of y-coordinates of radial edge ends.
        :param x_arcs: List of x-coordinates of arc segments for tree branches.
        :param y_arcs: List of y-coordinates of arc segments for tree branches.
        """
        x_right = self.node_radius[clade]
        y_right = self.node_ycoord[clade]

        X, Y = self._get_points_on_lines(line_type="radial", x_left=x_left, x_right=x_right, y_right=y_right)
        x_lines.extend(X)
        y_lines.extend(Y)

        if clade.clades:
            y_top = self.node_ycoord[clade.clades[0]]
            y_bottom = self.node_ycoord[clade.clades[-1]]

            X, Y = self._get_points_on_lines(line_type="angular", x_right=x_right, y_bottom=y_bottom, y_top=y_top)
            x_arcs.extend(X)
            y_arcs.extend(Y)

            for child in clade:
                self._get_line_lists(child, x_right, x_lines, y_lines, x_arcs, y_arcs)

    def _get_circular_tree_data(
        self,
    ) -> t.Tuple[t.List[float], t.List[float], t.List[float], t.List[float], t.List[float], t.List[float]]:
        """
        Define data needed to get the Plotly plot of a circular tree.

        :return: Tuple containing xnodes, ynodes, xlines, ylines, xarc, yarc.
        """
        self.node_radius = self._get_radius()
        self.node_ycoord = self._get_vertical_position()
        y_vals = self.node_ycoord.values()
        self.ymin, self.ymax = min(y_vals), max(y_vals)
        self.ymin -= self.dist

        x_lines = []
        y_lines = []
        x_arcs = []
        y_arcs = []
        self._get_line_lists(self.tree.root, 0, x_lines, y_lines, x_arcs, y_arcs)

        x_nodes = []
        y_nodes = []

        for clade in self.tree.find_clades(order=self.order):
            theta = self._ycoord_to_theta(self.node_ycoord[clade], self.ymin, self.ymax)
            x_nodes.append(self.node_radius[clade] * np.cos(theta))
            y_nodes.append(self.node_radius[clade] * np.sin(theta))

        return x_nodes, y_nodes, x_lines, y_lines, x_arcs, y_arcs

    @property
    @lru_cache(maxsize=None)
    def plotly_figure(self) -> go.Figure:
        # Plot points
        trace_nodes = dict(
            type="scatter",
            name="data_nodes",
            x=self.x_nodes,
            y=self.y_nodes,
            mode="markers+text",
            text=self.text_labels,
            textposition="top center",
            textfont=dict(size=10),  # Base font size
            marker=dict(
                color=self.score_list,
                size=self.size_list,
                colorscale=self.score_colorscale,
                colorbar=dict(
                    thickness=10,
                    dtick=10,
                    ticklen=4,
                    title=dict(
                        text="Evidence score",
                        side="top",
                        font=dict(family="Arial", size=14),
                    ),
                    tickvals=[0, 0.25, 0.5, 0.75, 1],
                    ticktext=["0", "0.25", "0.5", "0.75", "1"],
                    tickmode="array",
                    lenmode="fraction",
                    len=0.4,
                    orientation="h",  # Horizontal orientation
                    xanchor="right",
                    yanchor="bottom",
                    x=1.0,
                    y=-0.1,
                    tickfont=dict(family="Arial", size=12),
                ),
                cmin=0.0,
                cmax=1.0,
                opacity=1.0,
            ),
            hovertext=self.tooltip_string_list,
            hoverinfo="text",
            opacity=1,
            # uirevision is needed to maintain pan/zoom state.  It must be updated to trigger a refresh
            uirevision="true",
        )

        # Straight lines
        trace_radial_lines = dict(
            type="scatter",
            name="radial_lines",
            x=self.x_lines,
            y=self.y_lines,
            mode="lines",
            line=dict(color=self.linecolor, width=1),
            hoverinfo="none",
        )

        # Curved lines
        trace_arcs = dict(
            type="scatter",
            name="arc_lines",
            x=self.x_arcs,
            y=self.y_arcs,
            mode="lines",
            line=dict(color=self.linecolor, width=1, shape="spline"),
            hoverinfo="none",
        )

        title = dict(
            text=self.title,
            font=dict(family="Arial", size=1) if self.title is None else dict(family="Arial", size=20),
        )
        layout = dict(
            title=title,
            font=dict(family="Arial", size=14),
            autosize=True,
            showlegend=False,
            xaxis=dict(
                visible=False,
                scaleanchor="y",
                scaleratio=1,
            ),
            yaxis=dict(visible=False),
            hovermode="closest",
            dragmode="pan",
            plot_bgcolor="rgb(255,255,255)",
            margin=dict(l=0, r=25, t=50, b=0),
        )

        fig = go.Figure(data=[trace_radial_lines, trace_arcs, trace_nodes], layout=layout)

        return fig

    def get_clade_path_from_index(self, selected_cl_idx: int) -> t.Optional[t.List[str]]:
        """
        Returns the path from the root to the selected clade.  Returns None if the index isn't found.
        """
        node = list(self.tree.find_clades(order="preorder"))[selected_cl_idx]
        if node is None:
            return None

        return [self.tree.root.name] + [cl.name for cl in self.tree.root.get_path(target=node)]

    def update_selected_nodes(self, selected_cl_path: t.List[str] = []):
        if len(selected_cl_path) == 0:
            line_widths = 0
        else:
            # Find the starting node to start navigating down the tree
            node = self.tree.find_clades(name=selected_cl_path[0], order="preorder")
            selected_nodes = []
            if node is not None:
                node = next(node)
                selected_nodes.append(self.clade_to_index[node])

            for node_name in selected_cl_path[1:]:
                node = next(filter(lambda cl: cl.name == node_name, node.clades))
                if node is not None:
                    selected_nodes.append(self.clade_to_index[node])
                else:
                    break
            last_index = selected_nodes.pop()

            line_widths = [0] * self.num_nodes
            for idx in selected_nodes:
                line_widths[idx] = 1.5
            # The selected note should have a thicker line
            line_widths[last_index] = 3

        self.plotly_figure.update_traces(
            marker=dict(line=dict(color="DarkSlateGrey", width=line_widths)),
            selector=({"name": "data_nodes"}),
            uirevision=self.plotly_figure["layout"]["uirevision"],
        )
