import typing as t

# the 'cell' node
CL_CELL_ROOT_NODE = "CL:0000000"

# the 'eukaryotic cell' node
CL_EUKARYOTIC_CELL_ROOT_NODE = "CL:0000255"


class CellOntologyCache:
    """
    A class representing a cache for the Cell Ontology (CL), populated from a precomputed resource dict.

    Attributes:
        cl_names (list): Ordered list of Cell Ontology term IDs.
        cl_labels (list): Ordered list of Cell Ontology term labels (same order as cl_names).
        cl_names_to_labels_map (dict): Mapping from term ID to human-readable label.
        cl_names_to_idx_map (dict): Mapping from term ID to integer column index.
        children_dict (dict): Mapping from term ID to list of direct children term IDs.
    """

    def __init__(self, resource: t.Dict[str, t.Any]) -> None:
        """
        Initialize from a precomputed cell ontology resource dict (as returned by the CAS API).

        :param resource: Dict with keys: cl_names, cell_ontology_term_id_to_cell_type,
            children_dictionary, shortest_path_lengths_from_cell_root,
            longest_path_lengths_from_cell_root.
        """
        cl_names = resource["cl_names"]
        cl_names_to_labels_map = resource["cell_ontology_term_id_to_cell_type"]

        self.cl_names: t.List[str] = cl_names
        self.cl_labels: t.List[str] = [cl_names_to_labels_map[n] for n in cl_names]
        self.cl_names_to_labels_map: t.Dict[str, str] = cl_names_to_labels_map
        self.cl_names_to_idx_map: t.Dict[str, int] = {name: idx for idx, name in enumerate(cl_names)}
        self.children_dict: t.Dict[str, t.List[str]] = resource["children_dictionary"]
        self._shortest_path_lengths: t.Dict[str, int] = resource["shortest_path_lengths_from_cell_root"]
        self._longest_path_lengths: t.Dict[str, int] = resource["longest_path_lengths_from_cell_root"]

    def get_shortest_path_lengths_from_target(self, target: str) -> t.Dict[str, int]:
        if target != CL_CELL_ROOT_NODE:
            raise ValueError(
                f"Only '{CL_CELL_ROOT_NODE}' is supported as a target for precomputed path lengths, got '{target}'"
            )
        return self._shortest_path_lengths

    def get_longest_path_lengths_from_target(self, target: str) -> t.Dict[str, int]:
        if target != CL_CELL_ROOT_NODE:
            raise ValueError(
                f"Only '{CL_CELL_ROOT_NODE}' is supported as a target for precomputed path lengths, got '{target}'"
            )
        return self._longest_path_lengths
