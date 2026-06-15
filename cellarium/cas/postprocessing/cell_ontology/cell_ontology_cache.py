import typing as t

# the 'cell' node
CL_CELL_ROOT_NODE = "CL:0000000"

# the 'eukaryotic cell' node
CL_EUKARYOTIC_CELL_ROOT_NODE = "CL:0000255"


def _bfs_reachable(start: str, adjacency: t.Dict[str, t.List[str]]) -> t.FrozenSet[str]:
    """Return *start* plus all nodes reachable from it via *adjacency* (BFS, self-inclusive)."""
    visited: t.Set[str] = {start}
    queue = [start]
    while queue:
        node = queue.pop(0)
        for neighbor in adjacency.get(node, ()):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    return frozenset(visited)


class CellOntologyCache:
    """
    A class representing a cache for the Cell Ontology (CL), populated from a precomputed resource dict.

    Attributes:
        cl_names (list): Ordered list of Cell Ontology term IDs.
        cl_labels (list): Ordered list of Cell Ontology term labels (same order as cl_names).
        cl_names_to_labels_map (dict): Mapping from term ID to human-readable label.
        cl_names_to_idx_map (dict): Mapping from term ID to integer column index.
        children_dict (dict): Mapping from term ID to list of direct children term IDs.
        ancestors_dict (dict): Mapping from term ID to frozen set of ancestors (self-inclusive, root-exclusive).
        descendants_dict (dict): Mapping from term ID to frozen set of descendants (self-inclusive).
    """

    def __init__(self, resource: t.Dict[str, t.Any]) -> None:
        """
        Initialize from a precomputed cell ontology resource dict (as returned by the CAS API).

        :param resource: Dict with keys: cl_names, cell_ontology_term_id_to_cell_type,
            children_dictionary, ancestors_dictionary (optional), descendants_dictionary (optional),
            shortest_path_lengths_from_cell_root, longest_path_lengths_from_cell_root.
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

        # Build parents map for ancestor traversal
        parents_map: t.Dict[str, t.List[str]] = {}
        for parent, children in self.children_dict.items():
            for child in children:
                parents_map.setdefault(child, []).append(parent)

        raw_ancestors = resource.get("ancestors_dictionary")
        if raw_ancestors:
            # CAS API does not include self in ancestors_dictionary; add it for consistency.
            self.ancestors_dict: t.Dict[str, t.FrozenSet[str]] = {
                term: frozenset(ancs) | {term} for term, ancs in raw_ancestors.items()
            }
        else:
            self.ancestors_dict = {term: _bfs_reachable(term, parents_map) for term in cl_names}

        raw_descendants = resource.get("descendants_dictionary")
        if raw_descendants:
            self.descendants_dict: t.Dict[str, t.FrozenSet[str]] = {
                term: frozenset(descs) | {term} for term, descs in raw_descendants.items()
            }
        else:
            self.descendants_dict = {term: _bfs_reachable(term, self.children_dict) for term in cl_names}

    def get_ancestors(self, term: str, remove_root: bool = False, remove_self: bool = False) -> t.FrozenSet[str]:
        """
        Return the ancestor set for *term* (self-inclusive, root-inclusive by default).

        Falls back to a singleton ``{term}`` for terms not in the ontology.

        :param remove_root: If ``True``, exclude the ontology root (``CL:0000000``).
        :param remove_self: If ``True``, exclude *term* itself.
        """
        result = self.ancestors_dict.get(term, frozenset({term}))
        if remove_root:
            result = result - {CL_CELL_ROOT_NODE}
        if remove_self:
            result = result - {term}
        return result

    def get_descendants(self, term: str, remove_root: bool = False, remove_self: bool = False) -> t.FrozenSet[str]:
        """
        Return the descendant set for *term* (self-inclusive by default).

        Falls back to a singleton ``{term}`` for terms not in the ontology.

        :param remove_root: If ``True``, exclude the ontology root (``CL:0000000``).
        :param remove_self: If ``True``, exclude *term* itself.
        """
        result = self.descendants_dict.get(term, frozenset({term}))
        if remove_root:
            result = result - {CL_CELL_ROOT_NODE}
        if remove_self:
            result = result - {term}
        return result

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
