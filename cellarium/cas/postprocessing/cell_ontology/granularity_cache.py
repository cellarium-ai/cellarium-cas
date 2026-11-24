import typing as t

# Curated granularity levels, from broad (L1) to most granular (L4).
GRANULARITY_LEVELS = ("l1", "l2", "l3", "l4")

_UNDEFINED_LEVEL_VALUES = frozenset({"", "unknown", "none", "null"})


def _normalize_level(value: t.Any) -> t.Optional[str]:
    """Normalize a granularity level value, mapping placeholders to ``None``."""
    if value is None:
        return None
    if isinstance(value, str) and value.strip().lower() in _UNDEFINED_LEVEL_VALUES:
        return None
    return value


class CellGranularityCache:
    """
    A cache for the curated L1-L4 cell type granularity lookup, populated from a preloaded resource dict.

    The resource maps each Cell Ontology term ID to its curated granularity levels::

        {
            "CL:0000084": {"l1": "CL:0000988", "l2": "CL:0000084", "l3": None, "l4": None},
            ...
        }

    Curated depth is the highest defined level (1-4), or ``0`` when no level is defined. A mapping that
    repeats the same CL ID at L3 and L4 intentionally has curated depth ``4``.

    Attributes:
        granularity_dict (dict): Mapping from CL term ID to normalized ``{"l1", "l2", "l3", "l4"}`` levels.
        depth_dict (dict): Mapping from CL term ID to curated depth (0-4).
    """

    def __init__(self, resource: t.Dict[str, t.Dict[str, t.Any]]) -> None:
        """
        Initialize from a preloaded granularity lookup resource.

        :param resource: Mapping from CL term ID to a dict containing ``l1``, ``l2``, ``l3``, and ``l4``.
        """
        self.granularity_dict: t.Dict[str, t.Dict[str, t.Optional[str]]] = {}
        self.depth_dict: t.Dict[str, int] = {}

        for cl_name, raw_levels in resource.items():
            levels = {level: _normalize_level(raw_levels.get(level)) for level in GRANULARITY_LEVELS}
            self.granularity_dict[cl_name] = levels
            self.depth_dict[cl_name] = max(
                (i for i in range(1, 5) if levels[f"l{i}"] is not None),
                default=0,
            )

    def get_depth(self, cl_name: str) -> int:
        """
        Return the curated granularity depth of *cl_name*, or ``0`` when it has no curated mapping.

        :param cl_name: Cell Ontology term ID.
        """
        return self.depth_dict.get(cl_name, 0)

    def get_levels(self, cl_name: str) -> t.Dict[str, t.Optional[str]]:
        """
        Return a copy of the curated L1-L4 levels for *cl_name*.

        Missing terms return all levels as ``None``.

        :param cl_name: Cell Ontology term ID.
        """
        levels = self.granularity_dict.get(cl_name)
        if levels is None:
            return {level: None for level in GRANULARITY_LEVELS}
        return dict(levels)
