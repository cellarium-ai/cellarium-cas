#!/usr/bin/env python3
"""
Build a CAS-compatible ``ontology_response.json`` from Azimuth cell type annotations.

Each cell's Azimuth annotations (one per granularity level) are mapped to CL ontology
terms via the HRA crosswalk and then propagated up the ontology ancestor chain.
Propagation follows these rules:

- Levels are processed **most granular first**.
- For each level, the mapped CL term and all its ancestors receive
  ``score = max(current_score, level_score)``.
- When a less-granular level's CL term is not a confirmed ancestor of the more-granular
  level's CL term (crosswalk linearity violation), a warning is emitted once per unique
  pair and propagation continues independently for both.

The resulting ``CellTypeOntologyAwareResults`` mimics the CAS API response format so
that ``cellarium-cas benchmark ontology-aware`` can evaluate Azimuth on equal footing.
Synthetic fields not produced by Azimuth are set as: ``total_neighbors=1``,
``total_neighbors_unrecognized=0``, ``total_weight=<most_granular_level_score>``.

Usage
-----
::

    python helpers/build_ontology_response.py \\
        --azimuth-csv            /data/azimuth_output.csv \\
        --h5ad-path              /data/my_dataset.h5ad \\
        --output-path            /data/azimuth_annotate_dir/ontology_response.json \\
        --crosswalk-csv          crosswalk.csv \\
        --crosswalk-azimuth-col  tool_cell_label \\
        --crosswalk-cl-id-col    cl_id \\
        --ontology-resource-path /path/to/ontology_resource.json \\
        --azimuth-ref-name       pbmcref \\
        --level predicted.celltype.l3:predicted.celltype.l3.score \\
        --level predicted.celltype.l2:predicted.celltype.l2.score \\
        --level predicted.celltype.l1:predicted.celltype.l1.score

``--level`` arguments must be ordered **most granular first**.  Each value is
``<azimuth_label_column>:<azimuth_score_column>`` as they appear in the Azimuth
output CSV.

The ``--ontology-resource-path`` is the ``ontology_resource.json`` saved by
``cellarium-cas annotate --save-ontology-resource``.
"""

import json
import typing as t
import warnings
from pathlib import Path

import click
from cellarium.cas.models import CellTypeOntologyAwareResults
from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache



def _build_parents_map(children_dict: t.Dict[str, t.List[str]]) -> t.Dict[str, t.List[str]]:
    """Invert ``children_dict`` (child → parents) for ancestor traversal."""
    parents_map: t.Dict[str, t.List[str]] = {}
    for parent, children in children_dict.items():
        for child in children:
            parents_map.setdefault(child, []).append(parent)
    return parents_map


def _get_node_and_ancestors(cl_id: str, parents_map: t.Dict[str, t.List[str]]) -> t.Set[str]:
    """Return *cl_id* and all its ontology ancestors via BFS."""
    visited: t.Set[str] = {cl_id}
    queue = [cl_id]
    while queue:
        node = queue.pop(0)
        for parent in parents_map.get(node, []):
            if parent not in visited:
                visited.add(parent)
                queue.append(parent)
    return visited


def build_ontology_response(
    azimuth_csv_path: str,
    h5ad_path: str,
    crosswalk_csv_path: str,
    crosswalk_azimuth_col: str,
    crosswalk_cl_id_col: str,
    level_specs: t.List[t.Tuple[str, str]],
    cell_ontology_cache: CellOntologyCache,
    azimuth_ref_name: str,
    output_path: str,
) -> CellTypeOntologyAwareResults:
    """
    Build a CAS-compatible ``CellTypeOntologyAwareResults`` from Azimuth annotations.

    :param azimuth_csv_path: Path to the Azimuth metadata CSV (row index = barcodes).
    :param h5ad_path: Path to the source ``.h5ad`` file; used for authoritative cell ordering.
    :param crosswalk_csv_path: Path to the HRA crosswalk CSV mapping Azimuth labels to CL IDs.
    :param crosswalk_azimuth_col: Column in the crosswalk containing Azimuth cell type labels.
    :param crosswalk_cl_id_col: Column in the crosswalk containing CL ontology term IDs.
    :param level_specs: List of ``(azimuth_label_col, azimuth_score_col)`` tuples ordered
        **most granular first**.
    :param cell_ontology_cache: :class:`~cellarium.cas.postprocessing.cell_ontology.CellOntologyCache`
        instance for ancestor traversal and label lookup.
    :param azimuth_ref_name: Azimuth reference name (e.g. ``"pbmcref"``), written into
        ``model_name`` as ``azimuth_<ref_name>``.
    :param output_path: Path to write ``ontology_response.json``.

    :returns: The constructed :class:`~cellarium.cas.models.CellTypeOntologyAwareResults`.
    :raises ValueError: If barcodes in the h5ad are missing from the Azimuth CSV, or if
        required columns are absent from either input CSV.
    """
    import anndata
    import pandas as pd

    # --- load inputs ---
    adata = anndata.read_h5ad(h5ad_path)
    obs_names = list(adata.obs_names)

    azimuth_df = pd.read_csv(azimuth_csv_path, index_col=0)
    missing_barcodes = set(obs_names) - set(azimuth_df.index)
    if missing_barcodes:
        sample = sorted(missing_barcodes)[:5]
        raise ValueError(
            f"{len(missing_barcodes)} barcodes from h5ad not found in Azimuth CSV "
            f"(first {len(sample)}): {sample}"
        )
    azimuth_df = azimuth_df.reindex(obs_names)

    crosswalk_df = pd.read_csv(crosswalk_csv_path)
    for col in (crosswalk_azimuth_col, crosswalk_cl_id_col):
        if col not in crosswalk_df.columns:
            raise ValueError(
                f"Column '{col}' not found in crosswalk CSV. "
                f"Available columns: {list(crosswalk_df.columns)}"
            )
    crosswalk_map: t.Dict[str, str] = dict(
        zip(crosswalk_df[crosswalk_azimuth_col].astype(str), crosswalk_df[crosswalk_cl_id_col].astype(str))
    )

    # --- pre-build ancestor traversal structures ---
    parents_map = _build_parents_map(cell_ontology_cache.children_dict)
    known_cl_ids = set(cell_ontology_cache.cl_names_to_idx_map.keys())
    ancestors_cache: t.Dict[str, t.Set[str]] = {}

    def get_ancestors_cached(cl_id: str) -> t.Set[str]:
        if cl_id not in ancestors_cache:
            ancestors_cache[cl_id] = _get_node_and_ancestors(cl_id, parents_map)
        return ancestors_cache[cl_id]

    # --- per-cell processing ---
    warned_missing_labels: t.Set[str] = set()
    warned_missing_cl_ids: t.Set[str] = set()
    warned_nonlinear_pairs: t.Set[t.Tuple[str, str]] = set()

    annotations: t.List[CellTypeOntologyAwareResults.OntologyAwareAnnotation] = []

    for barcode in obs_names:
        # Resolve (cl_id, score) for each level; None cl_id means unmapped
        level_cl_ids: t.List[t.Tuple[t.Optional[str], float]] = []
        for label_col, score_col in level_specs:
            azimuth_label = azimuth_df.at[barcode, label_col] if label_col in azimuth_df.columns else None
            azimuth_score = azimuth_df.at[barcode, score_col] if score_col in azimuth_df.columns else None

            score = float(azimuth_score) if azimuth_score is not None and azimuth_score == azimuth_score else 0.0
            cl_id: t.Optional[str] = None

            if azimuth_label is not None and azimuth_label == azimuth_label:  # not NaN
                label_str = str(azimuth_label)
                cl_id = crosswalk_map.get(label_str)
                if cl_id is None:
                    if label_str not in warned_missing_labels:
                        warnings.warn(
                            f"Azimuth label '{label_str}' (column '{label_col}') not found in crosswalk — skipped."
                        )
                        warned_missing_labels.add(label_str)
                elif cl_id not in known_cl_ids:
                    if cl_id not in warned_missing_cl_ids:
                        warnings.warn(
                            f"CL ID '{cl_id}' from crosswalk is not present in the cell ontology cache — skipped."
                        )
                        warned_missing_cl_ids.add(cl_id)
                    cl_id = None

            level_cl_ids.append((cl_id, score))

        # Linearity check (warn once per unique pair, A+C strategy)
        valid_ids = [(cl_id, score) for cl_id, score in level_cl_ids if cl_id is not None]
        for i in range(len(valid_ids) - 1):
            more_granular_cl = valid_ids[i][0]
            less_granular_cl = valid_ids[i + 1][0]
            pair = (more_granular_cl, less_granular_cl)
            if pair not in warned_nonlinear_pairs:
                ancestors_of_more_granular = get_ancestors_cached(more_granular_cl)
                if less_granular_cl not in ancestors_of_more_granular:
                    warnings.warn(
                        f"Crosswalk hierarchy violation: '{less_granular_cl}' is not an ancestor of "
                        f"'{more_granular_cl}' in the cell ontology — propagating independently."
                    )
                    warned_nonlinear_pairs.add(pair)

        # Score propagation: most-granular first, max-merge at shared ancestors
        node_scores: t.Dict[str, float] = {}
        for cl_id, score in level_cl_ids:
            if cl_id is None:
                continue
            for ancestor in get_ancestors_cached(cl_id):
                if ancestor in known_cl_ids:
                    node_scores[ancestor] = max(node_scores.get(ancestor, 0.0), score)

        matches = [
            CellTypeOntologyAwareResults.Match(
                cell_type_ontology_term_id=cl_id,
                cell_type=cell_ontology_cache.cl_names_to_labels_map.get(cl_id, cl_id),
                score=score,
            )
            for cl_id, score in node_scores.items()
        ]

        # total_weight: score of the most-granular non-null level
        most_granular_score = next((score for cl_id, score in level_cl_ids if cl_id is not None), 0.0)

        annotations.append(
            CellTypeOntologyAwareResults.OntologyAwareAnnotation(
                query_cell_id=barcode,
                matches=matches,
                total_weight=most_granular_score,
                total_neighbors=1,
                total_neighbors_unrecognized=0,
            )
        )

    result = CellTypeOntologyAwareResults(
        data=annotations,
        model_name=f"azimuth_{azimuth_ref_name}",
    )

    output_path_obj = Path(output_path)
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path_obj, "w") as f:
        f.write(result.model_dump_json())

    return result


def _parse_level(value: str, param: t.Any, ctx: t.Any) -> t.Tuple[str, str]:
    parts = value.split(":", 1)
    if len(parts) != 2:
        raise click.BadParameter(f"must be in the form LABEL_COL:SCORE_COL, got: '{value}'")
    return parts[0].strip(), parts[1].strip()


@click.command("build-ontology-response")
@click.option("--azimuth-csv", required=True, type=click.Path(exists=True), help="Path to Azimuth metadata CSV.")
@click.option("--h5ad-path", required=True, type=click.Path(exists=True), help="Path to source .h5ad file.")
@click.option("--output-path", required=True, type=click.Path(), help="Path to write ontology_response.json.")
@click.option("--crosswalk-csv", required=True, type=click.Path(exists=True), help="Path to HRA crosswalk CSV.")
@click.option(
    "--crosswalk-azimuth-col",
    required=True,
    help="Column in crosswalk containing Azimuth cell type labels.",
)
@click.option(
    "--crosswalk-cl-id-col",
    required=True,
    help="Column in crosswalk containing CL ontology term IDs.",
)
@click.option(
    "--ontology-resource-path",
    required=True,
    type=click.Path(exists=True),
    help=(
        "Path to ontology_resource.json (saved by 'cellarium-cas annotate "
        "--save-ontology-resource'). Used to construct the CellOntologyCache."
    ),
)
@click.option(
    "--azimuth-ref-name",
    required=True,
    help="Azimuth reference name (e.g. pbmcref). Written into model_name as azimuth_<ref_name>.",
)
@click.option(
    "--level",
    "levels",
    multiple=True,
    required=True,
    metavar="LABEL_COL:SCORE_COL",
    callback=lambda ctx, param, values: [_parse_level(v, param, ctx) for v in values],
    help=(
        "Azimuth label and score column pair for one annotation level, as "
        "LABEL_COL:SCORE_COL. Repeat for each level, most granular first."
    ),
)
def main(
    azimuth_csv: str,
    h5ad_path: str,
    output_path: str,
    crosswalk_csv: str,
    crosswalk_azimuth_col: str,
    crosswalk_cl_id_col: str,
    ontology_resource_path: str,
    azimuth_ref_name: str,
    levels: t.Tuple[t.Tuple[str, str], ...],
) -> None:
    """Build a CAS-compatible ontology_response.json from Azimuth annotations."""
    with open(ontology_resource_path) as f:
        ontology_resource = json.load(f)
    cell_ontology_cache = CellOntologyCache(ontology_resource)

    try:
        build_ontology_response(
            azimuth_csv_path=azimuth_csv,
            h5ad_path=h5ad_path,
            crosswalk_csv_path=crosswalk_csv,
            crosswalk_azimuth_col=crosswalk_azimuth_col,
            crosswalk_cl_id_col=crosswalk_cl_id_col,
            level_specs=list(levels),
            cell_ontology_cache=cell_ontology_cache,
            azimuth_ref_name=azimuth_ref_name,
            output_path=output_path,
        )
    except ValueError as e:
        raise click.UsageError(str(e)) from e

    click.echo(f"Saved ontology response → {output_path}")


if __name__ == "__main__":
    main()
