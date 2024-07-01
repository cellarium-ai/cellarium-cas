# classes
from .ontology_aware import (
    AggregatedCellOntologyScores,
    CellOntologyScoresAggregationOp,
    CellOntologyScoresAggregationDomain,
)

# methods
from .ontology_aware import (
    convert_cas_ontology_aware_response_to_score_matrix,
    insert_cas_ontology_aware_response_into_adata,
    get_aggregated_cas_ontology_aware_scores,
    convert_aggregated_cell_ontology_scores_to_rooted_tree,
    generate_phyloxml_from_scored_cell_ontology_tree,
    get_most_granular_top_k_calls,
    compute_most_granular_top_k_calls_single,
    compute_most_granular_top_k_calls_cluster,
)

from .cell_type_summary_statistics import (
    reduce_cas_query_result_by_majority_vote,
    reduce_cas_cell_type_summary_statistics_response_by_min_distance,
    reduce_cas_cell_type_summary_statistics_response_by_majority_vote_per_cluster,
    reduce_cas_cell_type_summary_statistics_response_by_wnn,
    get_interpolated_cell_type_colors,
)

from .common import get_obs_indices_for_cluster

# constants
from .ontology_aware import (
    CAS_CL_SCORES_ANNDATA_OBSM_KEY,
    CAS_METADATA_ANNDATA_UNS_KEY,
    CAS_SCORE_PROPERTY_REF,
    CAS_FRACTION_PROPERTY_REF,
    CAS_CL_LABEL_PROPERTY_REF,
)


__all__ = [
    "AggregatedCellOntologyScores",
    "CellOntologyScoresAggregationOp",
    "CellOntologyScoresAggregationDomain",
    "convert_cas_ontology_aware_response_to_score_matrix",
    "insert_cas_ontology_aware_response_into_adata",
    "get_aggregated_cas_ontology_aware_scores",
    "convert_aggregated_cell_ontology_scores_to_rooted_tree",
    "generate_phyloxml_from_scored_cell_ontology_tree",
    "get_obs_indices_for_cluster",
    # ontology aware
    "get_most_granular_top_k_calls",
    "compute_most_granular_top_k_calls_single",
    "compute_most_granular_top_k_calls_cluster",
    # cell type summary statistics
    "reduce_cas_query_result_by_majority_vote",
    "reduce_cas_cell_type_summary_statistics_response_by_min_distance",
    "reduce_cas_cell_type_summary_statistics_response_by_majority_vote_per_cluster",
    "reduce_cas_cell_type_summary_statistics_response_by_wnn",
    "get_interpolated_cell_type_colors",
    # constants
    "CAS_CL_SCORES_ANNDATA_OBSM_KEY",
    "CAS_METADATA_ANNDATA_UNS_KEY",
    "CAS_SCORE_PROPERTY_REF",
    "CAS_FRACTION_PROPERTY_REF",
    "CAS_CL_LABEL_PROPERTY_REF",
]
