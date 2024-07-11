_API_BASE_ENDPOINT = "api"
_API_CELL_OPERATIONS_ENDPOINT_PREFIX = f"{_API_BASE_ENDPOINT}/cellarium-cell-operations"
_API_GENERAL_ENDPOINT_PREFIX = f"{_API_BASE_ENDPOINT}/cellarium-general"
# General
VALIDATE_TOKEN = f"{_API_GENERAL_ENDPOINT_PREFIX}/validate-token"
APPLICATION_INFO = f"{_API_GENERAL_ENDPOINT_PREFIX}/application-info"
FEEDBACK_OPT_OUT = f"{_API_GENERAL_ENDPOINT_PREFIX}/feedback/opt-out"
FEEDBACK_ANSWER = f"{_API_GENERAL_ENDPOINT_PREFIX}/feedback/answer?client_session_id={{client_session_id}}&client_action_id={{client_action_id}}"
GET_FEATURE_SCHEMAS = f"{_API_GENERAL_ENDPOINT_PREFIX}/feature-schemas"
GET_SCHEMA_BY_NAME = f"{_API_GENERAL_ENDPOINT_PREFIX}/feature-schema/{{schema_name}}"
LIST_MODELS = f"{_API_GENERAL_ENDPOINT_PREFIX}/list-models"
GET_USER_QUOTA = f"{_API_GENERAL_ENDPOINT_PREFIX}/quota"
VALIDATE_VERSION = f"{_API_GENERAL_ENDPOINT_PREFIX}/validate-client-version"
# Cell Analysis
ANNOTATE_CELL_TYPE_SUMMARY_STATS_STRATEGY = (
    f"{_API_CELL_OPERATIONS_ENDPOINT_PREFIX}/annotate-cell-type-summary-statistics-strategy"
)
ANNOTATE_CELL_TYPE_ONTOLOGY_AWARE_STRATEGY = (
    f"{_API_CELL_OPERATIONS_ENDPOINT_PREFIX}/annotate-cell-type-ontology-aware-strategy"
)
NEAREST_NEIGHBOR_SEARCH = f"{_API_CELL_OPERATIONS_ENDPOINT_PREFIX}/nearest-neighbor-search"
QUERY_CELLS_BY_IDS = f"{_API_CELL_OPERATIONS_ENDPOINT_PREFIX}/query-cells-by-ids"
