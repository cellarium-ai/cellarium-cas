_API_BASE_ENDPOINT = "api"
_API_CAS_ENDPOINT = f"{_API_BASE_ENDPOINT}/cellarium-cas"
_API_GENERAL_ENDPOINT = f"{_API_BASE_ENDPOINT}/cellarium-general"
# General
VALIDATE_TOKEN = f"{_API_GENERAL_ENDPOINT}/validate-token"
APPLICATION_INFO = f"{_API_GENERAL_ENDPOINT}/application-info"
GET_FEATURE_SCHEMAS = f"{_API_GENERAL_ENDPOINT}/feature-schemas"
GET_SCHEMA_BY_NAME = f"{_API_GENERAL_ENDPOINT}/feature-schema/{{schema_name}}"
LIST_MODELS = f"{_API_GENERAL_ENDPOINT}/list-models"
# Cell Analysis
ANNOTATE = f"{_API_CAS_ENDPOINT}/annotate"
NEAREST_NEIGHBOR_SEARCH = f"{_API_CAS_ENDPOINT}/nearest-neighbor-search"
QUERY_CELLS_BY_IDS = f"{_API_CAS_ENDPOINT}/query-cells-by-ids"
