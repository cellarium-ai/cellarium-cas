try:
    from .flat import compute_flat_metrics, extract_predictions_from_adata  # noqa
    from .ontology_aware import compute_ontology_aware_metrics  # noqa

    __all__ = [
        "compute_flat_metrics",
        "compute_ontology_aware_metrics",
        "extract_predictions_from_adata",
    ]
except ImportError as e:
    raise ImportError(
        "Benchmarking dependencies are not installed. " "Install them with: pip install cellarium-cas[benchmark]"
    ) from e
