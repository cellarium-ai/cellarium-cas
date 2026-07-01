try:
    from .confusion_matrix import (  # noqa
        aggregate_confusion_matrices,
        build_confusion_matrix,
        load_confusion_matrix,
        save_confusion_matrix,
    )
    from .f_measure import compute_f_measure_from_cm  # noqa
    from .hierarchical_f_measure import compute_hierarchical_f_measure_from_cm  # noqa

    __all__ = [
        "build_confusion_matrix",
        "save_confusion_matrix",
        "load_confusion_matrix",
        "aggregate_confusion_matrices",
        "compute_f_measure_from_cm",
        "compute_hierarchical_f_measure_from_cm",
    ]
except ImportError as e:
    raise ImportError(
        "Benchmarking dependencies are not installed. Install them with: pip install cellarium-cas[benchmark]"
    ) from e
