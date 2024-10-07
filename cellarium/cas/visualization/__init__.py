from cellarium.cas.logging import logger

try:
    from .circular_tree_plot_umap_dash_app.app import CASCircularTreePlotUMAPDashApp  # noqa
    from .ui_utils import find_and_kill_process  # noqa

    __all__ = ["CASCircularTreePlotUMAPDashApp", "find_and_kill_process"]
except ImportError:
    logger.warning(
        """
Visualization dependencies not installed.
To install the Cellarium CAS Client with visualation dependencies, please run:
pip install --force-reinstall cellarium-cas[vis]
"""
    )
    __all__ = []
