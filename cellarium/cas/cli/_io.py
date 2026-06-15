"""
Internal I/O helpers for the cellarium-cas CLI.

All functions use plain local file paths. Callers are responsible for any data staging
(e.g. KFP's /gcs/ mount makes GCS paths transparent).
"""

import json
import typing as t
from pathlib import Path

import pandas as pd

from cellarium.cas import models

ONTOLOGY_RESPONSE_FILENAME = "ontology_response.json"
ONTOLOGY_RESOURCE_FILENAME = "ontology_resource.json"
INFERRED_LABELS_FILENAME = "inferred_labels.csv"
METADATA_FILENAME = "metadata.json"

_REQUIRED_FILES = (ONTOLOGY_RESPONSE_FILENAME, INFERRED_LABELS_FILENAME, METADATA_FILENAME)


def save_ontology_response(response: "models.CellTypeOntologyAwareResults", output_dir: Path) -> Path:
    """Serialize ``response`` to ``ontology_response.json`` inside *output_dir*."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / ONTOLOGY_RESPONSE_FILENAME
    with open(path, "w") as f:
        f.write(response.model_dump_json())
    return path


def load_ontology_response(output_dir: Path) -> "models.CellTypeOntologyAwareResults":
    """Load and validate ``ontology_response.json`` from *output_dir*."""
    path = Path(output_dir) / ONTOLOGY_RESPONSE_FILENAME
    if not path.exists():
        raise FileNotFoundError(f"ontology_response.json not found in {output_dir}")
    with open(path) as f:
        return models.CellTypeOntologyAwareResults.model_validate_json(f.read())


def save_inferred_labels(df: pd.DataFrame, output_dir: Path) -> Path:
    """Save label predictions DataFrame to ``inferred_labels.csv`` inside *output_dir*."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / INFERRED_LABELS_FILENAME
    df.to_csv(path)
    return path


def save_ontology_resource(resource: t.Dict[str, t.Any], output_dir: Path) -> Path:
    """Save the raw ontology resource dict to ``ontology_resource.json`` inside *output_dir*."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / ONTOLOGY_RESOURCE_FILENAME
    with open(path, "w") as f:
        json.dump(resource, f)
    return path


def load_ontology_resource(output_dir: Path) -> t.Dict[str, t.Any]:
    """
    Load ``ontology_resource.json`` from *output_dir*.

    Raises a descriptive :class:`FileNotFoundError` if the file is absent, reminding the
    caller to re-run ``cellarium-cas annotate`` with ``--save-ontology-resource``.
    """
    path = Path(output_dir) / ONTOLOGY_RESOURCE_FILENAME
    if not path.exists():
        raise FileNotFoundError(
            f"ontology_resource.json not found in {output_dir}. "
            "Re-run 'cellarium-cas annotate' with --save-ontology-resource to generate it."
        )
    with open(path) as f:
        return json.load(f)


def save_metadata(metadata: t.Dict[str, t.Any], output_dir: Path) -> Path:
    """Save the metadata dict to ``metadata.json`` inside *output_dir*."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / METADATA_FILENAME
    with open(path, "w") as f:
        json.dump(metadata, f, indent=2)
    return path


def load_metadata(output_dir: Path) -> t.Dict[str, t.Any]:
    """
    Load ``metadata.json`` from *output_dir*.

    Raises a descriptive :class:`FileNotFoundError` if the file is absent, reminding the
    caller to re-run ``cellarium-cas annotate`` with ``--save-metadata``.
    """
    path = Path(output_dir) / METADATA_FILENAME
    if not path.exists():
        raise FileNotFoundError(
            f"metadata.json not found in {output_dir}. "
            "Re-run 'cellarium-cas annotate' with --save-metadata (enabled by default) to generate it."
        )
    with open(path) as f:
        return json.load(f)


def validate_annotate_output_dir(
    output_dir: Path,
    required_files: t.Sequence[str] = _REQUIRED_FILES,
    missing_hint: str = "Ensure 'cellarium-cas annotate' was run with --infer-labels and --save-metadata.",
) -> None:
    """
    Validate that *output_dir* contains all files required for benchmarking.

    Required files default to ``ontology_response.json``, ``inferred_labels.csv``, ``metadata.json``.

    :raises ValueError: If any required file is missing, listing all absent files.
    """
    output_dir = Path(output_dir)
    missing = [f for f in required_files if not (output_dir / f).exists()]
    if missing:
        raise ValueError(
            f"Annotate output directory '{output_dir}' is missing required files: {missing}. {missing_hint}"
        )


def collect_annotate_output_dirs(
    input_path: t.Union[str, Path],
    required_files: t.Sequence[str] = _REQUIRED_FILES,
    missing_hint: str = "Ensure 'cellarium-cas annotate' was run with --infer-labels and --save-metadata.",
) -> t.List[Path]:
    """
    Resolve a list of annotate output directories from *input_path*.

    *input_path* may be either:

    - A ``.txt`` file where each line is a path to an annotate output directory.
    - A directory whose immediate subdirectories are the annotate output directories.

    Each resolved directory is validated via :func:`validate_annotate_output_dir`.

    :param input_path: Path to a ``.txt`` file or a parent directory.
    :returns: List of validated :class:`~pathlib.Path` objects.
    :raises ValueError: If any resolved directory fails validation.
    """
    input_path = Path(input_path)
    if input_path.suffix == ".txt":
        lines = input_path.read_text().splitlines()
        dirs = [Path(line.strip()) for line in lines if line.strip()]
    elif input_path.is_dir():
        dirs = sorted(p for p in input_path.iterdir() if p.is_dir())
    else:
        raise ValueError(f"--annotate-dirs must be a .txt file or a directory, got: {input_path}")

    if not dirs:
        raise ValueError(f"No annotate output directories found in: {input_path}")

    for d in dirs:
        validate_annotate_output_dir(d, required_files=required_files, missing_hint=missing_hint)

    return dirs
