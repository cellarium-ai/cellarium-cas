#!/usr/bin/env python3
"""
Run scTab inference on one or more ``.h5ad`` files and write per-cell top-k
predictions as CSV.

This script is the scTab comparator's container entrypoint.  It loads a single
checkpoint once per invocation, aligns each input's gene space to the model's
19331-gene vocabulary, applies scTab's normalization (10,000-count scaling +
log1p), and writes ``predicted_label_k`` / ``predicted_score_k`` columns (softmax
probabilities, k = 1..--top-k) for every cell.

Output layout (one directory per input, named from ``--output-dir-template``):

- ``<template with {stem}>/sctab_predictions.csv``
- ``<template with {stem}>/sctab_stats.json``

The CSV index is the input's ``obs_names`` (barcodes); the first column header
is empty, so downstream readers load it with ``index_col=0``.

Following the repo's CLI pattern (see ``cellarium/cas/cli/annotate.py``), the core
logic lives in :func:`run_sc_tab`, which is importable and click-free, and the
:func:`run_sctab_command` click command is a thin wrapper around it.

Top-level imports stay light (argparse-adjacent: click / json / re / warnings /
pathlib / numpy / pandas / anndata) so this module imports in environments
without torch.  torch, cellnet, and yaml are imported lazily inside functions.
"""

import json
import re
import typing as t
import warnings
from pathlib import Path

import anndata
import click
import numpy as np
import pandas as pd

# Fallback architecture defaults, mirroring cellnet/models.py (scTab @ 5ede7f2).
_TABNET_DEFAULTS = {
    "n_d": 128,
    "n_a": 64,
    "n_steps": 3,
    "gamma": 1.3,
    "n_independent": 5,
    "n_shared": 3,
    "epsilon": 1e-15,
    "virtual_batch_size": 256,
    "momentum": 0.02,
    "mask_type": "entmax",
}
_MLP_DEFAULTS = {"hidden_size": 128, "n_hidden": 8, "dropout": 0.2}

_MODEL_TYPE_CHOICES = ("tabnet", "mlp", "linear")


def strip_version(gene_id: str) -> str:
    """Return *gene_id* without a trailing Ensembl-style ``.<digits>`` version suffix."""
    return re.sub(r"\.\d+$", "", str(gene_id))


def align_to_model_genes(
    counts,
    input_gene_ids,
    model_genes,
    min_overlap: float,
) -> tuple:
    """Align *counts* (cells x input genes) to the model's gene order.

    Returns ``(matrix, matched_count, expected_count)`` where *matrix* is a dense
    ``float32`` array with one column per model gene (missing genes zero-filled)
    and input genes are matched by version-stripped id, summing duplicate input
    ids into one column with a warning.

    :raises ValueError: If the fraction of matched model genes is below
        *min_overlap*.
    """
    model_index: dict = {}
    for i, gene in enumerate(model_genes):
        key = strip_version(gene)
        if key in model_index:
            raise ValueError(f"Duplicate model gene id '{key}' (columns {model_index[key]} and {i}).")
        model_index[key] = i

    n_cells = counts.shape[0]
    out = np.zeros((n_cells, len(model_genes)), dtype=np.float32)
    matched: set = set()
    written_counts: dict = {}
    for j, gene in enumerate(input_gene_ids):
        key = strip_version(gene)
        target = model_index.get(key)
        if target is None:
            continue
        col = counts[:, j]
        if hasattr(col, "toarray"):  # scipy sparse column
            col = col.toarray().ravel()
        out[:, target] += np.asarray(col, dtype=np.float32).ravel()
        matched.add(key)
        if key in written_counts:
            if written_counts[key] == 1:
                warnings.warn(
                    f"Duplicate input gene id '{key}' — summing counts into one model-gene column.",
                    UserWarning,
                    stacklevel=2,
                )
            written_counts[key] += 1
        else:
            written_counts[key] = 1

    expected_count = len(model_genes)
    matched_count = len(matched)
    if matched_count / expected_count < min_overlap:
        sample_input = [str(g) for g in input_gene_ids[:5]]
        sample_model = [str(g) for g in model_genes[:5]]
        raise ValueError(
            f"Gene overlap {matched_count}/{expected_count} = {matched_count / expected_count:.3f} "
            f"below --min-gene-overlap {min_overlap}. Sample input ids: {sample_input}; "
            f"sample model genes: {sample_model}."
        )
    return out, matched_count, expected_count


def _read_hparams(hparams_path: str) -> dict:
    """Load the ``hparams.yaml`` dict (yaml imported lazily)."""
    import yaml

    with open(hparams_path) as f:
        hp = yaml.full_load(f)
    if not isinstance(hp, dict):
        raise ValueError(f"hparams file {hparams_path} must contain a YAML mapping.")
    return hp


def _build_model(model_type: str, hp: dict, checkpoint: str):
    """Build the bare classifier submodule and load the checkpoint's weights.

    The checkpoint stores the classifier under a ``classifier.`` prefix; it is
    stripped here so the submodule loads with ``strict=True``.  torch and cellnet
    are imported lazily.
    """
    import torch

    gene_dim = hp["gene_dim"]
    type_dim = hp["type_dim"]

    if model_type == "tabnet":
        from cellnet.tabnet.tab_network import TabNet

        model = TabNet(
            input_dim=gene_dim,
            output_dim=type_dim,
            n_d=hp.get("n_d", _TABNET_DEFAULTS["n_d"]),
            n_a=hp.get("n_a", _TABNET_DEFAULTS["n_a"]),
            n_steps=hp.get("n_steps", _TABNET_DEFAULTS["n_steps"]),
            gamma=hp.get("gamma", _TABNET_DEFAULTS["gamma"]),
            n_independent=hp.get("n_independent", _TABNET_DEFAULTS["n_independent"]),
            n_shared=hp.get("n_shared", _TABNET_DEFAULTS["n_shared"]),
            epsilon=hp.get("epsilon", _TABNET_DEFAULTS["epsilon"]),
            virtual_batch_size=hp.get("virtual_batch_size", _TABNET_DEFAULTS["virtual_batch_size"]),
            momentum=hp.get("momentum", _TABNET_DEFAULTS["momentum"]),
            mask_type=hp.get("mask_type", _TABNET_DEFAULTS["mask_type"]),
        )
    elif model_type == "mlp":
        from cellnet.models import MLP

        model = MLP(
            input_dim=gene_dim,
            output_dim=type_dim,
            hidden_size=hp.get("hidden_size", _MLP_DEFAULTS["hidden_size"]),
            n_hidden=hp.get("n_hidden", _MLP_DEFAULTS["n_hidden"]),
            dropout=hp.get("dropout", _MLP_DEFAULTS["dropout"]),
        )
    elif model_type == "linear":
        model = torch.nn.Linear(gene_dim, type_dim)
    else:
        raise ValueError(f"Unknown --model-type '{model_type}' (expected tabnet, mlp, or linear).")

    state_dict = torch.load(checkpoint, map_location="cpu", weights_only=False)["state_dict"]
    classifier_weights = {k[len("classifier.") :]: v for k, v in state_dict.items() if k.startswith("classifier.")}
    if not classifier_weights:
        raise ValueError(f"No 'classifier.*' keys found in checkpoint {checkpoint}.")
    model.load_state_dict(classifier_weights, strict=True)
    model.eval()
    return model


def _model_genes(var_parquet_path: str) -> list:
    """Return the model gene ids in output order from ``var.parquet``."""
    var_df = pd.read_parquet(var_parquet_path)
    if "feature_id" in var_df.columns:
        return var_df["feature_id"].astype(str).tolist()
    if "feature_name" in var_df.columns:
        return var_df["feature_name"].astype(str).tolist()
    raise ValueError(
        f"var.parquet {var_parquet_path} has neither 'feature_id' nor 'feature_name'. "
        f"Available columns: {list(var_df.columns)}"
    )


def _class_labels(cell_type_parquet_path: str) -> list:
    """Return the class labels in output-class order from ``cell_type.parquet``."""
    return pd.read_parquet(cell_type_parquet_path)["label"].tolist()


def _predict_to_csv(
    model,
    matrix: np.ndarray,
    obs_names: list,
    class_labels: list,
    top_k: int,
    chunk_size: int,
) -> pd.DataFrame:
    """Run chunked forward passes and return a DataFrame of top-k labels/scores."""
    import torch

    x = torch.from_numpy(matrix)
    label_cols = []
    score_cols = []
    with torch.no_grad():
        for start in range(0, matrix.shape[0], chunk_size):
            x_chunk = x[start : start + chunk_size]
            counts_sum = x_chunk.sum(dim=1, keepdim=True)
            counts_sum += counts_sum == 0
            x_norm = torch.log1p(10000.0 / counts_sum * x_chunk)
            logits = model(x_norm)
            if isinstance(logits, tuple):  # TabNet returns (logits, m_loss)
                logits = logits[0]
            probs = torch.softmax(logits, dim=1)
            values, indices = torch.topk(probs, k=min(top_k, probs.shape[1]), dim=1)
            label_cols.extend([[class_labels[i] for i in row] for row in indices.tolist()])
            score_cols.extend([row for row in values.tolist()])
    columns = []
    for k in range(1, top_k + 1):
        columns.append(f"predicted_label_{k}")
        columns.append(f"predicted_score_{k}")
    data = []
    for labels, scores in zip(label_cols, score_cols):
        row = []
        for k in range(top_k):
            row.append(labels[k] if k < len(labels) else None)
            row.append(scores[k] if k < len(scores) else None)
        data.append(row)
    return pd.DataFrame(data, index=pd.Index(obs_names, name=None), columns=columns)


def run_sc_tab(
    model_type: str,
    checkpoint: str,
    hparams: str,
    var_parquet: str,
    cell_type_parquet: str,
    input_h5ad: t.Sequence[str],
    output_dir_template: str,
    chunk_size: int = 2048,
    min_gene_overlap: float = 0.5,
    top_k: int = 5,
) -> t.List[t.Dict[str, t.Any]]:
    """
    Run scTab inference over one or more inputs and write top-k prediction CSVs.

    This is the plain-Python entry point.  The CLI command :func:`run_sctab_command`
    is a thin wrapper around this function and passes all arguments through.

    :param model_type: scTab classifier architecture, one of ``tabnet``, ``mlp``,
        ``linear``.
    :param checkpoint: Path to a scTab checkpoint (loads its ``classifier.*`` weights).
    :param hparams: Path to the ``hparams.yaml`` used to build the classifier.
    :param var_parquet: Path to the model's ``var.parquet`` (defines the 19331-gene order).
    :param cell_type_parquet: Path to the model's ``cell_type.parquet`` (defines the
        class-label order).
    :param input_h5ad: One or more ``.h5ad`` input paths.
    :param output_dir_template: Directory template per input, containing ``{stem}``
        (replaced by the input filename without extension).
    :param chunk_size: Number of cells per forward-pass chunk.
    :param min_gene_overlap: Minimum fraction of model genes that must be present in an
        input, else ``ValueError``.
    :param top_k: Number of top predictions (labels + softmax scores) to write per cell.

    :returns: List of per-input stats dicts (``stem``, ``n_cells``, ``n_model_genes``,
        ``n_matched_genes``, ``overlap_ratio``, ``model_type``).
    :raises ValueError: If gene overlap falls below *min_gene_overlap*, or an input's
        gene ids / the model metadata use unexpected columns.
    """
    hp = _read_hparams(hparams)
    model = _build_model(model_type, hp, checkpoint)
    model_genes = _model_genes(var_parquet)
    class_labels = _class_labels(cell_type_parquet)

    stats_list: t.List[t.Dict[str, t.Any]] = []
    for input_path in input_h5ad:
        adata = anndata.read_h5ad(input_path)
        if "feature_id" in adata.var.columns:
            input_gene_ids = adata.var["feature_id"].astype(str).tolist()
        else:
            input_gene_ids = adata.var_names.astype(str).tolist()
        matrix, matched, expected = align_to_model_genes(adata.X, input_gene_ids, model_genes, min_gene_overlap)
        predictions_df = _predict_to_csv(
            model, matrix, adata.obs_names.astype(str).tolist(), class_labels, top_k, chunk_size
        )

        stem = Path(input_path).name
        if stem.endswith(".h5ad"):
            stem = stem[: -len(".h5ad")]
        output_dir = Path(output_dir_template.format(stem=stem))
        output_dir.mkdir(parents=True, exist_ok=True)
        predictions_df.to_csv(output_dir / "sctab_predictions.csv")

        stats = {
            "input_h5ad": str(Path(input_path).resolve()),
            "n_cells": int(matrix.shape[0]),
            "n_model_genes": expected,
            "n_matched_genes": matched,
            "overlap_ratio": float(matched / expected),
            "model_type": model_type,
        }
        with open(output_dir / "sctab_stats.json", "w") as f:
            json.dump(stats, f, indent=2)

        stats_list.append({"stem": stem, **stats})
    return stats_list


@click.command("run-sctab")
@click.option(
    "--model-type",
    required=True,
    type=click.Choice(_MODEL_TYPE_CHOICES),
    help="scTab classifier architecture.",
)
@click.option(
    "--checkpoint",
    required=True,
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="Path to a scTab checkpoint (loads its 'classifier.*' weights).",
)
@click.option(
    "--hparams",
    required=True,
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="Path to the hparams.yaml used to build the classifier.",
)
@click.option(
    "--var-parquet",
    required=True,
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="Path to the model's var.parquet (defines the 19331-gene order).",
)
@click.option(
    "--cell-type-parquet",
    required=True,
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="Path to the model's cell_type.parquet (defines the class-label order).",
)
@click.option(
    "--input-h5ad",
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="Input .h5ad file. May be repeated to process several inputs with one checkpoint load.",
)
@click.option(
    "--output-dir-template",
    required=True,
    help="Directory template per input, containing '{stem}' (replaced by the input filename without extension).",
)
@click.option(
    "--chunk-size",
    default=2048,
    show_default=True,
    type=click.IntRange(min=1),
    help="Number of cells per forward-pass chunk.",
)
@click.option(
    "--min-gene-overlap",
    default=0.5,
    show_default=True,
    type=click.FloatRange(min=0.0, max=1.0),
    help="Minimum fraction of model genes that must be present in an input.",
)
@click.option(
    "--top-k",
    default=5,
    show_default=True,
    type=click.IntRange(min=1),
    help="Number of top predictions (labels + softmax scores) to write per cell.",
)
def run_sctab_command(
    model_type: str,
    checkpoint: str,
    hparams: str,
    var_parquet: str,
    cell_type_parquet: str,
    input_h5ad: t.Sequence[str],
    output_dir_template: str,
    chunk_size: int,
    min_gene_overlap: float,
    top_k: int,
) -> None:
    """Run scTab inference over one or more .h5ad inputs, writing top-k prediction CSVs."""
    click.echo(f"Loading scTab model ({model_type}) from {checkpoint} ...")

    stats_list = run_sc_tab(
        model_type=model_type,
        checkpoint=checkpoint,
        hparams=hparams,
        var_parquet=var_parquet,
        cell_type_parquet=cell_type_parquet,
        input_h5ad=input_h5ad,
        output_dir_template=output_dir_template,
        chunk_size=chunk_size,
        min_gene_overlap=min_gene_overlap,
        top_k=top_k,
    )

    for stats in stats_list:
        click.echo(
            f"[sctab] {stats['stem']}: cells={stats['n_cells']} "
            f"matched_genes={stats['n_matched_genes']}/{stats['n_model_genes']} "
            f"model={stats['model_type']}"
        )


if __name__ == "__main__":
    run_sctab_command()
