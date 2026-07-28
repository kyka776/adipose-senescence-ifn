"""Cell-level program scoring and animal-level aggregation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def load_expression(path: Path) -> pd.DataFrame:
    if path.suffix == ".h5ad":
        try:
            import anndata as ad
        except ImportError as exc:
            raise RuntimeError("Install the `h5ad` extra to read .h5ad files") from exc
        obj = ad.read_h5ad(path)
        matrix = obj.X.toarray() if hasattr(obj.X, "toarray") else np.asarray(obj.X)
        frame = pd.DataFrame(matrix, index=obj.obs_names, columns=obj.var_names)
    else:
        frame = pd.read_csv(path, sep="\t", index_col=0)
    if frame.index.duplicated().any() or frame.columns.duplicated().any():
        raise ValueError("Expression cell and gene identifiers must be unique")
    if not all(np.issubdtype(dtype, np.number) for dtype in frame.dtypes):
        raise ValueError("Expression matrix contains non-numeric columns")
    values = frame.to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Expression matrix must contain finite non-negative values")
    return frame


def _zscore(matrix: pd.DataFrame) -> pd.DataFrame:
    means = matrix.mean(axis=0)
    stds = matrix.std(axis=0, ddof=0)
    usable = stds > 0
    return (matrix.loc[:, usable] - means[usable]) / stds[usable]


def score_programs(
    expression: pd.DataFrame,
    signature_path: Path,
    min_binary_genes: int = 10,
    min_binary_fraction: float = 0.20,
    min_sencat_genes: int = 500,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    payload = json.loads(signature_path.read_text(encoding="utf-8"))
    variants: dict[str, dict[str, float]] = payload["sencat_variants"]
    axes: dict[str, list[str]] = payload["binary_axes"]
    needed = set().union(
        *(set(weights) for weights in variants.values()),
        *(set(genes) for genes in axes.values()),
    )
    detected = [gene for gene in expression.columns if gene in needed]
    if not detected:
        raise ValueError("No signature genes detected in expression matrix")

    totals = expression.sum(axis=1)
    if (totals <= 0).any():
        raise ValueError("Cells with zero total counts must be removed before scoring")
    log_cpm = np.log1p(expression.loc[:, detected].div(totals, axis=0) * 10_000)
    standardized = _zscore(log_cpm)
    scores = pd.DataFrame(index=expression.index)
    coverage_rows: list[dict] = []

    for name, weights in variants.items():
        present = [gene for gene in weights if gene in standardized.columns]
        if len(present) < min_sencat_genes:
            raise ValueError(
                f"{name} coverage too low: {len(present)} detected; "
                f"minimum is {min_sencat_genes}"
            )
        coefficients = np.array([weights[gene] for gene in present], dtype=float)
        denominator = float(np.abs(coefficients).sum())
        scores[name] = standardized.loc[:, present].to_numpy() @ coefficients / denominator
        coverage_rows.append(
            {
                "score": name,
                "source_genes": len(weights),
                "detected_nonconstant_genes": len(present),
                "coverage_fraction": len(present) / len(weights),
            }
        )

    for name, genes in axes.items():
        present = [gene for gene in genes if gene in standardized.columns]
        coverage = len(present) / len(genes) if genes else 0.0
        if len(present) < min_binary_genes or coverage < min_binary_fraction:
            raise ValueError(
                f"{name} coverage too low: {len(present)}/{len(genes)} "
                f"({coverage:.1%})"
            )
        scores[name] = standardized.loc[:, present].mean(axis=1)
        coverage_rows.append(
            {
                "score": name,
                "source_genes": len(genes),
                "detected_nonconstant_genes": len(present),
                "coverage_fraction": coverage,
            }
        )
    return scores, pd.DataFrame(coverage_rows)


def aggregate_animal_scores(
    scores: pd.DataFrame, metadata: pd.DataFrame, statistic: str = "mean"
) -> pd.DataFrame:
    if not scores.index.equals(metadata.index):
        metadata = metadata.set_index("cell_id").loc[scores.index]
    joined = metadata[
        ["animal_id", "age_months", "depot", "sex", "cell_type", "subtype"]
    ].join(scores)
    group_columns = ["animal_id", "age_months", "depot", "sex", "cell_type", "subtype"]
    score_columns = list(scores.columns)
    grouped = joined.groupby(group_columns, observed=True)[score_columns]
    if statistic == "mean":
        aggregate = grouped.mean()
    elif statistic == "median":
        aggregate = grouped.median()
    else:
        raise ValueError("statistic must be 'mean' or 'median'")
    counts = grouped.size().rename("n_cells")
    return aggregate.join(counts).reset_index()

