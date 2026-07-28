"""End-to-end post-access analysis orchestration."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .figures import plot_animal_results
from .gate import audit_metadata, write_gate_report
from .scoring import aggregate_animal_scores, load_expression, score_programs
from .statistics import run_inference


def run_analysis(
    metadata_path: Path,
    expression_path: Path,
    signature_path: Path,
    output_dir: Path,
    n_permutations: int = 9999,
    n_bootstrap: int = 2000,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    gate, metadata = audit_metadata(metadata_path)
    write_gate_report(gate, output_dir)
    if not gate.passed or metadata is None:
        raise RuntimeError(
            "Gate 0 failed. Biological scoring was not run; inspect gate-report.md."
        )
    if not expression_path.exists():
        raise RuntimeError(f"Expression file is missing: {expression_path}")
    if not signature_path.exists():
        raise RuntimeError(
            f"Private signature bundle is missing: {signature_path}. Run `make predata`."
        )

    eligible = metadata[metadata["qc_pass"] == True].copy()  # noqa: E712
    vascular = eligible["cell_type"].astype(str).str.contains(
        "endothelial|vascular", case=False, regex=True
    )
    eligible = eligible.loc[vascular].set_index("cell_id")
    expression = load_expression(expression_path)
    missing_expression = eligible.index.difference(expression.index)
    if len(missing_expression):
        raise RuntimeError(
            f"{len(missing_expression):,} eligible metadata cells are absent from expression"
        )
    expression = expression.loc[eligible.index]

    scores, coverage = score_programs(expression, signature_path)
    coverage.to_csv(output_dir / "signature-coverage.tsv", sep="\t", index=False)
    aggregates = aggregate_animal_scores(scores, eligible, statistic="mean")
    medians = aggregate_animal_scores(scores, eligible, statistic="median")
    aggregates.to_csv(output_dir / "animal-aggregates-mean.tsv", sep="\t", index=False)
    medians.to_csv(output_dir / "animal-aggregates-median.tsv", sep="\t", index=False)

    score_columns = list(scores.columns)
    inference = run_inference(
        aggregates,
        score_columns,
        n_permutations=n_permutations,
        n_bootstrap=n_bootstrap,
    )
    inference.to_csv(output_dir / "inference.tsv", sep="\t", index=False)
    median_inference = run_inference(
        medians,
        score_columns,
        n_permutations=n_permutations,
        n_bootstrap=n_bootstrap,
    )
    median_inference.to_csv(
        output_dir / "inference-median-sensitivity.tsv", sep="\t", index=False
    )

    figure_dir = output_dir / "figures"
    for score in score_columns:
        plot_animal_results(
            aggregates,
            inference,
            score,
            figure_dir / f"animal-{score}",
        )

    summary = {
        "gate_passed": True,
        "cells_scored": len(scores),
        "animals": int(aggregates["animal_id"].nunique()),
        "tests": len(inference),
        "interpretation": (
            "Results require comparison against the preregistered falsification "
            "matrix; no causal claim follows from transcriptomic association."
        ),
    }
    (output_dir / "analysis-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary

