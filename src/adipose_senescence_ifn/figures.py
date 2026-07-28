"""Deterministic project figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


COLORS = {
    "sencat": "#6D4C9F",
    "ifng_response": "#E66100",
    "mhcii_antigen_presentation": "#1B9E77",
    "cell_cycle_arrest": "#7570B3",
    "residual": "#2C3E50",
    "muted": "#667085",
    "background": "#F7F8FA",
}

# Stable element IDs plus fixed SVG metadata make source-backed figures
# byte-for-byte reproducible across identical builds.
plt.rcParams["svg.hashsalt"] = "adipose-senescence-ifn-0.1"


def _save(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        stem.with_suffix(".svg"),
        bbox_inches="tight",
        facecolor="white",
        metadata={"Date": "2026-07-29"},
    )
    fig.savefig(
        stem.with_suffix(".png"),
        dpi=180,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def plot_design_schematic(output_stem: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 5.8))
    fig.patch.set_facecolor("white")
    ax.set_facecolor(COLORS["background"])
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 6)
    ax.axis("off")

    ax.text(
        0.4,
        5.55,
        "PRE-DATA DESIGN — NOT A BIOLOGICAL RESULT",
        fontsize=11,
        weight="bold",
        color="#B42318",
    )
    ax.text(
        0.4,
        5.10,
        "Initiative 07: separate broad senescence from overlapping inflammatory and arrest programs",
        fontsize=14,
        weight="bold",
        color="#101828",
    )

    boxes = [
        (0.5, 2.0, 2.1, 1.55, "Matched atlas", "cells → animal IDs\nSAT/VAT + age + subtype", "#344054"),
        (3.2, 3.65, 2.1, 1.15, "Full SenCat", "signed weights", COLORS["sencat"]),
        (3.2, 2.10, 2.1, 1.15, "IFN-γ / MHC-II", "Reactome v97", COLORS["ifng_response"]),
        (3.2, 0.55, 2.1, 1.15, "Cell-cycle arrest", "TP53 G1/G2", COLORS["cell_cycle_arrest"]),
        (6.0, 2.0, 2.25, 1.55, "Overlap depletion", "full and four\npredeclared variants", COLORS["residual"]),
        (8.9, 2.0, 1.7, 1.55, "Inference", "animal × depot\n× subtype", "#027A48"),
    ]
    for x, y, width, height, title, body, color in boxes:
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.08,rounding_size=0.08",
            linewidth=1.8,
            edgecolor=color,
            facecolor="white",
        )
        ax.add_patch(patch)
        ax.text(x + 0.14, y + height - 0.35, title, fontsize=11, weight="bold", color=color)
        ax.text(x + 0.14, y + 0.25, body, fontsize=9.5, color="#344054", va="bottom")

    arrows = [
        ((2.62, 2.77), (3.12, 4.20)),
        ((2.62, 2.77), (3.12, 2.67)),
        ((2.62, 2.77), (3.12, 1.12)),
        ((5.32, 4.20), (5.92, 3.12)),
        ((5.32, 2.67), (5.92, 2.77)),
        ((5.32, 1.12), (5.92, 2.40)),
        ((8.28, 2.77), (8.82, 2.77)),
    ]
    for start, end in arrows:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.4,
                color="#98A2B3",
            )
        )

    ax.text(
        5.5,
        0.18,
        "Gate 0: no expression scoring until animal-level metadata and replication pass",
        ha="center",
        fontsize=10.5,
        color="#B42318",
        weight="bold",
    )
    _save(fig, output_stem)


def plot_overlap(summary: pd.DataFrame, output_stem: Path) -> None:
    labels = {
        "ifng_response": "IFN-γ signaling",
        "mhcii_antigen_presentation": "MHC-II presentation",
        "cell_cycle_arrest": "TP53 G1/G2 arrest",
    }
    ordered = summary.copy()
    ordered["label"] = ordered["axis"].map(labels).fillna(ordered["axis"])
    ordered = ordered.sort_values("overlap_with_sencat")
    colors = [COLORS.get(axis, COLORS["muted"]) for axis in ordered["axis"]]

    fig, (ax, table_ax) = plt.subplots(
        1,
        2,
        figsize=(13, 4.8),
        gridspec_kw={"width_ratios": [1.35, 1.15]},
    )
    bars = ax.barh(
        ordered["label"],
        ordered["overlap_with_sencat"],
        color=colors,
        edgecolor="white",
        height=0.58,
    )
    ax.bar_label(bars, padding=5, fontsize=10, weight="bold")
    ax.set_xlabel("Mapped mouse genes shared with SenCat")
    ax.set_title("Pre-data signature overlap audit", loc="left", weight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#E4E7EC", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_xlim(0, max(1, ordered["overlap_with_sencat"].max() * 1.25))
    ax.text(
        0,
        -0.28,
        "Aggregate source overlap only — not an expression result",
        transform=ax.transAxes,
        fontsize=9,
        color="#B42318",
    )

    table_ax.axis("off")
    table_data = [
        [
            row["label"],
            f"{int(row['mapped_mouse_genes']):,}",
            f"{row['fraction_of_axis_overlapping_sencat']:.1%}",
        ]
        for _, row in ordered.iterrows()
    ]
    table = table_ax.table(
        cellText=table_data,
        colLabels=["Axis", "Mapped", "% in SenCat"],
        loc="center",
        cellLoc="left",
        colLoc="left",
        colWidths=[0.50, 0.24, 0.31],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.8)
    table.scale(1.0, 1.55)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor("#EAECF0")
        cell.set_facecolor("#F9FAFB" if row == 0 else "white")
        if row == 0:
            cell.set_text_props(weight="bold")
    table_ax.set_title(
        "One-to-one human→mouse mapping",
        loc="left",
        fontsize=11,
        weight="bold",
        pad=14,
    )
    fig.suptitle(
        "Frozen axes before matched-atlas access",
        x=0.08,
        ha="left",
        fontsize=15,
        weight="bold",
    )
    fig.tight_layout()
    _save(fig, output_stem)


def plot_animal_results(
    aggregates: pd.DataFrame,
    inference: pd.DataFrame,
    score: str,
    output_stem: Path,
) -> None:
    subsets = list(aggregates.groupby(["depot", "subtype"]))
    if not subsets:
        raise ValueError("No animal aggregates to plot")
    ncols = min(3, len(subsets))
    nrows = int(np.ceil(len(subsets) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 3.4 * nrows),
        squeeze=False,
        sharex=False,
    )
    for ax, ((depot, subtype), frame) in zip(axes.flat, subsets):
        ax.scatter(
            frame["age_months"],
            frame[score],
            color=COLORS["sencat"] if score.startswith("sencat") else COLORS["ifng_response"],
            alpha=0.85,
            edgecolor="white",
            linewidth=0.6,
            s=38,
        )
        if frame["age_months"].nunique() > 1:
            coefficients = np.polyfit(frame["age_months"], frame[score], 1)
            x = np.linspace(frame["age_months"].min(), frame["age_months"].max(), 100)
            ax.plot(x, coefficients[0] * x + coefficients[1], color="#344054")
        row = inference[
            (inference["depot"] == depot)
            & (inference["subtype"] == subtype)
            & (inference["score"] == score)
        ]
        subtitle = ""
        if len(row) == 1:
            result = row.iloc[0]
            subtitle = (
                f"\nβ={result['slope_per_month']:.3g}, "
                f"95% CI [{result['slope_ci_low']:.3g}, {result['slope_ci_high']:.3g}], "
                f"q={result['q_value']:.3g}"
            )
        ax.set_title(f"{depot} · {subtype}{subtitle}", fontsize=9.5, loc="left")
        ax.set_xlabel("Age (months)")
        ax.set_ylabel(score)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(color="#EAECF0", linewidth=0.7)
    for ax in axes.flat[len(subsets) :]:
        ax.axis("off")
    fig.suptitle(
        f"Animal-level {score} by age",
        x=0.02,
        ha="left",
        fontsize=14,
        weight="bold",
    )
    fig.tight_layout()
    _save(fig, output_stem)
