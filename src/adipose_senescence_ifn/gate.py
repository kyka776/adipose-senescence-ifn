"""Study data Gate 0 checks."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = [
    "cell_id",
    "animal_id",
    "age_months",
    "depot",
    "sex",
    "cell_type",
    "subtype",
    "sample_id",
    "qc_pass",
]


@dataclass
class GateReport:
    passed: bool
    checks: list[dict] = field(default_factory=list)
    summary: dict = field(default_factory=dict)

    def add(self, name: str, passed: bool, detail: str) -> None:
        self.checks.append({"name": name, "passed": bool(passed), "detail": detail})
        self.passed = self.passed and bool(passed)

    def as_dict(self) -> dict:
        return {
            "schema_version": "1.0",
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "passed": self.passed,
            "checks": self.checks,
            "summary": self.summary,
        }


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_metadata(
    metadata_path: Path, minimum_animals_per_age_depot: int = 3
) -> tuple[GateReport, pd.DataFrame | None]:
    report = GateReport(passed=True)
    if not metadata_path.exists():
        report.add("metadata_exists", False, f"Missing: {metadata_path}")
        return report, None

    try:
        metadata = pd.read_csv(metadata_path, sep="\t", dtype={"cell_id": str})
    except Exception as exc:
        report.add("metadata_readable", False, f"{type(exc).__name__}: {exc}")
        return report, None

    report.add("metadata_readable", True, f"{len(metadata):,} rows")
    missing_columns = [c for c in REQUIRED_COLUMNS if c not in metadata.columns]
    report.add(
        "required_columns",
        not missing_columns,
        "all present" if not missing_columns else f"missing: {', '.join(missing_columns)}",
    )
    if missing_columns:
        return report, metadata

    missing_values = metadata[REQUIRED_COLUMNS].isna().sum()
    missing_total = int(missing_values.sum())
    report.add(
        "required_values_complete",
        missing_total == 0,
        "no missing values"
        if missing_total == 0
        else "; ".join(f"{k}={v}" for k, v in missing_values.items() if v),
    )

    duplicate_cells = int(metadata["cell_id"].duplicated().sum())
    report.add(
        "cell_ids_unique",
        duplicate_cells == 0,
        f"{duplicate_cells} duplicate cell IDs",
    )

    numeric_age = pd.to_numeric(metadata["age_months"], errors="coerce")
    age_valid = numeric_age.notna().all() and bool((numeric_age > 0).all())
    report.add("age_months_valid", age_valid, f"ages={sorted(numeric_age.dropna().unique())}")
    metadata = metadata.assign(age_months=numeric_age)

    normalized_depot = metadata["depot"].astype(str).str.upper()
    depot_valid = set(normalized_depot.unique()).issubset({"SAT", "VAT"})
    report.add(
        "depot_values",
        depot_valid,
        f"observed={sorted(normalized_depot.unique())}",
    )
    metadata = metadata.assign(depot=normalized_depot)

    for column in ("age_months", "sex"):
        consistency = metadata.groupby("animal_id", dropna=False)[column].nunique()
        inconsistent = consistency[consistency != 1]
        report.add(
            f"animal_{column}_consistent",
            inconsistent.empty,
            "constant within animal"
            if inconsistent.empty
            else f"inconsistent animals={list(inconsistent.index[:10])}",
        )

    qc = metadata["qc_pass"]
    if qc.dtype != bool:
        normalized_qc = qc.astype(str).str.lower().map(
            {"true": True, "false": False, "1": True, "0": False}
        )
    else:
        normalized_qc = qc
    qc_valid = normalized_qc.notna().all()
    report.add("qc_pass_boolean", qc_valid, f"valid={int(normalized_qc.notna().sum()):,}")
    metadata = metadata.assign(qc_pass=normalized_qc)
    eligible = metadata[metadata["qc_pass"] == True].copy()  # noqa: E712

    animal_counts = (
        eligible.groupby(["age_months", "depot"])["animal_id"]
        .nunique()
        .sort_index()
    )
    replication_ok = (
        len(animal_counts) > 0
        and bool((animal_counts >= minimum_animals_per_age_depot).all())
    )
    count_detail = ", ".join(
        f"{age:g}m/{depot}={count}"
        for (age, depot), count in animal_counts.items()
    )
    report.add(
        "animal_level_replication",
        replication_ok,
        count_detail or "no eligible age × depot groups",
    )

    cell_types = eligible["cell_type"].astype(str)
    vascular = cell_types.str.contains("endothelial|vascular", case=False, regex=True)
    report.add(
        "vascular_cells_identifiable",
        bool(vascular.any()),
        f"{int(vascular.sum()):,} matching cells",
    )

    animals_per_subtype = (
        eligible.loc[vascular]
        .groupby(["depot", "subtype"])["animal_id"]
        .nunique()
        .sort_values(ascending=False)
    )
    report.summary = {
        "metadata_sha256": _hash_file(metadata_path),
        "cells_total": len(metadata),
        "cells_qc_pass": len(eligible),
        "animals": int(eligible["animal_id"].nunique()),
        "ages": sorted(float(x) for x in eligible["age_months"].dropna().unique()),
        "animals_per_age_depot": {
            f"{age:g}m|{depot}": int(count)
            for (age, depot), count in animal_counts.items()
        },
        "animals_per_vascular_subtype": {
            f"{depot}|{subtype}": int(count)
            for (depot, subtype), count in animals_per_subtype.items()
        },
    }
    return report, metadata


def write_gate_report(report: GateReport, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = report.as_dict()
    (output_dir / "gate-report.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    lines = [
        "# Gate 0 report",
        "",
        f"Decision: **{'GO' if report.passed else 'STOP'}**",
        "",
        "| Check | Result | Detail |",
        "|---|---:|---|",
    ]
    for check in report.checks:
        detail = check["detail"].replace("|", "\\|")
        lines.append(
            f"| `{check['name']}` | {'PASS' if check['passed'] else 'FAIL'} | {detail} |"
        )
    lines.extend(["", "```json", json.dumps(report.summary, indent=2), "```", ""])
    (output_dir / "gate-report.md").write_text("\n".join(lines), encoding="utf-8")

