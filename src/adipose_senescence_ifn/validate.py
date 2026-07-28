"""Repository-level validation and public-release safety checks."""

from __future__ import annotations

import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


REQUIRED_PATHS = [
    "README.md",
    "availability-check.md",
    "analysis-plan.md",
    "data-request.md",
    "decision-log.md",
    "signatures/provenance.md",
    "signatures/spec.json",
    "outreach/author_email_draft.md",
    ".github/workflows/ci.yml",
    "notebooks/00_signature_overlap.ipynb",
    "results/predata/signature_overlap_summary.tsv",
    "results/predata/sencat_variant_counts.tsv",
    "results/predata/figure-01-design.svg",
    "results/predata/figure-01-design.png",
    "results/predata/figure-02-signature-overlap.svg",
    "results/predata/figure-02-signature-overlap.png",
    "results/verification.md",
]

SECRET_PATTERNS = [
    re.compile(r"AKIA[0-9A-Z]{16}"),
    re.compile(r"ghp_[A-Za-z0-9]{30,}"),
    re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    re.compile(r"(?i)(?:api[_-]?key|secret|token)\s*[:=]\s*['\"][^'\"]{12,}"),
]


def _tracked_files(root: Path) -> list[Path]:
    try:
        result = subprocess.run(
            ["git", "ls-files"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return [p for p in root.rglob("*") if p.is_file() and ".git" not in p.parts]
    return [root / line for line in result.stdout.splitlines() if line]


def validate_repository(root: Path) -> tuple[bool, list[dict]]:
    checks: list[dict] = []

    missing = [path for path in REQUIRED_PATHS if not (root / path).exists()]
    checks.append(
        {
            "name": "required_artifacts",
            "passed": not missing,
            "detail": "all present" if not missing else f"missing: {', '.join(missing)}",
        }
    )

    tracked = _tracked_files(root)
    prohibited = [
        str(path.relative_to(root))
        for path in tracked
        if str(path.relative_to(root)).startswith(
            ("data/raw/", "data/external/", "results/private/")
        )
        and path.name != ".gitkeep"
    ]
    checks.append(
        {
            "name": "no_private_or_raw_files_tracked",
            "passed": not prohibited,
            "detail": "none" if not prohibited else ", ".join(prohibited[:10]),
        }
    )

    secret_hits: list[str] = []
    for path in tracked:
        if not path.exists() or path.suffix.lower() in {".png", ".pdf", ".gz"}:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if any(pattern.search(text) for pattern in SECRET_PATTERNS):
            secret_hits.append(str(path.relative_to(root)))
    checks.append(
        {
            "name": "secret_scan",
            "passed": not secret_hits,
            "detail": "no candidate secrets" if not secret_hits else ", ".join(secret_hits),
        }
    )

    notebook_path = root / "notebooks/00_signature_overlap.ipynb"
    notebook_executed = False
    if notebook_path.exists():
        notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
        code_cells = [cell for cell in notebook.get("cells", []) if cell["cell_type"] == "code"]
        notebook_executed = bool(code_cells) and all(
            cell.get("execution_count") is not None for cell in code_cells
        )
    checks.append(
        {
            "name": "notebook_executed_top_to_bottom",
            "passed": notebook_executed,
            "detail": "all code cells executed" if notebook_executed else "unexecuted code cells",
        }
    )

    overlap_path = root / "results/predata/signature_overlap_summary.tsv"
    table_valid = False
    detail = "missing"
    if overlap_path.exists():
        try:
            table = pd.read_csv(overlap_path, sep="\t")
            expected_axes = {
                "ifng_response",
                "mhcii_antigen_presentation",
                "cell_cycle_arrest",
            }
            table_valid = set(table["axis"]) == expected_axes and bool(
                (table["overlap_with_sencat"] <= table["mapped_mouse_genes"]).all()
            )
            detail = f"{len(table)} axes; arithmetic {'valid' if table_valid else 'invalid'}"
        except Exception as exc:
            detail = f"{type(exc).__name__}: {exc}"
    checks.append(
        {"name": "overlap_summary_arithmetic", "passed": table_valid, "detail": detail}
    )

    manifest_path = root / "results/predata/signature_build_manifest.json"
    licensing_ok = False
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        licensing_ok = (
            manifest.get("redistribution", {}).get("detailed_sencat_material") is False
            and manifest.get("sencat_source_sha256")
            == "699dd16ee3b956b4b9442c22a36b03a82b483d80b93036d54fae1ea3589cb093"
        )
    checks.append(
        {
            "name": "unlicensed_sencat_not_redistributed",
            "passed": licensing_ok,
            "detail": "aggregate-only manifest and pinned source checksum"
            if licensing_ok
            else "manifest or redistribution flag invalid",
        }
    )

    passed = all(check["passed"] for check in checks)
    lines = [
        "# Validation report",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        f"Overall: **{'PASS' if passed else 'FAIL'}**",
        "",
        "| Check | Result | Detail |",
        "|---|---:|---|",
    ]
    for check in checks:
        escaped_detail = check["detail"].replace("|", "\\|")
        lines.append(
            f"| `{check['name']}` | {'PASS' if check['passed'] else 'FAIL'} | "
            f"{escaped_detail} |"
        )
    lines.extend(
        [
            "",
            "This report validates software/repository integrity and the pre-data "
            "signature audit. It does not validate a biological result because the "
            "matched atlas has not passed Gate 0.",
            "",
        ]
    )
    report_dir = root / "results"
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "validation-report.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    return passed, checks
