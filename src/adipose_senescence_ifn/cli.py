"""Command-line interface."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from .analysis import run_analysis
from .figures import plot_design_schematic, plot_overlap
from .gate import audit_metadata, write_gate_report
from .signatures import build_signatures
from .validate import validate_repository


def _build_signatures(args: argparse.Namespace) -> int:
    root = Path(__file__).resolve().parents[2]
    public_out = Path(args.public_out)
    result = build_signatures(
        spec_path=root / "signatures/spec.json",
        cache_dir=Path(args.cache),
        private_out=Path(args.private_out),
        public_out=public_out,
    )
    summary = pd.DataFrame(result["rows"])
    plot_design_schematic(public_out / "figure-01-design")
    plot_overlap(summary, public_out / "figure-02-signature-overlap")
    print(json.dumps(result["manifest"]["counts"], indent=2))
    return 0


def _audit(args: argparse.Namespace) -> int:
    report, _ = audit_metadata(
        Path(args.metadata),
        minimum_animals_per_age_depot=args.minimum_animals,
    )
    write_gate_report(report, Path(args.out))
    print("GO" if report.passed else "STOP")
    return 0 if report.passed else 2


def _analyze(args: argparse.Namespace) -> int:
    try:
        summary = run_analysis(
            Path(args.metadata),
            Path(args.expression),
            Path(args.signatures),
            Path(args.out),
            n_permutations=args.permutations,
            n_bootstrap=args.bootstrap,
        )
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    print(json.dumps(summary, indent=2))
    return 0


def _validate(args: argparse.Namespace) -> int:
    passed, checks = validate_repository(Path(args.root).resolve())
    for check in checks:
        print(f"{'PASS' if check['passed'] else 'FAIL'} {check['name']}: {check['detail']}")
    return 0 if passed else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="adipose-senescence-ifn",
        description="Animal-aware adipose senescence/IFN decomposition",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    signature_parser = subparsers.add_parser(
        "build-signatures", help="fetch pinned sources and build pre-data artifacts"
    )
    signature_parser.add_argument("--cache", required=True)
    signature_parser.add_argument("--private-out", required=True)
    signature_parser.add_argument("--public-out", required=True)
    signature_parser.set_defaults(func=_build_signatures)

    audit_parser = subparsers.add_parser(
        "audit-data", help="run Gate 0 without reading expression values"
    )
    audit_parser.add_argument("--metadata", required=True)
    audit_parser.add_argument("--out", required=True)
    audit_parser.add_argument("--minimum-animals", type=int, default=3)
    audit_parser.set_defaults(func=_audit)

    analysis_parser = subparsers.add_parser(
        "analyze", help="run the post-access pipeline after Gate 0"
    )
    analysis_parser.add_argument("--metadata", required=True)
    analysis_parser.add_argument("--expression", required=True)
    analysis_parser.add_argument("--signatures", required=True)
    analysis_parser.add_argument("--out", required=True)
    analysis_parser.add_argument("--permutations", type=int, default=9999)
    analysis_parser.add_argument("--bootstrap", type=int, default=2000)
    analysis_parser.set_defaults(func=_analyze)

    validate_parser = subparsers.add_parser(
        "validate-repository", help="check release safety and required artifacts"
    )
    validate_parser.add_argument("--root", default=".")
    validate_parser.set_defaults(func=_validate)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())

