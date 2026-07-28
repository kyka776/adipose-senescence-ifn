"""Fetch, validate, map, and summarize the frozen program signatures."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import tempfile
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


USER_AGENT = "adipose-senescence-ifn/0.1 (+https://github.com/kyka776/adipose-senescence-ifn)"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _download(url: str, destination: Path, expected_sha256: str | None = None) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json, text/plain, */*"},
    )
    fd, temporary_name = tempfile.mkstemp(prefix=".download-", dir=destination.parent)
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as out:
            while block := response.read(1024 * 1024):
                out.write(block)
        observed = sha256_file(temporary)
        if expected_sha256 and observed != expected_sha256:
            raise ValueError(
                f"Checksum mismatch for {url}: expected {expected_sha256}, observed {observed}"
            )
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return {
        "url": url,
        "path": str(destination),
        "sha256": sha256_file(destination),
        "bytes": destination.stat().st_size,
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }


def _ensure_download(
    url: str, destination: Path, expected_sha256: str | None = None
) -> dict:
    if destination.exists():
        observed = sha256_file(destination)
        if expected_sha256 and observed != expected_sha256:
            raise ValueError(
                f"Cached checksum mismatch for {destination}: "
                f"expected {expected_sha256}, observed {observed}"
            )
        return {
            "url": url,
            "path": str(destination),
            "sha256": observed,
            "bytes": destination.stat().st_size,
            "fetched_at": None,
            "cache_hit": True,
        }
    return _download(url, destination, expected_sha256)


def read_sencat(path: Path) -> dict[str, float]:
    weights: dict[str, float] = {}
    source_rows = 0
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["converted_alias", "gene", "coef"]:
            raise ValueError(f"Unexpected SenCat columns: {reader.fieldnames}")
        for row in reader:
            source_rows += 1
            symbol = row["gene"].strip().upper()
            coefficient = float(row["coef"])
            # The pinned source uses the literal "0" for 21 Ensembl features
            # that have no converted gene alias. They cannot be mapped by
            # one-to-one symbol homology and are explicitly excluded.
            if symbol in {"", "0", "NA", "N/A"}:
                continue
            if not symbol:
                raise ValueError("SenCat contains an empty gene symbol")
            if symbol in weights:
                raise ValueError(f"Duplicate SenCat gene symbol: {symbol}")
            weights[symbol] = coefficient
    if source_rows != 5000:
        raise ValueError(f"Expected 5,000 SenCat rows, observed {source_rows:,}")
    if not any(value > 0 for value in weights.values()) or not any(
        value < 0 for value in weights.values()
    ):
        raise ValueError("SenCat coefficients must contain both signs")
    return weights


def read_mgi_one_to_one(path: Path) -> tuple[dict[str, str], dict]:
    groups: dict[str, dict[str, list[str]]] = defaultdict(
        lambda: {"human": [], "mouse": []}
    )
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"DB Class Key", "NCBI Taxon ID", "Symbol"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Unexpected MGI homology report columns")
        for row in reader:
            taxon = row["NCBI Taxon ID"].strip()
            if taxon == "9606":
                groups[row["DB Class Key"]]["human"].append(row["Symbol"].upper())
            elif taxon == "10090":
                groups[row["DB Class Key"]]["mouse"].append(row["Symbol"])

    mapping: dict[str, str] = {}
    ambiguous = 0
    for members in groups.values():
        human = sorted(set(members["human"]))
        mouse = sorted(set(members["mouse"]))
        if len(human) == 1 and len(mouse) == 1:
            mapping[human[0]] = mouse[0]
        elif human and mouse:
            ambiguous += 1
    return mapping, {
        "homology_classes": len(groups),
        "one_to_one_pairs": len(mapping),
        "ambiguous_cross_species_classes": ambiguous,
    }


def _reactome_symbols(payload: object) -> set[str]:
    if not isinstance(payload, list):
        raise ValueError("Reactome participant response must be a list")
    symbols: set[str] = set()
    for participant in payload:
        if not isinstance(participant, dict):
            continue
        for entity in participant.get("refEntities", []):
            if entity.get("schemaClass") not in {
                "ReferenceGeneProduct",
                "ReferenceIsoform",
                "ReferenceDNASequence",
                "ReferenceRNASequence",
            }:
                continue
            display = str(entity.get("displayName", "")).strip()
            if " " not in display:
                continue
            symbol = display.rsplit(" ", 1)[-1].strip().upper()
            if symbol and symbol.replace("-", "").replace(".", "").isalnum():
                symbols.add(symbol)
    if not symbols:
        raise ValueError("No gene symbols parsed from Reactome participant response")
    return symbols


def _map_weights(
    weights: dict[str, float], homology: dict[str, str]
) -> tuple[dict[str, float], list[str]]:
    mapped: dict[str, float] = {}
    missing: list[str] = []
    for human_symbol, coefficient in weights.items():
        mouse_symbol = homology.get(human_symbol)
        if mouse_symbol is None:
            missing.append(human_symbol)
            continue
        if mouse_symbol in mapped:
            raise ValueError(f"Non-unique mouse target after one-to-one mapping: {mouse_symbol}")
        mapped[mouse_symbol] = coefficient
    return mapped, missing


def _map_set(
    genes: Iterable[str], homology: dict[str, str]
) -> tuple[set[str], list[str]]:
    mapped: set[str] = set()
    missing: list[str] = []
    for human_symbol in sorted(set(genes)):
        mouse_symbol = homology.get(human_symbol)
        if mouse_symbol is None:
            missing.append(human_symbol)
        else:
            mapped.add(mouse_symbol)
    return mapped, missing


def build_signatures(
    spec_path: Path,
    cache_dir: Path,
    private_out: Path,
    public_out: Path,
) -> dict:
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    private_out.mkdir(parents=True, exist_ok=True)
    public_out.mkdir(parents=True, exist_ok=True)

    sources: list[dict] = []
    sencat = spec["sencat"]
    sencat_url = (
        f"{sencat['repository'].replace('github.com', 'raw.githubusercontent.com')}/"
        f"{sencat['commit']}/{sencat['path']}"
    )
    sencat_path = cache_dir / "sencat-transcriptomic.csv"
    sources.append(
        _ensure_download(sencat_url, sencat_path, expected_sha256=sencat["sha256"])
    )

    mgi = spec["homology"]
    mgi_path = cache_dir / "HOM_MouseHumanSequence.rpt"
    sources.append(
        _ensure_download(mgi["url"], mgi_path, expected_sha256=mgi["sha256"])
    )

    release_url = "https://reactome.org/ContentService/data/database/version"
    release_path = cache_dir / "reactome-version.txt"
    sources.append(_download(release_url, release_path))
    observed_release = int(release_path.read_text(encoding="utf-8").strip())
    expected_release = int(spec["reactome"]["release"])
    if observed_release != expected_release:
        raise ValueError(
            f"Reactome release drift: spec freezes {expected_release}, API reports "
            f"{observed_release}. Update only with a documented review."
        )

    human_axes: dict[str, set[str]] = {}
    for axis, pathway_ids in spec["reactome"]["pathways"].items():
        axis_symbols: set[str] = set()
        for pathway_id in pathway_ids:
            url = f"https://reactome.org/ContentService/data/participants/{pathway_id}"
            path = cache_dir / f"reactome-{pathway_id}.json"
            sources.append(_download(url, path))
            axis_symbols.update(
                _reactome_symbols(json.loads(path.read_text(encoding="utf-8")))
            )
        human_axes[axis] = axis_symbols

    human_sencat = read_sencat(sencat_path)
    homology, homology_stats = read_mgi_one_to_one(mgi_path)
    mouse_sencat, sencat_missing = _map_weights(human_sencat, homology)
    mouse_axes: dict[str, set[str]] = {}
    missing_by_axis: dict[str, list[str]] = {}
    for axis, human_genes in human_axes.items():
        mapped, missing = _map_set(human_genes, homology)
        mouse_axes[axis] = mapped
        missing_by_axis[axis] = missing

    depletion = spec["depletion"]
    variants: dict[str, dict[str, float]] = {"sencat_full": mouse_sencat}
    for variant_name, axes in depletion.items():
        excluded = set().union(*(mouse_axes[axis] for axis in axes))
        variants[variant_name] = {
            gene: coefficient
            for gene, coefficient in mouse_sencat.items()
            if gene not in excluded
        }

    private_payload = {
        "schema_version": "1.0",
        "target_species": spec["target_species"],
        "sencat_weights": mouse_sencat,
        "binary_axes": {name: sorted(genes) for name, genes in mouse_axes.items()},
        "sencat_variants": variants,
        "source_manifest": [
            {
                "url": source["url"],
                "sha256": source["sha256"],
                "bytes": source["bytes"],
            }
            for source in sources
        ],
    }
    private_path = cache_dir / "mouse_signatures.json"
    private_path.write_text(
        json.dumps(private_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    detailed_path = private_out / "signature_membership.tsv"
    all_genes = sorted(set(mouse_sencat).union(*mouse_axes.values()))
    with detailed_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["mouse_gene", "in_sencat", *[f"in_{axis}" for axis in mouse_axes]]
        )
        for gene in all_genes:
            writer.writerow(
                [
                    gene,
                    int(gene in mouse_sencat),
                    *[int(gene in mouse_axes[axis]) for axis in mouse_axes],
                ]
            )

    rows: list[dict] = []
    for axis, genes in mouse_axes.items():
        overlap = genes & set(mouse_sencat)
        rows.append(
            {
                "axis": axis,
                "source_human_genes": len(human_axes[axis]),
                "mapped_mouse_genes": len(genes),
                "unmapped_human_genes": len(missing_by_axis[axis]),
                "overlap_with_sencat": len(overlap),
                "fraction_of_axis_overlapping_sencat": (
                    len(overlap) / len(genes) if genes else 0.0
                ),
            }
        )

    summary_tsv = public_out / "signature_overlap_summary.tsv"
    with summary_tsv.open("w", newline="", encoding="utf-8") as handle:
        fields = list(rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    variants_tsv = public_out / "sencat_variant_counts.tsv"
    with variants_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["variant", "mapped_mouse_genes", "genes_removed"])
        full_count = len(mouse_sencat)
        for name, genes in variants.items():
            writer.writerow([name, len(genes), full_count - len(genes)])

    public_manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "reactome_release": observed_release,
        "sencat_commit": sencat["commit"],
        "sencat_source_sha256": sencat["sha256"],
        "mgi_source_sha256": sha256_file(mgi_path),
        "source_sha256": {
            source["url"]: source["sha256"]
            for source in sources
        },
        "mouse_signature_bundle_sha256": sha256_file(private_path),
        "counts": {
            "sencat_human": len(human_sencat),
            "sencat_mouse_one_to_one": len(mouse_sencat),
            "sencat_unmapped": len(sencat_missing),
            "axes": {row["axis"]: row for row in rows},
            "variants": {name: len(genes) for name, genes in variants.items()},
            "homology": homology_stats,
        },
        "redistribution": {
            "detailed_sencat_material": False,
            "public_outputs_are_aggregate_counts_only": True,
        },
    }
    manifest_path = public_out / "signature_build_manifest.json"
    manifest_path.write_text(
        json.dumps(public_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return {
        "manifest": public_manifest,
        "rows": rows,
        "variants": {name: len(genes) for name, genes in variants.items()},
        "private_signature_path": private_path,
    }
