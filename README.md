# Adipose Senescence–IFN Decomposition

Status: **waiting for matched atlas access (`CONTACT_REQUIRED`)**.

This is an independent, reproducible implementation of Initiative 07 from
GeroScout Phase 2. It asks whether vascular MHC-II/IFN-γ activity in aging
mouse adipose tissue reflects a broader senescence program or a separable
inflammatory axis. The target dataset is the animal-resolved single-nucleus
atlas described by Xie et al. (2026), not a substitute atlas.

No biological age effect is reported here yet. The matched expression matrix,
cell metadata, and animal identifiers were not publicly downloadable when the
repository was prepared. The pre-data results are limited to signature
provenance, overlap/depletion diagnostics, and pipeline validation.

## What is complete

- a dated availability audit and explicit stop/go decision;
- frozen definitions for SenCat, IFN-γ, MHC-II, and cell-cycle-arrest axes;
- legal/provenance controls that keep unlicensed SenCat material out of Git;
- an animal-aware data gate, scoring, aggregation, permutation inference,
  falsification logic, and deterministic figures;
- unit and synthetic integration tests;
- an executed pre-data notebook and validation report;
- a precise data request and author email.

## Verified pre-data result

The pinned SenCat source has 4,979 usable human gene aliases; 4,630 map through
one-to-one MGI homology. Of the mapped Reactome axes, 31/72 IFN-γ genes,
59/108 MHC-II genes, and 16/29 arrest genes overlap SenCat. Removing their
union leaves a 4,524-gene residual SenCat variant.

These are source-overlap counts, not expression or age effects. See the
[executed notebook](notebooks/00_signature_overlap.ipynb) and
[figure](results/predata/figure-02-signature-overlap.svg).

## Reproduce the verified pre-data state

```bash
conda env create -f environment.yml
conda activate adipose-senescence-ifn
make all
```

`make predata` downloads SenCat only from a pinned Git commit and verifies its
SHA-256. It also enforces the frozen Reactome release and MGI report checksum.
The downloaded tables and all detailed derived gene lists remain in ignored
directories. Aggregate overlap counts are safe to publish.

## Run after access

Place the author-provided files under `data/raw/` without committing them:

```bash
make analysis
```

The command first enforces Gate 0. It refuses to score or test anything if
animal IDs, age, depot, vascular subtype, or animal-level replication are
missing. See [data-request.md](data-request.md),
[docs/data-contract.md](docs/data-contract.md), and
[analysis-plan.md](analysis-plan.md).

## Interpretation boundary

Cells are the measurement layer; animals are the inference unit. A cell-level
P value is never produced. Any eventual result must show individual animal
points and both full and overlap-depleted SenCat scores. Until the matched
atlas passes Gate 0, the honest conclusion is “not tested.”

## Sources and licensing

Our code is MIT licensed. Reactome pathway data are used under CC0 with
attribution; MGI homology data are used under CC BY 4.0. The SenCat repository
has no license as of the pinned commit, so neither its code nor marker table is
redistributed. The Xie et al. preprint is all-rights-reserved and is linked,
not copied. Details and checksums are in
[signatures/provenance.md](signatures/provenance.md).
