# Verification evidence

Verified on 2026-07-29 before publication.

## Software

- 10/10 unit and synthetic integration tests passed.
- The synthetic end-to-end test exercises Gate 0, score coverage, cell scoring,
  animal aggregation, animal-label permutation, bootstrap uncertainty,
  multiplicity correction, and figure export.
- Python source compiled without syntax errors.
- `git diff --check` reported no whitespace errors.

## Source-backed pre-data build

- SenCat source checksum matched
  `699dd16ee3b956b4b9442c22a36b03a82b483d80b93036d54fae1ea3589cb093`.
- MGI report checksum matched
  `4823cbf8ea7961d17419c891c56a2474b643515333f99492b9d2ae3107f335ea`.
- Reactome Content Service reported frozen release 97.
- Two consecutive builds produced byte-identical TSV tables, private signature
  bundle, PNG figures, and—after fixed SVG metadata—SVG figures.
- The notebook executed top-to-bottom with saved outputs.
- Both PNG figures were visually inspected; labels and values are legible and
  the pre-data/non-biological boundary is explicit.

## Release safety

- No raw study data are present.
- Downloaded SenCat, MGI, and Reactome files are Git-ignored.
- Detailed mapped signature membership is Git-ignored.
- The public SenCat-derived artifacts contain only aggregate counts and
  checksums.
- Repository scan found no candidate secrets.
- The all-rights-reserved preprint is linked and paraphrased, not redistributed.

## Required negative control

With no matched metadata under `data/raw/`, `make analysis` exits non-zero at
Gate 0 and writes only a private gate report. No expression scoring or
biological result is produced.

