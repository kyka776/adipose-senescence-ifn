# Decision log

## 2026-07-29 — Project boundary

- Created a separate project under `GeroScout/initiatives/`; no Phase 1
  collector, snapshot, queue, report, or finalist brief is modified.
- The Xie et al. single-nucleus atlas is the only eligible primary dataset.
  An unrelated public adipose atlas will not be substituted.

## 2026-07-29 — Gate 0 status

- Exact DOI/title searches found no downloadable matched atlas or code package
  in the paper page, GEO/SRA, GitHub, CELLxGENE, Single Cell Portal, Zenodo, or
  Figshare.
- The paper describes animal-level replication and a sample-aware pseudobulk,
  so the design is eligible in principle.
- Status is `CONTACT_REQUIRED`, not permanent no-go. Analysis remains blocked
  until files and animal identifiers arrive and pass `audit-data`.

## 2026-07-29 — Frozen axes

- Use the pinned SenCat transcriptomic coefficients.
- Use Reactome v97 pathways for IFN-γ, MHC-II, and TP53-mediated G1/G2 arrest.
- Convert human symbols with one-to-one MGI mouse–human homology only.
- Publish aggregate overlap diagnostics but not the unlicensed SenCat table or
  derived detailed SenCat gene lists.

## 2026-07-29 — Inference unit

- Cells are measurements; animals are biological replicates.
- All confirmatory tests use animal × depot × vascular subtype aggregates.
- Age labels are permuted at the animal level and kept aligned across depots
  and subtypes.

