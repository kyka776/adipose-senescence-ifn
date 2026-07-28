# Signature provenance

Frozen on 2026-07-29. Human source symbols are mapped to mouse only through
one-to-one MGI homology classes. Ambiguous, missing, and many-to-many mappings
are excluded and counted.

| Axis | Frozen source | Definition | Redistribution |
|---|---|---|---|
| SenCat | maragkakislab repository, commit `6c16d47d…`, marker SHA-256 `699dd16e…` | 5,000 source rows; 4,979 usable gene aliases with signed transcriptomic coefficients | marker table and derived detailed lists stay ignored because no license was found |
| IFN-γ | Reactome v97, `R-HSA-877300` | unique human gene products participating in Interferon gamma signaling | aggregate/frozen derived set allowed; Reactome annotations are CC0 |
| MHC-II | Reactome v97, `R-HSA-2132295` | unique human gene products participating in MHC class II antigen presentation | aggregate/frozen derived set allowed; Reactome annotations are CC0 |
| Arrest | Reactome v97, union of `R-HSA-6804116` and `R-HSA-6804114` | TP53-regulated G1 and G2 cell-cycle-arrest gene products | aggregate/frozen derived set allowed; Reactome annotations are CC0 |
| Homology | MGI `HOM_MouseHumanSequence.rpt`, fetched 2026-07-29 | one human + one mouse member per homology class | MGI data are CC BY 4.0 |

The build records exact response hashes and the MGI database timestamp in a
private manifest. The public summary contains only source identifiers,
checksums, counts, and aggregate overlaps.

Twenty-one SenCat rows contain the literal alias `0`. They remain represented
by the pinned source checksum but are excluded before homology mapping because
`0` is not a gene symbol.

## Scoring contract

1. Use raw counts, normalize each cell to 10,000 counts, then `log1p`.
2. Standardize each included gene across all eligible cells after the primary
   QC and subtype filter.
3. Binary axes are the mean standardized expression of detected genes.
4. SenCat is the signed weighted mean, divided by the sum of absolute weights.
5. Produce full SenCat plus each predeclared overlap-depleted variant.
6. Require at least 10 detected genes and 20% coverage for every binary axis;
   require at least 500 mapped and detected SenCat genes.
7. Aggregate cell scores to animal × depot × vascular subtype before inference.

These choices are frozen before access to the target matrix. Any later change
requires a dated decision-log entry and sensitivity analysis.

## License links

- Reactome license: https://reactome.org/license
- MGI copyright and license: https://www.informatics.jax.org/mgihome/other/copyright.shtml
- SenCat source: https://github.com/maragkakislab/wf-ml-markers-senescence
