# Data request

The minimum sufficient package is the longitudinal ScaleBio single-nucleus
cohort from Xie et al., DOI `10.64898/2026.05.13.724992`.

## Required

1. Raw or unlogged count matrix, preferably sparse `.h5ad`, `.h5`, or Matrix
   Market files.
2. One metadata row per nucleus with:
   - stable cell/nucleus ID matching the matrix;
   - de-identified animal ID;
   - age in months or exact age-group label;
   - depot (`SAT`/inguinal subcutaneous or `VAT`/perigonadal visceral);
   - sex;
   - broad cell type and final cell subtype;
   - sample/library/batch and final QC-pass flag.
3. The vascular/endothelial subtype annotation used for the manuscript
   figures.
4. A sample sheet mapping library/sample IDs to animal, age, depot, and batch.
5. Gene identifiers and reference build used for the delivered matrix.

Animal IDs may be arbitrary de-identified codes. No personal or sensitive data
are requested.

## Helpful but not required

- the authors' pseudobulk counts;
- the final Scanpy object with embeddings;
- annotation scripts or exact marker lists;
- the list of excluded doublets/low-quality nuclei;
- the sample-to-well mapping for the ScaleBio split-pool workflow.

## Why aggregate output is insufficient

The planned inference treats animals—not nuclei—as replicates and compares
overlap-depleted gene programs. Group means, plotted curves, or a matrix pooled
across animals cannot support valid uncertainty estimates or separate the
SenCat, IFN-γ, MHC-II, and arrest axes.

## Secure transfer

A private, expiring institutional download or author-chosen secure transfer is
preferred. Raw files will remain local and excluded from Git. Only code,
source/checksum metadata, aggregate statistics, and original figures will be
made public, subject to any data-use terms the authors specify.

