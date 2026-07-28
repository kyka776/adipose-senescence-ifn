# Input data contract

## Metadata TSV

UTF-8, tab-separated, one row per cell, unique `cell_id`.

| Column | Type | Allowed/example |
|---|---|---|
| `cell_id` | string | unique and identical to expression row index |
| `animal_id` | string | de-identified biological replicate |
| `age_months` | number | positive; constant within animal |
| `depot` | category | `SAT` or `VAT` after explicit mapping |
| `sex` | category | `male`, `female`, or supplied controlled value |
| `cell_type` | string | broad author annotation |
| `subtype` | string | final author subtype |
| `sample_id` | string | library/sample identifier |
| `qc_pass` | boolean | final author QC decision |

Extra columns are preserved. Missing required values fail Gate 0.

## Expression

The portable test interface accepts a tab-separated matrix (optionally gzip):
rows are `cell_id`, columns are gene symbols, values are non-negative raw
counts. For the full atlas use the optional `.h5ad` adapter to avoid a dense
file; row and gene identifiers must be unique.

The loader rejects negative values, duplicate identifiers, mismatched cells,
and matrices that appear centered or z-scored. Count-like non-integers require
an explicit decision-log entry.

## Gate report

`audit-data` writes machine-readable JSON and Markdown with:

- file hashes and dimensions;
- missing/duplicate IDs;
- animals per age × depot;
- consistency of animal-level age and sex;
- cells and animals per subtype;
- pass/fail for every gate clause.

`analyze` consumes the report and refuses to continue unless `passed=true`.

