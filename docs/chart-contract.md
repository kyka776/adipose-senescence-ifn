# Figure contracts

## Figure 1 — Pre-data decomposition

- Question: which program components are explicitly separated before data
  access?
- Takeaway: the primary test compares full SenCat with overlap-depleted
  variants at the animal level.
- Data: protocol definitions only; no biological measurements.
- Form: original schematic SVG.
- Mandatory label: “PRE-DATA DESIGN — NOT A BIOLOGICAL RESULT”.

## Figure 2 — Signature overlap audit

- Question: how many mapped mouse genes are shared between SenCat and each
  inflammatory/arrest axis?
- Takeaway: overlap is measurable and motivates depleted variants.
- Data: pinned SenCat coefficients, Reactome v97 participants, one-to-one MGI
  homology.
- Form: horizontal bar chart plus coverage table.
- Output: SVG and PNG, deterministic dimensions, values printed directly.
- QA: counts must equal the public TSV; axes start at zero; no inferred
  biological effect language.

## Post-access animal plots

- Question: how does each animal aggregate change with age within depot and
  subtype?
- Form: small multiples with one point per animal, a fitted slope and 95%
  animal bootstrap interval.
- Encoding: age on x; score on y; depot separated into panels, not only color;
  score variants use a colorblind-safe palette and direct labels.
- Mandatory: N animals per group, effect and interval, permutation P/Q,
  signature coverage, and a visible “exploratory” label where applicable.
- Prohibited: violin plots of cells presented as inferential evidence.

