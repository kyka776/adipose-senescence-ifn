# Matched-atlas availability check

Checked: **2026-07-29 (Europe/Moscow)**

Decision: **CONTACT_REQUIRED — proceed with author request; stop biological analysis**

Target:

- Qijing Xie et al., “Senescent cells induce vascular MHC II to recruit CD4+
  T cells and drive inflammation in aging adipose tissue”
- bioRxiv DOI: https://doi.org/10.64898/2026.05.13.724992
- version 1, posted 2026-05-16

## What the paper supports

The methods describe a longitudinal ScaleBio single-nucleus atlas from paired
inguinal subcutaneous and perigonadal visceral adipose tissue in male
C57BL/6J mice. The longitudinal cohort contains five age groups and four mice
per age group. Raw counts were summed by cell type and sample for the authors'
pseudobulk analysis. This is sufficient animal-level replication in principle.

The paper analyzes p16/p21 expression, IFN-γ response, and MHC-II biology, but
does not report SenCat scoring or the planned full-versus-overlap-depleted
SenCat decomposition. Therefore the proposed test is not a duplicate of the
reported analysis.

## Dated search log

| Location | Query/check | Result |
|---|---|---|
| bioRxiv article page and JATS metadata | supplementary material, data/code availability, accession strings | no downloadable expression/metadata package and no Data/Code Availability section |
| bioRxiv publication-status API | DOI lookup | no journal publication linked |
| NCBI GEO, SRA, PubMed | exact title and DOI; author/title combinations | no matched record/accession |
| GitHub | exact DOI, exact title, `Qijing Xie adipose`; author account repositories | no matched analysis repository |
| CELLxGENE and Single Cell Portal | exact title/DOI | no matched collection |
| Zenodo and Figshare | exact title/DOI | no matched deposit |
| Calico publications page | title/DOI | no downloadable dataset entry |

Negative search results are not proof that no private or newly deposited copy
exists. They establish that a reproducible public download was not locatable
on the check date.

## Gate 0

| Gate condition from the brief | Current evidence | Decision |
|---|---|---|
| no animal-level replication | paper states four mice per age group | not triggered |
| matched data cannot be obtained | no public download; authors not yet asked | unresolved |
| authors already performed the same decomposition | no SenCat/depletion analysis found | not triggered |
| only aggregate output, unable to separate axes | cell-level matrix described but not released | unresolved |

**Operational decision:** send the precise request in `data-request.md`. Keep the
analysis gate closed until an author response supplies a matrix plus
animal-resolved metadata. If authors decline, do not respond after two
documented follow-ups, or can provide only aggregate curves without animal IDs,
close the initiative as `NO_GO_DATA_UNAVAILABLE`.

## Prohibited substitution

No unrelated adipose or vascular atlas may be used to claim the target result.
Public atlases may only be used in a separately labelled methods smoke test
after licensing review; synthetic fixtures are used for current software tests.
