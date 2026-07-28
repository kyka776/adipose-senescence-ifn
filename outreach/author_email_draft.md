# Author email draft

To: Qijing Xie (`chichi@calicolabs.com`)

Cc: Cynthia Kenyon (`cynthia@calicolabs.com`)

From: `kyka776@gmail.com`

Subject: Request for animal-resolved snRNA-seq data from the aging adipose atlas

Dear Dr. Xie and Dr. Kenyon,

I am building a small, open, independently reproducible analysis inspired by
your bioRxiv preprint, “Senescent cells induce vascular MHC II to recruit CD4+
T cells and drive inflammation in aging adipose tissue”
(https://doi.org/10.64898/2026.05.13.724992).

The specific question is whether age-associated vascular MHC-II/IFN-γ activity
overlaps with, or is separable from, a broader transcriptomic senescence
signature. The analysis is preregistered to compare full and
overlap-depleted SenCat scores, with inference performed on animal × depot ×
vascular-subtype aggregates—not on nuclei as independent replicates.

Would you be willing to share the longitudinal ScaleBio count matrix and
de-identified cell/sample metadata? The minimum fields needed are animal ID,
age, SAT/VAT depot, sex, sample/library, broad cell type, and final vascular
subtype. An `.h5ad`, sparse Matrix Market package, or your preferred equivalent
would work. A precise field list is available here:
https://github.com/kyka776/adipose-senescence-ifn/blob/main/data-request.md

Raw data would remain private and excluded from Git. I would publish only the
analysis code, provenance/checksums, aggregate animal-level results, and
original figures, and I will follow any data-use or embargo terms you specify.
I would also be glad to share the result and manuscript-ready figures with you
before public interpretation.

Thank you for considering the request, and for the detailed experimental work.

Best regards,

Kyka
https://github.com/kyka776/adipose-senescence-ifn
