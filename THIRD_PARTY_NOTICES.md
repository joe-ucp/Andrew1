# Dataset and third-party attribution

## DCLP3, REPLACE-BG and CLVer

The Jaeb Center provides [DCLP3 release 3](https://public.jaeb.org/dataset/573), [REPLACE-BG](https://public.jaeb.org/dataset/546), and [CLVer](https://public.jaeb.org/dataset/591). Acquisition requires registration and acceptance of the provider's terms. Preserve the ReadMe and study-specific attribution supplied with each archive. This repository does not redistribute the original archives or participant-level derivatives and does not infer a blanket redistribution grant from public availability.

DCLP3 and REPLACE-BG analyses, content and conclusions are solely the responsibility of UCP Technology LLC and have not been reviewed or approved by the Jaeb Center, iDCL Trial Research Group or T1D Exchange study group.

Required CLVer statement: “The source of the data is the CLVer: Hybrid Closed Loop Therapy and Verapamil for Beta Cell Preservation in New Onset Type 1 Diabetes: A Proof-of-Concept Study, but the analyses, content and conclusions presented herein are solely the responsibility of the authors and have not been reviewed or approved by the Jaeb Center for Health Research.”

## TN10 expression and metadata

Lemaitre P, Mathieu C, Gysemans C, Type 1 Diabetes TrialNet. *Whole Blood RNAsequencing-Teplizumab for Prevention of Type 1 Diabetes in Relatives At-Risk (TN10)*, version 1 (2025). [Zenodo record and DOI](https://doi.org/10.5281/zenodo.17063348). Provider license: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/), confirmed in the record's license metadata.

The analyses filter baseline libraries, transform counts into CPM and log expression, construct scores and aggregate outcomes. Generated TN10 files are adaptations or analytical outputs described in the manuscript; the provider's attribution and share-alike requirements must accompany adapted source material where applicable. The original CSVs are acquired separately.

## Fixed AbATE gene signature

Sassi et al. *Neutrophil-enriched gene signature correlates with teplizumab therapy resistance in different stages of type 1 diabetes*. J Clin Invest. 2025;135:e176403. [Article](https://www.jci.org/articles/view/176403), [supplementary Table 6](https://www.jci.org/articles/view/176403/sd/6). Copyright © 2025 Sassi et al.; [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

The scientific implementation extracts the fixed 359-gene negative-correlation list and retains 312 exactly matching TN10 symbols. The included `evidence/card002/results/primary_genes.txt` is this transformed symbol list. Extraction, matching, score construction and interpretation are by UCP Technology LLC; no endorsement by source authors is implied.

## Dependencies and original work

Python packages in [requirements-lock.txt](requirements-lock.txt) retain their respective licenses. This repository includes the project's scientific source; installed dependency packages are obtained separately. No proprietary fonts or external document-rendering binaries are bundled. Original code and manuscript permissions are recorded separately in [LICENSE_STATUS.md](LICENSE_STATUS.md).
