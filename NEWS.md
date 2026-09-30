# scDonorAudit 0.99.2

- Added the maintainer ORCID and aligned the minimum `testthat` version with
  the stable mocking interface used by the tests.
- Excluded `data-raw` from source builds while retaining its tracked files.
- Added a clear vignette error when the Crowell19 selection rule yields fewer
  than two genes, and removed an unsupported fixed effect-direction claim.

# scDonorAudit 0.99.1

- Preserved backend warnings when the same donor refit later fails, and kept
  each deletion's original design diagnosis when the baseline is unavailable.
- Added sample-aligned fit records and package/contrast provenance, plus consistency checks
  for pseudobulk objects whose columns or audit metadata were modified.
- Clarified the output dictionary, related work, missing effects, observed
  ranges, and fixed-versus-re-estimated analysis choices in the vignette.

# scDonorAudit 0.99.0

- Fixed factor contrast coding so user R settings cannot reverse or rescale
  the named condition effect; added regression coverage for full and
  donor-deletion fits.
- Expanded the Crowell19 vignette with a reproducible interpretation of
  per-donor effect changes and their limitations.
- Added donor-aware pseudobulk preparation, design audits, leave-one-donor-out
  edgeR quasi-likelihood fits, summaries, and plots for two-condition studies.
- Added explicit handling of incomplete pairs, missing cell types, low cell
  counts, unavailable fits, and genes that become all-zero after deletion.
- Added simulated and Crowell19 Astrocyte pseudobulk walkthroughs, with
  runnable examples for the exported APIs.
- Added Bioconductor installation guidance and a structured vignette.
