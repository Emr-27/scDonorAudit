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
