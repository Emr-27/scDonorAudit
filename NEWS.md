# scDonorAudit 0.99.5

- Reject logical and pattern Matrix counts with an INVALID_COUNTS error,
  before accessing numeric slots. Numeric sparse counts stay sparse.
- Return the same 16 typed summary columns when no genes are retained.
- Split design planning and fit execution into internal helpers while keeping
  the public interfaces, fixed filtering family and numerical backend intact.
- Add package help, a complete output/status dictionary, synthetic data
  provenance, and a runnable mixed coverage/design-failure example.
- Clarify the in-memory count requirement and preaggregated entry point,
  method references, functional positioning and AI-assisted development.

# scDonorAudit 0.99.4

- Revalidate prepared counts and gene IDs at both public model entries.
  Invalid edits fail before auditing or fitting; valid count edits, row
  subsets, and column reordering remain supported.
- Use sample IDs for coverage plot positions, keeping donor/sample display
  labels separate so distinct samples cannot overlap due to label collisions.
- Add dense and sparse prepared-object mutation regression tests.

# scDonorAudit 0.99.3

- Reject reserved covariate names before they can overwrite model columns.
- Require syntactic covariate names and reject the special formula term `.`.
- Use unique run IDs in effect plots, including for donors named `baseline`,
  and preserve positions for unavailable effects.
- Add regression tests for both public model entries and plotting boundaries.

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
