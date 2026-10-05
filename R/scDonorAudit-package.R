#' Audit donor omission in two-condition pseudobulk analyses
#'
#' scDonorAudit describes how a specified edgeR quasi-likelihood analysis
#' changes after omitting each biological donor. It supports independent
#' donors or paired donors with two conditions and optional sample covariates.
#' The output does not classify donors as invalid, recommend exclusion, or
#' provide a new false discovery rate guarantee.
#'
#' @section Workflow:
#' [preparePseudobulk()] aggregates cells or validates preaggregated counts.
#' [auditDesign()] records coverage, cohort choices and design viability.
#' [assessDonorInfluence()] fits the baseline and planned donor omissions.
#' [summarizeInfluence()] reports observed changes and valid-effect coverage.
#' [plotInfluence()] shows coverage, effects and changes with missingness.
#'
#' @section Inputs and interoperability:
#' Inputs use `SingleCellExperiment` or `SummarizedExperiment` and in-memory
#' numeric base or `Matrix` count matrices. Numeric sparse cell-level counts
#' stay sparse during aggregation. Logical/pattern and delayed or disk-backed
#' assays are unsupported. Already aggregated raw counts can be supplied via
#' `SummarizedExperiment` with a sample registry and cell counts. Results use
#' `SimpleList`, `DataFrame`, `SummarizedExperiment` and `CharacterList`.
#'
#' @section Output and interpretation:
#' The baseline gene-filtering family and contrast are fixed within each cell
#' type. Library normalization and model parameters are re-estimated for every
#' completed deletion. Failed or skipped omissions remain in `n_planned`.
#' Missing effects are not zero changes; observed effect ranges are not
#' confidence intervals. Empty summaries retain the same 16 typed columns.
#' The installed output dictionary lists all fields, types and status codes:
#' `system.file('OUTPUT_SCHEMA.md', package = 'scDonorAudit')`.
#'
#' @section Data and methods:
#' Help examples use synthetic negative-binomial counts (seed 1103); see
#' `extdata/README-example.md` and `scripts/create_synthetic_example.R`.
#' The vignette also uses attributed Crowell19 Astrocyte pseudobulk counts;
#' see `extdata/README-crowell19.md`. Method and data citations are in the
#' vignette. For edgeR use [edgeR::edgeR] and `citation('edgeR')`.
#'
#' @seealso
#' - [preparePseudobulk()]
#' - [auditDesign()]
#' - [assessDonorInfluence()]
#' - [summarizeInfluence()]
#' - [plotInfluence()]
#' @name scDonorAudit
#' @aliases scDonorAudit-package
"_PACKAGE"
