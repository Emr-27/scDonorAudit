# scDonorAudit

This is an early prototype intended for Bioconductor submission. The analysis reports
descriptive sensitivity to donor omission; it does not identify invalid donors
or provide a new false discovery rate guarantee.

The implemented entry points are `preparePseudobulk()`, `auditDesign()`,
`assessDonorInfluence()`, `summarizeInfluence()` and `plotInfluence()`.

## Minimal local example

Use R with `edgeR`, `SingleCellExperiment`, `SummarizedExperiment` and the
other dependencies in `DESCRIPTION` installed. From the package directory,
install the local prototype with `install.packages(".", repos = NULL,
type = "source")`. The example uses already aggregated counts; each column
is one observed sample and cell type.

```r
set.seed(1103)
sample_id <- paste0("s", 1:8)
donor_id <- paste0("d", 1:8)
condition <- rep(c("ctrl", "stim"), each = 4)
counts <- matrix(rnbinom(250 * 8, mu = 25, size = 5), nrow = 250,
                 dimnames = list(paste0("g", 1:250), sample_id))
se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(sample_id, donor_id, condition,
                                   cell_type = "T", n_cells = 30L))
registry <- data.frame(sample_id, donor_id, condition)

pb <- scDonorAudit::preparePseudobulk(
    se, "sample_id", "donor_id", "condition", "cell_type",
    sample_table = registry, n_cells = "n_cells")
contrast <- c(numerator = "stim", denominator = "ctrl")
audit <- scDonorAudit::auditDesign(pb, "independent", contrast)
result <- scDonorAudit::assessDonorInfluence(pb, "independent", contrast)
fit_ledger <- result[["fits"]]
gene_summary <- scDonorAudit::summarizeInfluence(result,
                                                 effect_threshold = 0.5)
scDonorAudit::plotInfluence(result, "effect", cell_type = "T", gene_id = "g1")
```

Inspect `fit_ledger` and `result[["sample_membership"]]` before interpreting
gene effects. `n_planned` counts all eligible donor deletions, including failed
ones. `max_abs_delta_observed` describes only deletion fits with valid effects;
missing fits remain missing. Adjusted p values use one fixed gene family per
cell type and do not control error across cell types or deletion runs.

The work was developed with AI assistance. Substantial generated code and
its provenance will be disclosed in a future Bioconductor submission as
required by its contribution policy. The maintainer must review and validate
all code before such a submission.
