# scDonorAudit

scDonorAudit reports how a two-condition single-cell pseudobulk analysis
changes when each biological donor is omitted in turn. The results describe
sensitivity to donor omission; they do not identify invalid donors or provide
a new false discovery rate guarantee.

The implemented entry points are `preparePseudobulk()`, `auditDesign()`,
`assessDonorInfluence()`, `summarizeInfluence()` and `plotInfluence()`.

Packages such as [scuttle](https://bioconductor.org/packages/scuttle)
aggregate cells, while [muscat](https://bioconductor.org/packages/muscat)
and [dreamlet](https://bioconductor.org/packages/dreamlet) provide broader
pseudobulk differential-expression workflows. scDonorAudit uses edgeR for a
specified two-condition analysis and adds a donor-deletion sensitivity ledger:
which fits were planned, which ran, why others did not, and how each retained
gene's effect changed. It is an audit layer, not a replacement DE test or a
method for automatically excluding donors.

The vignette also runs an eight-mouse Astrocyte pseudobulk example derived
from the public Crowell19 single-nucleus dataset. The bundled counts retain
all 11,076 genes and are 0.47 MB; source, transformation, and CC BY 4.0
attribution are documented in the
[example data notes](inst/extdata/README-crowell19.md).

## Installation

The package is currently available from its source repository. From the
package directory, with the dependencies in `DESCRIPTION` installed, run:

```r
install.packages(".", repos = NULL, type = "source")
```

After the package is accepted and released by Bioconductor, the standard
installation command will be:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scDonorAudit")
```

## Minimal local example

The example uses already aggregated counts; each column is one observed
sample and cell type.

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

AI assistance was used to develop code and documentation. The
[development provenance](inst/CODE_PROVENANCE.md) records its scope and the
available validation evidence.
