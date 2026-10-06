# scDonorAudit

scDonorAudit reports how a specified two-condition edgeR pseudobulk analysis
changes when each biological donor is omitted in turn. It links planned
omissions to design checks, fit outcomes, and gene-effect changes.

The implemented entry points are `preparePseudobulk()`, `auditDesign()`,
`assessDonorInfluence()`, `summarizeInfluence()` and `plotInfluence()`.

See `?scDonorAudit`, the [runnable vignette](vignettes/scDonorAudit.Rmd) and the installed
[output dictionary](inst/OUTPUT_SCHEMA.md) for fields, types and statuses.
After installation, open the evaluated tutorial with
`vignette("scDonorAudit", package = "scDonorAudit")`.

## Installation

The package is currently available from its source repository. Download or
clone this repository, then start R in the directory containing `DESCRIPTION`.
With the dependencies in that file installed, run:

```r
install.packages(".", repos = NULL, type = "source")
```

After acceptance and a successful Bioconductor build, the package can be
installed from the Bioconductor version that contains it: initially devel,
then the applicable release. Set up the matching R/Bioconductor environment
using the [Bioconductor installation guide](https://bioconductor.org/install/),
then run:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scDonorAudit")
```

## Minimal local example

The example uses already aggregated counts; each column is one observed
sample and cell type. Inputs must contain raw counts in in-memory numeric
matrices; preaggregated counts in a `SummarizedExperiment` are supported.
Delayed/disk-backed assays and logical/pattern matrices are unsupported.

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

## Interpret the results

Inspect `fit_ledger` and `result[["sample_membership"]]` before interpreting
gene effects. `n_planned` counts all eligible donor deletions, including failed
or skipped runs. Read `max_abs_delta_observed` alongside the number of valid
effects: it describes only deletion fits with comparable effects. Missing
effects remain missing, rather than becoming zero changes or non-significant
tests. The minimum and maximum deletion effects form an observed range, not
a confidence interval.

Genes retained by the baseline filter form one fixed comparison family per
cell type. Each successful deletion fit re-estimates normalization and the
edgeR model on its remaining samples. Adjusted p values apply within that
family and do not control error across cell types or deletion runs. A large
observed change can reflect donor heterogeneity or limited precision; it
does not establish that a donor is invalid or should be excluded.

The fit ledger keeps each run's original design diagnosis in
`design_reason_code` even when a later baseline or backend problem prevents
execution. See the vignette for the output dictionary and interpretation rules.

The vignette demonstrates missing coverage and skipped design failures alongside
valid effects. It also runs an eight-mouse Astrocyte pseudobulk example derived
from the public Crowell19 single-nucleus dataset. The bundled counts retain
all 11,076 genes and are 0.47 MB; source, transformation, and CC BY 4.0
attribution are documented in the
[example data notes](inst/extdata/README-crowell19.md).

## Related tools

Packages such as [scuttle](https://bioconductor.org/packages/scuttle)
aggregate cells, while [muscat](https://bioconductor.org/packages/muscat)
and [dreamlet](https://bioconductor.org/packages/dreamlet) provide broader
pseudobulk differential-expression workflows. scDonorAudit uses edgeR for a
specified two-condition independent or paired design and reports the planned
donor omissions, their execution outcomes, and the observed effect changes.

| Task | Existing tool | scDonorAudit output |
|:--|:--|:--|
| Cell aggregation | scuttle | Registered coverage, including absent combinations |
| Cell-type DE | muscat | Fixed gene family across whole-donor omission fits |
| Complex designs | dreamlet | Specified two-condition edgeR deletion plan |
| Model fitting | edgeR | Linked design/execution ledger and per-run sample records |
| Error reporting | dreamlet assay/gene errors | Original design reasons plus warnings and backend errors |
| Interpretation | Custom DE diagnostics | Observed change alongside planned/valid coverage |
