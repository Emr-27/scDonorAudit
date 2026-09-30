#!/usr/bin/env Rscript

# AI-assisted source (OpenAI Codex); see inst/CODE_PROVENANCE.md.

# Create an attribution-ready Astrocyte pseudobulk from EH3297.
# Usage: Rscript create_crowell19_pseudobulk.R <EH3297-Rda> <output-directory>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("supply source Rda and output directory")
source_file <- args[[1L]]
output_dir <- args[[2L]]
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

suppressPackageStartupMessages(library(SingleCellExperiment))
loaded <- new.env(parent = emptyenv())
names <- load(source_file, envir = loaded)
if (!identical(names, "sce")) stop("expected one sce object")
sce <- loaded$sce
if (!identical(dim(sce), c(11076L, 25224L))) {
    stop("unexpected EH3297 dimensions")
}

keep <- !is.na(sce$cluster_id) & sce$cluster_id == "Astrocytes"
selected <- sce[, keep]
samples <- sort(unique(as.character(selected$sample_id)), method = "radix")
if (length(samples) != 8L) stop("expected eight Astrocyte samples")
groups <- vapply(samples, function(id) {
    observed <- unique(as.character(selected$group_id[
        as.character(selected$sample_id) == id]))
    if (length(observed) != 1L) stop("ambiguous group label")
    observed
}, character(1))
if (!identical(sort(as.integer(table(groups))), c(4L, 4L))) {
    stop("expected four samples per condition")
}

raw <- SummarizedExperiment::assay(selected, "counts")
aggregated <- vapply(samples, function(id) {
    Matrix::rowSums(raw[, as.character(selected$sample_id) == id,
                        drop = FALSE])
}, numeric(nrow(selected)))
gene_ids <- as.character(SummarizedExperiment::rowData(selected)$ENSEMBL)
if (anyNA(gene_ids) || anyDuplicated(gene_ids)) {
    stop("gene IDs must be unique and complete")
}
row_order <- order(gene_ids, method = "radix")
small <- aggregated[row_order, , drop = FALSE]
storage.mode(small) <- "integer"
rownames(small) <- gene_ids[row_order]
colnames(small) <- samples
registry <- data.frame(sample_id = samples, donor_id = samples,
    condition = groups, cell_type = "Astrocytes",
    n_cells = vapply(samples, function(id) {
        sum(as.character(selected$sample_id) == id)
    }, integer(1)))
if (any(registry$n_cells < 20L)) stop("one sample has fewer than 20 cells")

utils::write.csv(small, file.path(output_dir, "crowell19_astrocyte_counts.csv"))
utils::write.csv(registry,
    file.path(output_dir, "crowell19_astrocyte_samples.csv"),
    row.names = FALSE)
cat("Created", nrow(small), "gene x 8-sample example from",
    ncol(selected), "Astrocyte nuclei\n")
