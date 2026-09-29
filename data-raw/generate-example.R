#!/usr/bin/env Rscript

# Rebuild the small, synthetic files used by the package help examples.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("supply an output directory")
dir.create(args[[1L]], recursive = TRUE, showWarnings = FALSE)

set.seed(1103)
samples <- data.frame(sample_id = paste0("s", seq_len(8)),
                      donor_id = paste0("d", seq_len(8)),
                      condition = rep(c("ctrl", "stim"), each = 4))
counts <- matrix(stats::rnbinom(250 * nrow(samples), mu = 25, size = 5),
                 nrow = 250,
                 dimnames = list(paste0("g", seq_len(250)),
                                 samples$sample_id))
counts[seq_len(10), samples$condition == "stim"] <-
    counts[seq_len(10), samples$condition == "stim"] * 2L
utils::write.csv(data.frame(gene_id = rownames(counts), counts,
                            check.names = FALSE),
                 file.path(args[[1L]], "example_counts.csv"),
                 row.names = FALSE)
utils::write.csv(samples,
                 file.path(args[[1L]], "example_samples.csv"),
                 row.names = FALSE)
