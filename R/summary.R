#' Summarize observed donor influence on each gene
#'
#' A missing or failed deletion remains in the planned denominator. Returned
#' minimum and maximum effects are observed ranges, not confidence intervals
#' or significance tests. Interpret a change only alongside `n_planned` and
#' `n_effect_valid`; near-zero sign changes need the material threshold.
#'
#' @param result Output of [assessDonorInfluence()].
#' @param effect_threshold Optional positive absolute log2 fold-change threshold
#'   for a material sign reversal.
#' @return A `DataFrame` with one row per retained gene and cell type.
#' @examples
#' count_file <- system.file("extdata", "example_counts.csv",
#'                           package = "scDonorAudit")
#' sample_file <- system.file("extdata", "example_samples.csv",
#'                            package = "scDonorAudit")
#' counts <- as.matrix(utils::read.csv(count_file, row.names = 1,
#'                                     check.names = FALSE))
#' samples <- utils::read.csv(sample_file)
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = counts),
#'     colData = S4Vectors::DataFrame(samples, cell_type = "T"))
#' pb <- preparePseudobulk(se, "sample_id", "donor_id", "condition",
#'                         "cell_type", sample_table = samples)
#' result <- assessDonorInfluence(pb, "independent",
#'                                c(numerator = "stim", denominator = "ctrl"))
#' influence <- summarizeInfluence(result, effect_threshold = 0.5)
#' head(as.data.frame(influence)[, c("gene_id", "n_effect_valid",
#'                                   "max_abs_delta_observed")])
#' @export
summarizeInfluence <- function(result, effect_threshold = NULL) {
    if (!methods::is(result, "SimpleList") ||
        !all(c("results", "spec") %in% names(result))) {
        .scd_stop("INVALID_RESULT", "use assessDonorInfluence() first")
    }
    if (!is.null(effect_threshold) &&
        (!is.numeric(effect_threshold) || length(effect_threshold) != 1L ||
         !is.finite(effect_threshold) || effect_threshold <= 0)) {
        .scd_stop("INVALID_ARGUMENT", "effect_threshold must be positive")
    }
    output <- list()
    spec <- as.data.frame(result[["spec"]])
    for (cell in names(result[["results"]])) {
        se <- result[["results"]][[cell]]
        if (nrow(se) == 0L) next
        effects <- SummarizedExperiment::assay(se, "logFC")
        p <- SummarizedExperiment::assay(se, "p_value")
        valid_effect <- SummarizedExperiment::assay(se, "effect_valid")
        valid_test <- SummarizedExperiment::assay(se, "test_valid")
        omitted <- as.character(SummarizedExperiment::colData(se)$omitted_donor)[-1L]
        planned <- spec$n_planned[match(cell, spec$cell_type)]
        rows <- vector("list", nrow(se))
        ties <- vector("list", nrow(se))
        for (i in seq_len(nrow(se))) {
            baseline <- effects[i, 1L]
            baseline_valid <- valid_effect[i, 1L]
            baseline_test <- valid_test[i, 1L]
            deletion <- if (ncol(se) > 1L) effects[i, -1L] else numeric()
            evalid <- if (ncol(se) > 1L) valid_effect[i, -1L] else logical()
            tvalid <- if (ncol(se) > 1L) valid_test[i, -1L] else logical()
            evalid <- evalid & baseline_valid
            tvalid <- tvalid & baseline_test
            delta <- deletion[evalid] - baseline
            n_effect <- sum(evalid)
            n_test <- sum(tvalid)
            max_delta <- if (n_effect) max(abs(delta)) else NA_real_
            winners <- character()
            if (n_effect) {
                tolerance <- 1e-8 * max(1, max_delta)
                winners <- omitted[evalid][abs(abs(delta) - max_delta) <= tolerance]
            }
            ties[[i]] <- winners
            reversal <- if (n_effect) sum(baseline * deletion[evalid] < 0) else NA_integer_
            material <- if (is.null(effect_threshold)) NA_integer_ else
                if (n_effect) sum((baseline > effect_threshold &
                    deletion[evalid] < -effect_threshold) |
                    (baseline < -effect_threshold &
                    deletion[evalid] > effect_threshold)) else NA_integer_
            rows[[i]] <- data.frame(cell_type = cell, gene_id = rownames(se)[i],
                baseline_logFC = if (baseline_valid) baseline else NA_real_,
                baseline_p_value = if (baseline_test) p[i, 1L] else NA_real_,
                n_planned = planned, n_effect_valid = n_effect,
                n_test_valid = n_test,
                max_abs_delta_observed = max_delta,
                min_logFC_observed = if (n_effect) min(deletion[evalid]) else NA_real_,
                max_logFC_observed = if (n_effect) max(deletion[evalid]) else NA_real_,
                n_sign_reversal = reversal,
                n_material_reversal = material,
                baseline_near_zero = if (is.null(effect_threshold) ||
                    !baseline_valid) NA else abs(baseline) <= effect_threshold,
                complete_effect_coverage = isTRUE(baseline_valid) &&
                    !is.na(planned) && planned > 0L && n_effect == planned,
                complete_test_coverage = isTRUE(baseline_test) &&
                    !is.na(planned) && planned > 0L && n_test == planned,
                stringsAsFactors = FALSE)
        }
        table <- S4Vectors::DataFrame(do.call(rbind, rows))
        table$max_influence_donors <- IRanges::CharacterList(ties)
        output[[cell]] <- table
    }
    if (!length(output)) return(S4Vectors::DataFrame())
    do.call(rbind, output)
}
