.scd_plot_theme <- function() {
    ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = "white", color = NA),
            panel.background = ggplot2::element_rect(fill = "white", color = NA),
            text = ggplot2::element_text(color = "#222222"))
}

#' Plot donor coverage or gene-level influence
#'
#' @param result Output of [assessDonorInfluence()].
#' @param type One of `"coverage"`, `"influence"`, or `"effect"`.
#' @param cell_type Required for gene-level plots.
#' @param gene_id Explicit gene IDs; a nonempty vector for `"influence"`
#'   and one ID for `"effect"`.
#' @return A `ggplot` object.
#' @details Grey influence tiles mark unavailable effects, which may arise
#'   from design, backend, or gene-level status. Inspect `fits` and the gene
#'   status assay for the specific reason.
#'   Effect plots use unique run IDs for positions and label deletions as
#'   `omit: <donor ID>`. Unavailable effects retain their position without a point.
#' @importFrom rlang .data
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
#' plotInfluence(result, "effect", cell_type = "T", gene_id = "g1")
#' @export
plotInfluence <- function(result, type, cell_type = NULL, gene_id = NULL) {
    if (!methods::is(result, "SimpleList") ||
        !all(c("results", "sample_membership") %in% names(result))) {
        .scd_stop("INVALID_RESULT", "use assessDonorInfluence() first")
    }
    if (!is.character(type) || length(type) != 1L ||
        !type %in% c("coverage", "influence", "effect")) {
        .scd_stop("INVALID_ARGUMENT", "unsupported plot type")
    }
    if (type == "coverage") {
        table <- as.data.frame(result[["sample_membership"]])
        table$status <- ifelse(table$in_baseline, "included",
            tolower(table$reason_code))
        table$sample <- paste(table$donor_id, table$sample_id, sep = " / ")
        return(ggplot2::ggplot(table,
            ggplot2::aes(x = .data$cell_type,
                y = .data$sample, fill = .data$status)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::labs(x = "Cell type", y = "Donor / sample",
                          fill = "Status") +
            .scd_plot_theme())
    }
    .scd_scalar_string(cell_type, "cell_type")
    if (!cell_type %in% names(result[["results"]])) {
        .scd_stop("UNKNOWN_CELL_TYPE", cell_type)
    }
    se <- result[["results"]][[cell_type]]
    if (!is.character(gene_id) || !length(gene_id) || anyNA(gene_id) ||
        any(!nzchar(gene_id)) || anyDuplicated(gene_id) ||
        !all(gene_id %in% rownames(se))) {
        .scd_stop("UNKNOWN_GENE_ID", "supply retained gene IDs explicitly")
    }
    if (type == "effect" && length(gene_id) != 1L) {
        .scd_stop("INVALID_ARGUMENT", "effect plot needs one gene ID")
    }
    effects <- SummarizedExperiment::assay(se, "logFC")
    valid <- SummarizedExperiment::assay(se, "effect_valid")
    donors <- as.character(SummarizedExperiment::colData(se)$omitted_donor)
    if (type == "influence") {
        if (ncol(se) < 2L) .scd_stop("NO_DELETIONS", "no donor deletion was planned")
        values <- effects[gene_id, -1L, drop = FALSE] -
            effects[gene_id, 1L]
        usable <- valid[gene_id, -1L, drop = FALSE] &
            valid[gene_id, 1L]
        values[!usable] <- NA_real_
        table <- expand.grid(gene_id = gene_id, donor_id = donors[-1L],
                             stringsAsFactors = FALSE)
        table$delta_logFC <- as.vector(values)
        return(ggplot2::ggplot(table,
            ggplot2::aes(x = .data$donor_id,
                y = .data$gene_id,
                fill = .data$delta_logFC)) +
            ggplot2::geom_tile(color = "grey85") +
            ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white",
                high = "#b2182b", midpoint = 0, na.value = "grey70") +
            ggplot2::labs(x = "Omitted donor", y = "Gene ID",
                fill = "Change in log2FC",
                title = paste("Donor influence:", cell_type),
                subtitle = paste("Grey: effect unavailable; see fit and gene",
                                 "status for the reason")) +
            .scd_plot_theme() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(
                angle = 60, hjust = 1, vjust = 1, size = 8)))
    }
    runs <- as.character(SummarizedExperiment::colData(se)$run_id)
    labels <- c("baseline", paste0("omit: ", donors[-1L]))
    table <- data.frame(run = factor(runs, levels = runs),
        logFC = as.numeric(effects[gene_id, ]),
        valid = as.logical(valid[gene_id, ]))
    table$logFC[!table$valid] <- NA_real_
    ggplot2::ggplot(table,
        ggplot2::aes(x = .data$run, y = .data$logFC)) +
        ggplot2::geom_hline(yintercept = table$logFC[1L],
            linetype = "dashed", color = "grey50", na.rm = TRUE) +
        ggplot2::geom_point(size = 2, na.rm = TRUE) +
        ggplot2::scale_x_discrete(limits = runs, labels = labels, drop = FALSE) +
        ggplot2::labs(x = "Omitted donor", y = "log2 fold-change",
            title = paste(gene_id, "in", cell_type),
            subtitle = "Baseline is the first position; missing refits have no point") +
        .scd_plot_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 60, hjust = 1, vjust = 1, size = 8))
}
