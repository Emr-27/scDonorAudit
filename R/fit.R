# AI-assisted source (OpenAI Codex); see inst/CODE_PROVENANCE.md.

.scd_fit_settings <- function(filter_args, backend_args) {
    defaults <- list(min.count = 10, min.total.count = 15,
                     large.n = 10, min.prop = 0.7)
    if (!is.list(filter_args) || is.null(names(filter_args)) &&
        length(filter_args) || anyDuplicated(names(filter_args)) ||
        !all(names(filter_args) %in% names(defaults))) {
        .scd_stop("INVALID_FILTER_ARGS", "unsupported filter_args")
    }
    for (name in names(filter_args)) {
        value <- filter_args[[name]]
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 0 || (name == "min.prop" && value > 1)) {
            .scd_stop("INVALID_FILTER_ARGS", paste0("invalid ", name))
        }
        defaults[[name]] <- value
    }
    backend <- list(robust = TRUE)
    if (!is.list(backend_args) || is.null(names(backend_args)) &&
        length(backend_args) || anyDuplicated(names(backend_args)) ||
        !all(names(backend_args) %in% names(backend))) {
        .scd_stop("INVALID_BACKEND_ARGS", "unsupported backend_args")
    }
    for (name in names(backend_args)) backend[[name]] <- backend_args[[name]]
    if (!is.logical(backend$robust) || length(backend$robust) != 1L ||
        is.na(backend$robust)) {
        .scd_stop("INVALID_BACKEND_ARGS", "invalid robust")
    }
    list(filter = defaults, backend = backend)
}

.scd_empty_result <- function(genes, runs, cell) {
    nr <- length(genes)
    nc <- nrow(runs)
    numeric_assay <- matrix(NA_real_, nr, nc,
        dimnames = list(genes, runs$run_id))
    logical_assay <- matrix(FALSE, nr, nc,
        dimnames = dimnames(numeric_assay))
    status_assay <- matrix(5L, nr, nc,
        dimnames = dimnames(numeric_assay))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(logFC = numeric_assay,
            p_value = numeric_assay,
            padj_within_cell_type = numeric_assay,
            effect_valid = logical_assay,
            test_valid = logical_assay,
            gene_status_code = status_assay),
        colData = S4Vectors::DataFrame(runs[, c("run_id", "run_type",
            "omitted_donor"), drop = FALSE]),
        metadata = list(cell_type = cell,
            status_codes = c(OK = 0L, ALL_ZERO_AFTER_DELETION = 1L,
                NONFINITE_EFFECT = 2L, INVALID_PVALUE = 3L,
                ONE_CONDITION_ZERO = 4L, RUN_UNAVAILABLE = 5L)))
    se
}

.scd_adjust_p <- function(p, effect_valid) {
    valid <- effect_valid & is.finite(p) & p >= 0 & p <= 1
    adjusted <- rep(NA_real_, length(p))
    if (any(valid)) {
        adjusted[valid] <- stats::p.adjust(p[valid], method = "BH",
                                           n = length(p))
    }
    list(test_valid = valid, adjusted = adjusted)
}

.scd_result_from_table <- function(tab, counts, data, active, family_size) {
    all_zero <- Matrix::rowSums(counts) == 0
    effect <- rep(NA_real_, family_size)
    p <- padj <- rep(NA_real_, family_size)
    effect[active] <- tab$logFC
    p[active] <- tab$PValue
    effect_valid <- is.finite(effect)
    adjustment <- .scd_adjust_p(p, effect_valid)
    test_valid <- adjustment$test_valid
    padj <- adjustment$adjusted
    condition_zero <- Matrix::rowSums(counts[, data$condition ==
        unique(data$condition)[1L], drop = FALSE]) == 0 |
        Matrix::rowSums(counts[, data$condition ==
        unique(data$condition)[2L], drop = FALSE]) == 0
    code <- rep(0L, family_size)
    code[condition_zero & !all_zero] <- 4L
    code[!effect_valid & !all_zero] <- 2L
    code[effect_valid & !test_valid] <- 3L
    code[all_zero] <- 1L
    list(effect = effect, p = p, padj = padj,
         effect_valid = effect_valid, test_valid = test_valid,
         code = code)
}

.scd_fit_run <- function(counts, data, design, settings, family_size) {
    active <- which(Matrix::rowSums(counts) != 0)
    if (!length(active)) .scd_stop("ALL_ZERO_RUN", "all genes are zero")
    use_counts <- as.matrix(counts[active, , drop = FALSE])
    y <- edgeR::DGEList(counts = use_counts)
    y <- edgeR::normLibSizes(y, method = "TMM")
    fitted <- edgeR::glmQLFit(y, design = design$matrix,
        dispersion = NULL, abundance.trend = TRUE,
        robust = settings$backend$robust, legacy = FALSE,
        top.proportion = NULL)
    tested <- edgeR::glmQLFTest(fitted, contrast = design$contrast)
    result <- .scd_result_from_table(tested$table, counts, data,
        active, family_size)
    c(result, list(library_size = y$samples$lib.size,
        norm_factors = y$samples$norm.factors,
        design_columns = colnames(design$matrix)))
}

#' Assess the influence of omitting each donor
#'
#' Fits an edgeR quasi-likelihood model to a fixed baseline gene family in
#' each cell type, then removes one donor at a time and refits. This is a
#' sensitivity analysis, not an automatic rule for excluding donors.
#'
#' @inheritParams auditDesign
#' @param filter_args Named list of `filterByExpr` options: `min.count`,
#'   `min.total.count`, `large.n`, and `min.prop`.
#' @param backend_args Named list with `robust`. Other backend options,
#'   including `prior.count`, use the defaults of the installed edgeR version.
#' @return A `SimpleList` containing per-cell-type results, fit ledger,
#'   sample membership, gene filter, coverage, issues, and configuration.
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
#' head(as.data.frame(result[["fits"]]))
#' @export
assessDonorInfluence <- function(pb, design, contrast, min_cells = NULL,
                                 covariates = character(),
                                 pair_policy = "strict",
                                 min_replicates_warn = 3L,
                                 filter_args = list(), backend_args = list()) {
    options <- .scd_options(design, contrast, min_cells, covariates,
                            pair_policy, min_replicates_warn)
    settings <- .scd_fit_settings(filter_args, backend_args)
    audited <- .scd_audit(pb, options)
    full_counts <- SummarizedExperiment::assay(pb, "counts")
    genes <- rownames(pb)
    results <- list()
    fit_tables <- list()
    filters <- list()
    problems <- list()
    for (cell in names(audited$internal)) {
        info <- audited$internal[[cell]]
        runs <- audited$runs[audited$runs$cell_type == cell, , drop = FALSE]
        runs$library_sizes <- NA_character_
        runs$norm_factors <- NA_character_
        runs$design_columns <- NA_character_
        kept <- rep(FALSE, length(genes))
        filter_status <- "not_run"
        if (runs$reason_code[1L] == "OK") {
            raw <- full_counts[, info$baseline_indices, drop = FALSE]
            filter_status <- "filtered"
            try_filter <- tryCatch(do.call(edgeR::filterByExpr,
                c(list(y = as.matrix(raw),
                       group = info$baseline_data$condition),
                  settings$filter)), error = function(e) e)
            if (inherits(try_filter, "error")) {
                runs$execution_status[1L] <- "failed"
                runs$reason_code[1L] <- "FILTER_ERROR"
                runs$stage[1L] <- "filter"
                filter_status <- "failed"
                problems[[length(problems) + 1L]] <- data.frame(
                    cell_type = cell, run_id = "baseline", stage = "filter",
                    reason_code = "FILTER_ERROR",
                    message = conditionMessage(try_filter))
            } else {
                kept <- as.logical(try_filter)
                if (!any(kept)) {
                    runs$execution_status[1L] <- "skipped"
                    runs$reason_code[1L] <- "NO_GENES_AFTER_FILTER"
                    runs$stage[1L] <- "filter"
                }
            }
        }
        filters[[cell]] <- data.frame(cell_type = cell, gene_id = genes,
            kept = kept, filter_status = filter_status,
            stringsAsFactors = FALSE)
        selected <- genes[kept]
        se <- .scd_empty_result(selected, runs, cell)
        if (runs$reason_code[1L] != "OK") {
            if (nrow(runs) > 1L) {
                runs$execution_status[-1L] <- "skipped"
                runs$reason_code[-1L] <- "BASELINE_UNAVAILABLE"
                runs$stage[-1L] <- "baseline"
            }
        } else {
            for (j in seq_len(nrow(runs))) {
                if (j > 1L && runs$reason_code[1L] != "OK") {
                    runs$execution_status[j] <- "skipped"
                    runs$reason_code[j] <- "BASELINE_UNAVAILABLE"
                    runs$stage[j] <- "baseline"
                    next
                }
                if (runs$reason_code[j] != "OK") next
                use <- if (j == 1L) rep(TRUE, nrow(info$baseline_data)) else
                    info$baseline_data$donor_id != runs$omitted_donor[j]
                data <- info$baseline_data[use, , drop = FALSE]
                local_counts <- full_counts[kept, info$baseline_indices[use],
                                            drop = FALSE]
                design_info <- .scd_make_design(data, options, info$cov_spec)
                warnings <- new.env(parent = emptyenv())
                warnings$messages <- character()
                fitted <- tryCatch(withCallingHandlers(
                    .scd_fit_run(local_counts, data, design_info,
                                 settings, length(selected)),
                    warning = function(w) {
                        warnings$messages <- c(warnings$messages,
                                               conditionMessage(w))
                        invokeRestart("muffleWarning")
                    }), error = function(e) e)
                if (inherits(fitted, "error")) {
                    runs$execution_status[j] <- "failed"
                    runs$reason_code[j] <- "BACKEND_ERROR"
                    runs$stage[j] <- "backend"
                    problems[[length(problems) + 1L]] <- data.frame(
                        cell_type = cell, run_id = runs$run_id[j],
                        stage = "backend", reason_code = "BACKEND_ERROR",
                        message = conditionMessage(fitted))
                    next
                }
                for (field in c("logFC", "p_value", "padj_within_cell_type",
                                "effect_valid", "test_valid", "gene_status_code")) {
                    value <- switch(field, logFC = fitted$effect,
                        p_value = fitted$p,
                        padj_within_cell_type = fitted$padj,
                        effect_valid = fitted$effect_valid,
                        test_valid = fitted$test_valid,
                        gene_status_code = fitted$code)
                    current <- SummarizedExperiment::assay(se, field)
                    current[, j] <- value
                    SummarizedExperiment::assay(se, field) <- current
                }
                runs$execution_status[j] <- "completed"
                runs$reason_code[j] <- "OK"
                runs$stage[j] <- "backend"
                runs$library_sizes[j] <- paste(fitted$library_size,
                    collapse = ",")
                runs$norm_factors[j] <- paste(fitted$norm_factors,
                    collapse = ",")
                runs$design_columns[j] <- paste(fitted$design_columns,
                    collapse = ",")
                if (length(warnings$messages)) {
                    problems[[length(problems) + 1L]] <- data.frame(
                        cell_type = cell, run_id = runs$run_id[j],
                        stage = "backend", reason_code = "BACKEND_WARNING",
                        message = paste(unique(warnings$messages),
                                        collapse = " | "))
                }
            }
        }
        results[[cell]] <- se
        fit_tables[[cell]] <- runs
    }
    fits <- do.call(rbind, fit_tables)
    rownames(fits) <- NULL
    gene_filter <- do.call(rbind, filters)
    rownames(gene_filter) <- NULL
    issues <- if (length(problems)) do.call(rbind, problems) else
        data.frame(cell_type = character(), run_id = character(),
                   stage = character(), reason_code = character(),
                   message = character())
    failed_design <- fits[fits$reason_code != "OK", c("cell_type", "run_id",
        "stage", "reason_code"), drop = FALSE]
    if (nrow(failed_design)) {
        failed_design$message <- ""
        issues <- rbind(issues, failed_design)
    }
    S4Vectors::SimpleList(
        results = S4Vectors::SimpleList(results),
        fits = S4Vectors::DataFrame(fits),
        spec = S4Vectors::DataFrame(audited$spec),
        sample_membership = S4Vectors::DataFrame(audited$samples),
        gene_filter = S4Vectors::DataFrame(gene_filter),
        coverage = S4Vectors::DataFrame(audited$coverage),
        issues = S4Vectors::DataFrame(issues),
        config = c(options, settings,
            list(adjustment = "BH_within_cell_type",
                 schema_version = "0.2")),
        provenance = list(input = audited$preparation,
            edgeR_version = as.character(utils::packageVersion("edgeR")),
            R_version = as.character(getRversion())))
}
