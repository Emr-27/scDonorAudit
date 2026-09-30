# AI-assisted source (OpenAI Codex); see inst/CODE_PROVENANCE.md.

.scd_validate_registry <- function(registry, options) {
    observed_conditions <- sort(unique(registry$condition), method = "radix")
    expected <- sort(unname(options$contrast), method = "radix")
    if (!identical(observed_conditions, expected)) {
        .scd_stop("UNSUPPORTED_CONDITIONS",
            "sample registry must contain exactly the two contrast conditions")
    }
    count_by_donor <- split(registry$condition, registry$donor_id)
    if (options$design == "paired" &&
        any(vapply(count_by_donor, length, integer(1)) != 2L)) {
        .scd_stop("INCOMPLETE_REGISTRY_PAIR",
            "paired design needs both conditions for every registered donor")
    }
    if (options$design == "independent" &&
        any(vapply(count_by_donor, length, integer(1)) != 1L)) {
        .scd_stop("INVALID_INDEPENDENT_DESIGN",
            "an independent donor may belong to only one condition")
    }
    if (!all(options$covariates %in% names(registry))) {
        .scd_stop("MISSING_COVARIATE", paste(
            setdiff(options$covariates, names(registry)), collapse = ", "))
    }
    for (field in options$covariates) {
        value <- registry[[field]]
        if (!(is.numeric(value) || is.factor(value) || is.character(value)) ||
            anyNA(value) || (is.numeric(value) && any(!is.finite(value)))) {
            .scd_stop("INVALID_COVARIATE", paste0(field,
                " must be a finite numeric or complete categorical sample variable"))
        }
    }
    invisible(TRUE)
}

.scd_covariate_spec <- function(data, covariates) {
    output <- list()
    for (field in covariates) {
        value <- data[[field]]
        if (is.numeric(value)) {
            center <- mean(value)
            scale <- stats::sd(value)
            if (!is.finite(scale) || scale == 0) scale <- 1
            output[[field]] <- list(type = "numeric", center = center,
                                    scale = scale)
        } else {
            output[[field]] <- list(type = "factor",
                                    levels = sort(unique(as.character(value)),
                                                  method = "radix"))
        }
    }
    output
}

.scd_make_design <- function(data, options, cov_spec) {
    if (nrow(data) == 0L) {
        return(list(status = "NO_ELIGIBLE_SAMPLES", matrix = NULL,
                    contrast = NULL, rank = 0L, residual_df = 0L,
                    target_estimable = FALSE))
    }
    if (!all(unname(options$contrast) %in% data$condition)) {
        return(list(status = "MISSING_CONDITION", matrix = NULL,
                    contrast = NULL, rank = NA_integer_, residual_df = NA_integer_,
                    target_estimable = FALSE))
    }
    input <- data.frame(
        donor_id = factor(data$donor_id,
            levels = sort(unique(data$donor_id), method = "radix")),
        condition = factor(data$condition,
            levels = c(options$contrast[["denominator"]],
                       options$contrast[["numerator"]])))
    for (field in options$covariates) {
        value <- data[[field]]
        spec <- cov_spec[[field]]
        input[[field]] <- if (spec$type == "numeric") {
            (value - spec$center) / spec$scale
        } else factor(as.character(value), levels = spec$levels)
    }
    predictors <- c(if (options$design == "paired") "donor_id",
                    "condition", options$covariates)
    formula <- stats::reformulate(predictors)
    matrix <- tryCatch(stats::model.matrix(formula, data = input),
                       error = function(e) NULL)
    if (is.null(matrix)) {
        return(list(status = "DESIGN_MATRIX_ERROR", matrix = NULL,
                    contrast = NULL, rank = NA_integer_, residual_df = NA_integer_,
                    target_estimable = NA))
    }
    terms <- attr(stats::terms(formula), "term.labels")
    condition_term <- match("condition", terms)
    condition_columns <- which(attr(matrix, "assign") == condition_term)
    if (length(condition_columns) != 1L) {
        return(list(status = "CONTRAST_NOT_ESTIMABLE", matrix = matrix,
                    contrast = NULL, rank = qr(matrix)$rank,
                    residual_df = nrow(matrix) - qr(matrix)$rank,
                    target_estimable = FALSE))
    }
    contrast_vector <- rep(0, ncol(matrix))
    contrast_vector[condition_columns] <- 1
    names(contrast_vector) <- colnames(matrix)
    decomposition <- svd(matrix)
    threshold <- max(dim(matrix)) * max(decomposition$d) * .Machine$double.eps
    rank <- sum(decomposition$d > threshold)
    if (rank > 0L) {
        basis <- decomposition$v[, seq_len(rank), drop = FALSE]
        projected <- as.vector(basis %*% crossprod(basis, contrast_vector))
    } else projected <- rep(0, length(contrast_vector))
    estimable <- sqrt(sum((contrast_vector - projected)^2)) <=
        1e-7 * max(1, sqrt(sum(contrast_vector^2)))
    residual_df <- nrow(matrix) - rank
    status <- if (rank < ncol(matrix)) "DESIGN_RANK_DEFICIENT" else
        if (residual_df <= 0L) "NO_RESIDUAL_DF" else
            if (!estimable) "CONTRAST_NOT_ESTIMABLE" else "OK"
    list(status = status, matrix = matrix, contrast = contrast_vector,
         rank = rank, residual_df = residual_df,
         target_estimable = estimable)
}

.scd_audit <- function(pb, options) {
    .scd_check_metadata(pb)
    meta <- S4Vectors::metadata(pb)$scdonoraudit
    registry <- as.data.frame(meta$sample_table)
    .scd_validate_registry(registry, options)
    coverage <- as.data.frame(meta$coverage)
    if (!is.null(options$min_cells) &&
        any(coverage$observed & is.na(coverage$n_cells))) {
        .scd_stop("CELL_COUNT_UNKNOWN",
            "min_cells requires n_cells for every observed pseudobulk")
    }
    pb_columns <- as.data.frame(SummarizedExperiment::colData(pb))
    cell_types <- sort(unique(coverage$cell_type), method = "radix")
    all_specs <- list()
    all_samples <- list()
    all_runs <- list()
    internal <- list()
    counts <- SummarizedExperiment::assay(pb, "counts")
    for (cell in cell_types) {
        local <- coverage[coverage$cell_type == cell, , drop = FALSE]
        col_index <- which(pb_columns$cell_type == cell)
        local$pb_index <- col_index[match(local$sample_id,
                                         pb_columns$sample_id[col_index])]
        below_minimum <- rep(FALSE, nrow(local))
        if (!is.null(options$min_cells)) {
            below_minimum <- local$observed & !is.na(local$n_cells) &
                local$n_cells < options$min_cells
        }
        local$reason <- ifelse(!local$observed, "MISSING_CELL_TYPE",
            ifelse(below_minimum, "BELOW_MIN_CELLS", "OK"))
        local$in_baseline <- local$reason == "OK"
        pair_blocked <- FALSE
        if (options$design == "paired") {
            group_count <- tapply(local$in_baseline, local$donor_id, sum)
            incomplete <- names(group_count)[group_count == 1L]
            if (length(incomplete) && options$pair_policy == "strict") {
                pair_blocked <- TRUE
                local$in_baseline[] <- FALSE
                local$reason[local$reason == "OK"] <- "STRICT_PAIR_BLOCK"
            } else if (length(incomplete)) {
                displaced <- local$donor_id %in% incomplete & local$in_baseline
                local$in_baseline[displaced] <- FALSE
                local$reason[displaced] <- "PAIR_MEMBER_EXCLUDED"
            }
        }
        baseline_samples <- local$sample_id[local$in_baseline]
        baseline_indices <- local$pb_index[local$in_baseline]
        data <- registry[match(baseline_samples, registry$sample_id),,
                         drop = FALSE]
        cov_spec <- .scd_covariate_spec(data, options$covariates)
        donors <- sort(unique(data$donor_id), method = "radix")
        created <- !pair_blocked
        run_ids <- c("baseline", sprintf("omit_%03d", seq_along(donors)))
        omissions <- c(NA_character_, donors)
        if (!created) {
            run_ids <- "baseline"
            omissions <- NA_character_
        }
        rows <- vector("list", length(run_ids))
        for (j in seq_along(run_ids)) {
            if (pair_blocked) {
                checked <- list(status = "INCOMPLETE_PAIR", rank = NA_integer_,
                    residual_df = NA_integer_, target_estimable = NA)
                run_indices <- integer(0)
            } else {
                keep <- if (is.na(omissions[j])) rep(TRUE, nrow(data)) else
                    data$donor_id != omissions[j]
                run_data <- data[keep, , drop = FALSE]
                checked <- .scd_make_design(run_data,
                                            options, cov_spec)
                run_indices <- baseline_indices[keep]
                if (checked$status == "OK" &&
                    any(Matrix::colSums(counts[, run_indices,
                                              drop = FALSE]) == 0)) {
                    checked$status <- "ZERO_LIBRARY"
                }
            }
            replicate_count <- if (pair_blocked) NA_integer_ else
                if (options$design == "paired")
                    length(unique(run_data$donor_id)) else
                    min(tabulate(match(run_data$condition,
                                       unname(options$contrast)), nbins = 2L))
            rows[[j]] <- data.frame(cell_type = cell, run_id = run_ids[j],
                run_type = if (j == 1L) "baseline" else "omission",
                omitted_donor = omissions[j],
                execution_status = if (checked$status == "OK") "ready" else
                    "skipped", reason_code = checked$status,
                stage = "design", n_samples = length(run_indices),
                n_donors = if (pair_blocked) 0L else
                    length(unique(data$donor_id[if (is.na(omissions[j]))
                        rep(TRUE, nrow(data)) else data$donor_id != omissions[j]])),
                design_rank = checked$rank,
                residual_df = checked$residual_df,
                target_estimable = checked$target_estimable,
                low_replication = if (is.na(replicate_count)) NA else
                    replicate_count < options$min_replicates_warn,
                stringsAsFactors = FALSE)
        }
        runs <- do.call(rbind, rows)
        all_runs[[cell]] <- runs
        all_samples[[cell]] <- data.frame(
            cell_type = cell, sample_id = local$sample_id,
            donor_id = local$donor_id, condition = local$condition,
            observed = local$observed, n_cells = local$n_cells,
            in_baseline = local$in_baseline, reason_code = local$reason,
            stringsAsFactors = FALSE)
        all_specs[[cell]] <- data.frame(cell_type = cell,
            plan_status = if (created) "created" else "not_created",
            n_planned = if (created) length(donors) else NA_integer_,
            n_samples = nrow(data), n_donors = length(donors),
            baseline_status = runs$reason_code[1L],
            stringsAsFactors = FALSE)
        internal[[cell]] <- list(baseline_indices = baseline_indices,
            baseline_data = data, cov_spec = cov_spec,
            donors = donors, run_ids = run_ids)
    }
    spec <- do.call(rbind, all_specs)
    samples <- do.call(rbind, all_samples)
    runs <- do.call(rbind, all_runs)
    rownames(spec) <- rownames(samples) <- rownames(runs) <- NULL
    list(spec = spec, samples = samples, runs = runs,
         issues = runs[runs$reason_code != "OK", , drop = FALSE],
         config = options, internal = internal,
         coverage = coverage,
         preparation = meta$preparation)
}

#' Audit a two-condition donor design
#'
#' Reports which samples enter each cell type's analysis and whether a model
#' can be fit before and after omitting each included donor. Low replication is
#' recorded as a reminder; structural model checks determine whether fitting
#' can proceed.
#'
#' @param pb Output of [preparePseudobulk()].
#' @param design Either `"independent"` or `"paired"`.
#' @param contrast Named character vector such as
#'   `c(numerator="stim", denominator="ctrl")`.
#' @param min_cells Optional minimum cell count per sample/cell type.
#' @param covariates Names of sample-level covariates retained during preparation.
#' @param pair_policy `"strict"` blocks a cell type with an incomplete eligible
#'   pair; `"complete_pairs"` drops both samples of that donor from that cell
#'   type's baseline analysis.
#' @param min_replicates_warn Threshold for a descriptive low-replication flag.
#' @return A `SimpleList` of `DataFrame` tables and resolved settings.
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
#' audit <- auditDesign(pb, "independent",
#'                      c(numerator = "stim", denominator = "ctrl"))
#' head(as.data.frame(audit[["runs"]]))
#' @export
auditDesign <- function(pb, design, contrast, min_cells = NULL,
                        covariates = character(), pair_policy = "strict",
                        min_replicates_warn = 3L) {
    options <- .scd_options(design, contrast, min_cells, covariates,
                            pair_policy, min_replicates_warn)
    audit <- .scd_audit(pb, options)
    S4Vectors::SimpleList(
        spec = S4Vectors::DataFrame(audit$spec),
        sample_membership = S4Vectors::DataFrame(audit$samples),
        runs = S4Vectors::DataFrame(audit$runs),
        issues = S4Vectors::DataFrame(audit$issues),
        config = audit$config,
        provenance = list(audit_universe =
            audit$preparation$audit_universe))
}
