# Assisted-by: OpenAI Codex and DeepSeek V4.1 Flash.
# See inst/CODE_PROVENANCE.md for scope and maintenance responsibility.

.scd_stop <- function(code, detail) {
    stop(sprintf("[%s] %s", code, detail), call. = FALSE)
}

.scd_scalar_string <- function(x, what) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
        .scd_stop("INVALID_ARGUMENT", paste0(what,
            " must be one non-empty string"))
    }
    x
}

.scd_field <- function(df, name, what) {
    .scd_scalar_string(name, what)
    if (!name %in% colnames(df)) {
        .scd_stop("MISSING_FIELD", paste0(what, " column '", name,
            "' is missing"))
    }
    df[[name]]
}

.scd_ids <- function(x, what) {
    x <- as.character(x)
    if (anyNA(x) || any(!nzchar(x))) {
        .scd_stop("MISSING_ID", paste0(what,
            " contains a missing or empty value"))
    }
    x
}

.scd_counts <- function(x, assay_name) {
    if (!(is.matrix(x) || methods::is(x, "Matrix"))) {
        .scd_stop("UNSUPPORTED_COUNTS", paste0(assay_name,
            " must be a matrix or an in-memory Matrix object"))
    }
    if (methods::is(x, "lMatrix") || methods::is(x, "nMatrix")) {
        .scd_stop("INVALID_COUNTS", paste0(assay_name,
            " must contain finite non-negative integer-valued counts"))
    }
    values <- if (methods::is(x, "sparseMatrix"))
        x@x else as.vector(x)
    if (!is.numeric(values) || any(!is.finite(values)) || any(values <
        0) || any(values != floor(values))) {
        .scd_stop("INVALID_COUNTS", paste0(assay_name,
            " must contain finite non-negative integer-valued counts"))
    }
    invisible(x)
}

.scd_gene_ids <- function(x) {
    ids <- rownames(x)
    if (is.null(ids))
        .scd_stop("MISSING_GENE_ID", "input row names must contain gene IDs")
    ids <- .scd_ids(ids, "gene IDs")
    if (anyDuplicated(ids))
        .scd_stop("DUPLICATE_GENE_ID", "gene IDs must be unique")
    ids
}

.scd_one_value_per_sample <- function(values, samples, what) {
    if (length(values) != length(samples)) {
        .scd_stop("INVALID_METADATA", paste0(what, " has incorrect length"))
    }
    for (id in unique(samples)) {
        current <- values[samples == id]
        if (anyNA(current) || length(unique(as.character(current))) !=
            1L) {
            .scd_stop("METADATA_CONFLICT", paste0(what,
                " is missing or inconsistent within sample '",
                id, "'"))
        }
    }
    invisible(TRUE)
}

.scd_sample_table <- function(source_data, sample_id, donor_id, condition,
    sample_vars, sample_table) {
    observed <- data.frame(sample_id = .scd_ids(.scd_field(source_data,
        sample_id, "sample_id"), "sample_id"), stringsAsFactors = FALSE)
    for (field in sample_vars) {
        if (field %in% c("sample_id", "donor_id", "condition", "cell_type",
            "n_cells", "pb_id")) {
            .scd_stop("INVALID_ARGUMENT", paste0("sample variable '", field,
                "' uses a reserved name"))
        }
    }
    if (is.null(sample_table)) {
        observed$donor_id <- .scd_ids(.scd_field(source_data, donor_id,
            "donor_id"), "donor_id")
        observed$condition <- .scd_ids(.scd_field(source_data, condition,
            "condition"), "condition")
        for (field in sample_vars) observed[[field]] <- .scd_field(source_data,
            field, field)
        universe <- "observed_only"
    } else {
        if (!is.data.frame(sample_table)) {
            .scd_stop("INVALID_SAMPLE_TABLE",
                "sample_table must be a data.frame")
        }
        needed <- c("sample_id", "donor_id", "condition", sample_vars)
        if (!all(needed %in% names(sample_table))) {
            .scd_stop("INVALID_SAMPLE_TABLE", paste0("sample_table needs: ",
                paste(setdiff(needed, names(sample_table)), collapse = ", ")))
        }
        table <- as.data.frame(sample_table, stringsAsFactors = FALSE)
        for (field in c("sample_id", "donor_id", "condition")) {
            table[[field]] <- .scd_ids(table[[field]], paste0("sample_table$",
                field))
        }
        if (anyDuplicated(table$sample_id)) {
            .scd_stop("INVALID_SAMPLE_TABLE",
                "sample_table has duplicate sample_id")
        }
        if (anyNA(table[, sample_vars, drop = FALSE])) {
            .scd_stop("INVALID_SAMPLE_TABLE",
                "sample_table has missing sample variables")
        }
        matched <- match(observed$sample_id, table$sample_id)
        if (anyNA(matched)) {
            .scd_stop("SAMPLE_NOT_REGISTERED",
                "input contains a sample absent from sample_table")
        }
        for (field in c("donor_id", "condition")) {
            source_field <- if (field == "donor_id")
                donor_id else condition
            if (source_field %in% colnames(source_data) && any(as.character(
                source_data[[source_field]]) !=
                as.character(table[[field]][matched]) | is.na(
                    source_data[[source_field]]))) {
                .scd_stop("METADATA_CONFLICT", paste0(source_field,
                    " disagrees with sample_table"))
            }
        }
        for (field in sample_vars) {
            if (field %in% colnames(source_data) && any(as.character(
                source_data[[field]]) !=
                as.character(table[[field]][matched]) | is.na(
                    source_data[[field]]))) {
                .scd_stop("METADATA_CONFLICT", paste0(field,
                    " disagrees with sample_table"))
            }
        }
        observed$donor_id <- table$donor_id[matched]
        observed$condition <- table$condition[matched]
        for (field in sample_vars) observed[[field]] <- table[[field]][matched]
        universe <- "registered_samples"
    }
    for (field in setdiff(colnames(observed), "sample_id")) {
        .scd_one_value_per_sample(observed[[field]], observed$sample_id,
            field)
    }
    table <- if (is.null(sample_table)) {
        observed[!duplicated(observed$sample_id), , drop = FALSE]
    } else table[, c("sample_id", "donor_id", "condition", sample_vars),
        drop = FALSE]
    donor_levels <- sort(unique(table$donor_id), method = "radix")
    condition_levels <- sort(unique(table$condition), method = "radix")
    pair_key <- (match(table$donor_id, donor_levels) - 1L) * length(
        condition_levels) +
        match(table$condition, condition_levels)
    if (anyDuplicated(pair_key)) {
        .scd_stop("UNSUPPORTED_REPEATED_MEASURES",
            "each donor and condition must have at most one biological sample")
    }
    table <- table[order(table$sample_id, method = "radix"), , drop = FALSE]
    rownames(table) <- NULL
    list(table = table, observed = observed, universe = universe)
}

.scd_check_metadata <- function(pb) {
    if (!methods::is(pb, "SummarizedExperiment") || is.null(
        S4Vectors::metadata(pb)$scdonoraudit)) {
        .scd_stop("INVALID_PSEUDOBULK", "use preparePseudobulk() first")
    }
    if (!"counts" %in% SummarizedExperiment::assayNames(pb)) {
        .scd_stop("MISSING_ASSAY", "assay 'counts' is missing")
    }
    if (nrow(pb) == 0L || ncol(pb) == 0L) {
        .scd_stop("EMPTY_INPUT", "input needs at least one gene and one column")
    }
    .scd_gene_ids(pb)
    .scd_counts(SummarizedExperiment::assay(pb, "counts"), "counts")
    meta <- S4Vectors::metadata(pb)$scdonoraudit
    columns <- as.data.frame(SummarizedExperiment::colData(pb))
    registry <- meta$sample_table
    coverage <- meta$coverage
    core <- c("sample_id", "donor_id", "condition")
    if (!is.data.frame(registry) || !is.data.frame(coverage) || !all(c(core,
        "cell_type", "observed", "n_cells") %in% names(coverage)) || !all(c(
            core,
        "cell_type", "pb_id", "n_cells") %in% names(columns)) || !all(core %in%
        names(registry)) || !is.logical(coverage$observed) || anyNA(
            coverage$observed) ||
        anyNA(columns[, c(core, "cell_type", "pb_id"), drop = FALSE]) ||
        anyNA(coverage[, c(core, "cell_type"), drop = FALSE])) {
        .scd_stop("INVALID_PSEUDOBULK", "audit metadata is incomplete")
    }
    key <- function(x) paste0(nchar(x$sample_id), ":", x$sample_id, x$cell_type)
    column_key <- key(columns)
    coverage_key <- key(coverage)
    observed <- coverage[coverage$observed, , drop = FALSE]
    observed_key <- key(observed)
    sample_index <- match(coverage$sample_id, registry$sample_id)
    column_index <- match(column_key, observed_key)
    expected_rows <- nrow(registry) * length(unique(coverage$cell_type))
    if (anyDuplicated(registry$sample_id) || anyDuplicated(columns$pb_id) ||
        anyDuplicated(column_key) || anyDuplicated(coverage_key) || nrow(
            coverage) !=
        expected_rows || nrow(columns) != nrow(observed) || anyNA(
            sample_index) ||
        anyNA(column_index) || !identical(colnames(pb), as.character(
            columns$pb_id))) {
        .scd_stop("STALE_PSEUDOBULK",
            "registry or coverage differs from columns")
    }
    for (field in core[-1L]) {
        if (!identical(as.character(coverage[[field]]), as.character(
            registry[[field]][sample_index])) ||
            !identical(as.character(columns[[field]]), as.character(
                observed[[field]][column_index]))) {
            .scd_stop("STALE_PSEUDOBULK", paste0(field,
                " differs between registry, coverage, and columns"))
        }
    }
    if (!identical(as.character(columns$n_cells), as.character(
        observed$n_cells[column_index])) ||
        any(!coverage$observed & !is.na(coverage$n_cells))) {
        .scd_stop("STALE_PSEUDOBULK", "cell counts differ from coverage")
    }
    for (field in meta$preparation$sample_vars) {
        if (!field %in% names(registry) || !field %in% names(columns) ||
            !identical(as.character(columns[[field]]), as.character(
                registry[[field]][match(columns$sample_id,
                registry$sample_id)]))) {
            .scd_stop("STALE_PSEUDOBULK", paste0(field,
                " differs between registry and columns"))
        }
    }
    invisible(TRUE)
}

.scd_options <- function(design, contrast, min_cells, covariates, pair_policy,
    min_replicates_warn) {
    if (!identical(design, "independent") && !identical(design, "paired")) {
        .scd_stop("INVALID_DESIGN", "design must be independent or paired")
    }
    if (!is.character(contrast) || length(contrast) != 2L || !identical(names(
        contrast),
        c("numerator", "denominator")) || anyNA(contrast) || any(!nzchar(
            contrast)) ||
        identical(contrast[[1L]], contrast[[2L]])) {
        .scd_stop("INVALID_CONTRAST",
            "contrast must be c(numerator='...', denominator='...')")
    }
    if (!is.null(min_cells) && (!is.numeric(min_cells) || length(min_cells) !=
        1L || !is.finite(min_cells) || min_cells < 1 || min_cells != floor(
            min_cells))) {
        .scd_stop("INVALID_ARGUMENT",
            "min_cells must be NULL or a positive integer")
    }
    if (!identical(pair_policy, "strict") && !identical(pair_policy,
        "complete_pairs")) {
        .scd_stop("INVALID_ARGUMENT",
            "pair_policy must be strict or complete_pairs")
    }
    if (!is.numeric(min_replicates_warn) || length(min_replicates_warn) !=
        1L || !is.finite(min_replicates_warn) || min_replicates_warn <
        1 || min_replicates_warn != floor(min_replicates_warn)) {
        .scd_stop("INVALID_ARGUMENT", "min_replicates_warn must be positive")
    }
    if (!is.character(covariates) || anyNA(covariates) || any(!nzchar(
        covariates)) ||
        anyDuplicated(covariates)) {
        .scd_stop("INVALID_ARGUMENT", "covariates must be unique column names")
    }
    reserved <- c("sample_id", "donor_id", "condition", "cell_type", "n_cells",
        "pb_id")
    if (any(covariates %in% reserved)) {
        .scd_stop("INVALID_ARGUMENT", "covariates must not use a reserved name")
    }
    if (any(covariates == "." | make.names(covariates) != covariates)) {
        .scd_stop("INVALID_ARGUMENT",
            "covariates must be syntactic column names other than '.'")
    }
    list(design = design, contrast = contrast, min_cells = min_cells,
        covariates = covariates,
        pair_policy = pair_policy, min_replicates_warn = min_replicates_warn)
}
