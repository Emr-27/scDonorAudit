#' Prepare donor-aware pseudobulk counts
#'
#' Aggregates an in-memory `SingleCellExperiment` by sample and cell type, or
#' checks an already aggregated `SummarizedExperiment`. Missing sample/cell-type
#' combinations are recorded in `metadata(pb)$scdonoraudit$coverage`, not
#' inserted as zero-expression columns. A sample registry is needed to see
#' samples that are entirely absent from the input object.
#'
#' @param x A `SingleCellExperiment` of cells or a `SummarizedExperiment` with
#'   one column for each observed sample/cell-type combination.
#' @param sample_id,donor_id,condition,cell_type Names of metadata columns.
#'   If `sample_table` is supplied, `donor_id` and `condition` may be absent
#'   from the object; the table supplies them by `sample_id`.
#' @param assay Name of the raw count assay.
#' @param sample_vars Names of additional sample-level metadata columns to keep.
#' @param sample_table Optional `data.frame` with canonical columns `sample_id`,
#'   `donor_id`, `condition`, and selected `sample_vars`.
#' @param n_cells For already aggregated input, name of the metadata column
#'   containing the number of cells. It is computed for single-cell input.
#' @return A `SummarizedExperiment` with raw pseudobulk counts and audit metadata.
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
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
#' pb
#' @export
preparePseudobulk <- function(x, sample_id, donor_id, condition, cell_type,
                              assay = "counts", sample_vars = character(),
                              sample_table = NULL, n_cells = NULL) {
    if (!methods::is(x, "SummarizedExperiment")) {
        .scd_stop("INVALID_INPUT", "x must inherit from SummarizedExperiment")
    }
    is_cells <- methods::is(x, "SingleCellExperiment")
    for (arg in list(sample_id = sample_id, donor_id = donor_id,
                     condition = condition, cell_type = cell_type,
                     assay = assay)) {
        .scd_scalar_string(arg, "metadata mapping")
    }
    if (!is.character(sample_vars) || anyNA(sample_vars) ||
        any(!nzchar(sample_vars)) || anyDuplicated(sample_vars)) {
        .scd_stop("INVALID_ARGUMENT", "sample_vars must be unique column names")
    }
    if (!assay %in% SummarizedExperiment::assayNames(x)) {
        .scd_stop("MISSING_ASSAY", paste0("assay '", assay, "' is missing"))
    }
    if (ncol(x) == 0L || nrow(x) == 0L) {
        .scd_stop("EMPTY_INPUT", "input needs at least one gene and one column")
    }
    gene_ids <- .scd_gene_ids(x)
    row_order <- order(gene_ids, method = "radix")
    gene_ids <- gene_ids[row_order]
    counts <- SummarizedExperiment::assay(x, assay, withDimnames = TRUE)
    .scd_counts(counts, assay)
    counts <- counts[row_order, , drop = FALSE]
    source <- as.data.frame(SummarizedExperiment::colData(x))
    parsed <- .scd_sample_table(source, sample_id, donor_id, condition,
                                sample_vars, sample_table)
    registry <- parsed$table
    cell_types <- .scd_ids(.scd_field(source, cell_type, "cell_type"),
                           "cell_type")
    cell_levels <- sort(unique(cell_types), method = "radix")
    sample_index <- match(parsed$observed$sample_id, registry$sample_id)
    type_index <- match(cell_types, cell_levels)
    n_types <- length(cell_levels)
    source_keys <- (as.double(sample_index) - 1) * n_types + type_index
    observed_keys <- sort(unique(source_keys))
    if (is_cells) {
        group_index <- match(source_keys, observed_keys)
        membership <- Matrix::sparseMatrix(
            i = seq_len(ncol(x)), j = group_index,
            x = rep(1, ncol(x)),
            dims = c(ncol(x), length(observed_keys)))
        pb_counts <- counts %*% membership
        cell_counts <- tabulate(group_index, nbins = length(observed_keys))
    } else {
        if (anyDuplicated(source_keys)) {
            .scd_stop("DUPLICATE_PSEUDOBULK",
                "already aggregated input has duplicate sample/cell-type columns")
        }
        column_order <- order(source_keys)
        pb_counts <- counts[, column_order, drop = FALSE]
        if (is.null(n_cells)) {
            cell_counts <- rep(NA_integer_, ncol(x))
        } else {
            cell_counts <- .scd_field(source, n_cells, "n_cells")
            if (!is.numeric(cell_counts) || anyNA(cell_counts) ||
                any(!is.finite(cell_counts)) || any(cell_counts < 1) ||
                any(cell_counts != floor(cell_counts))) {
                .scd_stop("INVALID_CELL_COUNT",
                    "n_cells must be positive integer-valued for observed columns")
            }
            cell_counts <- as.integer(cell_counts)[column_order]
        }
    }
    observed_sample_index <- ((observed_keys - 1) %/% n_types) + 1
    observed_type_index <- ((observed_keys - 1) %% n_types) + 1
    pb_ids <- sprintf("pb_%05d", seq_along(observed_keys))
    rownames(pb_counts) <- gene_ids
    colnames(pb_counts) <- pb_ids
    column_data <- data.frame(
        pb_id = pb_ids,
        sample_id = registry$sample_id[observed_sample_index],
        donor_id = registry$donor_id[observed_sample_index],
        condition = registry$condition[observed_sample_index],
        cell_type = cell_levels[observed_type_index],
        n_cells = cell_counts,
        stringsAsFactors = FALSE)
    for (field in sample_vars) {
        column_data[[field]] <- registry[[field]][observed_sample_index]
    }
    full_sample_index <- rep(seq_len(nrow(registry)), each = n_types)
    full_type_index <- rep(seq_len(n_types), times = nrow(registry))
    full_keys <- (as.double(full_sample_index) - 1) * n_types + full_type_index
    matched <- match(full_keys, observed_keys)
    coverage <- data.frame(
        sample_id = registry$sample_id[full_sample_index],
        donor_id = registry$donor_id[full_sample_index],
        condition = registry$condition[full_sample_index],
        cell_type = cell_levels[full_type_index],
        observed = !is.na(matched),
        n_cells = cell_counts[matched],
        stringsAsFactors = FALSE)
    coverage$n_cells[!coverage$observed] <- NA_integer_
    output <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = pb_counts),
        rowData = SummarizedExperiment::rowData(x)[row_order, , drop = FALSE],
        colData = S4Vectors::DataFrame(column_data),
        metadata = list(scdonoraudit = list(
            schema_version = "0.2",
            sample_table = registry,
            coverage = coverage,
            preparation = list(input_type = if (is_cells) "single_cell" else
                "pseudobulk", source_assay = assay,
                audit_universe = parsed$universe,
                sample_vars = sample_vars))))
    output
}
