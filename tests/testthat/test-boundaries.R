boundary_fixture <- function(paired = FALSE, collision_labels = FALSE) {
    set.seed(1103)
    donor <- if (paired) rep(paste0("d", 1:6), each = 2) else paste0("d", 1:12)
    donor[donor == "d1"] <- "baseline"
    condition <- if (paired) rep(c("ctrl", "stim"), 6) else
        rep(c("ctrl", "stim"), each = 6)
    sample_ids <- paste0("s", seq_along(donor))
    if (collision_labels) {
        donor[1:2] <- c("a / b", "a")
        sample_ids[1:2] <- c("c", "b / c")
    }
    samples <- data.frame(sample_id = sample_ids,
        donor_id = donor, condition = condition,
        age = rep(c(21, 35, 28, 42, 33, 26), 2),
        batch = rep(c("a", "b", "a"), 4))
    counts <- matrix(rnbinom(100 * nrow(samples), mu = 25, size = 5),
        nrow = 100, dimnames = list(paste0("g", 1:100), samples$sample_id))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(samples, cell_type = "T"))
    preparePseudobulk(se, "sample_id", "donor_id", "condition", "cell_type",
        sample_table = samples, sample_vars = c("age", "batch"))
}

test_that("both public model entries reject reserved and unsafe covariates", {
    contrast <- c(numerator = "ctrl", denominator = "stim")
    for (paired in c(FALSE, TRUE)) {
        pb <- boundary_fixture(paired)
        design <- if (paired) "paired" else "independent"
        for (entry in list(auditDesign, assessDonorInfluence)) {
            for (field in c("condition", "donor_id", "sample_id", "cell_type",
                            "n_cells", "pb_id")) {
                expect_error(entry(pb, design, contrast, covariates = field),
                    "INVALID_ARGUMENT.*reserved name")
            }
            for (field in c("age years", "age/batch", ".", "if", "age+batch")) {
                expect_error(entry(pb, design, contrast, covariates = field),
                    "INVALID_ARGUMENT.*syntactic column names")
            }
        }
    }
})

test_that("prepared objects are rechecked at both public model entries", {
    contrast <- c(numerator = "stim", denominator = "ctrl")
    for (paired in c(FALSE, TRUE)) {
        pb <- boundary_fixture(paired)
        design <- if (paired) "paired" else "independent"
        for (entry in list(auditDesign, assessDonorInfluence)) {
            for (sparse in c(FALSE, TRUE)) {
                for (value in c(0.5, -1, NA_real_, NaN, Inf)) {
                    bad <- pb
                    counts <- SummarizedExperiment::assay(bad, "counts")
                    if (sparse) counts <- Matrix::Matrix(counts, sparse = TRUE)
                    counts[1, 1] <- value
                    SummarizedExperiment::assay(bad, "counts") <- counts
                    expect_error(entry(bad, design, contrast), "INVALID_COUNTS")
                }
            }
            bad <- pb
            rownames(bad)[2] <- rownames(bad)[1]
            expect_error(entry(bad, design, contrast), "DUPLICATE_GENE_ID")
            rownames(bad) <- NULL
            expect_error(entry(bad, design, contrast), "MISSING_GENE_ID")
            rownames(bad) <- c("", rownames(pb)[-1])
            expect_error(entry(bad, design, contrast), "MISSING_ID")
            expect_error(entry(pb[0, ], design, contrast), "EMPTY_INPUT")
            expect_error(entry(pb[, 0], design, contrast), "EMPTY_INPUT")
            bad <- pb
            SummarizedExperiment::assayNames(bad) <- "renamed_counts"
            expect_error(entry(bad, design, contrast), "MISSING_ASSAY")
        }
    }
})

test_that("legal prepared-object edits and sparse counts remain usable", {
    contrast <- c(numerator = "stim", denominator = "ctrl")
    for (paired in c(FALSE, TRUE)) {
        pb <- boundary_fixture(paired)
        design <- if (paired) "paired" else "independent"
        pb <- pb[seq_len(50), rev(seq_len(ncol(pb)))]
        counts <- SummarizedExperiment::assay(pb, "counts")
        counts[1, 1] <- counts[1, 1] + 1
        for (sparse in c(FALSE, TRUE)) {
            changed <- pb
            SummarizedExperiment::assay(changed, "counts") <-
                if (sparse) Matrix::Matrix(counts, sparse = TRUE) else counts
            expect_no_error(auditDesign(changed, design, contrast))
            result <- assessDonorInfluence(changed, design, contrast)
            expect_true(all(as.data.frame(result[["fits"]])$
                execution_status == "completed"))
            expect_identical(rownames(result[["results"]][["T"]]), rownames(changed))
        }
    }
})

test_that("coverage positions use sample IDs despite colliding display labels", {
    result <- assessDonorInfluence(boundary_fixture(collision_labels = TRUE),
        "independent", c(numerator = "stim", denominator = "ctrl"))
    original <- result
    p <- plotInfluence(result, "coverage")
    built <- ggplot2::ggplot_build(p)
    rows <- which(p$data$sample == "a / b / c")
    expect_length(rows, 2)
    expect_equal(length(unique(built$data[[1]]$y[rows])), 2)
    expect_identical(as.character(p$data$sample_key), p$data$sample_id)
    expect_equal(length(built$layout$panel_params[[1]]$y$get_limits()), 12)
    expect_identical(result, original)
})

test_that("legal covariates retain the requested reverse contrast", {
    expect_no_error(scDonorAudit:::.scd_options("independent",
        c(numerator = "stim", denominator = "ctrl"), NULL,
        "age.years", "strict", 3L))
    pb <- boundary_fixture()
    forward <- assessDonorInfluence(pb, "independent",
        c(numerator = "stim", denominator = "ctrl"), covariates = c("age", "batch"))
    reverse <- assessDonorInfluence(pb, "independent",
        c(numerator = "ctrl", denominator = "stim"), covariates = c("age", "batch"))
    expect_true(all(as.data.frame(forward[["fits"]])$execution_status == "completed"))
    expect_equal(SummarizedExperiment::assay(reverse[["results"]][["T"]], "logFC"),
        -SummarizedExperiment::assay(forward[["results"]][["T"]], "logFC"),
        tolerance = 1e-8)
    paired <- assessDonorInfluence(boundary_fixture(TRUE), "paired",
        c(numerator = "ctrl", denominator = "stim"),
        covariates = c("age", "batch"))
    expect_true(all(as.data.frame(paired[["fits"]])$execution_status == "completed"))
})

test_that("effect plot separates baseline from a donor named baseline", {
    result <- assessDonorInfluence(boundary_fixture(), "independent",
        c(numerator = "stim", denominator = "ctrl"))
    se <- result[["results"]][["T"]]
    original <- SummarizedExperiment::assays(se)
    p <- plotInfluence(result, "effect", "T", "g1")
    expect_no_error(ggplot2::ggplot_build(p))
    expect_identical(levels(p$data$run), colnames(se))
    expect_equal(length(unique(p$data$run)), ncol(se))
    expect_identical(SummarizedExperiment::assays(result[["results"]][["T"]]), original)
    SummarizedExperiment::assay(se, "effect_valid")[1, 2] <- FALSE
    result[["results"]][["T"]] <- se
    p <- plotInfluence(result, "effect", "T", "g1")
    built <- ggplot2::ggplot_build(p)
    expect_true(is.na(p$data$logFC[2]))
    expect_equal(length(built$layout$panel_params[[1]]$x$get_limits()), ncol(se))
    expect_true("omit: baseline" %in% built$layout$panel_params[[1]]$x$get_labels())
})

counts_input_fixture <- function() {
    set.seed(2024)
    samples <- data.frame(sample_id = paste0("s", 1:8), donor_id = paste0("d",
        1:8), condition = rep(c("ctrl", "stim"), each = 4), stringsAsFactors = FALSE)
    counts <- matrix(rnbinom(40 * nrow(samples), mu = 30, size = 5), nrow = 40,
        dimnames = list(paste0("g", 1:40), samples$sample_id))
    SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts),
        colData = S4Vectors::DataFrame(samples, cell_type = "T"))
}

influence_fixture <- function(cell_cells = c(T = 100L, B = 1L)) {
    set.seed(5150)
    samples <- data.frame(sample_id = paste0("s", 1:8), donor_id = paste0("d",
        1:8), condition = rep(c("ctrl", "stim"), each = 4), stringsAsFactors = FALSE)
    genes <- paste0("g", 1:100)
    counts <- do.call(cbind, lapply(names(cell_cells), function(cell) matrix(rnbinom(length(genes) *
        nrow(samples), mu = 40, size = 5), nrow = length(genes), dimnames = list(genes,
        paste0(cell, "_", samples$sample_id)))))
    column_data <- do.call(rbind, lapply(names(cell_cells), function(cell) data.frame(sample_id = samples$sample_id,
        cell_type = cell, n_cells = cell_cells[[cell]], stringsAsFactors = FALSE)))
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts),
        colData = S4Vectors::DataFrame(column_data))
    preparePseudobulk(se, "sample_id", "donor_id", "condition", "cell_type",
        sample_table = samples, n_cells = "n_cells")
}

logical_sparse_counts <- function(counts) {
    dims <- dim(counts)
    Matrix::sparseMatrix(i = c(1L, 2L, dims[1L]), j = c(1L, 2L, dims[2L]),
        dims = dims, x = TRUE, dimnames = dimnames(counts))
}

pattern_sparse_counts <- function(counts) {
    methods::as(logical_sparse_counts(counts), "nMatrix")
}

influence_columns <- c("cell_type", "gene_id", "baseline_logFC", "baseline_p_value",
    "n_planned", "n_effect_valid", "n_test_valid", "max_abs_delta_observed",
    "min_logFC_observed", "max_logFC_observed", "n_sign_reversal", "n_material_reversal",
    "baseline_near_zero", "complete_effect_coverage", "complete_test_coverage",
    "max_influence_donors")

expect_empty_influence_schema <- function(x) {
    expect_s4_class(x, "DFrame")
    expect_identical(nrow(x), 0L)
    expect_identical(names(x), influence_columns)
    expect_identical(class(x$cell_type), "character")
    expect_identical(class(x$gene_id), "character")
    for (field in c("baseline_logFC", "baseline_p_value", "max_abs_delta_observed",
        "min_logFC_observed", "max_logFC_observed")) {
        expect_identical(class(x[[field]]), "numeric")
    }
    for (field in c("n_planned", "n_effect_valid", "n_test_valid", "n_sign_reversal",
        "n_material_reversal")) {
        expect_identical(class(x[[field]]), "integer")
    }
    for (field in c("baseline_near_zero", "complete_effect_coverage", "complete_test_coverage")) {
        expect_identical(class(x[[field]]), "logical")
    }
    expect_true(methods::is(x$max_influence_donors, "CharacterList"))
    expect_identical(length(x$max_influence_donors), 0L)
}

test_that("summarizeInfluence keeps its typed schema when no gene is retained",
    {
        result <- assessDonorInfluence(boundary_fixture(), "independent",
            c(numerator = "stim", denominator = "ctrl"), filter_args = list(min.count = 1000))
        expect_identical(nrow(result[["results"]][["T"]]), 0L)
        expect_empty_influence_schema(summarizeInfluence(result, effect_threshold = 0.5))
    })

test_that("summarizeInfluence keeps its typed schema when every cell type is blocked",
    {
        result <- assessDonorInfluence(influence_fixture(c(B = 1L)), "independent",
            c(numerator = "stim", denominator = "ctrl"), min_cells = 10L)
        expect_identical(nrow(result[["results"]][["B"]]), 0L)
        expect_empty_influence_schema(summarizeInfluence(result))
    })

test_that("summarizeInfluence keeps typed columns for mixed valid and blocked cell types",
    {
        influence <- summarizeInfluence(assessDonorInfluence(influence_fixture(),
            "independent", c(numerator = "stim", denominator = "ctrl"),
            min_cells = 10L), effect_threshold = 0.5)
        expect_true(nrow(influence) > 0L)
        expect_identical(unique(influence$cell_type), "T")
        expect_identical(names(influence), influence_columns)
        expect_true(methods::is(influence$max_influence_donors, "CharacterList"))
        expect_true(all(lengths(influence$max_influence_donors) >= 1L))
        expect_identical(class(influence$n_planned), "integer")
        expect_identical(class(influence$n_effect_valid), "integer")
        expect_identical(class(influence$max_abs_delta_observed), "numeric")
        expect_identical(class(influence$complete_effect_coverage), "logical")
        expect_true(all(influence$n_planned == 8L))
        blocked <- summarizeInfluence(assessDonorInfluence(influence_fixture(c(B = 1L)),
            "independent", c(numerator = "stim", denominator = "ctrl"),
            min_cells = 10L))
        expect_identical(names(influence), names(blocked))
        expect_identical(vapply(as.list(influence), function(col) class(col)[1],
            character(1)), vapply(as.list(blocked), function(col) class(col)[1],
            character(1)))
    })

test_that("numeric sparse counts are supported while logical and pattern counts are rejected",
    {
        contrast <- c(numerator = "stim", denominator = "ctrl")
        input <- counts_input_fixture()
        input_counts <- SummarizedExperiment::assay(input, "counts")
        sparse_input <- input
        SummarizedExperiment::assay(sparse_input, "counts") <- Matrix::Matrix(input_counts,
            sparse = TRUE)
        pb <- preparePseudobulk(sparse_input, "sample_id", "donor_id",
            "condition", "cell_type")
        expect_true(methods::is(SummarizedExperiment::assay(pb, "counts"),
            "sparseMatrix"))
        expect_no_error(auditDesign(pb, "independent", contrast))

        for (bad in list(logical_sparse_counts(input_counts), pattern_sparse_counts(input_counts),
            input_counts > 0)) {
            broken <- input
            SummarizedExperiment::assay(broken, "counts") <- bad
            expect_error(preparePseudobulk(broken, "sample_id", "donor_id",
                "condition", "cell_type"), "INVALID_COUNTS")
        }

        prepared <- boundary_fixture()
        prepared_counts <- SummarizedExperiment::assay(prepared, "counts")
        sparse_prepared <- prepared
        SummarizedExperiment::assay(sparse_prepared, "counts") <- Matrix::Matrix(prepared_counts,
            sparse = TRUE)
        expect_no_error(auditDesign(sparse_prepared, "independent", contrast))
        expect_true(all(as.data.frame(assessDonorInfluence(sparse_prepared,
            "independent", contrast)[["fits"]])$execution_status == "completed"))

        for (bad in list(logical_sparse_counts(prepared_counts), pattern_sparse_counts(prepared_counts))) {
            broken <- prepared
            SummarizedExperiment::assay(broken, "counts") <- bad
            expect_error(auditDesign(broken, "independent", contrast),
                "INVALID_COUNTS")
            expect_error(assessDonorInfluence(broken, "independent", contrast),
                "INVALID_COUNTS")
        }
    })
