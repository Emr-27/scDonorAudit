storage_counts <- function(value = 24) {
    counts <- outer(seq_len(8), seq_len(8), function(i, j) 10 + 2 * i + j)
    counts[(row(counts) + col(counts)) %% 3L == 0L] <- 0
    counts[1L, 1L] <- value
    compressed <- methods::as(Matrix::Matrix(counts, sparse = TRUE),
        "generalMatrix")
    list(compressed = compressed,
        row_compressed = methods::as(compressed, "RsparseMatrix"),
        triplet = methods::as(compressed, "TsparseMatrix"),
        symmetric = Matrix::forceSymmetric(compressed),
        triangular = Matrix::triu(compressed),
        diagonal = Matrix::Diagonal(8, x = c(value, rep(24, 7))))
}

storage_input <- function(counts, paired = FALSE, cells = FALSE) {
    dimnames(counts) <- rep(list(paste0("g", seq_len(nrow(counts)))), 2)
    if (cells) {
        samples <- rep(paste0("s", seq_len(4)), each = 2)
        donors <- rep(paste0("d", seq_len(4)), each = 2)
        condition <- rep(c("ctrl", "stim"), each = 4)
        cell_type <- rep(c("T", "T", "T", "B"), 2)
    } else {
        samples <- paste0("s", seq_len(8))
        donors <- if (paired) rep(paste0("d", seq_len(4)), each = 2) else
            paste0("d", seq_len(8))
        condition <- if (paired) rep(c("ctrl", "stim"), 4) else
            rep(c("ctrl", "stim"), each = 4)
        cell_type <- rep("T", 8)
    }
    metadata <- S4Vectors::DataFrame(sample_id = samples, donor_id = donors,
        condition = condition, cell_type = cell_type)
    constructor <- if (cells) SingleCellExperiment::SingleCellExperiment else
        SummarizedExperiment::SummarizedExperiment
    constructor(assays = list(counts = counts), colData = metadata)
}

storage_prepare <- function(x) {
    preparePseudobulk(x, "sample_id", "donor_id", "condition", "cell_type")
}

test_that("numeric Matrix storage preserves aggregation and missing coverage", {
    cases <- storage_counts()
    cases$unit_diagonal <- Matrix::Diagonal(8)
    cases$unit_triangular <- methods::as(Matrix::Diagonal(8), "CsparseMatrix")
    cases$duplicate_triplet <- Matrix::sparseMatrix(
        i = c(seq_len(8), 1L, 1L, 8L),
        j = c(seq_len(8), 2L, 2L, 1L),
        x = c(rep(24, 8), 11, 13, 0), dims = c(8L, 8L), repr = "T")
    for (counts in cases) {
        sparse <- storage_prepare(storage_input(counts, cells = TRUE))
        dense <- storage_prepare(storage_input(as.matrix(counts), cells = TRUE))
        expect_equal(as.matrix(SummarizedExperiment::assay(sparse, "counts")),
            as.matrix(SummarizedExperiment::assay(dense, "counts")))
        expect_equal(SummarizedExperiment::colData(sparse),
            SummarizedExperiment::colData(dense))
        expect_identical(S4Vectors::metadata(sparse), S4Vectors::metadata(dense))
        expect_true(methods::is(SummarizedExperiment::assay(sparse, "counts"),
            "sparseMatrix"))
    }
})

test_that("Matrix storage remains usable at both public model entries", {
    contrast <- c(numerator = "stim", denominator = "ctrl")
    cases <- storage_counts()
    cases$unit_diagonal <- Matrix::Diagonal(8)
    cases$unit_triangular <- methods::as(Matrix::Diagonal(8), "CsparseMatrix")
    for (paired in c(FALSE, TRUE)) {
        design <- if (paired) "paired" else "independent"
        for (counts in cases) {
            pb <- storage_prepare(storage_input(as.matrix(counts), paired))
            SummarizedExperiment::assay(pb, "counts", withDimnames = FALSE) <-
                counts
            expect_true(all(as.data.frame(auditDesign(pb, design, contrast)[[
                "runs"]])$reason_code == "OK"))
            result <- assessDonorInfluence(pb, design, contrast)
            expect_identical(as.data.frame(result[["fits"]])$design_reason_code,
                rep("OK", if (paired) 5L else 9L))
        }
    }
})

test_that("invalid stored counts retain public errors across Matrix classes", {
    contrast <- c(numerator = "stim", denominator = "ctrl")
    prepared <- lapply(c(FALSE, TRUE), function(paired) {
        storage_prepare(storage_input(as.matrix(storage_counts()[[1L]]), paired))
    })
    for (value in c(NA_real_, NaN, Inf, -1, 0.5)) {
        for (counts in storage_counts(value)) {
            expect_error(storage_prepare(storage_input(counts, cells = TRUE)),
                "INVALID_COUNTS")
            for (i in seq_along(prepared)) {
                pb <- prepared[[i]]
                SummarizedExperiment::assay(pb, "counts", withDimnames = FALSE) <-
                    counts
                design <- c("independent", "paired")[[i]]
                expect_error(auditDesign(pb, design, contrast), "INVALID_COUNTS")
                expect_error(assessDonorInfluence(pb, design, contrast),
                    "INVALID_COUNTS")
            }
        }
    }
})

test_that("triplet duplicates do not conceal invalid stored counts", {
    for (values in list(c(-1, 2), c(0.5, 0.5))) {
        counts <- Matrix::sparseMatrix(i = c(1L, 1L), j = c(1L, 1L),
            x = values, dims = c(8L, 8L), repr = "T")
        expect_equal(as.matrix(counts)[1L, 1L], 1)
        expect_error(storage_prepare(storage_input(counts)), "INVALID_COUNTS")
    }
})
