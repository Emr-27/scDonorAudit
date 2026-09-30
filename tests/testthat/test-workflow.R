make_fixture <- function(paired = FALSE, missing_pair = FALSE) {
    set.seed(1103)
    donor <- if (paired) rep(paste0("d", 1:4), each = 2) else paste0("d", 1:8)
    condition <- if (paired) rep(c("ctrl", "stim"), 4) else
        rep(c("ctrl", "stim"), each = 4)
    samples <- paste0("s", seq_along(donor))
    sample_table <- data.frame(sample_id = samples, donor_id = donor,
                               condition = condition)
    keep <- if (missing_pair) seq_along(samples) != 8L else
        rep(TRUE, length(samples))
    values <- matrix(rnbinom(250 * sum(keep), mu = 25, size = 5),
        nrow = 250, dimnames = list(paste0("g", seq_len(250)), samples[keep]))
    col_data <- S4Vectors::DataFrame(sample_id = samples[keep],
        donor_id = donor[keep], condition = condition[keep],
        cell_type = rep("T", sum(keep)), n_cells = rep(30L, sum(keep)))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = values), colData = col_data)
    preparePseudobulk(se, "sample_id", "donor_id", "condition",
        "cell_type", sample_table = sample_table, n_cells = "n_cells")
}

test_that("independent workflow returns comparable donor refits", {
    pb <- make_fixture()
    contrast <- c(numerator = "stim", denominator = "ctrl")
    audit <- auditDesign(pb, "independent", contrast)
    expect_true(all(as.data.frame(audit[["runs"]])$reason_code == "OK"))
    result <- assessDonorInfluence(pb, "independent", contrast)
    fits <- as.data.frame(result[["fits"]])
    expect_equal(nrow(fits), 9L)
    expect_true(all(fits$execution_status == "completed"))
    observed <- result[["results"]][["T"]]
    expect_equal(ncol(observed), 9L)
    expect_equal(nrow(observed), 250L)
    p <- SummarizedExperiment::assay(observed, "p_value")[, 1L]
    padj <- SummarizedExperiment::assay(observed,
        "padj_within_cell_type")[, 1L]
    expect_equal(padj, stats::p.adjust(p, "BH", n = nrow(observed)))
    summary <- summarizeInfluence(result, effect_threshold = 0.5)
    expect_equal(nrow(summary), nrow(observed))
    expect_true(all(summary$n_planned == 8L))
    expect_true(all(summary$n_effect_valid == 8L))
})

test_that("pair policy records a blocked and a reduced cohort", {
    pb <- make_fixture(paired = TRUE, missing_pair = TRUE)
    contrast <- c(numerator = "stim", denominator = "ctrl")
    strict <- assessDonorInfluence(pb, "paired", contrast)
    expect_identical(as.data.frame(strict[["spec"]])$plan_status,
                     "not_created")
    expect_true(is.na(as.data.frame(strict[["spec"]])$n_planned))
    expect_equal(nrow(strict[["results"]][["T"]]), 0L)
    expect_equal(ncol(strict[["results"]][["T"]]), 1L)
    expect_identical(as.data.frame(strict[["fits"]])$reason_code,
                     "INCOMPLETE_PAIR")
    reduced <- assessDonorInfluence(pb, "paired", contrast,
                                    pair_policy = "complete_pairs")
    fits <- as.data.frame(reduced[["fits"]])
    expect_equal(nrow(fits), 4L)
    expect_true(all(fits$execution_status == "completed"))
    membership <- as.data.frame(reduced[["sample_membership"]])
    expect_true(any(membership$reason_code == "PAIR_MEMBER_EXCLUDED"))
    expect_true(any(membership$reason_code == "MISSING_CELL_TYPE"))
    expect_identical(S4Vectors::metadata(pb)$scdonoraudit$preparation$
        audit_universe, "registered_samples")
})

test_that("invalid raw counts are rejected", {
    pb <- make_fixture()
    original <- SummarizedExperiment::assay(pb, "counts")
    for (value in c(-1, Inf, NaN, 0.2)) {
        current <- original
        current[1L, 1L] <- value
        SummarizedExperiment::assay(pb, "counts") <- current
        expect_error(preparePseudobulk(pb, "sample_id", "donor_id",
            "condition", "cell_type"), "INVALID_COUNTS")
    }
    SummarizedExperiment::assay(pb, "counts") <- original
    rownames(pb)[2L] <- rownames(pb)[1L]
    expect_error(preparePseudobulk(pb, "sample_id", "donor_id",
        "condition", "cell_type"), "DUPLICATE")
})

test_that("sparse cells aggregate without inventing missing columns", {
    counts <- Matrix::Matrix(matrix(c(1, 0, 2, 3, 4, 0, 0, 5),
        nrow = 2, dimnames = list(c("g1", "g2"), NULL)), sparse = TRUE)
    cd <- S4Vectors::DataFrame(sample_id = c("s1", "s1", "s2", "s2"),
        donor_id = c("d1", "d1", "d2", "d2"),
        condition = c("ctrl", "ctrl", "stim", "stim"),
        cell_type = c("T", "T", "T", "B"))
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts), colData = cd)
    pb <- preparePseudobulk(sce, "sample_id", "donor_id", "condition",
                            "cell_type")
    expect_equal(ncol(pb), 3L)
    expect_equal(as.matrix(SummarizedExperiment::assay(pb, "counts")),
        matrix(c(3, 3, 0, 5, 4, 0), nrow = 2,
               dimnames = list(c("g1", "g2"),
                               paste0("pb_0000", 1:3))))
    cover <- S4Vectors::metadata(pb)$scdonoraudit$coverage
    expect_equal(sum(!cover$observed), 1L)
    expect_identical(S4Vectors::metadata(pb)$scdonoraudit$preparation$
        audit_universe, "observed_only")
})

test_that("summary retains ties and the planned denominator", {
    pb <- make_fixture()
    result <- assessDonorInfluence(pb, "independent",
        c(numerator = "stim", denominator = "ctrl"))
    se <- result[["results"]][["T"]]
    for (field in c("logFC", "effect_valid", "test_valid")) {
        matrix <- SummarizedExperiment::assay(se, field)
        matrix[1L, ] <- if (field == "logFC")
            c(0.8, 0.9, 0.7, -0.6, rep(NA_real_, 5L)) else
            c(TRUE, TRUE, TRUE, TRUE, rep(FALSE, 5L))
        SummarizedExperiment::assay(se, field) <- matrix
    }
    result[["results"]][["T"]] <- se
    summary <- summarizeInfluence(result, effect_threshold = 0.5)
    first <- which(summary$gene_id == rownames(se)[1L])
    expect_equal(summary$n_planned[first], 8L)
    expect_equal(summary$n_effect_valid[first], 3L)
    expect_equal(summary$max_abs_delta_observed[first], 1.4)
    expect_equal(summary$n_material_reversal[first], 1L)
    expect_false(summary$complete_effect_coverage[first])
    expect_equal(length(summary$max_influence_donors[[first]]), 1L)
    effects <- SummarizedExperiment::assay(se, "logFC")
    valid <- SummarizedExperiment::assay(se, "effect_valid")
    effects[2L, ] <- NA_real_
    valid[2L, ] <- FALSE
    effects[3L, ] <- c(0.01, -0.01, rep(NA_real_, 7L))
    valid[3L, ] <- c(TRUE, TRUE, rep(FALSE, 7L))
    effects[1L, 4L] <- -0.6
    effects[1L, 5L] <- -0.6
    valid[1L, 5L] <- TRUE
    SummarizedExperiment::assay(se, "logFC") <- effects
    SummarizedExperiment::assay(se, "effect_valid") <- valid
    result[["results"]][["T"]] <- se
    summary <- summarizeInfluence(result, effect_threshold = 0.5)
    first <- which(summary$gene_id == rownames(se)[1L])
    second <- which(summary$gene_id == rownames(se)[2L])
    third <- which(summary$gene_id == rownames(se)[3L])
    expect_equal(summary$n_effect_valid[first], 4L)
    expect_identical(as.character(summary$max_influence_donors[[first]]),
                     c("d3", "d4"))
    expect_equal(summary$n_effect_valid[second], 0L)
    expect_true(is.na(summary$max_abs_delta_observed[second]))
    expect_length(summary$max_influence_donors[[second]], 0L)
    expect_true(summary$baseline_near_zero[third])
    expect_equal(summary$n_sign_reversal[third], 1L)
    expect_equal(summary$n_material_reversal[third], 0L)
})

test_that("low cell count excludes both sides of an incomplete pair", {
    pb <- make_fixture(paired = TRUE)
    cd <- SummarizedExperiment::colData(pb)
    cd$n_cells[1L] <- 19L
    SummarizedExperiment::colData(pb) <- cd
    source <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = SummarizedExperiment::assay(pb, "counts")),
        colData = cd)
    registry <- S4Vectors::metadata(pb)$scdonoraudit$sample_table
    prepared <- preparePseudobulk(source, "sample_id", "donor_id",
        "condition", "cell_type", sample_table = registry,
        n_cells = "n_cells")
    contrast <- c(numerator = "stim", denominator = "ctrl")
    strict <- auditDesign(prepared, "paired", contrast, min_cells = 20L)
    expect_identical(as.data.frame(strict[["spec"]])$plan_status,
                     "not_created")
    reduced <- auditDesign(prepared, "paired", contrast, min_cells = 20L,
                           pair_policy = "complete_pairs")
    members <- as.data.frame(reduced[["sample_membership"]])
    expect_identical(members$reason_code[members$sample_id == "s1"],
                     "BELOW_MIN_CELLS")
    expect_identical(members$reason_code[members$sample_id == "s2"],
                     "PAIR_MEMBER_EXCLUDED")
    expect_equal(as.data.frame(reduced[["spec"]])$n_planned, 3L)
})

test_that("one blocked cell type does not suppress a valid cell type", {
    pb <- make_fixture(paired = TRUE)
    values <- SummarizedExperiment::assay(pb, "counts")
    combined_values <- cbind(values, values[, 1:7])
    colnames(combined_values) <- paste0("column", seq_len(ncol(combined_values)))
    cd <- as.data.frame(SummarizedExperiment::colData(pb))
    b_cd <- cd[1:7, , drop = FALSE]
    b_cd$cell_type <- "B"
    combined_cd <- rbind(cd, b_cd)
    rownames(combined_cd) <- NULL
    registry <- S4Vectors::metadata(pb)$scdonoraudit$sample_table
    source <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = combined_values),
        colData = S4Vectors::DataFrame(combined_cd))
    prepared <- preparePseudobulk(source, "sample_id", "donor_id",
        "condition", "cell_type", sample_table = registry,
        n_cells = "n_cells")
    result <- assessDonorInfluence(prepared, "paired",
        c(numerator = "stim", denominator = "ctrl"))
    fits <- as.data.frame(result[["fits"]])
    expect_equal(nrow(fits[fits$cell_type == "B", ]), 1L)
    expect_equal(nrow(fits[fits$cell_type == "T", ]), 5L)
    expect_identical(fits$reason_code[fits$cell_type == "B"],
                     "INCOMPLETE_PAIR")
    expect_true(all(fits$execution_status[fits$cell_type == "T"] ==
                    "completed"))
})

test_that("baseline and donor refit agree with a direct edgeR analysis", {
    pb <- make_fixture()
    result <- assessDonorInfluence(pb, "independent",
        c(numerator = "stim", denominator = "ctrl"))
    observed <- result[["results"]][["T"]]
    counts <- SummarizedExperiment::assay(pb, "counts")
    info <- as.data.frame(SummarizedExperiment::colData(pb))
    kept <- edgeR::filterByExpr(as.matrix(counts), group = info$condition,
        min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7)
    expect_identical(rownames(observed), rownames(pb)[kept])
    omitted <- as.character(SummarizedExperiment::colData(observed)$
        omitted_donor)
    for (j in seq_len(ncol(observed))) {
        use <- if (j == 1L) rep(TRUE, nrow(info)) else
            info$donor_id != omitted[j]
        y <- edgeR::DGEList(counts = as.matrix(counts[kept, use, drop = FALSE]))
        y <- edgeR::normLibSizes(y, method = "TMM")
        group <- factor(info$condition[use], levels = c("ctrl", "stim"))
        design <- stats::model.matrix(~ group)
        fit <- edgeR::glmQLFit(y, design = design, dispersion = NULL,
            abundance.trend = TRUE, robust = TRUE, legacy = FALSE,
            top.proportion = NULL)
        reference <- edgeR::glmQLFTest(fit, contrast = c(0, 1))$table
        expect_equal(as.numeric(SummarizedExperiment::assay(observed,
            "logFC")[, j]), reference$logFC, tolerance = 1e-8)
        expect_equal(as.numeric(SummarizedExperiment::assay(observed,
            "p_value")[, j]), reference$PValue, tolerance = 1e-8)
    }
})

test_that("all-zero-after-deletion keeps a gene row with missing values", {
    pb <- make_fixture()
    counts <- SummarizedExperiment::assay(pb, "counts")
    counts[1L, ] <- 0
    counts[1L, 1L] <- 1000
    SummarizedExperiment::assay(pb, "counts") <- counts
    result <- assessDonorInfluence(pb, "independent",
        c(numerator = "stim", denominator = "ctrl"),
        filter_args = list(min.count = 0, min.total.count = 1))
    se <- result[["results"]][["T"]]
    expect_true(rownames(pb)[1L] %in% rownames(se))
    donor <- as.character(SummarizedExperiment::colData(se)$omitted_donor)
    column <- which(donor == "d1")
    gene <- match(rownames(pb)[1L], rownames(se))
    expect_equal(SummarizedExperiment::assay(se,
        "gene_status_code")[gene, column], 1L)
    expect_false(SummarizedExperiment::assay(se,
        "effect_valid")[gene, column])
    expect_true(is.na(SummarizedExperiment::assay(se,
        "logFC")[gene, column]))
})

test_that("paired refit after removing the reference donor matches edgeR", {
    pb <- make_fixture(paired = TRUE)
    result <- assessDonorInfluence(pb, "paired",
        c(numerator = "stim", denominator = "ctrl"))
    observed <- result[["results"]][["T"]]
    counts <- SummarizedExperiment::assay(pb, "counts")
    info <- as.data.frame(SummarizedExperiment::colData(pb))
    kept <- edgeR::filterByExpr(as.matrix(counts), group = info$condition,
        min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7)
    omitted <- as.character(SummarizedExperiment::colData(observed)$
        omitted_donor)
    for (j in seq_len(ncol(observed))) {
        use <- if (j == 1L) rep(TRUE, nrow(info)) else
            info$donor_id != omitted[j]
        donor_id <- factor(info$donor_id[use])
        condition <- factor(info$condition[use],
                            levels = c("ctrl", "stim"))
        design <- stats::model.matrix(~ donor_id + condition)
        y <- edgeR::DGEList(counts = as.matrix(counts[kept, use, drop = FALSE]))
        y <- edgeR::normLibSizes(y, method = "TMM")
        fit <- edgeR::glmQLFit(y, design = design, dispersion = NULL,
            abundance.trend = TRUE, robust = TRUE, legacy = FALSE,
            top.proportion = NULL)
        contrast <- as.numeric(colnames(design) == "conditionstim")
        reference <- edgeR::glmQLFTest(fit, contrast = contrast)$table
        expect_equal(as.numeric(SummarizedExperiment::assay(observed,
            "logFC")[, j]), reference$logFC, tolerance = 1e-8)
        expect_equal(as.numeric(SummarizedExperiment::assay(observed,
            "p_value")[, j]), reference$PValue, tolerance = 1e-8)
    }
})

test_that("rank deficient nuisance term is reported without dropping it", {
    pb <- make_fixture()
    cd <- SummarizedExperiment::colData(pb)
    cd$constant <- rep(1, ncol(pb))
    SummarizedExperiment::colData(pb) <- cd
    registry <- S4Vectors::metadata(pb)$scdonoraudit$sample_table
    registry$constant <- 1
    source <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = SummarizedExperiment::assay(pb, "counts")),
        colData = cd)
    prepared <- preparePseudobulk(source, "sample_id", "donor_id",
        "condition", "cell_type", sample_vars = "constant",
        sample_table = registry, n_cells = "n_cells")
    audit <- auditDesign(prepared, "independent",
        c(numerator = "stim", denominator = "ctrl"),
        covariates = "constant")
    runs <- as.data.frame(audit[["runs"]])
    expect_true(all(runs$reason_code == "DESIGN_RANK_DEFICIENT"))
    expect_true(all(runs$target_estimable))
    result <- assessDonorInfluence(prepared, "independent",
        c(numerator = "stim", denominator = "ctrl"),
        covariates = "constant")
    expect_true(all(as.data.frame(result[["fits"]])$execution_status ==
                    "skipped"))
})

test_that("design audit distinguishes confounding, no df, missing group and zero library", {
    pb <- make_fixture()
    contrast <- c(numerator = "stim", denominator = "ctrl")
    registry <- S4Vectors::metadata(pb)$scdonoraudit$sample_table
    registry$group_alias <- as.integer(registry$condition == "stim")
    source <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = SummarizedExperiment::assay(pb, "counts")),
        colData = SummarizedExperiment::colData(pb))
    confounded <- preparePseudobulk(source, "sample_id", "donor_id",
        "condition", "cell_type", sample_vars = "group_alias",
        sample_table = registry, n_cells = "n_cells")
    runs <- as.data.frame(auditDesign(confounded, "independent", contrast,
        covariates = "group_alias")[["runs"]])
    expect_identical(runs$reason_code[1L], "DESIGN_RANK_DEFICIENT")
    expect_false(runs$target_estimable[1L])

    selected <- c(1L, 5L)
    small <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = SummarizedExperiment::assay(pb, "counts")[,
            selected, drop = FALSE]),
        colData = SummarizedExperiment::colData(pb)[selected, ])
    two_samples <- preparePseudobulk(small, "sample_id", "donor_id",
        "condition", "cell_type", sample_table = registry[selected,
            c("sample_id", "donor_id", "condition")], n_cells = "n_cells")
    runs <- as.data.frame(auditDesign(two_samples, "independent",
        contrast)[["runs"]])
    expect_identical(runs$reason_code[1L], "NO_RESIDUAL_DF")
    expect_true(all(runs$reason_code[-1L] == "MISSING_CONDITION"))

    zero <- SummarizedExperiment::assay(pb, "counts")
    zero[, 1L] <- 0
    SummarizedExperiment::assay(pb, "counts") <- zero
    runs <- as.data.frame(auditDesign(pb, "independent",
        contrast)[["runs"]])
    expect_identical(runs$reason_code[1L], "ZERO_LIBRARY")
})

test_that("cell threshold requires observed cell counts", {
    pb <- make_fixture()
    source <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = SummarizedExperiment::assay(pb, "counts")),
        colData = SummarizedExperiment::colData(pb))
    registry <- S4Vectors::metadata(pb)$scdonoraudit$sample_table
    unknown <- preparePseudobulk(source, "sample_id", "donor_id",
        "condition", "cell_type", sample_table = registry)
    expect_error(auditDesign(unknown, "independent",
        c(numerator = "stim", denominator = "ctrl"), min_cells = 20),
        "CELL_COUNT_UNKNOWN")
})

test_that("numeric covariates are finite and thresholds cannot be infinite", {
    pb <- make_fixture()
    contrast <- c(numerator = "stim", denominator = "ctrl")
    registry <- S4Vectors::metadata(pb)$scdonoraudit$sample_table
    registry$age <- seq_len(nrow(registry)) + 30
    source <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = SummarizedExperiment::assay(pb, "counts")),
        colData = SummarizedExperiment::colData(pb))
    prepared <- preparePseudobulk(source, "sample_id", "donor_id",
        "condition", "cell_type", sample_vars = "age",
        sample_table = registry, n_cells = "n_cells")
    expect_true(all(as.data.frame(auditDesign(prepared, "independent",
        contrast, covariates = "age")[["runs"]])$reason_code == "OK"))
    registry$age[1L] <- Inf
    invalid <- preparePseudobulk(source, "sample_id", "donor_id",
        "condition", "cell_type", sample_vars = "age",
        sample_table = registry, n_cells = "n_cells")
    expect_error(auditDesign(invalid, "independent", contrast,
        covariates = "age"), "INVALID_COVARIATE")
    expect_error(auditDesign(pb, "independent", contrast,
        min_cells = Inf), "INVALID_ARGUMENT")
    expect_error(auditDesign(pb, "independent", contrast,
        min_replicates_warn = Inf), "INVALID_ARGUMENT")
})

test_that("one-condition-zero genes retain an explicit diagnostic code", {
    pb <- make_fixture()
    counts <- SummarizedExperiment::assay(pb, "counts")
    condition <- as.character(SummarizedExperiment::colData(pb)$condition)
    counts[1L, condition == "ctrl"] <- 100L
    counts[1L, condition == "stim"] <- 0L
    SummarizedExperiment::assay(pb, "counts") <- counts
    result <- assessDonorInfluence(pb, "independent",
        c(numerator = "stim", denominator = "ctrl"))
    observed <- result[["results"]][["T"]]
    expect_true(all(SummarizedExperiment::assay(observed,
        "gene_status_code")["g1", ] == 4L))
    expect_true(all(SummarizedExperiment::assay(observed,
        "effect_valid")["g1", ]))
})

test_that("unavailable baseline retains planned deletions and missing effects", {
    pb <- make_fixture()
    counts <- SummarizedExperiment::assay(pb, "counts")
    counts[,] <- 0
    counts[1L, ] <- 1
    SummarizedExperiment::assay(pb, "counts") <- counts
    result <- assessDonorInfluence(pb, "independent",
        c(numerator = "stim", denominator = "ctrl"))
    fits <- as.data.frame(result[["fits"]])
    expect_identical(fits$reason_code[1L], "NO_GENES_AFTER_FILTER")
    expect_true(all(fits$reason_code[-1L] == "BASELINE_UNAVAILABLE"))
    expect_equal(ncol(result[["results"]][["T"]]), 9L)
    expect_equal(nrow(result[["results"]][["T"]]), 0L)
    expect_equal(as.data.frame(result[["spec"]])$n_planned, 8L)
})

test_that("partial missing p values retain the fixed BH family", {
    adjusted <- scDonorAudit:::.scd_adjust_p(
        c(0.01, 0.04, NA_real_), c(TRUE, TRUE, TRUE))
    expect_identical(adjusted$test_valid, c(TRUE, TRUE, FALSE))
    expect_equal(adjusted$adjusted, c(0.03, 0.06, NA_real_))
    no_effect <- scDonorAudit:::.scd_adjust_p(c(0.01, 0.04, 0.5),
                                              c(FALSE, TRUE, TRUE))
    expect_false(no_effect$test_valid[1L])
    expect_true(is.na(no_effect$adjusted[1L]))
})

test_that("a partial invalid backend p value preserves effect and status", {
    counts <- matrix(c(3, 2, 0, 4, 2, 0, 5, 2, 0, 6, 2, 0),
        nrow = 3L)
    data <- data.frame(condition = c("ctrl", "ctrl", "stim", "stim"))
    tab <- data.frame(logFC = c(1, -1), PValue = c(0.01, NA_real_))
    output <- scDonorAudit:::.scd_result_from_table(tab, counts, data,
        active = c(1L, 2L), family_size = 3L)
    expect_identical(output$effect_valid, c(TRUE, TRUE, FALSE))
    expect_identical(output$test_valid, c(TRUE, FALSE, FALSE))
    expect_identical(output$code, c(0L, 3L, 1L))
    expect_equal(output$padj, c(0.03, NA_real_, NA_real_))
    expect_equal(output$effect[2L], -1)
})

test_that("input order is canonical and contrast reversal changes effect sign", {
    pb <- make_fixture()
    contrast <- c(numerator = "stim", denominator = "ctrl")
    original <- assessDonorInfluence(pb, "independent", contrast)
    set.seed(809)
    shuffled <- pb[sample(seq_len(nrow(pb))), sample(seq_len(ncol(pb)))]
    registry <- S4Vectors::metadata(pb)$scdonoraudit$sample_table
    reordered <- preparePseudobulk(shuffled, "sample_id", "donor_id",
        "condition", "cell_type", sample_table = registry,
        n_cells = "n_cells")
    expect_equal(SummarizedExperiment::assay(reordered, "counts"),
                 SummarizedExperiment::assay(pb, "counts"))
    replay <- assessDonorInfluence(reordered, "independent", contrast)
    expect_equal(SummarizedExperiment::assay(
        original[["results"]][["T"]], "logFC"),
        SummarizedExperiment::assay(replay[["results"]][["T"]], "logFC"))
    reversed <- assessDonorInfluence(pb, "independent",
        c(numerator = "ctrl", denominator = "stim"))
    expect_equal(SummarizedExperiment::assay(
        original[["results"]][["T"]], "logFC"),
        -SummarizedExperiment::assay(reversed[["results"]][["T"]], "logFC"),
        tolerance = 1e-8)
    expect_equal(SummarizedExperiment::assay(
        original[["results"]][["T"]], "p_value"),
        SummarizedExperiment::assay(reversed[["results"]][["T"]], "p_value"),
        tolerance = 1e-8)
    forward <- SummarizedExperiment::assay(
        original[["results"]][["T"]], "logFC")
    backward <- SummarizedExperiment::assay(
        reversed[["results"]][["T"]], "logFC")
    expect_equal(sweep(forward[, -1L], 1L, forward[, 1L], "-"),
        -sweep(backward[, -1L], 1L, backward[, 1L], "-"),
        tolerance = 1e-8)
})

test_that("named contrast has a fixed meaning under user factor settings", {
    old <- options("contrasts")
    on.exit(options(old), add = TRUE)
    for (paired in c(FALSE, TRUE)) {
        data <- data.frame(
            donor_id = if (paired) rep(paste0("d", 1:4), each = 2) else
                paste0("d", 1:8),
            condition = if (paired) rep(c("ctrl", "stim"), 4) else
                rep(c("ctrl", "stim"), each = 4),
            batch = if (paired) c("A", "A", "A", "B", "B", "A", "B", "B")
                else rep(c("A", "B"), 4))
        spec <- scDonorAudit:::.scd_covariate_spec(data, "batch")
        config <- list(design = if (paired) "paired" else "independent",
            contrast = c(numerator = "stim", denominator = "ctrl"),
            covariates = "batch")
        response <- 2 * (data$condition == "stim") +
            3 * (data$batch == "B") +
            if (paired) rep(seq_len(4), each = 2) else 0
        designs <- lapply(c("contr.treatment", "contr.sum"), function(code) {
            options(contrasts = c(code, "contr.poly"))
            before <- getOption("contrasts")
            design <- scDonorAudit:::.scd_make_design(data, config, spec)
            expect_identical(getOption("contrasts"), before)
            design
        })
        expect_identical(designs[[1L]]$status, "OK")
        expect_equal(designs[[1L]]$matrix, designs[[2L]]$matrix)
        for (design in designs) {
            coefficients <- stats::lm.fit(design$matrix, response)$coefficients
            expect_equal(sum(design$contrast * coefficients), 2,
                         tolerance = 1e-10)
        }
        config$contrast <- c(numerator = "ctrl", denominator = "stim")
        reversed <- scDonorAudit:::.scd_make_design(data, config, spec)
        coefficients <- stats::lm.fit(reversed$matrix, response)$coefficients
        expect_equal(sum(reversed$contrast * coefficients), -2,
                     tolerance = 1e-10)
    }
})

test_that("baseline and donor deletions ignore global contrast options", {
    old <- options("contrasts")
    on.exit(options(old), add = TRUE)
    contrast <- c(numerator = "stim", denominator = "ctrl")
    for (paired in c(FALSE, TRUE)) {
        pb <- make_fixture(paired = paired)
        covariates <- character()
        if (!paired) {
            registry <- S4Vectors::metadata(pb)$scdonoraudit$sample_table
            registry$batch <- rep(c("A", "B"), 4)
            cd <- SummarizedExperiment::colData(pb)
            cd$batch <- registry$batch
            source <- SummarizedExperiment::SummarizedExperiment(
                assays = list(counts = SummarizedExperiment::assay(pb,
                    "counts")), colData = cd)
            pb <- preparePseudobulk(source, "sample_id", "donor_id",
                "condition", "cell_type", sample_vars = "batch",
                sample_table = registry, n_cells = "n_cells")
            covariates <- "batch"
        }
        design <- if (paired) "paired" else "independent"
        options(contrasts = c("contr.treatment", "contr.poly"))
        reference <- assessDonorInfluence(pb, design, contrast,
                                          covariates = covariates)
        options(contrasts = c("contr.sum", "contr.poly"))
        before <- getOption("contrasts")
        changed <- assessDonorInfluence(pb, design, contrast,
                                        covariates = covariates)
        expect_identical(getOption("contrasts"), before)
        for (field in c("logFC", "p_value")) {
            expect_equal(SummarizedExperiment::assay(
                reference[["results"]][["T"]], field),
                SummarizedExperiment::assay(
                    changed[["results"]][["T"]], field),
                tolerance = 1e-8)
        }
        expect_identical(as.data.frame(reference[["fits"]])$design_columns,
                         as.data.frame(changed[["fits"]])$design_columns)
        reversed <- assessDonorInfluence(pb, design,
            c(numerator = "ctrl", denominator = "stim"),
            covariates = covariates)
        expect_equal(SummarizedExperiment::assay(
            changed[["results"]][["T"]], "logFC"),
            -SummarizedExperiment::assay(
                reversed[["results"]][["T"]], "logFC"),
            tolerance = 1e-8)
    }
})

test_that("a planted one-donor shift is traceable through the fit ledger", {
    set.seed(672)
    donor <- rep(paste0("d", 1:4), each = 2)
    condition <- rep(c("ctrl", "stim"), 4)
    sample <- paste0("s", seq_along(donor))
    counts <- matrix(stats::rnbinom(250 * 8, mu = 25, size = 5),
        nrow = 250, dimnames = list(paste0("g", 1:250), sample))
    counts["g1", "s8"] <- counts["g1", "s8"] * 20
    registry <- data.frame(sample_id = sample, donor_id = donor,
                           condition = condition)
    source <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        colData = S4Vectors::DataFrame(registry, cell_type = "T",
                                        n_cells = 30L))
    pb <- preparePseudobulk(source, "sample_id", "donor_id", "condition",
        "cell_type", sample_table = registry, n_cells = "n_cells")
    result <- assessDonorInfluence(pb, "paired",
        c(numerator = "stim", denominator = "ctrl"))
    fits <- as.data.frame(result[["fits"]])
    expect_true(all(fits$execution_status == "completed"))
    expect_equal(nrow(fits), 5L)
    summary <- summarizeInfluence(result)
    row <- which(summary$gene_id == "g1")
    expect_length(row, 1L)
    expect_equal(summary$n_planned[row], 4L)
    expect_equal(summary$n_effect_valid[row], 4L)
    expect_identical(as.character(summary$max_influence_donors[[row]]),
                     "d4")
    effects <- SummarizedExperiment::assay(
        result[["results"]][["T"]], "logFC")["g1", ]
    expect_equal(summary$max_abs_delta_observed[row],
                 abs(effects["omit_004"] - effects["baseline"]),
                 tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("serialized results preserve state and plots build", {
    pb <- make_fixture()
    result <- assessDonorInfluence(pb, "independent",
        c(numerator = "stim", denominator = "ctrl"))
    file <- tempfile(fileext = ".rds")
    saveRDS(result, file)
    loaded <- readRDS(file)
    expect_identical(as.data.frame(result[["fits"]]),
                     as.data.frame(loaded[["fits"]]))
    expect_equal(SummarizedExperiment::assay(
        result[["results"]][["T"]], "logFC"),
        SummarizedExperiment::assay(
            loaded[["results"]][["T"]], "logFC"))
    for (type in c("coverage", "influence", "effect")) {
        plot <- if (type == "coverage") plotInfluence(loaded, type) else
            plotInfluence(loaded, type, cell_type = "T", gene_id = "g1")
        expect_s3_class(plot, "ggplot")
        expect_silent(ggplot2::ggplot_build(plot))
    }
})

test_that("one backend failure does not erase later donor runs", {
    pb <- make_fixture()
    original_fit <- scDonorAudit:::.scd_fit_run
    calls <- 0L
    testthat::local_mocked_bindings(
        .scd_fit_run = function(...) {
            calls <<- calls + 1L
            if (calls == 2L) stop("injected failure")
            original_fit(...)
        }, .package = "scDonorAudit")
    result <- assessDonorInfluence(pb, "independent",
        c(numerator = "stim", denominator = "ctrl"))
    fits <- as.data.frame(result[["fits"]])
    expect_identical(fits$execution_status[1:3],
                     c("completed", "failed", "completed"))
    expect_identical(fits$reason_code[2L], "BACKEND_ERROR")
    expect_equal(nrow(fits), 9L)
    observed <- result[["results"]][["T"]]
    expect_true(all(!SummarizedExperiment::assay(observed,
        "effect_valid")[, 2L]))
    expect_true(all(SummarizedExperiment::assay(observed,
        "effect_valid")[, 3L]))
    summary <- summarizeInfluence(result)
    expect_true(all(summary$n_planned == 8L))
    expect_true(all(summary$n_effect_valid == 7L))
    expect_false(any(summary$complete_effect_coverage))
    issues <- as.data.frame(result[["issues"]])
    expect_true(any(grepl("injected failure", issues$message,
                          fixed = TRUE)))
    plot <- plotInfluence(result, "influence", cell_type = "T",
        gene_id = "g1")
    expect_true(is.na(plot$data$delta_logFC[1L]))
    expect_true(any(is.finite(plot$data$delta_logFC[-1L])))
})

test_that("backend warnings are recorded while the fit completes", {
    pb <- make_fixture()
    original_fit <- scDonorAudit:::.scd_fit_run
    calls <- 0L
    testthat::local_mocked_bindings(
        .scd_fit_run = function(...) {
            calls <<- calls + 1L
            if (calls == 2L) warning("injected backend warning")
            original_fit(...)
        }, .package = "scDonorAudit")
    result <- assessDonorInfluence(pb, "independent",
        c(numerator = "stim", denominator = "ctrl"))
    fits <- as.data.frame(result[["fits"]])
    issues <- as.data.frame(result[["issues"]])
    expect_true(all(fits$execution_status == "completed"))
    expect_true(any(issues$reason_code == "BACKEND_WARNING" &
        grepl("injected backend warning", issues$message, fixed = TRUE)))
})

test_that("cell-level and preaggregated inputs agree on the same counts", {
    baseline <- make_fixture()
    sample_counts <- SummarizedExperiment::assay(baseline, "counts")
    registry <- S4Vectors::metadata(baseline)$scdonoraudit$sample_table
    set.seed(71)
    first <- matrix(stats::rbinom(length(sample_counts),
        size = as.vector(sample_counts), prob = 0.5),
        nrow = nrow(sample_counts))
    second <- sample_counts - first
    cell_counts <- cbind(first, second)
    rownames(cell_counts) <- rownames(sample_counts)
    colnames(cell_counts) <- paste0("cell", seq_len(ncol(cell_counts)))
    cd <- as.data.frame(SummarizedExperiment::colData(baseline))
    cells_cd <- rbind(cd, cd)
    rownames(cells_cd) <- NULL
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = Matrix::Matrix(cell_counts, sparse = TRUE)),
        colData = S4Vectors::DataFrame(cells_cd))
    from_cells <- preparePseudobulk(sce, "sample_id", "donor_id",
        "condition", "cell_type", sample_table = registry)
    source <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = sample_counts),
        colData = SummarizedExperiment::colData(baseline))
    SummarizedExperiment::colData(source)$n_cells <- 2L
    from_bulk <- preparePseudobulk(source, "sample_id", "donor_id",
        "condition", "cell_type", sample_table = registry,
        n_cells = "n_cells")
    expect_equal(as.matrix(SummarizedExperiment::assay(from_cells, "counts")),
                 SummarizedExperiment::assay(from_bulk, "counts"))
    expect_equal(as.data.frame(SummarizedExperiment::colData(from_cells)),
                 as.data.frame(SummarizedExperiment::colData(from_bulk)))
    contrast <- c(numerator = "stim", denominator = "ctrl")
    cell_result <- assessDonorInfluence(from_cells, "independent", contrast)
    bulk_result <- assessDonorInfluence(from_bulk, "independent", contrast)
    expect_equal(SummarizedExperiment::assay(cell_result[["results"]][["T"]],
        "logFC"), SummarizedExperiment::assay(bulk_result[["results"]][["T"]],
        "logFC"))
})

test_that("registry conflicts and repeated biological samples are rejected", {
    pb <- make_fixture()
    registry <- S4Vectors::metadata(pb)$scdonoraudit$sample_table
    source <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = SummarizedExperiment::assay(pb, "counts")),
        colData = SummarizedExperiment::colData(pb))
    bad_registry <- registry
    bad_registry$donor_id[1L] <- "wrong_donor"
    expect_error(preparePseudobulk(source, "sample_id", "donor_id",
        "condition", "cell_type", sample_table = bad_registry,
        n_cells = "n_cells"), "METADATA_CONFLICT")
    repeated <- registry
    repeated$donor_id[2L] <- repeated$donor_id[1L]
    repeated$condition[2L] <- repeated$condition[1L]
    SummarizedExperiment::colData(source)$donor_id[2L] <-
        repeated$donor_id[2L]
    SummarizedExperiment::colData(source)$condition[2L] <-
        repeated$condition[2L]
    expect_error(preparePseudobulk(source, "sample_id", "donor_id",
        "condition", "cell_type", sample_table = repeated,
        n_cells = "n_cells"), "UNSUPPORTED_REPEATED_MEASURES")
})
