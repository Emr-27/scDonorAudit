boundary_fixture <- function(paired = FALSE) {
    set.seed(1103)
    donor <- if (paired) rep(paste0("d", 1:6), each = 2) else paste0("d", 1:12)
    donor[donor == "d1"] <- "baseline"
    condition <- if (paired) rep(c("ctrl", "stim"), 6) else
        rep(c("ctrl", "stim"), each = 6)
    samples <- data.frame(sample_id = paste0("s", seq_along(donor)),
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
