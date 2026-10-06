# Output dictionary

Preparation metadata uses schema `0.2`; fitted results use schema `0.3`. Tables are
S4Vectors `DataFrame`s unless stated otherwise. Missing values are `NA`, never
zero-filled. IDs are strings, counts are integers unless noted below.

## Prepared pseudobulk

`preparePseudobulk()` returns a `SummarizedExperiment`. Its numeric `counts`
assay contains genes by observed sample/cell-type combinations. Sparse numeric
input is not densified at cell level. `rowData` is retained in gene-ID order.
`colData` contains character `pb_id`, `sample_id`, `donor_id`, `condition`,
`cell_type`; integer `n_cells` (unknown is NA); and retained `sample_vars` in
their original types. A `pb_id` identifies one observed combination.

`metadata(pb)$scdonoraudit` contains `schema_version` (character),
`sample_table` (base data.frame registry), `coverage` (base data.frame), and
`preparation` (list). Coverage has character `sample_id`, `donor_id`,
`condition`, `cell_type`; logical `observed`; integer `n_cells`. It spans the
registered samples and observed cell-type labels; a completely unobserved
cell-type label is not invented. Missing combinations have no count column.
Preparation has character `input_type` (`single_cell` or `pseudobulk`),
`source_assay`, `audit_universe` (`registered_samples` or `observed_only`),
and character vector `sample_vars`.

## Audit and fit tables

`auditDesign()` returns a `SimpleList` with `spec`, `sample_membership`, `runs`,
`issues`, `config`, `provenance`. `assessDonorInfluence()` returns a `SimpleList`
with `results`, `fits`, `fit_samples`, `spec`, `sample_membership`, `gene_filter`,
`coverage`, `issues`, `config`, `provenance`. Audit `issues` is the subset of
design run rows with non-OK reasons; fit `issues` is the event table below.

| Table | Fields and types |
|:--|:--|
| `spec` | character `cell_type`, `plan_status`, `baseline_status`; integer `n_planned`, `n_samples`, `n_donors` |
| `sample_membership` | character `cell_type`, `sample_id`, `donor_id`, `condition`, `reason_code`; logical `observed`, `in_baseline`; integer `n_cells` |
| `coverage` | Preparation coverage fields above, returned as DataFrame in fitted results |
| `runs` / `fits` | character `cell_type`, `run_id`, `run_type`, `omitted_donor`, `execution_status`, `reason_code`, `design_reason_code`, `stage`; integer `n_samples`, `n_donors`, `design_rank`, `residual_df`; logical `target_estimable`, `low_replication` |
| additional `fits` fields | character `library_sizes`, `norm_factors`, `design_columns`, comma-delimited and for completed fits only; prefer `fit_samples` for numeric use |
| `fit_samples` | character `cell_type`, `run_id`, `sample_id`; numeric `library_size`, `norm_factor`; one row per included sample per completed fit, in fit order |
| `gene_filter` | character `cell_type`, `gene_id`, `filter_status`; logical `kept`; includes every input gene in every cell type |
| fit `issues` | character `cell_type`, `run_id`, `stage`, `reason_code`, `message`; zero or more events per run |

`n_planned` counts omissions, excluding the baseline. An incomplete strict
pair yields `plan_status=not_created`, `n_planned=NA` and only a blocked baseline
row. Otherwise a plan is `created`, even if fitting is blocked later.
`baseline_status` in `spec` is the design audit status, not the eventual backend
status; use the baseline row in `fits` to inspect execution. `run_type` is
`baseline` or `omission`. `omitted_donor=NA` denotes the baseline, whose
`run_id` is `baseline`; deletions use `omit_001`, etc. within each cell type.

`execution_status` is `ready` (audit only), `completed`, `skipped` or `failed`.
`stage` is `design`, `filter`, `baseline` or `backend`. A structurally invalid
design is `skipped`, not a backend failure. `design_reason_code` remains the
original design diagnosis when `reason_code` subsequently changes.
`filter_status` is `not_run`, `filtered` or `failed`. Low replication is an
advisory flag, not a structural veto. A backend warning can coexist with a
completed fit or a subsequent error; consult every event in `issues`.

## Per-gene results

`results` is a named `SimpleList` of per-cell-type `SummarizedExperiment`s.
Rows are the fixed baseline-retained gene IDs. Columns include all planned
runs, including unavailable ones. `colData` has character `run_id`, `run_type`
and `omitted_donor`. `metadata` contains character `cell_type` and named
integer vector `status_codes`.

| Assay | Type | Meaning |
|:--|:--|:--|
| `logFC` | numeric matrix | log2 numerator/denominator effect |
| `p_value` | numeric matrix | edgeR test p value |
| `padj_within_cell_type` | numeric matrix | BH adjustment within one cell type and run; fixed baseline family size |
| `effect_valid` | logical matrix | effect is finite |
| `test_valid` | logical matrix | effect is valid and p value finite within [0,1] |
| `gene_status_code` | integer matrix | diagnostic code below |

| Code | Name | Meaning |
|--:|:--|:--|
| 0 | `OK` | finite effect and valid test without another flag |
| 1 | `ALL_ZERO_AFTER_DELETION` | gene is zero across remaining samples; no effect/test |
| 2 | `NONFINITE_EFFECT` | nonfinite effect; no comparable effect/test |
| 3 | `INVALID_PVALUE` | finite effect but unavailable/invalid test |
| 4 | `ONE_CONDITION_ZERO` | counts zero in one condition; a finite effect/test may still be valid |
| 5 | `RUN_UNAVAILABLE` | whole run unavailable; effects/tests are NA and validity flags FALSE |

Validity flags are authoritative for comparison. Code 4 does not automatically
invalidate a finite effect. BH does not control multiplicity across cell types
or repeated deletions.

## The 16 summary columns

`summarizeInfluence()` returns zero or more rows with exactly these columns
and types, including when no genes survive or all cell types are blocked.

| Field | Type | Meaning |
|:--|:--|:--|
| `cell_type` | character | cell-type ID |
| `gene_id` | character | retained stable gene ID |
| `baseline_logFC` | numeric | baseline effect |
| `baseline_p_value` | numeric | baseline unadjusted p value |
| `n_planned` | integer | planned omissions, including skipped/failed ones |
| `n_effect_valid` | integer | omissions with finite comparable effect AND valid baseline effect |
| `n_test_valid` | integer | omissions with valid test AND valid baseline test |
| `max_abs_delta_observed` | numeric | largest observed absolute effect change; NA if no valid comparison |
| `min_logFC_observed` | numeric | smallest valid deletion effect; observed range, not confidence bound |
| `max_logFC_observed` | numeric | largest valid deletion effect |
| `n_sign_reversal` | integer | comparable effects with opposite signs; zero does not reverse sign; NA if none are comparable |
| `n_material_reversal` | integer | opposite effects each strictly exceeding `effect_threshold` in magnitude; NA if threshold absent or no comparable effects |
| `baseline_near_zero` | logical | abs(baseline effect) at or below threshold; NA if threshold absent/unavailable |
| `complete_effect_coverage` | logical | valid baseline, positive planned count and comparable effect count equals planned count |
| `complete_test_coverage` | logical | valid baseline test, positive planned count and comparable test count equals planned count |
| `max_influence_donors` | IRanges CharacterList | donors tied within 1e-8 * max(1, largest change); empty if unavailable |

## Reasons and input errors

Membership reasons: `OK`, `MISSING_CELL_TYPE`, `BELOW_MIN_CELLS`,
`STRICT_PAIR_BLOCK`, `PAIR_MEMBER_EXCLUDED`.

Design/run reasons: `OK`, `NO_ELIGIBLE_SAMPLES`, `MISSING_CONDITION`,
`INCOMPLETE_PAIR`, `DESIGN_MATRIX_ERROR`, `DESIGN_RANK_DEFICIENT`,
`NO_RESIDUAL_DF`, `CONTRAST_NOT_ESTIMABLE`, `ZERO_LIBRARY`,
`NO_GENES_AFTER_FILTER`, `FILTER_ERROR`, `BASELINE_UNAVAILABLE`,
`BACKEND_ERROR`. `BACKEND_WARNING` is an event reason, not a failed-fit status.
Rank deficiency can coexist with an estimable target; fitting is still skipped.

Input errors stop with `[CODE] message`. These are not fit rows:
`INVALID_INPUT`, `INVALID_ARGUMENT`, `INVALID_COUNTS`, `UNSUPPORTED_COUNTS`,
`MISSING_ASSAY`, `EMPTY_INPUT`, `MISSING_GENE_ID`, `DUPLICATE_GENE_ID`,
`INVALID_SAMPLE_TABLE`, `SAMPLE_NOT_REGISTERED`, `METADATA_CONFLICT`,
`MISSING_ID`, `MISSING_FIELD`, `INVALID_METADATA`,
`UNSUPPORTED_REPEATED_MEASURES`,
`DUPLICATE_PSEUDOBULK`, `INVALID_CELL_COUNT`, `INVALID_PSEUDOBULK`,
`STALE_PSEUDOBULK`, `INVALID_DESIGN`, `INVALID_CONTRAST`,
`UNSUPPORTED_CONDITIONS`, `INCOMPLETE_REGISTRY_PAIR`,
`INVALID_INDEPENDENT_DESIGN`, `MISSING_COVARIATE`, `INVALID_COVARIATE`,
`CELL_COUNT_UNKNOWN`, `INVALID_FILTER_ARGS`, `INVALID_BACKEND_ARGS`,
`INVALID_RESULT`, `UNKNOWN_CELL_TYPE`, `UNKNOWN_GENE_ID`, `NO_DELETIONS`.
`ALL_ZERO_RUN` is captured as a backend
error if no active genes remain in a run. See individual help pages for the
conditions under which an argument is accepted.

## Configuration, provenance and plotting

`config` is a list: character `design`, named character `contrast`, optional
numeric `min_cells`, character `covariates`, character `pair_policy`, numeric
`min_replicates_warn`; fitted results also include a named list `filter`
with numeric values,
list `backend` with logical `robust`, character `adjustment` and
`schema_version`. Factor levels and numeric covariate scaling are frozen at
baseline; paired donor levels are adjusted to the remaining donors.

Audit `provenance` records character `audit_universe`. Fit `provenance`
records input preparation, character `scDonorAudit_version`, `edgeR_version`,
`R_version`, and the named character `contrast`.
`plotInfluence()` returns a ggplot object. Coverage displays membership;
effect and influence displays require explicit retained gene IDs and leave
unavailable comparisons visible.
