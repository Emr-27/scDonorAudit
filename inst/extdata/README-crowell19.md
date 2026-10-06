# Crowell19 Astrocyte pseudobulk example

`crowell19_astrocyte_counts.csv` contains summed raw counts for all 11,076 genes across eight mouse samples (four LPS, four Vehicle). `crowell19_astrocyte_samples.csv` records the 1,740 Astrocyte nuclei contributing to the eight sample columns. These files provide a compact, real-data example for preparing pseudobulk counts and interpreting donor-deletion results under a specified model.

Source: ExperimentHub resource `EH3297` (`muscData::Crowell19_4vs4()`), derived from Mark Robinson's [Figshare dataset](https://doi.org/10.6084/m9.figshare.8976473.v1), licensed **CC BY 4.0**. Cite that dataset and Crowell et al., *Nature Communications* 11, 6077 (2020), DOI [10.1038/s41467-020-19894-4](https://doi.org/10.1038/s41467-020-19894-4), when reusing the data. The source resource had SHA-256 `3193A33FA7ACB650FDB1E682822705D7C68057ABA93D84DAC63A364E40FDB1FE` when downloaded on 2026-09-29.

The bundled preparation script `inst/scripts/create_crowell19_pseudobulk.R` loads a saved `EH3297` Rda containing one `sce` object, selects `cluster_id == "Astrocytes"`, sums the raw `counts` assay by `sample_id`, uses `rowData$ENSEMBL` as gene IDs, and sorts genes and samples by ID. It does not filter or rank genes using treatment labels or model results. To reproduce the files, retrieve `EH3297` with `ExperimentHub`, save it as `sce` in an Rda file, and run the script with the Rda path and an output directory. The raw single-nucleus object is not bundled with the software package.

Derived file SHA-256 values:

- `crowell19_astrocyte_counts.csv`: `6AAC213B8DA8118DDACAAD7B023DB9761367551DDFCF2BE12A87806FBD157234`
- `crowell19_astrocyte_samples.csv`: `6E6520E600F1F0348163F13B6494F0FE128E811A37DC4CA66BE620AA17FAA159`
