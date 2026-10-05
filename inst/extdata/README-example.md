# Synthetic help-example data

`example_counts.csv` contains 250 genes by eight independent samples.
`example_samples.csv` assigns one donor per sample, with four `ctrl` and four
`stim` samples. Seed 1103 generates negative-binomial counts with mean 25 and
size 5. The first ten genes' generated `stim` counts are multiplied by two.
These artificial data demonstrate API use; they are not biological evidence
or a calibrated statistical performance benchmark. No human data are present.

Rebuild byte-equivalent CSVs with the installed script:

```text
Rscript create_synthetic_example.R <output-directory>
```

The script is at `system.file("scripts", "create_synthetic_example.R",
package = "scDonorAudit")`. The independently generated vignette examples
state their own seeds and cohort construction. Source files are distributed
under the package's Artistic-2.0 license; development provenance is recorded
in `CODE_PROVENANCE.md`.
