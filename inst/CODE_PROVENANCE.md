# Development provenance

OpenAI Codex assisted with the R implementation, regression tests, documentation,
and scripts for example-data preparation and validation. DeepSeek V4.1
Flash assisted with logical/pattern-count validation, empty-summary
handling, and related regression tests.

Automated validation includes direct edgeR reference calculations on
independent and paired synthetic fixtures. It compares baseline and
donor-deletion log2 fold-changes and p values, and runs package checks
across platforms.

Fuhao Jiang is the author and maintainer responsible for reviewing the
code, its reliability, scientific interpretation, and long-term
maintenance. Package code is distributed under Artistic-2.0.

The Crowell19 example data and its CC BY 4.0 license are documented
separately in [extdata/README-crowell19.md](extdata/README-crowell19.md).
