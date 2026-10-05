# Development provenance

OpenAI Codex assisted in drafting and revising the R implementation, tests,
documentation, and the script that prepares the bundled Crowell19 Astrocyte
pseudobulk example. This disclosure covers those components as a whole.
For version 0.99.5, OpenAI Codex (Sol) assisted with requirements, internal
refactoring, documentation, validation scripts and independent review.
DeepSeek V4.1 Flash assisted with the logical/pattern-count validation,
empty-summary implementation and their regression tests. Its bounded task
timed out before validation; Codex inspected the actual changes and ran the
full regression suite independently. The candidate's checks and numerical
comparisons are documented separately; AI assistance does not establish
scientific validity or transfer authorship or maintenance responsibility.
Fuhao Jiang is the package author and maintainer responsible for reviewing
and maintaining the code and its scientific interpretation.

The original Crowell19 dataset and its CC BY 4.0 license are attributed
separately in `extdata/README-crowell19.md`.
