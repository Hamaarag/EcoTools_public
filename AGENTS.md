# EcoTools Agent Instructions

## Package Workflow

- This is an R package (R >= 4.1) using testthat edition 3. Implement functions in `R/`; keep roxygen documentation next to the owning function. Regenerate `NAMESPACE` and `man/` from roxygen rather than editing generated files directly.
- Run tests from the package root with `testthat::test_local()`. For a focused test, use `testthat::test_local(filter = "<test-file-name-fragment>")`, for example `testthat::test_local(filter = "gllvm-policy")`.
- When changing exported functions or their documentation, run `roxygen2::roxygenise()` and then the relevant tests. Skip tests that require unavailable optional dependencies only through the established testthat skip conventions.
- Add regression tests under `tests/testthat/`, close to the behavior being changed. Preserve coverage for model serialization, cache stability, validity statuses, and restoration of process-global settings when touching those areas.

## Ownership Map

- `R/fit_and_simulate.R`: GLMM/GAM fitting, diagnostics, and model bundles.
- `R/fit_gllvm_cached.R`: GLLVM fitting, caching, response screening, design checks, and reference-level stability.
- `R/get_analysis_fingerprint.R`, `R/get_gllvm_fingerprint.R`, and `R/canonicalize.R`: cache identity and deterministic serialization.
- Other focused behavior is generally owned by the correspondingly named file in `R/`; use that file and its neighboring test as the starting point.

## Modeling Contracts

- A converged fit is not automatically suitable for inference. Preserve standardized `model_status` behavior and use `assert_valid_model()` as the inference gate; do not bypass it in downstream inference workflows.
- If a fitting or validation option can change the resulting model, its interpretation, or its eligibility, check that the corresponding cache fingerprint includes it. Add or update a fingerprint regression test when changing this contract.
- Keep cache canonicalization deterministic across sessions. Extra arguments that cannot be represented stably (including functions and environments) are intentionally rejected; do not weaken that policy without tests.
- `with_gllvm_tmb_threads()` manages process-level TMB/OpenMP settings. Ensure settings are restored on both normal return and error paths; requested threads do not imply parallel R workers.
- Treat sparse-response checks and simulation-based diagnostics as explicit validation signals, not universal guarantees of model quality.

## Further Reading

- [README](README.md) covers installation and the package function reference.
- [Workflow vignette](vignettes/ecotools-workflow.Rmd) describes the intended fitting, validation, and inference sequence.
- [tests/testthat](tests/testthat/) contains focused examples of package contracts and regression-test patterns.
