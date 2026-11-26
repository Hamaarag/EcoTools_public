test_that("reference-level stability diagnostics are stored when factor predictors exist", {
  skip_if_not_installed("gllvm")

  set.seed(123)

  X <- data.frame(
    Treatment = factor(
      rep(
        c(
          "Far",
          "Near",
          "Agriculture"
        ),
        each = 20
      )
    )
  )

  y <- matrix(
    stats::rpois(
      60L,
      lambda = c(
        rep(1, 20),
        rep(5, 20),
        rep(10, 20)
      )
    ),
    ncol = 1L
  )

  colnames(y) <- "Guild"
  y[y == 0] <- 1 # set any zeros to one for this test

  fit <- fit_gllvm_cached(
    model_name = "test_reference_stability",
    y = y,
    X = X,
    studyDesign = X,
    formula = ~Treatment,
    family = "poisson",
    num.lv = 0,
    sd.errors = FALSE,
    force_refit = TRUE,
    reference_stability_check = "always",
    cache_dir = tempdir(),
    sparse_screen = FALSE,
    plot_diag = FALSE
  )

  reference_metadata <- fit$reference_stability
  reference_result <- reference_metadata$result

  expect_type(
    reference_metadata,
    "list"
  )

  expect_identical(
    reference_metadata$requested_policy,
    "always"
  )

  expect_true(
    reference_metadata$candidate_eligible
  )

  expect_true(
    reference_metadata$check_run
  )

  expect_null(
    reference_metadata$skip_reason
  )

  expect_identical(
    reference_metadata$detected_factor_predictors,
    "Treatment"
  )

  expect_type(
    reference_result,
    "list"
  )

  expect_true(
    reference_result$checked
  )

  expect_true(
    "Treatment" %in%
      reference_result$factor_predictors
  )

  expect_true(
    "Treatment" %in%
      names(
        reference_result$factor_results
      )
  )

  expect_setequal(
    names(
      reference_result$factor_results$Treatment
    ),
    levels(X$Treatment)
  )

  expect_identical(
    fit$diagnostics$reference_stability,
    fit$reference_stability
  )

  expect_identical(
    fit$reference_sensitive,
    reference_metadata$reference_sensitive
  )
})

test_that("reference-level stability records when no factor predictors exist", {
  skip_if_not_installed("gllvm")

  set.seed(123)

  X <- data.frame(
    temperature = stats::rnorm(40L)
  )

  y <- matrix(
    stats::rpois(
      40L,
      lambda = 3
    ),
    ncol = 1L
  )

  colnames(y) <- "Guild"

  fit <- fit_gllvm_cached(
    model_name = "test_reference_stability_no_factor",
    y = y,
    X = X,
    studyDesign = X,
    formula = ~temperature,
    family = "poisson",
    num.lv = 0,
    sd.errors = FALSE,
    force_refit = TRUE,
    reference_stability_check = "always",
    cache_dir = tempdir(),
    sparse_screen = FALSE,
    plot_diag = FALSE
  )

  reference_metadata <- fit$reference_stability

  expect_type(
    reference_metadata,
    "list"
  )

  expect_identical(
    reference_metadata$requested_policy,
    "always"
  )

  expect_true(
    reference_metadata$candidate_eligible
  )

  expect_false(
    reference_metadata$check_run
  )

  expect_identical(
    reference_metadata$skip_reason,
    "no_factor_predictors"
  )

  expect_identical(
    reference_metadata$detected_factor_predictors,
    character(0)
  )

  expect_null(
    reference_metadata$result
  )

  expect_false(
    reference_metadata$reference_sensitive
  )

  expect_false(
    fit$reference_sensitive
  )

  expect_identical(
    fit$diagnostics$reference_stability,
    fit$reference_stability
  )
})

test_that("character predictors are detected as factor candidates", {
  X <- data.frame(
    Treatment = c(
      "Far",
      "Near",
      "Agriculture"
    ),
    stringsAsFactors = FALSE
  )

  detected <- detect_formula_factors(
    X = X,
    formula = ~Treatment
  )

  expect_identical(
    detected,
    "Treatment"
  )
})

test_that("high-cardinality factors are detected before level-limit filtering", {
  X <- data.frame(
    Plot = factor(
      seq_len(20)
    )
  )

  detected <- detect_formula_factors(
    X = X,
    formula = ~Plot
  )

  expect_identical(
    detected,
    "Plot"
  )
})

test_that("high-cardinality factor filtering works", {
  X <- data.frame(
    Plot = factor(
      seq_len(20)
    )
  )

  factor_names <- "Plot"

  all_factor_names <- factor_names

  filtered_factor_names <- factor_names[
    vapply(
      factor_names,
      function(factor_name) {
        length(
          levels(
            factor(
              X[[factor_name]]
            )
          )
        ) <=
          5
      },
      logical(1)
    )
  ]

  expect_identical(
    filtered_factor_names,
    character(0)
  )
})

test_that("detect_formula_factors detects a factor predictor", {
  X <- data.frame(
    Treatment = factor(
      c("A", "B", "C", "A")
    ),
    continuous = c(1, 2, 3, 4)
  )

  result <- EcoTools:::detect_formula_factors(
    X = X,
    formula = ~Treatment
  )

  expect_identical(
    result,
    "Treatment"
  )
})

test_that("detect_formula_factors ignores unused columns", {
  X <- data.frame(
    Treatment = factor(
      c("A", "B", "C", "A")
    ),
    Site = factor(
      c("one", "one", "two", "two")
    )
  )

  result <- EcoTools:::detect_formula_factors(
    X = X,
    formula = ~Treatment
  )

  expect_identical(
    result,
    "Treatment"
  )

  expect_false(
    "Site" %in% result
  )
})

test_that("reference report distinguishes skipped gradient check from no factors", {
  test_model <- list(
    model_name = "test_reference_skip",
    model_status = "suspicious_optimization_incomplete",
    opt_converged = TRUE,
    gradient_assessed = TRUE,
    gradient_valid = FALSE,
    max_abs_gradient = 0.02,
    gradient_threshold = 0.01,
    sd_errors_requested = TRUE,
    se_valid = TRUE,
    sd = list(
      Xcoef = matrix(
        0.1,
        nrow = 1L,
        ncol = 1L
      )
    ),
    diagnostics = list(
      inference_quality = list(
        se_requested = TRUE,
        se_valid = TRUE
      ),
      coef_screen = list(
        n_extreme = 0L,
        coef_threshold = 5
      )
    ),
    reference_stability = list(
      requested_policy = "candidate",
      candidate_eligible = FALSE,
      check_run = FALSE,
      skip_reason = "gradient_failed",
      detected_factor_predictors = "Treatment",
      result = NULL,
      reference_sensitive = FALSE
    ),
    reference_sensitive = FALSE,
    sparse_screen = NULL,
    sparse_species = character(0),
    warning_msgs = NULL
  )

  captured_warning <- NULL

  withCallingHandlers(
    print_convergence_report(
      test_model,
      "test_reference_skip"
    ),
    warning = function(w) {
      captured_warning <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  expect_match(
    captured_warning,
    "SKIPPED: gradient check failed"
  )

  expect_match(
    captured_warning,
    "Treatment"
  )

  expect_false(
    grepl(
      "NO FACTOR PREDICTORS DETECTED",
      captured_warning,
      fixed = TRUE
    )
  )
})

test_that("convergence report distinguishes disabled sparse screening from failure", {
  test_model <- list(
    model_name = "test_sparse_not_run",
    model_status = "valid_fit_no_standard_errors",
    opt_converged = TRUE,
    gradient_assessed = TRUE,
    gradient_valid = TRUE,
    max_abs_gradient = 0.001,
    gradient_threshold = 0.01,
    sd_errors_requested = FALSE,
    se_valid = NA,
    sd = NULL,
    diagnostics = list(
      inference_quality = list(
        se_requested = FALSE,
        se_valid = NA
      ),
      coef_screen = list(
        n_extreme = 0L,
        coef_threshold = 5
      )
    ),
    reference_stability = list(
      requested_policy = "always",
      candidate_eligible = TRUE,
      check_run = FALSE,
      skip_reason = "no_factor_predictors",
      detected_factor_predictors = character(0),
      result = NULL,
      reference_sensitive = FALSE
    ),
    reference_sensitive = FALSE,
    sparse_screen = NULL,
    sparse_species = character(0),
    warning_msgs = NULL
  )

  captured_message <- paste(
    capture.output(
      print_convergence_report(
        test_model,
        "test_sparse_not_run"
      ),
      type = "message"
    ),
    collapse = "\n"
  )

  expect_match(
    captured_message,
    "Sparse Responses\\s+: NOT RUN / NOT AVAILABLE"
  )

  expect_false(
    grepl(
      "SCREEN FAILED / NOT ASSESSED",
      captured_message,
      fixed = TRUE
    )
  )
})

test_that("reference stability metadata distinguishes detection from execution", {
  metadata <- list(
    requested_policy = "candidate",
    candidate_eligible = FALSE,
    check_run = FALSE,
    skip_reason = "gradient_failed",
    detected_factor_predictors = "Treatment",
    result = NULL,
    reference_sensitive = FALSE
  )

  expect_identical(
    metadata$detected_factor_predictors,
    "Treatment"
  )

  expect_false(
    metadata$check_run
  )

  expect_identical(
    metadata$skip_reason,
    "gradient_failed"
  )

  expect_null(
    metadata$result
  )
})
