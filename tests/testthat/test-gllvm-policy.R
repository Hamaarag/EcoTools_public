test_that("GLLVM design policy changes fingerprint", {
  fp_error <- get_gllvm_fingerprint(
    y = matrix(1:6, ncol = 2),
    X = data.frame(x = 1:3),
    studyDesign = data.frame(site = 1:3),
    formula = ~x,
    family = "poisson",
    num.lv = 1,
    sd.errors = TRUE,
    method = "VA",
    control = list(),
    control.start = list(),
    gradient_threshold = 0.01,
    coef_threshold = 5,
    sparse_screen = TRUE,
    sparse_action = "proceed",
    abort_on = "all",
    sparse_min_events_per_parameter = 5,
    sparse_min_total_events = 5,
    sparse_min_nonzero_rows = 3,
    sparse_event_basis = "nonzero",
    sparse_count_like_only = TRUE,
    design_check = TRUE,
    design_action = "error",
    design_tolerance = 1e-8,
    dependency_versions = list(gllvm = "x", TMB = "x")
  )

  fp_warn <- get_gllvm_fingerprint(
    y = matrix(1:6, ncol = 2),
    X = data.frame(x = 1:3),
    studyDesign = data.frame(site = 1:3),
    formula = ~x,
    family = "poisson",
    num.lv = 1,
    sd.errors = TRUE,
    method = "VA",
    control = list(),
    control.start = list(),
    gradient_threshold = 0.01,
    coef_threshold = 5,
    sparse_screen = TRUE,
    sparse_action = "proceed",
    abort_on = "all",
    sparse_min_events_per_parameter = 5,
    sparse_min_total_events = 5,
    sparse_min_nonzero_rows = 3,
    sparse_event_basis = "nonzero",
    sparse_count_like_only = TRUE,
    design_check = TRUE,
    design_action = "warn",
    design_tolerance = 1e-8,
    dependency_versions = list(gllvm = "x", TMB = "x")
  )

  expect_false(identical(fp_error, fp_warn))
})

test_that("mixed GLLVM response families change the fingerprint", {
  fingerprint_for_family <- function(family) {
    get_gllvm_fingerprint(
      y = matrix(1:8, ncol = 2),
      X = data.frame(x = 1:4),
      studyDesign = data.frame(site = 1:4),
      formula = ~x,
      family = family,
      num.lv = 1,
      sd.errors = FALSE,
      method = "VA",
      control = list(),
      control.start = list(),
      gradient_threshold = 0.01,
      coef_threshold = 5,
      sparse_screen = FALSE,
      sparse_action = "proceed",
      abort_on = "all",
      sparse_min_events_per_parameter = 5,
      sparse_min_total_events = 5,
      sparse_min_nonzero_rows = 3,
      sparse_event_basis = "nonzero",
      sparse_count_like_only = TRUE,
      design_check = FALSE,
      dependency_versions = list(gllvm = "x", TMB = "x")
    )
  }

  mixed_family <- c("negative.binomial", "betaH")

  expect_type(fingerprint_for_family(mixed_family), "character")
  expect_false(
    identical(
      fingerprint_for_family(mixed_family),
      fingerprint_for_family(rev(mixed_family))
    )
  )
})

test_that("count-like detector rejects non-finite values", {
  expect_false(is_count_like_response(c(0, 1, Inf)))
  expect_false(is_count_like_response(c(0, 1, -Inf)))
  expect_false(is_count_like_response(c(0, 1, NaN)))
})

test_that("binary screening uses the minority outcome count", {
  skip_if_not_installed("glmmTMB")

  test_data <- data.frame(
    y = c(rep(1, 95), rep(0, 5)),
    x = seq_len(100)
  )

  fitted_model <- glmmTMB::glmmTMB(
    y ~ 1,
    family = stats::binomial(),
    data = test_data
  )

  result <- screen_sparse_count_model(
    model = fitted_model,
    data = test_data,
    response = "y",
    estimated_parameters = 2,
    min_events_per_parameter = 3,
    event_basis = "nonzero"
  )

  expect_true(result$details$is_binary_response)
  expect_identical(result$details$minority_class_n, 5L)
  expect_equal(result$details$informative_event_n, 5)
  expect_equal(
    result$details$informative_events_per_parameter,
    2.5
  )
  expect_false(result$pass)
  expect_true(
    "events_per_parameter" %in% result$warning_types
  )
})

test_that("balanced binary data use the smaller class count", {
  skip_if_not_installed("glmmTMB")

  test_data <- data.frame(
    y = rep(c(0, 1), each = 20),
    x = seq_len(40)
  )

  fitted_model <- glmmTMB::glmmTMB(
    y ~ 1,
    family = stats::binomial(),
    data = test_data
  )

  result <- screen_sparse_count_model(
    model = fitted_model,
    data = test_data,
    response = "y",
    estimated_parameters = 2,
    min_events_per_parameter = 5
  )

  expect_identical(result$details$minority_class_n, 20L)
  expect_equal(
    result$details$informative_events_per_parameter,
    10
  )
  expect_true(result$pass)
})

test_that("all-missing GLLVM response cannot pass", {
  y <- data.frame(
    species_a = c(NA, NA, NA),
    species_b = c(0, 1, 2)
  )

  result <- screen_sparse_gllvm(
    y = y,
    X = data.frame(x = 1:3),
    formula = ~x
  )

  expect_true(
    result$details$sparse_flag[
      result$details$response == "species_a"
    ]
  )
})

test_that("GLLVM sparse screening counts fixed-effect columns correctly", {
  response_matrix <- data.frame(
    species_a = c(0, 1, 2, 0, 3, 1)
  )

  predictor_data <- data.frame(
    treatment = factor(c("A", "A", "B", "B", "C", "C")),
    gradient = seq_len(6)
  )

  result <- screen_sparse_gllvm(
    y = response_matrix,
    X = predictor_data,
    formula = species_a ~ treatment + gradient,
    min_events_per_parameter = 0,
    min_total_events = 0,
    min_nonzero_rows = 0
  )

  expected_parameters <- ncol(
    stats::model.matrix(
      ~ treatment + gradient,
      data = predictor_data
    )
  )

  expect_identical(
    result$estimated_parameters,
    as.integer(expected_parameters)
  )
})

test_that("GLLVM fingerprint changes with reference diagnostics settings", {
  common_arguments <- list(
    y = matrix(
      1:8,
      nrow = 4L,
      ncol = 2L
    ),
    X = data.frame(
      treatment = factor(
        c(
          "A",
          "B",
          "A",
          "B"
        )
      )
    ),
    studyDesign = data.frame(
      site = seq_len(4L)
    ),
    formula = ~treatment,
    TR = NULL,
    row.eff = FALSE,
    family = rep(
      "poisson",
      2L
    ),
    num.lv = 1L,
    sd.errors = FALSE,
    method = "VA",
    control = list(),
    control.start = list(),
    gradient_threshold = 0.01,
    coef_threshold = 5,
    sparse_screen = FALSE,
    sparse_action = "proceed",
    abort_on = "all",
    sparse_min_events_per_parameter = 5,
    sparse_min_total_events = 5,
    sparse_min_nonzero_rows = 3,
    sparse_event_basis = "nonzero",
    sparse_count_like_only = TRUE,
    design_check = TRUE,
    design_action = "error",
    design_tolerance = sqrt(
      .Machine$double.eps
    ),
    reference_stability_check = "candidate",
    reference_stability_max_levels = 5L,
    tmb_threads = 1L,
    dependency_versions = list(
      gllvm = "test",
      TMB = "test"
    ),
    extra_args = list()
  )

  default_fingerprint <- do.call(
    EcoTools:::get_gllvm_fingerprint,
    common_arguments
  )

  changed_tolerance_arguments <- common_arguments
  changed_tolerance_arguments$reference_logLik_tolerance <- 1e-5

  changed_seed_arguments <- common_arguments
  changed_seed_arguments$reference_stability_seed <- 2L

  changed_tolerance_fingerprint <- do.call(
    EcoTools:::get_gllvm_fingerprint,
    changed_tolerance_arguments
  )

  changed_seed_fingerprint <- do.call(
    EcoTools:::get_gllvm_fingerprint,
    changed_seed_arguments
  )

  expect_false(
    identical(
      default_fingerprint,
      changed_tolerance_fingerprint
    )
  )

  expect_false(
    identical(
      default_fingerprint,
      changed_seed_fingerprint
    )
  )
})

test_that("GLLVM fingerprint validates reference diagnostics settings", {
  base_arguments <- list(
    y = matrix(
      1:6,
      nrow = 3L
    ),
    X = data.frame(
      x = 1:3
    ),
    studyDesign = data.frame(
      site = 1:3
    ),
    formula = ~x,
    family = rep(
      "poisson",
      2L
    ),
    num.lv = 1L,
    sd.errors = FALSE,
    method = "VA",
    control = list(),
    control.start = list(),
    gradient_threshold = 0.01,
    coef_threshold = 5,
    sparse_screen = FALSE,
    sparse_action = "proceed",
    abort_on = "all",
    sparse_min_events_per_parameter = 5,
    sparse_min_total_events = 5,
    sparse_min_nonzero_rows = 3,
    sparse_event_basis = "nonzero",
    sparse_count_like_only = TRUE,
    dependency_versions = list(
      gllvm = "test",
      TMB = "test"
    )
  )

  invalid_tolerance_arguments <- base_arguments
  invalid_tolerance_arguments$reference_logLik_tolerance <- -1

  expect_error(
    do.call(
      EcoTools:::get_gllvm_fingerprint,
      invalid_tolerance_arguments
    ),
    "non-negative"
  )

  invalid_seed_arguments <- base_arguments
  invalid_seed_arguments$reference_stability_seed <- 1.5

  expect_error(
    do.call(
      EcoTools:::get_gllvm_fingerprint,
      invalid_seed_arguments
    ),
    "finite integer"
  )

  null_seed_arguments <- base_arguments
  null_seed_arguments$reference_stability_seed <- NULL

  expect_type(
    do.call(
      EcoTools:::get_gllvm_fingerprint,
      null_seed_arguments
    ),
    "character"
  )
})

test_that("non-finite SE component names are extracted safely", {
  se_component_summary <- data.table::data.table(
    component = c(
      "Xcoef",
      "phi",
      "sigma"
    ),
    n_values = c(
      10L,
      3L,
      2L
    ),
    n_finite = c(
      10L,
      1L,
      0L
    ),
    n_nonfinite = c(
      0L,
      2L,
      2L
    ),
    n_zero = c(
      0L,
      0L,
      0L
    ),
    min_value = c(
      0.1,
      0.2,
      NA_real_
    ),
    max_value = c(
      1,
      0.5,
      NA_real_
    )
  )

  affected_components <- as.character(
    se_component_summary$component[
      se_component_summary$n_nonfinite > 0L
    ]
  )

  expect_identical(
    affected_components,
    c(
      "phi",
      "sigma"
    )
  )
})
