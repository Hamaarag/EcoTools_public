test_that("mgcv GAMs produce DHARMa residuals when mgcViz is installed", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("mgcViz")
  skip_if_not_installed("DHARMa")

  set.seed(1)

  test_data <- data.frame(
    y = stats::rpois(
      80L,
      lambda = 3
    ),
    x = seq(
      0,
      1,
      length.out = 80L
    )
  )

  result <- fit_and_simulate(
    model_name = "test_mgcv_dharma",
    package = "mgcv",
    model_formula = y ~ s(x, k = 5),
    family = "poisson",
    data = test_data,
    output_dir = tempdir(),
    force_refit = TRUE,
    dharma_n_sim = 50L,
    dharma_seed = 123L
  )

  expect_s3_class(
    result$model,
    "gam"
  )

  expect_false(
    is.null(result$residuals)
  )

  expect_null(
    result$simulation_error
  )
})

test_that("fitted models remain usable when data are local to a function", {
  set.seed(1)
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")

  fit_locally <- function() {
    d_local <- data.frame(
      y = stats::rpois(
        60L,
        lambda = 2
      ),
      x = stats::rnorm(60L),
      id = factor(
        rep(
          seq_len(12L),
          each = 5L
        )
      )
    )

    fit_and_simulate(
      model_name = "test_local_data",
      package = "glmmTMB",
      model_formula = y ~ x + (1 | id),
      family = "poisson",
      data = d_local,
      output_dir = tempdir(),
      force_refit = TRUE,
      num_cores = 1L,
      dharma_n_sim = 50L,
      dharma_seed = 123L
    )
  }

  result <- fit_locally()

  expect_s3_class(
    result$model,
    "glmmTMB"
  )

  expect_false(
    is.null(result$residuals)
  )

  expect_false(
    identical(
      result$model_status,
      "valid_fit_no_residual_diagnostics"
    )
  )
})

test_that("a locally fitted cached model remains usable after serialization", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")

  set.seed(123)

  cache_dir <- tempfile(
    pattern = "ecotools-local-serialization-"
  )

  dir.create(
    cache_dir,
    recursive = TRUE
  )

  withr::defer(
    unlink(
      cache_dir,
      recursive = TRUE,
      force = TRUE
    )
  )

  fit_locally <- function() {
    d_local <- data.frame(
      y = stats::rpois(
        60L,
        lambda = 2
      ),
      x = stats::rnorm(60L),
      id = factor(
        rep(
          seq_len(12L),
          each = 5L
        )
      )
    )

    fit_and_simulate(
      model_name = "test_local_serialized",
      package = "glmmTMB",
      model_formula = y ~ x + (1 | id),
      family = "poisson",
      data = d_local,
      output_dir = cache_dir,
      force_refit = TRUE,
      num_cores = 1L,
      dharma_n_sim = 50L,
      dharma_seed = 123L
    )
  }

  original_result <- fit_locally()

  cache_file <- file.path(
    cache_dir,
    "test_local_serialized.rds"
  )

  expect_true(
    file.exists(cache_file)
  )

  reloaded_result <- readRDS(
    cache_file
  )

  expect_s3_class(
    reloaded_result,
    "ecotools_model_bundle"
  )

  expect_s3_class(
    reloaded_result$model,
    "glmmTMB"
  )

  stored_formula <- stats::formula(
    reloaded_result$model
  )

  stored_formula_environment <- environment(
    stored_formula
  )

  expect_true(
    exists(
      "d_local",
      envir = stored_formula_environment,
      inherits = FALSE
    )
  )

  expect_s3_class(
    get(
      "d_local",
      envir = stored_formula_environment,
      inherits = FALSE
    ),
    "data.frame"
  )

  expect_false(
    is.null(reloaded_result$residuals)
  )

  expect_equal(
    stats::predict(
      original_result$model,
      type = "response"
    ),
    stats::predict(
      reloaded_result$model,
      type = "response"
    )
  )

  reloaded_residuals <- DHARMa::simulateResiduals(
    fittedModel = reloaded_result$model,
    n = 50L,
    seed = 123L
  )

  expect_s3_class(
    reloaded_residuals,
    "DHARMa"
  )
})

test_that("failed fits retain fitting-error metadata", {
  skip_if_not_installed("glmmTMB")

  set.seed(123)

  bad_data <- data.frame(
    y = stats::rpois(
      20L,
      lambda = 2
    ),
    x = factor(
      rep(
        "only_level",
        20L
      )
    )
  )

  cache_dir <- tempfile(
    pattern = "ecotools-failed-fit-"
  )

  dir.create(
    cache_dir,
    recursive = TRUE
  )

  withr::defer(
    unlink(
      cache_dir,
      recursive = TRUE,
      force = TRUE
    )
  )

  result <- NULL

  expect_warning(
    result <- fit_and_simulate(
      model_name = "test_failed_fit",
      package = "glmmTMB",
      model_formula = y ~ x,
      family = "poisson",
      data = bad_data,
      output_dir = cache_dir,
      force_refit = TRUE,
      design_check = FALSE,
      sparse_screen = FALSE,
      num_cores = 1L,
      dharma_n_sim = 50L,
      dharma_seed = 123L
    ),
    regexp = "FATAL ERROR fitting"
  )

  expect_s3_class(
    result,
    "ecotools_model_bundle"
  )

  expect_null(
    result$model
  )

  expect_identical(
    result$model_status,
    "invalid"
  )

  expect_false(
    result$converged
  )

  expect_false(
    result$numerically_converged
  )

  expect_false(
    result$valid_for_inference
  )

  expect_false(
    result$valid_for_model_selection
  )

  expect_type(
    result$fit_error,
    "list"
  )

  expect_true(
    is.character(
      result$fit_error$message
    )
  )

  expect_length(
    result$fit_error$message,
    1L
  )

  expect_true(
    nzchar(
      result$fit_error$message
    )
  )

  expect_null(
    result$simulation_error
  )

  expect_true(
    file.exists(
      file.path(
        cache_dir,
        "test_failed_fit.rds"
      )
    )
  )
})

test_that("glmmTMB stored call uses a compact named data reference", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")

  set.seed(123)

  test_data <- data.frame(
    response = rpois(80, lambda = 3),
    treatment = factor(
      rep(c("A", "B"), each = 40)
    ),
    x = rnorm(80)
  )

  bundle <- fit_and_simulate(
    model_name = "test_compact_named_data_call",
    package = "glmmTMB",
    model_formula = response ~ treatment + x,
    family = stats::poisson(),
    data = test_data,
    output_dir = tempdir(),
    force_refit = TRUE,
    sparse_screen = FALSE,
    dharma_n_sim = 50L,
    dharma_seed = 123L
  )

  expect_s3_class(
    bundle$model,
    "glmmTMB"
  )

  expect_true(
    is.symbol(bundle$model$call$data)
  )

  expect_identical(
    as.character(bundle$model$call$data),
    "test_data"
  )

  expect_identical(
    bundle$data_call_name,
    "test_data"
  )

  summary_text <- paste(
    capture.output(
      print(summary(bundle$model))
    ),
    collapse = "\n"
  )

  expect_match(
    summary_text,
    "Data:\\s+test_data"
  )

  expect_false(
    grepl(
      "Data:\\s+structure\\(",
      summary_text
    )
  )
})

test_that("glmmTMB model with compact data call survives serialization", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")
  set.seed(123)

  test_data <- data.frame(
    response = rpois(80, lambda = 3),
    treatment = factor(
      rep(c("A", "B"), each = 40)
    ),
    x = rnorm(80)
  )

  bundle <- fit_and_simulate(
    model_name = "test_serialized_compact_data_call",
    package = "glmmTMB",
    model_formula = response ~ treatment + x,
    family = stats::poisson(),
    data = test_data,
    output_dir = tempdir(),
    force_refit = TRUE,
    sparse_screen = FALSE,
    dharma_n_sim = 50L,
    dharma_seed = 123L
  )

  model_file <- tempfile(
    fileext = ".rds"
  )

  saveRDS(
    bundle$model,
    model_file
  )

  restored_model <- readRDS(
    model_file
  )

  expect_true(
    is.symbol(restored_model$call$data)
  )

  expect_identical(
    as.character(restored_model$call$data),
    "test_data"
  )

  expect_silent(
    summary(restored_model)
  )

  restored_predictions <- stats::predict(
    restored_model,
    type = "response"
  )

  expect_length(
    restored_predictions,
    nrow(test_data)
  )

  expect_true(
    all(is.finite(restored_predictions))
  )
})

test_that("DHARMa_diag supports a directly fitted glmmTMB model", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")

  set.seed(123)

  test_data <- data.frame(
    y = stats::rpois(60, lambda = 2),
    x = stats::rnorm(60)
  )

  fitted_model <- glmmTMB::glmmTMB(
    y ~ x,
    family = stats::poisson(),
    data = test_data
  )

  result <- DHARMa_diag(
    res_object = fitted_model,
    data = test_data,
    print_plots = FALSE,
    export_plots = FALSE,
    n_sim = 100L,
    seed = 123L
  )

  expect_type(result, "list")

  expect_true(
    all(
      c(
        "model_title",
        "dispersion_test",
        "dispersion_ratio",
        "dispersion_p_value",
        "zero_inflation_test",
        "zero_inflation_ratio",
        "zero_inflation_p_value",
        "response_has_zero",
        "predictors_tested",
        "offset_variables_tested",
        "offset_expressions"
      ) %in%
        names(result)
    )
  )

  expect_identical(
    result$model_title,
    "DHARMa Diagnostics"
  )
})
