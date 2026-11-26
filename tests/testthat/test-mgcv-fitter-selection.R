test_that("fit_and_simulate uses mgcv gam by default", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("mgcViz")
  skip_if_not_installed("DHARMa")

  set.seed(123)

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
    model_name = "test_mgcv_gam_selection",
    package = "mgcv",
    model_formula = y ~ s(x, k = 5),
    family = "poisson",
    data = test_data,
    output_dir = tempdir(),
    force_refit = TRUE,
    sparse_screen = FALSE,
    num_cores = 1L,
    dharma_n_sim = 50L,
    dharma_seed = 123L,
    residual_alpha = 1e-8
  )

  expect_s3_class(
    result$model,
    "gam"
  )

  expect_identical(
    result$mgcv_fitter,
    "gam"
  )

  expect_identical(
    result$model$call[[1L]],
    quote(mgcv::gam)
  )
})

test_that("fit_and_simulate can use mgcv bam", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("mgcViz")
  skip_if_not_installed("DHARMa")

  set.seed(456)

  test_data <- data.frame(
    x = seq(
      0,
      1,
      length.out = 400L
    )
  )

  test_data$y <-
    1 +
    sin(
      2 * pi * test_data$x
    ) +
    stats::rnorm(
      nrow(test_data),
      mean = 0,
      sd = 0.15
    )

  result <- fit_and_simulate(
    model_name = "test_mgcv_bam_selection",
    package = "mgcv",
    mgcv_fitter = "bam",
    model_formula = y ~ s(x, k = 5),
    family = stats::gaussian(),
    data = test_data,
    output_dir = tempdir(),
    force_refit = TRUE,
    sparse_screen = FALSE,
    num_cores = 1L,
    discrete = TRUE,
    method = "fREML",
    dharma_n_sim = 50L,
    dharma_seed = 123L,
    residual_alpha = 1e-8
  )

  expect_s3_class(
    result$model,
    "gam"
  )

  expect_identical(
    result$mgcv_fitter,
    "bam"
  )

  expect_identical(
    result$model$call[[1L]],
    quote(mgcv::bam)
  )

  expect_silent(
    stats::predict(
      result$model,
      type = "response"
    )
  )
})

test_that("mgcv fitter changes the analysis fingerprint", {
  set.seed(123)

  test_data <- data.frame(
    y = stats::rpois(
      20L,
      lambda = 2
    ),
    x = seq_len(20L)
  )

  common_arguments <- list(
    data = test_data,
    model_formula = y ~ x,
    zi_formula = ~0,
    disp_formula = ~1,
    family = stats::poisson(),
    weights = NULL,
    offset = NULL,
    package = "mgcv",
    dependency_versions = list(
      mgcv = "test",
      mgcViz = "test",
      DHARMa = "test"
    )
  )

  gam_arguments <- common_arguments
  gam_arguments$mgcv_fitter <- "gam"

  bam_arguments <- common_arguments
  bam_arguments$mgcv_fitter <- "bam"

  gam_fingerprint <- do.call(
    EcoTools:::get_analysis_fingerprint,
    gam_arguments
  )

  bam_fingerprint <- do.call(
    EcoTools:::get_analysis_fingerprint,
    bam_arguments
  )

  expect_false(
    identical(
      gam_fingerprint,
      bam_fingerprint
    )
  )
})

test_that("bam threading arguments cannot conflict with wrapper control", {
  skip_if_not_installed("mgcv")

  test_data <- data.frame(
    y = stats::rpois(
      20L,
      lambda = 2
    ),
    x = seq_len(20L)
  )

  expect_error(
    fit_and_simulate(
      model_name = "test_conflicting_bam_threads",
      package = "mgcv",
      mgcv_fitter = "bam",
      model_formula = y ~ x,
      family = "poisson",
      data = test_data,
      output_dir = tempdir(),
      force_refit = TRUE,
      num_cores = 1L,
      nthreads = 2L
    ),
    "Do not pass"
  )
})

test_that("DHARMa quantile warnings do not escape unclassified", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("mgcViz")
  skip_if_not_installed("DHARMa")

  set.seed(456)

  test_data <- data.frame(
    x = seq(
      0,
      1,
      length.out = 400L
    )
  )

  test_data$y <-
    1 +
    sin(
      2 * pi * test_data$x
    ) +
    stats::rnorm(
      nrow(test_data),
      mean = 0,
      sd = 0.15
    )

  result <- expect_no_warning(
    fit_and_simulate(
      model_name = "test_bam_quantile_warning_capture",
      package = "mgcv",
      mgcv_fitter = "bam",
      model_formula = y ~ s(
        x,
        k = 8
      ),
      family = stats::gaussian(),
      data = test_data,
      output_dir = tempdir(),
      force_refit = TRUE,
      sparse_screen = FALSE,
      num_cores = 1L,
      discrete = TRUE,
      method = "fREML",
      dharma_n_sim = 50L,
      dharma_seed = 456L,
      residual_alpha = 1e-8
    )
  )

  expect_s3_class(
    result$model,
    "gam"
  )

  if (
    !is.null(
      result$residual_diagnostics$predicted_quantile_warning
    )
  ) {
    expect_identical(
      result$model_status,
      "valid_fit_no_residual_diagnostics"
    )
  }
})
