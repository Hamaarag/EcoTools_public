test_that("DHARMa_diag titles and exports model-specific plots", {
  skip_if_not_installed("DHARMa")
  skip_if_not_installed("glmmTMB")

  set.seed(123)

  test_data <- data.frame(
    y = stats::rpois(
      80L,
      lambda = 3
    ),
    x = stats::runif(
      80L
    )
  )

  fitted_model <- glmmTMB::glmmTMB(
    y ~ x,
    family = stats::poisson(),
    data = test_data
  )

  residual_object <- DHARMa::simulateResiduals(
    fittedModel = fitted_model,
    n = 100L,
    seed = 123L,
    plot = FALSE
  )

  output_directory <- tempfile(
    pattern = "dharma-title-test-"
  )

  dir.create(
    output_directory,
    recursive = TRUE
  )

  on.exit(
    unlink(
      output_directory,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  diagnostic_result <- DHARMa_diag(
    res_object = list(
      model = fitted_model,
      residuals = residual_object,
      model_status = "valid",
      analysis_data = test_data
    ),
    model_title = "title_test_model",
    print_plots = FALSE,
    export_plots = TRUE,
    output_dir = output_directory
  )

  expected_files <- c(
    "title_test_model_DHARMa_overview.png",
    "title_test_model_residuals_vs_x.png",
    "title_test_model_dispersion_test.png"
  )

  if (isTRUE(diagnostic_result$response_has_zero)) {
    expected_files <- c(
      expected_files,
      "title_test_model_zero_inflation_test.png"
    )
  }

  expect_true(
    all(
      file.exists(
        file.path(
          output_directory,
          expected_files
        )
      )
    )
  )

  expect_identical(
    diagnostic_result$model_title,
    "title_test_model"
  )
})
