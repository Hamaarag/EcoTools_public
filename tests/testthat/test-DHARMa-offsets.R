test_that("DHARMa_diag separates offset variables from fixed predictors", {
  skip_if_not_installed(
    "DHARMa"
  )

  set.seed(1L)

  test_data <- data.frame(
    count = stats::rpois(100L, lambda = 4),
    habitat = factor(rep(c("A", "B"), each = 50L)),
    n_traps = rep(c(2, 4), 50L)
  )

  model <- stats::glm(
    count ~ habitat + offset(log(n_traps)),
    family = stats::poisson(),
    data = test_data
  )

  simulation_matrix <- replicate(
    100L,
    stats::rpois(
      nrow(test_data),
      lambda = stats::fitted(
        model
      )
    )
  )

  residual_object <- DHARMa::createDHARMa(
    simulatedResponse = simulation_matrix,
    observedResponse = test_data$count,
    fittedPredictedResponse = stats::fitted(
      model
    ),
    integerResponse = TRUE
  )

  diagnostic_bundle <- list(
    model = model,
    residuals = residual_object,
    model_status = "valid",
    analysis_data = test_data
  )

  diagnostics <- DHARMa_diag(
    res_object = diagnostic_bundle,
    data = test_data,
    print_plots = FALSE,
    export_plots = FALSE,
    plot_offset_diagnostics = TRUE
  )

  expect_type(
    diagnostics,
    "list"
  )

  expect_true(
    "habitat" %in%
      diagnostics$predictors_tested
  )

  expect_false(
    "n_traps" %in%
      diagnostics$predictors_tested
  )

  expect_true(
    "n_traps" %in%
      diagnostics$offset_variables_tested
  )

  normalized_offset_expressions <- gsub(
    "\\s+",
    "",
    diagnostics$offset_expressions
  )

  expect_true(
    "offset(log(n_traps))" %in%
      normalized_offset_expressions
  )
})

test_that("DHARMa_diag accepts explicitly named offset variables", {
  skip_if_not_installed("DHARMa")

  set.seed(2L)

  test_data <- data.frame(
    count = stats::rpois(60L, lambda = 3),
    x = stats::rnorm(60L),
    effort = rep(c(1, 2, 3), each = 20L)
  )

  model <- stats::glm(
    count ~ x,
    family = stats::poisson(),
    offset = log(test_data$effort),
    data = test_data
  )

  simulation_matrix <- replicate(
    100L,
    stats::rpois(
      nrow(test_data),
      lambda = stats::fitted(
        model
      )
    )
  )

  residual_object <- DHARMa::createDHARMa(
    simulatedResponse = simulation_matrix,
    observedResponse = test_data$count,
    fittedPredictedResponse = stats::fitted(model),
    integerResponse = TRUE
  )

  diagnostic_bundle <- list(
    model = model,
    residuals = residual_object,
    model_status = "valid",
    analysis_data = test_data
  )

  diagnostics <- DHARMa_diag(
    res_object = diagnostic_bundle,
    data = test_data,
    print_plots = FALSE,
    export_plots = FALSE,
    offset_variables = "effort"
  )

  expect_type(
    diagnostics,
    "list"
  )

  expect_true(
    "effort" %in%
      diagnostics$offset_variables_tested
  )

  expect_false(
    "effort" %in%
      diagnostics$predictors_tested
  )

  expect_length(
    diagnostics$offset_expressions,
    0L
  )
})
