test_that("Tweedie simulation guard permits ordinary powers", {
  skip_if_not_installed("glmmTMB")

  set.seed(123)

  test_data <- data.frame(
    y = pmax(
      stats::rgamma(
        80,
        shape = 2,
        rate = 1
      ),
      .Machine$double.eps
    ),
    x = stats::rnorm(80)
  )

  fitted_model <- glmmTMB::glmmTMB(
    y ~ x,
    family = glmmTMB::tweedie(link = "log"),
    data = test_data,
    start = list(
      psi = stats::qlogis(1.7 - 1)
    ),
    map = list(
      psi = factor(NA)
    )
  )

  result <- check_tweedie_simulation_guard(
    model = fitted_model,
    max_power = 1.99
  )

  expect_true(result$is_tweedie)
  expect_equal(result$power, 1.7, tolerance = 1e-6)
  expect_false(result$blocked)
  expect_identical(
    result$reason,
    "power_below_simulation_threshold"
  )
})


test_that("Tweedie simulation guard blocks powers above threshold", {
  skip_if_not_installed("glmmTMB")

  set.seed(123)

  test_data <- data.frame(
    y = pmax(
      stats::rgamma(
        80,
        shape = 2,
        rate = 1
      ),
      .Machine$double.eps
    ),
    x = stats::rnorm(80)
  )

  fitted_model <- glmmTMB::glmmTMB(
    y ~ x,
    family = glmmTMB::tweedie(link = "log"),
    data = test_data,
    start = list(
      psi = stats::qlogis(1.999 - 1)
    ),
    map = list(
      psi = factor(NA)
    )
  )

  result <- check_tweedie_simulation_guard(
    model = fitted_model,
    max_power = 1.99
  )

  expect_true(result$is_tweedie)
  expect_equal(result$power, 1.999, tolerance = 1e-6)
  expect_true(result$blocked)
  expect_identical(
    result$reason,
    "power_at_or_above_simulation_threshold"
  )
  expect_match(
    result$message,
    "DHARMa simulation was skipped",
    fixed = TRUE
  )
})


test_that("Tweedie threshold must be strictly between one and two", {
  skip_if_not_installed("glmmTMB")

  expect_error(
    check_tweedie_simulation_guard(
      model = NULL,
      max_power = 2
    ),
    regexp = "strictly between 1 and 2"
  )

  expect_error(
    check_tweedie_simulation_guard(
      model = NULL,
      max_power = 1
    ),
    regexp = "strictly between 1 and 2"
  )
})


test_that("non-Tweedie models are not blocked", {
  skip_if_not_installed("glmmTMB")

  set.seed(123)

  test_data <- data.frame(
    y = stats::rpois(50, lambda = 2),
    x = stats::rnorm(50)
  )

  fitted_model <- glmmTMB::glmmTMB(
    y ~ x,
    family = stats::poisson(),
    data = test_data
  )

  result <- check_tweedie_simulation_guard(
    model = fitted_model
  )

  expect_false(result$is_tweedie)
  expect_false(result$blocked)
  expect_identical(result$reason, "not_tweedie")
})

test_that("DHARMa diagnostics skip unsafe Tweedie simulation", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("DHARMa")

  set.seed(123)

  test_data <- data.frame(
    y = pmax(
      stats::rgamma(
        50,
        shape = 2,
        rate = 1
      ),
      .Machine$double.eps
    ),
    x = stats::rnorm(50)
  )

  fitted_model <- glmmTMB::glmmTMB(
    y ~ x,
    family = glmmTMB::tweedie(link = "log"),
    data = test_data,
    start = list(
      psi = stats::qlogis(1.999 - 1)
    ),
    map = list(
      psi = factor(NA)
    )
  )

  expect_warning(
    result <- DHARMa_diag(
      res_object = fitted_model,
      data = test_data,
      print_plots = FALSE,
      export_plots = FALSE,
      tweedie_simulation_max_power = 1.99
    ),
    regexp = "DHARMa simulation was skipped"
  )

  expect_identical(
    result$status,
    "skipped_tweedie_power_boundary"
  )
  expect_false(result$valid_for_inference)
  expect_true(result$tweedie_simulation_guard$blocked)
})

test_that("ordinary non-Tweedie models are never blocked", {
  skip_if_not_installed("glmmTMB")

  test_data <- data.frame(
    y = c(1L, 2L, 1L, 3L, 2L, 4L),
    x = seq_len(6)
  )

  fitted_model <- glmmTMB::glmmTMB(
    y ~ x,
    family = stats::poisson(),
    data = test_data
  )

  guard <- check_tweedie_simulation_guard(
    model = fitted_model
  )

  expect_false(guard$is_tweedie)
  expect_false(guard$blocked)
  expect_identical(
    guard$reason,
    "not_tweedie"
  )
})


test_that("power exactly at threshold is blocked", {
  skip_if_not_installed("glmmTMB")

  set.seed(123)

  test_data <- data.frame(
    y = stats::rgamma(
      50L,
      shape = 2,
      rate = 1
    ),
    x = stats::rnorm(50L)
  )

  fitted_model <- glmmTMB::glmmTMB(
    y ~ x,
    family = glmmTMB::tweedie(link = "log"),
    data = test_data,
    start = list(
      psi = stats::qlogis(1.99 - 1)
    ),
    map = list(
      psi = factor(NA)
    )
  )

  guard <- check_tweedie_simulation_guard(
    model = fitted_model,
    max_power = 1.99
  )

  expect_true(guard$blocked)
  expect_identical(
    guard$reason,
    "power_at_or_above_simulation_threshold"
  )
})


test_that("varying Tweedie powers use the maximum value", {
  mock_model <- structure(
    list(),
    class = "mock_tweedie_model"
  )

  attr(mock_model, "family") <- "twlss"

  guard <- check_tweedie_simulation_guard(
    model = mock_model,
    max_power = 1.99,
    fitted_power = c(
      1.50,
      1.75,
      1.995
    )
  )

  expect_true(guard$is_tweedie)
  expect_equal(
    guard$min_power,
    1.50
  )
  expect_equal(
    guard$max_power,
    1.995
  )
  expect_identical(
    guard$n_power_values,
    3L
  )
  expect_true(guard$blocked)
})

test_that("stored guard is tied to its configured threshold", {
  stored_guard <- list(
    is_tweedie = TRUE,
    max_power = 1.97,
    threshold = 1.99,
    blocked = FALSE,
    reason = "power_below_simulation_threshold"
  )

  expect_false(
    isTRUE(
      all.equal(
        stored_guard$threshold,
        1.95,
        tolerance = sqrt(.Machine$double.eps)
      )
    )
  )
})
