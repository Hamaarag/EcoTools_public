test_that("GLLVM response families normalize to one value per response", {
  expect_identical(
    normalize_gllvm_response_family("poisson", 3L),
    rep("poisson", 3L)
  )
  expect_identical(
    normalize_gllvm_response_family(
      c("negative.binomial", "betaH"),
      2L
    ),
    c("negative.binomial", "betaH")
  )
})

test_that("betaH summary preparation trims only a temporary model copy", {
  model <- structure(
    list(
      y = matrix(
        1:6,
        nrow = 2L,
        dimnames = list(
          NULL,
          c(
            "cover_a",
            "count",
            "cover_b"
          )
        )
      ),
      family = c(
        "betaH",
        "negative.binomial",
        "betaH",
        "betaH",
        "betaH"
      ),
      response_family = c(
        "betaH",
        "negative.binomial",
        "betaH"
      ),
      params = list(
        beta0 = seq_len(5L),
        theta = matrix(
          seq_len(10L),
          nrow = 5L
        )
      ),
      sd = list(
        beta0 = seq_len(5L),
        theta = matrix(
          seq_len(10L),
          nrow = 5L
        )
      )
    ),
    class = "gllvm"
  )

  original_model <- model

  summary_model <- prepare_gllvm_summary_object(
    model
  )

  # The temporary summary copy has one entry per observed response.
  expect_length(
    summary_model$params$beta0,
    3L
  )

  expect_equal(
    dim(summary_model$params$theta),
    c(
      3L,
      2L
    )
  )

  expect_length(
    summary_model$sd$beta0,
    3L
  )

  expect_equal(
    dim(summary_model$sd$theta),
    c(
      3L,
      2L
    )
  )

  expect_identical(
    summary_model$family,
    c(
      "betaH",
      "negative.binomial",
      "betaH"
    )
  )

  # The original fitted object retains gllvm's component-level representation.
  expect_identical(
    model$family,
    original_model$family
  )

  expect_identical(
    model$response_family,
    original_model$response_family
  )

  expect_identical(
    model$params$beta0,
    original_model$params$beta0
  )

  expect_identical(
    model$params$theta,
    original_model$params$theta
  )

  expect_identical(
    model$sd$beta0,
    original_model$sd$beta0
  )

  expect_identical(
    model$sd$theta,
    original_model$sd$theta
  )
})

test_that("betaH predictions combine hurdle and conditional means", {
  component_predictions <- cbind(
    conditional_mean_1 = c(0.2, 0.8),
    conditional_mean_2 = c(0.5, 0.4),
    positive_probability_1 = c(0.7, 0.25),
    positive_probability_2 = c(0.9, 0.6)
  )

  response_predictions <- collapse_gllvm_family_predictions(
    predictions = component_predictions,
    family = "betaH",
    response_values = matrix(0.1, nrow = 2L, ncol = 2L)
  )

  expect_equal(
    response_predictions,
    cbind(
      c(0.14, 0.2),
      c(0.45, 0.24)
    )
  )

  expect_error(
    collapse_gllvm_family_predictions(
      predictions = component_predictions[, 1:2],
      family = "betaH",
      response_values = matrix(0.1, nrow = 2L, ncol = 2L)
    ),
    "Unexpected prediction dimensions"
  )
})

test_that("ordinal predictions reduce to expected category scores", {
  category_probabilities <- array(
    0,
    dim = c(3L, 2L, 2L)
  )
  category_probabilities[, 1L, 1L] <- c(0.2, 0.5, 0.3)
  category_probabilities[, 2L, 1L] <- c(0.1, 0.4, 0.5)
  category_probabilities[, 1L, 2L] <- c(0.3, 0.4, 0.3)
  category_probabilities[, 2L, 2L] <- c(0.2, 0.3, 0.5)
  observed_responses <- cbind(
    ordinal_zero_based = c(0, 2),
    ordinal_one_based = c(1, 3)
  )

  expected_scores <- collapse_gllvm_family_predictions(
    predictions = category_probabilities,
    family = "ordinal",
    response_values = observed_responses
  )

  expect_equal(
    expected_scores,
    structure(
      cbind(
        c(1.1, 1.4),
        c(2, 2.3)
      ),
      dimnames = dimnames(observed_responses)
    )
  )
})

test_that("orderedBeta predictions reduce to expected response means", {
  category_probabilities <- array(
    0,
    dim = c(3L, 2L, 2L)
  )
  category_probabilities[, 1L, 1L] <- c(0.2, 0.5, 0.3)
  category_probabilities[, 2L, 1L] <- c(0.1, 0.4, 0.5)
  category_probabilities[, 1L, 2L] <- c(0.3, 0.4, 0.3)
  category_probabilities[, 2L, 2L] <- c(0.2, 0.3, 0.5)
  conditional_mean <- cbind(
    c(0.4, 0.6),
    c(0.5, 0.7)
  )

  expected_values <- collapse_gllvm_family_predictions(
    predictions = category_probabilities,
    family = "orderedBeta",
    response_values = matrix(0.5, nrow = 2L, ncol = 2L),
    conditional_mean = conditional_mean
  )

  expect_equal(
    expected_values,
    cbind(
      c(0.5, 0.74),
      c(0.5, 0.71)
    )
  )
})

test_that("binomial-family predictions include trials and inflation", {
  base_probability <- matrix(c(0.2, 0.8), ncol = 1L)
  response_values <- matrix(c(2, 8), ncol = 1L)
  n_trials <- matrix(10, nrow = 2L, ncol = 1L)

  expect_equal(
    collapse_gllvm_family_predictions(
      predictions = base_probability,
      family = "binomial",
      response_values = response_values,
      n_trials = n_trials
    ),
    matrix(c(2, 8), ncol = 1L)
  )
  expect_equal(
    collapse_gllvm_family_predictions(
      predictions = base_probability,
      family = "ZIB",
      response_values = response_values,
      n_trials = n_trials,
      inflation_probability = matrix(0.25, nrow = 2L, ncol = 1L)
    ),
    matrix(c(1.5, 6), ncol = 1L)
  )
  expect_equal(
    collapse_gllvm_family_predictions(
      predictions = base_probability,
      family = "ZNIB",
      response_values = response_values,
      n_trials = n_trials,
      inflation_probability = matrix(0.1, nrow = 2L, ncol = 1L),
      upper_inflation_probability = matrix(0.2, nrow = 2L, ncol = 1L)
    ),
    matrix(c(3.4, 7.6), ncol = 1L)
  )
  expect_equal(
    collapse_gllvm_family_predictions(
      predictions = base_probability,
      family = "beta.binomial",
      response_values = response_values,
      n_trials = n_trials
    ),
    matrix(c(2, 8), ncol = 1L)
  )
})

test_that("single-component GLLVM predictions pass through by response", {
  supported_single_column_families <- c(
    "poisson",
    "negative.binomial",
    "negative.binomial1",
    "ZIP",
    "ZINB",
    "gaussian",
    "tweedie",
    "gamma",
    "exponential",
    "beta"
  )
  predictions <- matrix(c(0.2, 0.8), ncol = 1L)
  response_values <- matrix(c(0.1, 0.9), ncol = 1L)

  for (family in supported_single_column_families) {
    expect_identical(
      collapse_gllvm_family_predictions(
        predictions = predictions,
        family = family,
        response_values = response_values
      ),
      predictions,
      info = family
    )
  }
})

test_that("GLLVM response family counts are validated before fitting", {
  expect_error(
    fit_gllvm_cached(
      y = matrix(1:6, ncol = 2L),
      X = data.frame(x = 1:3),
      studyDesign = data.frame(site = 1:3),
      formula = ~x,
      model_name = "invalid_family_count",
      family = c("negative.binomial", "betaH", "poisson"),
      design_check = FALSE,
      plot_diag = FALSE
    ),
    "one family per response column"
  )

  expect_error(
    normalize_gllvm_response_family(
      c("poisson", NA_character_),
      2L
    ),
    "non-missing"
  )
})

test_that("betaH summary preparation accepts response-level arrays", {
  model <- structure(
    list(
      y = matrix(
        1:6,
        nrow = 2L,
        dimnames = list(
          NULL,
          c(
            "cover_a",
            "count",
            "cover_b"
          )
        )
      ),
      family = c(
        "betaH",
        "negative.binomial",
        "betaH"
      ),
      response_family = c(
        "betaH",
        "negative.binomial",
        "betaH"
      ),
      params = list(
        beta0 = seq_len(3L),
        theta = matrix(
          seq_len(6L),
          nrow = 3L
        )
      ),
      sd = list(
        beta0 = seq_len(3L),
        theta = matrix(
          seq_len(6L),
          nrow = 3L
        )
      )
    ),
    class = "gllvm"
  )

  summary_model <- prepare_gllvm_summary_object(
    model
  )

  expect_length(
    summary_model$params$beta0,
    3L
  )

  expect_equal(
    dim(summary_model$params$theta),
    c(
      3L,
      2L
    )
  )

  expect_identical(
    summary_model$family,
    model$response_family
  )
})

test_that("expanded betaH summaries require response-family metadata", {
  model <- structure(
    list(
      y = matrix(
        1:6,
        nrow = 2L,
        ncol = 3L
      ),
      family = c(
        "betaH",
        "negative.binomial",
        "betaH",
        "betaH",
        "betaH"
      ),
      params = list(
        beta0 = seq_len(5L),
        theta = matrix(
          seq_len(10L),
          nrow = 5L
        )
      )
    ),
    class = "gllvm"
  )

  expect_error(
    prepare_gllvm_summary_object(
      model
    ),
    "has no valid 'response_family' metadata"
  )
})
