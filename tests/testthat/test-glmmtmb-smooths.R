test_that("namespaced mgcv smooth calls are normalized for glmmTMB", {
  input_formula <-
    y ~ mgcv::s(x, k = 3) + treatment + (1 | site)

  normalized_formula <-
    EcoTools:::normalize_glmmtmb_smooth_formula(
      input_formula
    )

  normalized_text <- paste(
    deparse(
      normalized_formula
    ),
    collapse = " "
  )

  expect_match(
    normalized_text,
    "s\\(x, k = 3\\)"
  )

  expect_false(
    grepl(
      "mgcv::s",
      normalized_text,
      fixed = TRUE
    )
  )
})


test_that("parametric design formula excludes smooth terms", {
  input_formula <-
    y ~ s(x, by = treatment, k = 3) + treatment

  parametric_formula <-
    EcoTools:::strip_mgcv_smooth_terms(
      input_formula
    )

  parametric_text <- paste(
    deparse(
      parametric_formula
    ),
    collapse = " "
  )

  expect_match(
    parametric_text,
    "treatment"
  )

  expect_false(
    grepl(
      "s\\(",
      parametric_text
    )
  )
})


test_that("glmmTMB can fit an mgcv smooth through fit_and_simulate", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("mgcv")
  skip_if_not_installed("DHARMa")

  set.seed(123)

  n <- 120L

  test_data <- data.frame(
    x = seq(
      -1,
      1,
      length.out = n
    ),
    treatment = factor(
      rep(
        c(
          "A",
          "B"
        ),
        length.out = n
      )
    ),
    site = factor(
      rep(
        seq_len(24L),
        each = 5L
      )
    )
  )

  test_data$y <- stats::rpois(
    n,
    lambda = exp(
      1 +
        0.4 *
          sin(
            pi * test_data$x
          ) +
        0.2 *
          (test_data$treatment == "B")
    )
  )

  result <- fit_and_simulate(
    model_name = "test_glmmtmb_smooth",
    package = "glmmTMB",
    model_formula = y ~
      s(
        x,
        k = 3
      ) +
      treatment +
      (1 | site),
    family = stats::poisson(),
    data = test_data,
    output_dir = tempdir(),
    force_refit = TRUE,
    sparse_screen = FALSE,
    num_cores = 1L,
    dharma_n_sim = 50L,
    dharma_seed = 123L,
    REML = TRUE
  )

  expect_s3_class(
    result$model,
    "glmmTMB"
  )

  expect_false(
    anyDuplicated(
      names(result)
    ) >
      0L,
    info = paste0(
      "Duplicated output-bundle names: ",
      paste(
        unique(
          names(result)[
            duplicated(
              names(result)
            )
          ]
        ),
        collapse = ", "
      )
    )
  )

  expect_true(
    result$glmmTMB_smooth_terms
  )

  expect_true(
    result$design_check$smooth_terms_deferred_to_glmmTMB
  )

  expect_identical(
    result$dharma_status,
    "completed_unconditional_glmmTMB_smooth"
  )

  expect_silent(
    stats::predict(
      result$model,
      type = "response"
    )
  )
})

test_that("glmmTMB by-variable smooths are rejected clearly", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("mgcv")

  set.seed(123)

  test_data <- data.frame(
    y = stats::rpois(
      40L,
      lambda = 3
    ),
    x = seq(
      -1,
      1,
      length.out = 40L
    ),
    treatment = factor(
      rep(
        c(
          "A",
          "B"
        ),
        each = 20L
      )
    ),
    site = factor(
      rep(
        seq_len(8L),
        each = 5L
      )
    )
  )

  expect_error(
    fit_and_simulate(
      model_name = "test_unsupported_glmmtmb_by_smooth",
      package = "glmmTMB",
      model_formula = y ~
        s(
          x,
          by = treatment,
          k = 3
        ) +
        treatment +
        (1 | site),
      family = stats::poisson(),
      data = test_data,
      output_dir = tempdir(),
      force_refit = TRUE,
      sparse_screen = FALSE,
      num_cores = 1L
    ),
    regexp = "cannot handle by-variable smooths"
  )
})

test_that("glmmTMB can combine a smooth with the COM-Poisson family", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("mgcv")
  skip_if_not_installed("DHARMa")

  set.seed(456)

  n <- 100L

  test_data <- data.frame(
    x = seq(
      -1,
      1,
      length.out = n
    ),
    group = factor(
      rep(
        seq_len(20L),
        each = 5L
      )
    )
  )

  test_data$y <- stats::rpois(
    n,
    lambda = exp(
      1 +
        0.3 *
          sin(
            pi * test_data$x
          )
    )
  )

  result <- fit_and_simulate(
    model_name = "test_glmmtmb_compois_smooth",
    package = "glmmTMB",
    model_formula = y ~
      s(
        x,
        k = 3
      ) +
      (1 | group),
    family = glmmTMB::compois(),
    data = test_data,
    output_dir = tempdir(),
    force_refit = TRUE,
    sparse_screen = FALSE,
    num_cores = 1L,
    dharma_n_sim = 50L,
    dharma_seed = 456L,
    REML = TRUE
  )

  expect_s3_class(
    result$model,
    "glmmTMB"
  )

  expect_false(
    anyDuplicated(
      names(result)
    ) >
      0L,
    info = paste0(
      "Duplicated output-bundle names: ",
      paste(
        unique(
          names(result)[
            duplicated(
              names(result)
            )
          ]
        ),
        collapse = ", "
      )
    )
  )

  expect_true(
    result$glmmTMB_smooth_terms
  )

  expect_false(
    is.null(
      result$residuals
    )
  )

  expect_null(
    result$simulation_error
  )

  expect_identical(
    result$dharma_status,
    "completed_unconditional_glmmTMB_smooth"
  )

  expect_identical(
    result$model_call_data_storage,
    "embedded_evaluated_data"
  )

  expect_true(
    isTRUE(
      result$glmmTMB_smooth_reconstruction$pass
    )
  )

  expect_identical(
    result$glmmTMB_smooth_reconstruction$n_predictions,
    nrow(test_data)
  )

  expect_null(
    result$glmmTMB_smooth_reconstruction$error
  )

  expect_true(
    is.data.frame(
      result$model$call$data
    )
  )

  expect_equal(
    nrow(
      result$model$call$data
    ),
    nrow(test_data)
  )

  expect_identical(
    names(
      result$model$call$data
    ),
    names(test_data)
  )

  population_predictions <- stats::predict(
    result$model,
    type = "response",
    re.form = ~0
  )

  expect_length(
    population_predictions,
    nrow(test_data)
  )

  expect_true(
    all(
      is.finite(
        population_predictions
      )
    )
  )

  expect_match(
    result$model$modelInfo$family$family,
    "compois",
    ignore.case = TRUE
  )
})

test_that("serialized glmmTMB smooth model retains local data access", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("mgcv")
  skip_if_not_installed("DHARMa")

  set.seed(789)

  cache_dir <- tempfile(
    pattern = "ecotools-smooth-cache-"
  )

  dir.create(
    cache_dir,
    recursive = TRUE
  )

  on.exit(
    unlink(
      cache_dir,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  fit_result <- local({
    local_data <- data.frame(
      x = seq(
        -1,
        1,
        length.out = 80L
      ),
      group = factor(
        rep(
          seq_len(16L),
          each = 5L
        )
      )
    )

    local_data$y <- stats::rpois(
      nrow(local_data),
      lambda = exp(
        1 +
          0.25 *
            sin(
              pi * local_data$x
            )
      )
    )

    fit_and_simulate(
      model_name = "test_serialized_glmmtmb_smooth",
      package = "glmmTMB",
      model_formula = y ~
        s(
          x,
          k = 3
        ) +
        (1 | group),
      family = stats::poisson(),
      data = local_data,
      output_dir = cache_dir,
      force_refit = TRUE,
      sparse_screen = FALSE,
      num_cores = 1L,
      dharma_n_sim = 50L,
      dharma_seed = 789L,
      residual_alpha = 1e-8,
      REML = TRUE
    )
  })

  bundle_path <- file.path(
    cache_dir,
    "test_serialized_glmmtmb_smooth.rds"
  )

  expect_true(
    file.exists(bundle_path)
  )

  reloaded_result <- readRDS(
    bundle_path
  )

  expect_s3_class(
    reloaded_result$model,
    "glmmTMB"
  )

  expect_identical(
    reloaded_result$model_call_data_storage,
    "embedded_evaluated_data"
  )

  expect_true(
    isTRUE(
      reloaded_result$glmmTMB_smooth_reconstruction$pass
    )
  )

  expect_identical(
    reloaded_result$glmmTMB_smooth_reconstruction$n_predictions,
    nrow(
      reloaded_result$analysis_data
    )
  )

  expect_null(
    reloaded_result$glmmTMB_smooth_reconstruction$error
  )

  expect_true(
    is.data.frame(
      reloaded_result$model$call$data
    )
  )

  expect_equal(
    reloaded_result$model$call$data,
    reloaded_result$analysis_data
  )

  population_predictions <- stats::predict(
    reloaded_result$model,
    type = "response",
    re.form = ~0
  )

  expect_length(
    population_predictions,
    nrow(
      reloaded_result$analysis_data
    )
  )

  expect_true(
    all(
      is.finite(
        population_predictions
      )
    )
  )
})
