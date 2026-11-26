# ==============================================================================
# Offset handling through fit_and_simulate() and DHARMa_diag()
# ==============================================================================

make_offset_test_data <- function(n = 120L, seed = 1201L) {
  if (
    !is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      !is.finite(n) ||
      n < 8L ||
      n != as.integer(n) ||
      n %% 4L != 0L
  ) {
    stop(
      "'n' must be an integer of at least 8 and divisible by 4.",
      call. = FALSE
    )
  }

  n <- as.integer(
    n
  )

  set.seed(seed)

  habitat <- factor(
    rep(
      c("A", "B"),
      length.out = n
    )
  )

  plot_code <- factor(
    rep(
      seq_len(n / 4L),
      each = 4L
    )
  )

  n_traps <- sample(
    x = 1:4,
    size = n,
    replace = TRUE
  )

  x <- stats::rnorm(n)

  linear_predictor <- -0.5 + 0.35 * x + 0.25 * (habitat == "B") + log(n_traps)

  expected_count <- exp(linear_predictor)

  count <- stats::rpois(
    n = n,
    lambda = expected_count
  )

  data.frame(
    count = count,
    x = x,
    habitat = habitat,
    plot_code = plot_code,
    n_traps = n_traps
  )
}

test_that("fit_and_simulate preserves and diagnoses a glmmTMB offset", {
  skip_if_not_installed(
    "glmmTMB"
  )

  skip_if_not_installed(
    "DHARMa"
  )

  skip_if_not_installed(
    "insight"
  )

  test_data <- make_offset_test_data(
    n = 120L,
    seed = 1201L
  )

  model_name <- paste0(
    "test_glmmtmb_offset_",
    as.integer(
      Sys.getpid()
    )
  )

  cache_directory <- tempfile(
    pattern = "ecotools-glmmtmb-offset-"
  )

  dir.create(
    cache_directory,
    recursive = TRUE
  )

  on.exit(
    unlink(
      cache_directory,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  bundle <- fit_and_simulate(
    model_name = model_name,
    package = "glmmTMB",
    model_formula = count ~
      x +
      habitat +
      offset(
        log(
          n_traps
        )
      ) +
      (1 | plot_code),
    family = stats::poisson(),
    data = test_data,
    output_dir = cache_directory,
    force_refit = TRUE,
    num_cores = 1L,
    sparse_screen = FALSE,
    dharma_n_sim = 100L,
    dharma_seed = 1201L
  )

  expect_type(
    bundle,
    "list"
  )

  expect_s3_class(
    bundle$model,
    "glmmTMB"
  )

  expect_false(
    is.null(
      bundle$residuals
    )
  )

  fitted_formula <- stats::formula(
    bundle$model
  )

  fitted_formula_text <- paste(
    deparse(
      fitted_formula,
      width.cutoff = 500L
    ),
    collapse = ""
  )

  expect_match(
    fitted_formula_text,
    "offset"
  )

  expect_match(
    fitted_formula_text,
    "n_traps"
  )

  expect_false(
    is.null(
      bundle$model
    )
  )

  expect_false(
    identical(
      bundle$model_status,
      "invalid"
    )
  )

  expect_false(
    is.null(
      bundle$residuals
    ),
    info = if (
      is.null(
        bundle$simulation_error
      )
    ) {
      paste0(
        "DHARMa residual object is missing; dharma_status = ",
        bundle$dharma_status
      )
    } else {
      paste0(
        "DHARMa simulation failed: ",
        bundle$simulation_error$message
      )
    }
  )

  diagnostic_result <- DHARMa_diag(
    res_object = bundle,
    data = test_data,
    model_title = "glmmTMB offset integration test",
    print_plots = FALSE,
    export_plots = FALSE,
    plot_offset_diagnostics = TRUE
  )

  expect_type(
    diagnostic_result,
    "list"
  )

  expect_true(
    "x" %in%
      diagnostic_result$predictors_tested
  )

  expect_true(
    "habitat" %in%
      diagnostic_result$predictors_tested
  )

  expect_false(
    "n_traps" %in%
      diagnostic_result$predictors_tested
  )

  expect_true(
    "n_traps" %in%
      diagnostic_result$offset_variables_tested
  )

  expect_true(
    any(
      grepl(
        "offset(log(n_traps))",
        gsub(
          "\\s+",
          "",
          diagnostic_result$offset_expressions
        ),
        fixed = TRUE
      )
    )
  )
})

test_that("fit_and_simulate preserves and diagnoses an mgcv GAM offset", {
  skip_if_not_installed(
    "mgcv"
  )

  skip_if_not_installed(
    "mgcViz"
  )

  skip_if_not_installed(
    "DHARMa"
  )

  skip_if_not_installed(
    "insight"
  )

  test_data <- make_offset_test_data(
    n = 120L,
    seed = 1202L
  )

  model_name <- paste0(
    "test_mgcv_gam_offset_",
    as.integer(
      Sys.getpid()
    )
  )

  cache_directory <- tempfile(
    pattern = "ecotools-mgcv-gam-offset-"
  )

  dir.create(
    cache_directory,
    recursive = TRUE
  )

  on.exit(
    unlink(
      cache_directory,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  bundle <- fit_and_simulate(
    model_name = model_name,
    package = "mgcv",
    mgcv_fitter = "gam",
    model_formula = count ~
      s(
        x,
        k = 5
      ) +
      habitat +
      offset(
        log(
          n_traps
        )
      ) +
      s(
        plot_code,
        bs = "re"
      ),
    family = stats::poisson(),
    data = test_data,
    output_dir = cache_directory,
    force_refit = TRUE,
    num_cores = 1L,
    sparse_screen = FALSE,
    dharma_n_sim = 100L,
    dharma_seed = 1202L,
    method = "REML"
  )

  expect_type(
    bundle,
    "list"
  )

  expect_true(
    inherits(
      bundle$model,
      "gam"
    )
  )

  expect_false(
    is.null(
      bundle$residuals
    )
  )

  fitted_formula <- stats::formula(
    bundle$model
  )

  fitted_formula_text <- paste(
    deparse(
      fitted_formula,
      width.cutoff = 500L
    ),
    collapse = ""
  )

  expect_match(
    fitted_formula_text,
    "offset"
  )

  expect_match(
    fitted_formula_text,
    "n_traps"
  )

  expect_false(
    is.null(
      bundle$model
    )
  )

  expect_false(
    identical(
      bundle$model_status,
      "invalid"
    )
  )

  expect_false(
    is.null(
      bundle$residuals
    ),
    info = if (
      is.null(
        bundle$simulation_error
      )
    ) {
      paste0(
        "DHARMa residual object is missing; dharma_status = ",
        bundle$dharma_status
      )
    } else {
      paste0(
        "DHARMa simulation failed: ",
        bundle$simulation_error$message
      )
    }
  )

  diagnostic_result <- DHARMa_diag(
    res_object = bundle,
    data = test_data,
    model_title = "mgcv GAM offset integration test",
    print_plots = FALSE,
    export_plots = FALSE,
    plot_offset_diagnostics = TRUE
  )

  expect_type(
    diagnostic_result,
    "list"
  )

  expect_true(
    "habitat" %in%
      diagnostic_result$predictors_tested
  )

  expect_false(
    "n_traps" %in%
      diagnostic_result$predictors_tested
  )

  expect_true(
    "n_traps" %in%
      diagnostic_result$offset_variables_tested
  )

  expect_true(
    any(
      grepl(
        "offset(log(n_traps))",
        gsub(
          "\\s+",
          "",
          diagnostic_result$offset_expressions
        ),
        fixed = TRUE
      )
    )
  )
})

test_that("fit_and_simulate preserves and diagnoses an mgcv BAM offset", {
  skip_if_not_installed(
    "mgcv"
  )

  skip_if_not_installed(
    "mgcViz"
  )

  skip_if_not_installed(
    "DHARMa"
  )

  skip_if_not_installed(
    "insight"
  )

  test_data <- make_offset_test_data(
    n = 160L,
    seed = 1203L
  )

  model_name <- paste0(
    "test_mgcv_bam_offset_",
    as.integer(
      Sys.getpid()
    )
  )

  cache_directory <- tempfile(
    pattern = "ecotools-mgcv-bam-offset-"
  )

  dir.create(
    cache_directory,
    recursive = TRUE
  )

  on.exit(
    unlink(
      cache_directory,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  bundle <- fit_and_simulate(
    model_name = model_name,
    package = "mgcv",
    mgcv_fitter = "bam",
    model_formula = count ~
      s(
        x,
        k = 5
      ) +
      habitat +
      offset(
        log(
          n_traps
        )
      ),
    family = stats::poisson(),
    data = test_data,
    output_dir = cache_directory,
    force_refit = TRUE,
    num_cores = 1L,
    sparse_screen = FALSE,
    dharma_n_sim = 100L,
    dharma_seed = 1203L,
    method = "fREML",
    discrete = TRUE
  )

  expect_type(
    bundle,
    "list"
  )

  expect_true(
    inherits(
      bundle$model,
      "bam"
    ) ||
      inherits(
        bundle$model,
        "gam"
      )
  )

  expect_false(
    is.null(
      bundle$model
    )
  )

  expect_false(
    identical(
      bundle$model_status,
      "invalid"
    )
  )

  expect_false(
    is.null(
      bundle$residuals
    ),
    info = if (
      is.null(
        bundle$simulation_error
      )
    ) {
      paste0(
        "DHARMa residual object is missing; dharma_status = ",
        bundle$dharma_status
      )
    } else {
      paste0(
        "DHARMa simulation failed: ",
        bundle$simulation_error$message
      )
    }
  )

  diagnostic_result <- DHARMa_diag(
    res_object = bundle,
    data = test_data,
    model_title = "mgcv BAM offset integration test",
    print_plots = FALSE,
    export_plots = FALSE,
    plot_offset_diagnostics = TRUE
  )

  expect_false(
    "n_traps" %in%
      diagnostic_result$predictors_tested
  )

  expect_true(
    "n_traps" %in%
      diagnostic_result$offset_variables_tested
  )
})

test_that("explicit offset variables can be diagnosed through the wrapper", {
  skip_if_not_installed(
    "glmmTMB"
  )

  skip_if_not_installed(
    "DHARMa"
  )

  test_data <- make_offset_test_data(
    n = 100L,
    seed = 1204L
  )

  cache_directory <- tempfile(
    pattern = "ecotools-explicit-offset-"
  )

  dir.create(
    cache_directory,
    recursive = TRUE
  )

  on.exit(
    unlink(
      cache_directory,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  bundle <- fit_and_simulate(
    model_name = "test_explicit_offset_wrapper",
    package = "glmmTMB",
    model_formula = count ~ x + habitat,
    family = stats::poisson(),
    data = test_data,
    offset = log(test_data$n_traps),
    output_dir = cache_directory,
    force_refit = TRUE,
    num_cores = 1L,
    sparse_screen = FALSE,
    dharma_n_sim = 100L,
    dharma_seed = 1204L
  )

  expect_false(
    is.null(
      bundle$model
    )
  )

  expect_false(
    identical(
      bundle$model_status,
      "invalid"
    )
  )

  expect_false(
    is.null(
      bundle$residuals
    ),
    info = if (
      is.null(
        bundle$simulation_error
      )
    ) {
      paste0(
        "DHARMa residual object is missing; dharma_status = ",
        bundle$dharma_status
      )
    } else {
      paste0(
        "DHARMa simulation failed: ",
        bundle$simulation_error$message
      )
    }
  )

  fitted_formula_environment <- environment(
    stats::formula(
      bundle$model
    )
  )

  expect_true(
    exists(
      "..offset",
      envir = fitted_formula_environment,
      inherits = FALSE
    )
  )

  expect_true(
    exists(
      ".ecotools_offset",
      envir = fitted_formula_environment,
      inherits = FALSE
    )
  )

  expect_equal(
    get(
      ".ecotools_offset",
      envir = fitted_formula_environment,
      inherits = FALSE
    ),
    log(
      test_data$n_traps
    )
  )

  expect_identical(
    bundle$model$call$offset,
    as.name(
      ".ecotools_offset"
    )
  )

  regenerated_residuals <- DHARMa::simulateResiduals(
    fittedModel = bundle$model,
    n = 50L,
    seed = 1204L
  )

  expect_s3_class(
    regenerated_residuals,
    "DHARMa"
  )

  diagnostic_result <- DHARMa_diag(
    res_object = bundle,
    data = test_data,
    model_title = "Explicit offset wrapper test",
    print_plots = FALSE,
    export_plots = FALSE,
    plot_offset_diagnostics = TRUE,
    offset_variables = "n_traps"
  )

  expect_false(
    "n_traps" %in%
      diagnostic_result$predictors_tested
  )

  expect_true(
    "n_traps" %in%
      diagnostic_result$offset_variables_tested
  )
})

test_that("explicit glmmTMB offsets survive serialization", {
  skip_if_not_installed("glmmTMB")

  skip_if_not_installed("DHARMa")

  test_data <- make_offset_test_data(
    n = 100L,
    seed = 1205L
  )

  cache_directory <- tempfile(
    pattern = "ecotools-explicit-offset-serialization-"
  )

  dir.create(
    cache_directory,
    recursive = TRUE
  )

  on.exit(
    unlink(
      cache_directory,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  bundle <- fit_and_simulate(
    model_name = "test_explicit_offset_serialized",
    package = "glmmTMB",
    model_formula = count ~
      x +
      habitat,
    family = stats::poisson(),
    data = test_data,
    offset = log(
      test_data$n_traps
    ),
    output_dir = cache_directory,
    force_refit = TRUE,
    num_cores = 1L,
    sparse_screen = FALSE,
    dharma_n_sim = 50L,
    dharma_seed = 1205L
  )

  expect_false(
    is.null(
      bundle$residuals
    ),
    info = if (
      is.null(
        bundle$simulation_error
      )
    ) {
      paste0(
        "DHARMa residual object is missing; dharma_status = ",
        bundle$dharma_status
      )
    } else {
      paste0(
        "DHARMa simulation failed: ",
        bundle$simulation_error$message
      )
    }
  )

  cache_file <- file.path(
    cache_directory,
    "test_explicit_offset_serialized.rds"
  )

  expect_true(
    file.exists(
      cache_file
    )
  )

  restored_bundle <- readRDS(
    cache_file
  )

  expect_s3_class(
    restored_bundle$model,
    "glmmTMB"
  )

  restored_formula_environment <- environment(
    stats::formula(
      restored_bundle$model
    )
  )

  expect_true(
    exists(
      "..offset",
      envir = restored_formula_environment,
      inherits = FALSE
    )
  )

  expect_true(
    exists(
      ".ecotools_offset",
      envir = restored_formula_environment,
      inherits = FALSE
    )
  )

  expect_equal(
    get(
      ".ecotools_offset",
      envir = restored_formula_environment,
      inherits = FALSE
    ),
    log(
      test_data$n_traps
    )
  )

  restored_residuals <- DHARMa::simulateResiduals(
    fittedModel = restored_bundle$model,
    n = 50L,
    seed = 1205L
  )

  expect_s3_class(
    restored_residuals,
    "DHARMa"
  )

  restored_predictions <- stats::predict(
    restored_bundle$model,
    type = "response"
  )

  expect_length(
    restored_predictions,
    nrow(
      test_data
    )
  )

  expect_true(
    all(
      is.finite(
        restored_predictions
      )
    )
  )
})

test_that("fit_and_simulate rejects explicit offsets of incorrect length", {
  skip_if_not_installed(
    "glmmTMB"
  )

  test_data <- make_offset_test_data(
    n = 100L,
    seed = 1206L
  )

  expect_error(
    fit_and_simulate(
      model_name = "test_bad_offset_length",
      package = "glmmTMB",
      model_formula = count ~ x,
      family = stats::poisson(),
      data = test_data,
      offset = rep(
        0,
        99L
      ),
      output_dir = tempdir(),
      force_refit = TRUE,
      sparse_screen = FALSE
    ),
    regexp = "length 1 or match"
  )
})

test_that("fit_and_simulate rejects non-finite explicit offsets", {
  skip_if_not_installed(
    "glmmTMB"
  )

  test_data <- make_offset_test_data(
    n = 100L,
    seed = 1207L
  )

  invalid_offset <- log(
    test_data$n_traps
  )

  invalid_offset[5L] <- NA_real_

  expect_error(
    fit_and_simulate(
      model_name = "test_nonfinite_offset",
      package = "glmmTMB",
      model_formula = count ~ x,
      family = stats::poisson(),
      data = test_data,
      offset = invalid_offset,
      output_dir = tempdir(),
      force_refit = TRUE,
      sparse_screen = FALSE
    ),
    regexp = "missing or non-finite"
  )
})

test_that("fit_and_simulate rejects duplicate offset specifications", {
  skip_if_not_installed(
    "glmmTMB"
  )

  test_data <- make_offset_test_data(
    n = 100L,
    seed = 1208L
  )

  expect_error(
    fit_and_simulate(
      model_name = "test_duplicate_offset",
      package = "glmmTMB",
      model_formula = count ~
        x +
        offset(
          log(
            n_traps
          )
        ),
      family = stats::poisson(),
      data = test_data,
      offset = log(
        test_data$n_traps
      ),
      output_dir = tempdir(),
      force_refit = TRUE,
      sparse_screen = FALSE
    ),
    regexp = "but not both"
  )
})
