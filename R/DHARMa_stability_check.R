# ==============================================================================
# Repeated DHARMa diagnostic stability assessment
# ==============================================================================

#' Assess the Stability of DHARMa Diagnostics
#'
#' @description
#' Repeats DHARMa residual simulations using different random seeds and
#' summarizes how frequently global and predictor-specific diagnostic tests are
#' significant or borderline.
#'
#' The function is intended to distinguish persistent residual problems from
#' results that occur only under a small proportion of DHARMa simulations.
#'
#' Models with suspicious residual diagnostics are not rejected by this
#' function because repeated diagnostics are useful for assessing the
#' persistence of the problem. Models with status `"invalid"` or `"aborted"`
#' are rejected because residual diagnostics should not be used to rehabilitate
#' a failed numerical fit.
#'
#' @param model_bundle
#' A model bundle returned by [fit_and_simulate()] or a compatible directly
#' fitted model object.
#'
#' @param predictors
#' Optional data frame containing predictors to test. When `NULL`, the function
#' first uses `model_bundle$analysis_data` when available and otherwise attempts
#' to recover model data and fixed-effect predictors using the insight package.
#'
#' @param n_runs
#' Positive integer giving the number of repeated DHARMa simulation runs.
#'
#' @param n_sim
#' Positive integer giving the number of simulated responses generated during
#' each run.
#'
#' @param tweedie_simulation_max_power
#' Finite numeric threshold strictly between 1 and 2. The repeated simulation
#' procedure is skipped when a fitted Tweedie model has an estimated power
#' at or above this threshold. The default is `1.99`.
#'
#' @param alpha
#' Numeric significance threshold used to classify diagnostic results as
#' significant.
#'
#' @param borderline
#' Numeric threshold used to classify diagnostic results as borderline.
#' This value must be greater than or equal to `alpha`.
#'
#' @param seed
#' Optional integer seed controlling generation of the sequence of DHARMa
#' simulation seeds. Each run receives a different generated seed. When `NULL`,
#' the current random-number stream is used.
#'
#' @param min_category_n
#' Minimum number of complete observations required within a categorical level
#' before a within-level uniformity test is performed. Levels below this
#' threshold are recorded as excluded rather than treated as having passed.
#'
#' @param p_adjust_method
#' Method passed to [stats::p.adjust()] when combining the within-level
#' uniformity p-values for a categorical predictor. The default is
#' `"BH"`.
#'
#' @param verbose
#' Logical. If `TRUE`, print progress and the final stability summary.
#'
#' @param cache_dir
#' Optional character string specifying a directory in which the complete
#' stability-check result is cached as an RDS file. When `NULL`, caching is
#' disabled.
#'
#' @param cache_name
#' Optional character scalar used as the cached RDS file name, without the
#' extension. When `NULL`, the function attempts to derive a name from the
#' model bundle. This argument is ignored when `cache_dir = NULL`.
#'
#' @param force_rerun
#' Logical. If `TRUE`, ignore any matching cached result and repeat all DHARMa
#' simulations. If `FALSE`, reuse a cache only when its stored fingerprint
#' matches the current model and diagnostic settings.
#'
#' @return
#' When the Tweedie simulation guard is triggered, no simulation runs are
#' performed. The returned object records zero completed runs and has
#' `status = "skipped_tweedie_power_boundary"` and
#' `valid_for_inference = FALSE`.
#'
#' Invisibly returns a list containing:
#'
#' \describe{
#'   \item{results}{
#'     A data table with one row for every diagnostic attempted in every run.
#'   }
#'   \item{summary}{
#'     A data table summarizing valid-run counts, p-value quantiles, and the
#'     proportions of significant and borderline results.
#'   }
#'   \item{excluded_levels}{
#'     A data table listing categorical levels excluded because they contained
#'     fewer than `min_category_n` complete observations.
#'   }
#'   \item{seeds}{
#'     Integer vector containing the DHARMa seed used for every repeated run.
#'   }
#'   \item{predictors_used}{
#'     Character vector naming the predictors examined.
#'   }
#'   \item{settings}{
#'     List containing the diagnostic thresholds and simulation settings.
#'   }
#'   \item{analysis_fingerprint}{
#'     Character hash identifying the parent model and stability-check
#'     configuration.
#'   }
#'   \item{cache_file}{
#'     Path to the cached RDS file, or `NULL` when caching is disabled.
#'   }
#' }
#'
#' @details
#' The global diagnostics include:
#'
#' \itemize{
#'   \item residual uniformity;
#'   \item residual dispersion;
#'   \item outlier frequency;
#'   \item zero inflation when the observed response contains zeros;
#'   \item conditional residual quantiles against fitted predictions.
#' }
#'
#' Numeric predictors are assessed using [DHARMa::testQuantiles()].
#'
#' Categorical predictors are assessed using within-level uniformity tests and
#' a robust Fligner-Killeen homogeneity-of-variance test. Within-level
#' uniformity p-values are adjusted within each predictor and simulation run.
#'
#' The returned proportions describe diagnostic stability rather than creating
#' an automatic validity decision. A diagnostic that is significant in most
#' runs provides stronger evidence of persistent lack of fit than an isolated
#' significant result.
#'
#' @seealso
#' [fit_and_simulate()],
#' [DHARMa_diag()],
#' [DHARMa::simulateResiduals()],
#' [DHARMa::testUniformity()],
#' [DHARMa::testDispersion()],
#' [DHARMa::testQuantiles()]
#'
#' @export
DHARMa_stability_check <- function(
  model_bundle,
  predictors = NULL,
  n_runs = 50L,
  n_sim = 1000L,
  tweedie_simulation_max_power = 1.99,
  alpha = 0.05,
  borderline = 0.10,
  seed = NULL,
  min_category_n = 8L,
  p_adjust_method = "BH",
  verbose = TRUE,
  cache_dir = NULL,
  cache_name = NULL,
  force_rerun = FALSE
) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    stop(
      "Package 'DHARMa' is required.",
      call. = FALSE
    )
  }

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
      "Package 'data.table' is required.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 1. Validate arguments
  # ---------------------------------------------------------------------------
  validate_positive_integer <- function(x, argument_name) {
    if (
      !is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        !is.finite(x) ||
        x < 1 ||
        x != as.integer(x)
    ) {
      stop(
        paste0(
          "'",
          argument_name,
          "' must be one positive integer."
        ),
        call. = FALSE
      )
    }

    as.integer(x)
  }

  n_runs <- validate_positive_integer(
    n_runs,
    "n_runs"
  )

  n_sim <- validate_positive_integer(
    n_sim,
    "n_sim"
  )

  min_category_n <- validate_positive_integer(
    min_category_n,
    "min_category_n"
  )

  if (
    !is.numeric(tweedie_simulation_max_power) ||
      length(tweedie_simulation_max_power) != 1L ||
      is.na(tweedie_simulation_max_power) ||
      !is.finite(tweedie_simulation_max_power) ||
      tweedie_simulation_max_power <= 1 ||
      tweedie_simulation_max_power >= 2
  ) {
    stop(
      paste0(
        "'tweedie_simulation_max_power' must be one finite numeric ",
        "value strictly between 1 and 2."
      ),
      call. = FALSE
    )
  }

  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1
  ) {
    stop(
      "'alpha' must be one finite number strictly between zero and one.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(borderline) ||
      length(borderline) != 1L ||
      is.na(borderline) ||
      !is.finite(borderline) ||
      borderline < alpha ||
      borderline >= 1
  ) {
    stop(
      paste0(
        "'borderline' must be finite, greater than or equal to 'alpha', ",
        "and less than one."
      ),
      call. = FALSE
    )
  }

  if (
    !is.character(p_adjust_method) ||
      length(p_adjust_method) != 1L ||
      is.na(p_adjust_method) ||
      !p_adjust_method %in% stats::p.adjust.methods
  ) {
    stop(
      paste0(
        "'p_adjust_method' must be one of: ",
        paste(
          stats::p.adjust.methods,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  if (
    !is.logical(verbose) ||
      length(verbose) != 1L ||
      is.na(verbose)
  ) {
    stop(
      "'verbose' must be one non-missing logical value.",
      call. = FALSE
    )
  }

  if (
    !is.logical(force_rerun) ||
      length(force_rerun) != 1L ||
      is.na(force_rerun)
  ) {
    stop(
      "'force_rerun' must be one non-missing logical value.",
      call. = FALSE
    )
  }

  if (
    !is.null(cache_dir) &&
      (!is.character(cache_dir) ||
        length(cache_dir) != 1L ||
        is.na(cache_dir) ||
        !nzchar(cache_dir))
  ) {
    stop(
      "'cache_dir' must be NULL or one non-empty character value.",
      call. = FALSE
    )
  }

  if (
    !is.null(cache_name) &&
      (!is.character(cache_name) ||
        length(cache_name) != 1L ||
        is.na(cache_name) ||
        !nzchar(cache_name))
  ) {
    stop(
      "'cache_name' must be NULL or one non-empty character value.",
      call. = FALSE
    )
  }

  if (
    !is.null(cache_dir) &&
      is.null(seed)
  ) {
    stop(
      paste0(
        "A non-NULL 'seed' is required when stability-result caching is ",
        "enabled. This ensures that the cached simulation experiment is ",
        "fully reproducible."
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 2. Extract the fitted model and stored metadata
  # ---------------------------------------------------------------------------
  is_bundle <- is.list(model_bundle) &&
    "model" %in% names(model_bundle)

  fitted_model <- if (is_bundle) {
    model_bundle$model
  } else {
    model_bundle
  }

  if (is.null(fitted_model)) {
    stop(
      "No fitted model was found in 'model_bundle'.",
      call. = FALSE
    )
  }

  stored_tweedie_guard <- if (
    is_bundle &&
      !is.null(model_bundle$tweedie_simulation_guard)
  ) {
    model_bundle$tweedie_simulation_guard
  } else {
    NULL
  }

  stored_guard_matches_threshold <- !is.null(stored_tweedie_guard) &&
    isTRUE(stored_tweedie_guard$is_tweedie) &&
    !is.null(stored_tweedie_guard$threshold) &&
    length(stored_tweedie_guard$threshold) == 1L &&
    is.finite(stored_tweedie_guard$threshold) &&
    isTRUE(
      all.equal(
        as.numeric(stored_tweedie_guard$threshold),
        as.numeric(tweedie_simulation_max_power),
        tolerance = sqrt(.Machine$double.eps)
      )
    )

  tweedie_simulation_guard <- if (isTRUE(stored_guard_matches_threshold)) {
    stored_tweedie_guard
  } else {
    check_tweedie_simulation_guard(
      model = fitted_model,
      max_power = tweedie_simulation_max_power
    )
  }

  if (isTRUE(tweedie_simulation_guard$blocked)) {
    warning(
      tweedie_simulation_guard$message,
      call. = FALSE
    )

    skipped_result <- list(
      status = "skipped_tweedie_power_boundary",
      valid_for_inference = FALSE,
      tweedie_simulation_guard = tweedie_simulation_guard,
      runs_requested = as.integer(n_runs),
      runs_completed = 0L,
      results = data.table::data.table(),
      summary = data.table::data.table(
        diagnostic = "Tweedie simulation guard",
        status = "skipped",
        reason = tweedie_simulation_guard$reason,
        fitted_power = tweedie_simulation_guard$max_power,
        threshold = tweedie_simulation_guard$threshold
      ),
      excluded_levels = data.table::data.table(),
      seeds = integer(0),
      predictors_used = character(0),
      settings = list(
        n_runs = n_runs,
        n_sim = n_sim,
        tweedie_simulation_max_power = tweedie_simulation_max_power,
        alpha = alpha,
        borderline = borderline,
        seed = seed,
        min_category_n = min_category_n,
        p_adjust_method = p_adjust_method
      ),
      analysis_fingerprint = NULL,
      cache_file = NULL
    )

    return(
      invisible(skipped_result)
    )
  }

  model_status <- if (
    is_bundle &&
      !is.null(model_bundle$model_status) &&
      length(model_bundle$model_status) == 1L
  ) {
    as.character(model_bundle$model_status)
  } else {
    NA_character_
  }

  if (
    !is.na(model_status) &&
      model_status %in% c("invalid", "aborted")
  ) {
    stop(
      paste0(
        "DHARMa stability diagnostics are not available for a model with ",
        "status '",
        model_status,
        "'."
      ),
      call. = FALSE
    )
  }

  if (
    !is.na(model_status) &&
      model_status %in%
        c(
          "suspicious_extreme_coefficients",
          "suspicious_residual_diagnostics",
          "valid_fit_no_residual_diagnostics"
        )
  ) {
    warning(
      paste0(
        "Running repeated DHARMa diagnostics for a model with status '",
        model_status,
        "'. The output should be used to assess diagnostic persistence, ",
        "not to override a numerical fitting failure."
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 3. Generate reproducible run-specific seeds
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) {
    if (
      !is.numeric(seed) ||
        length(seed) != 1L ||
        is.na(seed) ||
        !is.finite(seed)
    ) {
      stop(
        "'seed' must be NULL or one finite numeric value.",
        call. = FALSE
      )
    }

    existing_seed <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )

    if (existing_seed) {
      old_seed <- get(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    }

    on.exit(
      {
        if (existing_seed) {
          assign(
            ".Random.seed",
            old_seed,
            envir = .GlobalEnv
          )
        } else if (
          exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
          )
        ) {
          rm(
            ".Random.seed",
            envir = .GlobalEnv
          )
        }
      },
      add = TRUE
    )

    set.seed(as.integer(seed))
  }

  dharma_seeds <- sample.int(
    .Machine$integer.max,
    n_runs,
    replace = FALSE
  )

  # ---------------------------------------------------------------------------
  # 4. Local result-extraction helpers
  # ---------------------------------------------------------------------------
  extract_p_value <- function(test_object) {
    if (
      is.null(test_object) ||
        is.null(test_object$p.value) ||
        length(test_object$p.value) == 0L
    ) {
      return(NA_real_)
    }

    p_value <- suppressWarnings(
      as.numeric(test_object$p.value[1L])
    )

    if (
      length(p_value) != 1L ||
        !is.finite(p_value) ||
        p_value < 0 ||
        p_value > 1
    ) {
      return(NA_real_)
    }

    p_value
  }

  append_result <- function(
    result_list,
    run,
    diagnostic,
    variable = NA_character_,
    p_value = NA_real_
  ) {
    result_list[[length(result_list) + 1L]] <-
      data.table::data.table(
        run = as.integer(run),
        diagnostic = as.character(diagnostic),
        variable = as.character(variable),
        p.value = as.numeric(p_value)
      )

    result_list
  }

  categorical_tests <- function(
    scaled_residuals,
    group,
    variable_name,
    run_number
  ) {
    result <- list(
      uniformity_p = NA_real_,
      variance_p = NA_real_,
      excluded = data.table::data.table(
        run = integer(0),
        variable = character(0),
        level = character(0),
        n = integer(0),
        reason = character(0)
      )
    )

    if (length(group) != length(scaled_residuals)) {
      return(result)
    }

    valid_rows <- !is.na(group) &
      is.finite(scaled_residuals)

    group_use <- droplevels(
      as.factor(group[valid_rows])
    )

    residuals_use <- scaled_residuals[
      valid_rows
    ]

    if (
      length(group_use) == 0L ||
        nlevels(group_use) < 2L
    ) {
      return(result)
    }

    group_sizes <- table(group_use)

    eligible_levels <- names(group_sizes)[
      group_sizes >= min_category_n
    ]

    excluded_levels <- names(group_sizes)[
      group_sizes < min_category_n
    ]

    if (length(excluded_levels) > 0L) {
      result$excluded <- data.table::data.table(
        run = as.integer(run_number),
        variable = variable_name,
        level = excluded_levels,
        n = as.integer(
          group_sizes[excluded_levels]
        ),
        reason = paste0(
          "fewer than ",
          min_category_n,
          " complete observations"
        )
      )
    }

    if (length(eligible_levels) > 0L) {
      uniformity_p_values <- vapply(
        eligible_levels,
        function(group_level) {
          values <- residuals_use[
            group_use == group_level
          ]

          values <- values[
            is.finite(values)
          ]

          if (length(values) < min_category_n) {
            return(NA_real_)
          }

          tryCatch(
            suppressWarnings(
              stats::ks.test(
                values,
                "punif"
              )$p.value
            ),
            error = function(e) {
              NA_real_
            }
          )
        },
        numeric(1)
      )

      valid_uniformity <- is.finite(
        uniformity_p_values
      )

      if (any(valid_uniformity)) {
        adjusted_uniformity <- stats::p.adjust(
          uniformity_p_values[
            valid_uniformity
          ],
          method = p_adjust_method
        )

        # The minimum adjusted within-level p-value represents the strongest
        # evidence that at least one eligible level is non-uniform.
        result$uniformity_p <- min(
          adjusted_uniformity,
          na.rm = TRUE
        )
      }
    }

    # Variance testing uses all levels with at least two complete observations.
    variance_levels <- names(group_sizes)[
      group_sizes >= 2L
    ]

    if (length(variance_levels) >= 2L) {
      keep_for_variance <- group_use %in%
        variance_levels

      variance_group <- droplevels(
        group_use[keep_for_variance]
      )

      variance_residuals <- residuals_use[
        keep_for_variance
      ]

      result$variance_p <- tryCatch(
        stats::fligner.test(
          variance_residuals ~ variance_group
        )$p.value,
        error = function(e) {
          NA_real_
        }
      )
    }

    result
  }

  # ---------------------------------------------------------------------------
  # 5a. Recover predictors
  # ---------------------------------------------------------------------------
  if (is.null(predictors)) {
    model_data <- if (
      is_bundle &&
        !is.null(model_bundle$analysis_data)
    ) {
      model_bundle$analysis_data
    } else {
      tryCatch(
        insight::get_data(fitted_model),
        error = function(e) {
          NULL
        }
      )
    }

    fixed_predictor_names <- tryCatch(
      {
        unique(
          unlist(
            insight::find_predictors(
              fitted_model,
              effects = "fixed",
              component = "conditional"
            )
          )
        )
      },
      error = function(e) {
        character(0)
      }
    )

    if (
      !is.null(model_data) &&
        length(fixed_predictor_names) > 0L
    ) {
      model_data_df <- as.data.frame(
        model_data
      )

      fixed_predictor_names <- fixed_predictor_names[
        fixed_predictor_names %in%
          names(model_data_df)
      ]

      if (length(fixed_predictor_names) > 0L) {
        predictors <- model_data_df[,
          fixed_predictor_names,
          drop = FALSE
        ]
      }
    }
  } else {
    predictors <- as.data.frame(
      predictors
    )
  }

  if (!is.null(predictors)) {
    predictors <- as.data.frame(
      predictors
    )

    # glmmTMB residual methods may not always return a directly comparable
    # vector, so the definitive length check also occurs within each run.
    if (
      nrow(predictors) == 0L ||
        ncol(predictors) == 0L
    ) {
      predictors <- NULL
    }
  }

  # ---------------------------------------------------------------------------
  # 5b. Stability-result fingerprint and cache evaluation
  # ---------------------------------------------------------------------------

  parent_model_fingerprint <- if (
    is_bundle &&
      !is.null(model_bundle$analysis_fingerprint)
  ) {
    as.character(
      model_bundle$analysis_fingerprint
    )[1L]
  } else {
    NULL
  }

  predictors_for_fingerprint <- if (!is.null(predictors)) {
    as.data.frame(
      predictors,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }

  if (
    !is.null(cache_dir) &&
      is.null(parent_model_fingerprint)
  ) {
    stop(
      paste0(
        "Caching requires a model bundle containing ",
        "'analysis_fingerprint'. Supply an object returned by ",
        "fit_and_simulate(), or set 'cache_dir = NULL'."
      ),
      call. = FALSE
    )
  }

  stability_fingerprint_input <- list(
    wrapper_version = get_wrapper_version(),
    parent_model_fingerprint = parent_model_fingerprint,
    model_class = class(fitted_model),
    predictors = predictors_for_fingerprint,
    n_runs = n_runs,
    n_sim = n_sim,
    tweedie_simulation_max_power = tweedie_simulation_max_power,
    alpha = alpha,
    borderline = borderline,
    seed = if (is.null(seed)) {
      NULL
    } else {
      as.integer(seed)
    },
    min_category_n = min_category_n,
    p_adjust_method = p_adjust_method,
    dharma_version = as.character(
      utils::packageVersion("DHARMa")
    )
  )

  stability_fingerprint <- digest::digest(
    stability_fingerprint_input,
    algo = "xxhash64",
    serialize = TRUE
  )

  cache_file <- NULL

  if (!is.null(cache_dir)) {
    cache_dir_full <- normalizePath(
      cache_dir,
      mustWork = FALSE
    )

    if (!dir.exists(cache_dir_full)) {
      dir.create(
        cache_dir_full,
        recursive = TRUE
      )
    }

    if (is.null(cache_name)) {
      inferred_model_name <- if (
        is_bundle &&
          !is.null(model_bundle$model_name)
      ) {
        as.character(
          model_bundle$model_name
        )[1L]
      } else {
        "model"
      }

      cache_name <- paste0(
        inferred_model_name,
        "_DHARMa_stability"
      )
    }

    sanitized_cache_name <- gsub(
      pattern = "[^A-Za-z0-9_.-]",
      replacement = "_",
      x = cache_name
    )

    cache_file <- file.path(
      cache_dir_full,
      paste0(
        sanitized_cache_name,
        ".rds"
      )
    )

    if (
      !isTRUE(force_rerun) &&
        file.exists(cache_file)
    ) {
      cached_result <- tryCatch(
        readRDS(cache_file),
        error = function(e) {
          warning(
            paste0(
              "The cached DHARMa stability result could not be read and ",
              "will be regenerated: ",
              conditionMessage(e)
            ),
            call. = FALSE
          )

          NULL
        }
      )

      cached_fingerprint <- if (
        is.list(cached_result) &&
          !is.null(cached_result$analysis_fingerprint)
      ) {
        cached_result$analysis_fingerprint
      } else {
        NULL
      }

      if (
        !is.null(cached_result) &&
          identical(
            cached_fingerprint,
            stability_fingerprint
          )
      ) {
        if (isTRUE(verbose)) {
          message(
            paste0(
              "Loading fingerprint-verified cached DHARMa stability ",
              "result: ",
              cache_file
            )
          )

          cat(
            "\n--- Cached DHARMa stability summary ---\n"
          )

          print(
            cached_result$summary
          )
        }

        return(
          invisible(cached_result)
        )
      }

      if (isTRUE(verbose)) {
        message(
          paste0(
            "DHARMa stability cache was invalidated because the model, ",
            "predictors, settings, or package version changed."
          )
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 6. Repeated DHARMa simulations
  # ---------------------------------------------------------------------------
  all_results <- vector(
    "list",
    n_runs
  )

  all_excluded_levels <- vector(
    "list",
    n_runs
  )

  observed_has_zero <- response_has_zero(
    model = fitted_model
  )

  for (run_index in seq_len(n_runs)) {
    if (isTRUE(verbose)) {
      message(
        sprintf(
          "DHARMa stability run %d / %d",
          run_index,
          n_runs
        )
      )
    }

    simulation_output <- tryCatch(
      DHARMa::simulateResiduals(
        fittedModel = fitted_model,
        n = n_sim,
        seed = dharma_seeds[run_index]
      ),
      error = function(e) {
        NULL
      }
    )

    if (is.null(simulation_output)) {
      all_results[[run_index]] <-
        data.table::data.table(
          run = run_index,
          diagnostic = "simulation_failed",
          variable = NA_character_,
          p.value = NA_real_,
          significant = NA,
          borderline = NA
        )

      next
    }

    run_result_list <- list()
    run_excluded_list <- list()

    uniformity_test <- tryCatch(
      DHARMa::testUniformity(
        simulation_output,
        plot = FALSE
      ),
      error = function(e) {
        NULL
      }
    )

    dispersion_test <- tryCatch(
      DHARMa::testDispersion(
        simulation_output,
        alternative = "two.sided",
        plot = FALSE
      ),
      error = function(e) {
        NULL
      }
    )

    outlier_test <- tryCatch(
      DHARMa::testOutliers(
        simulation_output,
        nBoot = min(
          100L,
          n_sim
        ),
        plot = FALSE
      ),
      error = function(e) {
        NULL
      }
    )

    zero_inflation_test <- if (isTRUE(observed_has_zero)) {
      tryCatch(
        DHARMa::testZeroInflation(
          simulation_output,
          plot = FALSE
        ),
        error = function(e) {
          NULL
        }
      )
    } else {
      NULL
    }

    run_result_list <- append_result(
      run_result_list,
      run = run_index,
      diagnostic = "uniformity",
      p_value = extract_p_value(
        uniformity_test
      )
    )

    run_result_list <- append_result(
      run_result_list,
      run = run_index,
      diagnostic = "dispersion",
      p_value = extract_p_value(
        dispersion_test
      )
    )

    run_result_list <- append_result(
      run_result_list,
      run = run_index,
      diagnostic = "outliers",
      p_value = extract_p_value(
        outlier_test
      )
    )

    if (isTRUE(observed_has_zero)) {
      run_result_list <- append_result(
        run_result_list,
        run = run_index,
        diagnostic = "zero_inflation",
        p_value = extract_p_value(
          zero_inflation_test
        )
      )
    }

    fitted_predictions <- tryCatch(
      simulation_output$fittedPredictedResponse,
      error = function(e) {
        NULL
      }
    )

    if (
      is.null(fitted_predictions) ||
        length(fitted_predictions) != length(simulation_output$scaledResiduals)
    ) {
      fitted_predictions <- tryCatch(
        stats::predict(
          fitted_model,
          type = "response"
        ),
        error = function(e) {
          NULL
        }
      )
    }

    predicted_quantile_test <- if (
      !is.null(fitted_predictions) &&
        length(fitted_predictions) == length(simulation_output$scaledResiduals)
    ) {
      tryCatch(
        DHARMa::testQuantiles(
          simulation_output,
          predictor = fitted_predictions,
          plot = FALSE
        ),
        error = function(e) {
          NULL
        }
      )
    } else {
      NULL
    }

    run_result_list <- append_result(
      run_result_list,
      run = run_index,
      diagnostic = "predicted_quantile",
      variable = "model_prediction",
      p_value = extract_p_value(
        predicted_quantile_test
      )
    )

    # -------------------------------------------------------------------------
    # Predictor-specific diagnostics
    # -------------------------------------------------------------------------
    if (
      !is.null(predictors) &&
        ncol(predictors) > 0L &&
        nrow(predictors) == length(simulation_output$scaledResiduals)
    ) {
      for (predictor_name in names(predictors)) {
        predictor_values <- predictors[[predictor_name]]

        if (
          is.factor(predictor_values) ||
            is.character(predictor_values)
        ) {
          categorical_result <- categorical_tests(
            scaled_residuals = simulation_output$scaledResiduals,
            group = predictor_values,
            variable_name = predictor_name,
            run_number = run_index
          )

          run_result_list <- append_result(
            run_result_list,
            run = run_index,
            diagnostic = "predictor_categorical_uniformity",
            variable = predictor_name,
            p_value = categorical_result$uniformity_p
          )

          run_result_list <- append_result(
            run_result_list,
            run = run_index,
            diagnostic = "predictor_categorical_variance",
            variable = predictor_name,
            p_value = categorical_result$variance_p
          )

          if (nrow(categorical_result$excluded) > 0L) {
            next_excluded_index <- length(run_excluded_list) + 1L
            run_excluded_list[[next_excluded_index]] <-
              categorical_result$excluded
          }
        } else {
          numeric_test <- tryCatch(
            DHARMa::testQuantiles(
              simulation_output,
              predictor = predictor_values,
              plot = FALSE
            ),
            error = function(e) {
              NULL
            }
          )

          run_result_list <- append_result(
            run_result_list,
            run = run_index,
            diagnostic = "predictor_quantile",
            variable = predictor_name,
            p_value = extract_p_value(
              numeric_test
            )
          )
        }
      }
    }

    run_results <- data.table::rbindlist(
      run_result_list,
      fill = TRUE,
      use.names = TRUE
    )

    significant_values <-
      is.finite(run_results[["p.value"]]) &
      run_results[["p.value"]] < alpha

    borderline_values <-
      is.finite(run_results[["p.value"]]) &
      run_results[["p.value"]] >= alpha &
      run_results[["p.value"]] < borderline

    data.table::set(
      run_results,
      j = "significant",
      value = significant_values
    )

    data.table::set(
      run_results,
      j = "borderline",
      value = borderline_values
    )

    all_results[[run_index]] <- run_results

    if (length(run_excluded_list) > 0L) {
      all_excluded_levels[[run_index]] <-
        data.table::rbindlist(
          run_excluded_list,
          fill = TRUE,
          use.names = TRUE
        )
    }
  }

  # ---------------------------------------------------------------------------
  # 7. Combine and summarize results
  # ---------------------------------------------------------------------------
  results_table <- data.table::rbindlist(
    all_results,
    fill = TRUE,
    use.names = TRUE
  )

  excluded_levels_table <- data.table::rbindlist(
    all_excluded_levels,
    fill = TRUE,
    use.names = TRUE
  )

  if (nrow(results_table) == 0L) {
    stop(
      "No DHARMa stability results were generated.",
      call. = FALSE
    )
  }

  results_data_frame <- as.data.frame(
    results_table,
    stringsAsFactors = FALSE
  )

  # Use an explicit placeholder so global tests with variable = NA are retained
  # rather than dropped by split().
  summary_group_key <- paste(
    results_data_frame[["diagnostic"]],
    ifelse(
      is.na(results_data_frame[["variable"]]),
      "<GLOBAL>",
      results_data_frame[["variable"]]
    ),
    sep = "\r"
  )

  summary_groups <- split(
    results_data_frame,
    summary_group_key,
    drop = TRUE
  )

  summary_rows <- lapply(
    summary_groups,
    function(group_data) {
      valid_p_values <- group_data[["p.value"]][
        is.finite(group_data[["p.value"]])
      ]

      valid_significant <- group_data[["significant"]][
        !is.na(group_data[["significant"]])
      ]

      valid_borderline <- group_data[["borderline"]][
        !is.na(group_data[["borderline"]])
      ]

      data.frame(
        diagnostic = as.character(
          group_data[["diagnostic"]][1L]
        ),
        variable = if (all(is.na(group_data[["variable"]]))) {
          NA_character_
        } else {
          as.character(
            group_data[["variable"]][
              which(!is.na(group_data[["variable"]]))[1L]
            ]
          )
        },
        n_runs_requested = n_runs,
        n_rows = nrow(group_data),
        n_valid_p = length(valid_p_values),
        prop_test_available = length(valid_p_values) / n_runs,
        min_p = if (length(valid_p_values) > 0L) {
          min(valid_p_values)
        } else {
          NA_real_
        },
        q10_p = if (length(valid_p_values) > 0L) {
          as.numeric(
            stats::quantile(
              valid_p_values,
              probs = 0.10,
              names = FALSE,
              type = 8
            )
          )
        } else {
          NA_real_
        },
        median_p = if (length(valid_p_values) > 0L) {
          stats::median(valid_p_values)
        } else {
          NA_real_
        },
        q90_p = if (length(valid_p_values) > 0L) {
          as.numeric(
            stats::quantile(
              valid_p_values,
              probs = 0.90,
              names = FALSE,
              type = 8
            )
          )
        } else {
          NA_real_
        },
        prop_significant = if (length(valid_significant) > 0L) {
          mean(valid_significant)
        } else {
          NA_real_
        },
        prop_borderline = if (length(valid_borderline) > 0L) {
          mean(valid_borderline)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  )

  summary_table <- data.table::rbindlist(
    summary_rows,
    fill = TRUE,
    use.names = TRUE
  )

  data.table::setorderv(
    summary_table,
    cols = c(
      "prop_significant",
      "prop_borderline",
      "min_p"
    ),
    order = c(
      -1L,
      -1L,
      1L
    ),
    na.last = TRUE
  )

  # ---------------------------------------------------------------------------
  # 8. Report results
  # ---------------------------------------------------------------------------
  if (isTRUE(verbose)) {
    cat(
      "\n--- DHARMa stability summary ---\n"
    )

    print(summary_table)

    if (nrow(excluded_levels_table) > 0L) {
      cat(
        "\nCategorical levels excluded from within-level ",
        "uniformity testing:\n",
        sep = ""
      )

      excluded_levels_to_print <- unique(
        as.data.frame(
          excluded_levels_table[,
            c(
              "variable",
              "level",
              "n",
              "reason"
            ),
            with = FALSE
          ]
        )
      )

      print(
        excluded_levels_to_print,
        row.names = FALSE
      )
    }

    simulation_failures <- sum(
      results_table[["diagnostic"]] == "simulation_failed",
      na.rm = TRUE
    )

    if (simulation_failures > 0L) {
      warning(
        paste0(
          simulation_failures,
          " of ",
          n_runs,
          " DHARMa simulations failed."
        ),
        call. = FALSE
      )
    }
  }

  stability_result <- list(
    results = results_table,
    summary = summary_table,
    excluded_levels = excluded_levels_table,
    seeds = dharma_seeds,
    predictors_used = if (!is.null(predictors)) {
      names(predictors)
    } else {
      character(0)
    },
    settings = list(
      n_runs = n_runs,
      n_sim = n_sim,
      tweedie_simulation_max_power = tweedie_simulation_max_power,
      alpha = alpha,
      borderline = borderline,
      seed = seed,
      min_category_n = min_category_n,
      p_adjust_method = p_adjust_method
    ),
    model_status = model_status,
    parent_model_fingerprint = parent_model_fingerprint,
    analysis_fingerprint = stability_fingerprint,
    cache_file = cache_file,
    dharma_version = as.character(
      utils::packageVersion("DHARMa")
    ),
    generated_at = Sys.time()
  )

  if (!is.null(cache_file)) {
    saveRDS(
      stability_result,
      file = cache_file
    )

    if (isTRUE(verbose)) {
      message(
        paste0(
          "Saved DHARMa stability result to: ",
          cache_file
        )
      )
    }
  }

  invisible(stability_result)
}
