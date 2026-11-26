check_reference_stability_gllvm <- function(
  y,
  X,
  factor_names,
  TR,
  studyDesign,
  formula,
  row.eff,
  family,
  num.lv,
  method,
  control,
  control.start,
  tmb_threads,
  gradient_threshold,
  reference_stability_max_levels = 5L,
  reference_logLik_tolerance = 1e-6,
  reference_stability_seed = 1L,
  ...
) {
  results <- list()
  factor_summaries <- list()
  sensitivity_flags <- logical()
  all_warning_messages <- character()

  if (
    !is.numeric(tmb_threads) ||
      length(tmb_threads) != 1L ||
      is.na(tmb_threads) ||
      !is.finite(tmb_threads) ||
      tmb_threads < 1L ||
      tmb_threads != as.integer(tmb_threads)
  ) {
    stop(
      "'tmb_threads' must be one positive integer.",
      call. = FALSE
    )
  }

  tmb_threads <- as.integer(
    tmb_threads
  )

  if (
    !is.numeric(reference_logLik_tolerance) ||
      length(reference_logLik_tolerance) != 1L ||
      is.na(reference_logLik_tolerance) ||
      !is.finite(reference_logLik_tolerance) ||
      reference_logLik_tolerance < 0
  ) {
    stop(
      "'reference_logLik_tolerance' must be one non-negative finite number.",
      call. = FALSE
    )
  }

  if (!is.null(reference_stability_seed)) {
    if (
      !is.numeric(reference_stability_seed) ||
        length(reference_stability_seed) != 1L ||
        is.na(reference_stability_seed) ||
        !is.finite(reference_stability_seed) ||
        reference_stability_seed != as.integer(reference_stability_seed)
    ) {
      stop(
        paste0(
          "'reference_stability_seed' must be NULL or one finite ",
          "integer."
        ),
        call. = FALSE
      )
    }

    reference_stability_seed <- as.integer(
      reference_stability_seed
    )
  }

  had_random_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )

  if (had_random_seed) {
    old_random_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  }

  on.exit(
    {
      if (had_random_seed) {
        assign(
          ".Random.seed",
          old_random_seed,
          envir = .GlobalEnv
        )
      } else if (
        exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )
      ) {
        remove(
          ".Random.seed",
          envir = .GlobalEnv
        )
      }
    },
    add = TRUE
  )

  all_factor_names <- factor_names

  factor_names <- factor_names[
    vapply(
      factor_names,
      function(factor_name) {
        length(
          levels(
            factor(
              X[[factor_name]]
            )
          )
        ) <=
          reference_stability_max_levels
      },
      logical(1)
    )
  ]

  excluded_factor_predictors <- setdiff(
    all_factor_names,
    factor_names
  )

  if (length(factor_names) == 0L) {
    return(
      list(
        checked = TRUE,
        completed = NA,
        completed_with_warnings = FALSE,
        completed_with_structural_warnings = FALSE,
        completed_with_data_warnings = FALSE,
        n_reference_fits = 0L,
        n_successful_reference_fits = 0L,
        n_failed_reference_fits = 0L,
        factor_predictors = factor_names,
        excluded_factor_predictors = excluded_factor_predictors,
        level_limit = reference_stability_max_levels,
        reference_logLik_tolerance = reference_logLik_tolerance,
        reference_stability_seed = reference_stability_seed,
        factor_results = results,
        factor_summaries = factor_summaries,
        warning_messages = character(0),
        warning_classification = list(
          vector_recycling = FALSE,
          vector_recycling_in_betaH_model = FALSE,
          all_zero_response_rows = FALSE,
          structural = character(0),
          data = character(0),
          unclassified = character(0)
        ),
        reference_sensitive = FALSE
      )
    )
  }

  for (factor_name in factor_names) {
    factor_levels <- levels(
      factor(
        X[[factor_name]]
      )
    )

    factor_results <- list()

    for (ref_level in factor_levels) {
      X_ref <- as.data.frame(
        X
      )

      X_ref[[factor_name]] <- stats::relevel(
        factor(
          X_ref[[factor_name]]
        ),
        ref = ref_level
      )

      reference_warning_messages <- character()
      reference_error_message <- NULL

      if (!is.null(reference_stability_seed)) {
        set.seed(
          reference_stability_seed
        )
      }

      fit_ref <- with_gllvm_tmb_threads(
        threads = tmb_threads,
        code = tryCatch(
          withCallingHandlers(
            {
              gllvm::gllvm(
                y = y,
                X = X_ref,
                TR = TR,
                studyDesign = studyDesign,
                formula = formula,
                row.eff = row.eff,
                family = family,
                num.lv = num.lv,
                sd.errors = FALSE,
                method = method,
                control = control,
                control.start = control.start,
                ...
              )
            },
            warning = function(w) {
              reference_warning_messages <<- c(
                reference_warning_messages,
                conditionMessage(w)
              )

              invokeRestart(
                "muffleWarning"
              )
            }
          ),
          error = function(e) {
            reference_error_message <<-
              conditionMessage(e)

            NULL
          }
        )
      )

      reference_warning_messages <- unique(
        trimws(
          reference_warning_messages
        )
      )

      reference_warning_messages <-
        reference_warning_messages[
          !is.na(reference_warning_messages) &
            nzchar(reference_warning_messages)
        ]

      all_warning_messages <- c(
        all_warning_messages,
        reference_warning_messages
      )

      if (is.null(fit_ref)) {
        factor_results[[ref_level]] <- list(
          fit_success = FALSE,
          optimizer_converged = FALSE,
          max_abs_gradient = NA_real_,
          gradient_valid = NA,
          logLik = NA_real_,
          warning_messages = reference_warning_messages,
          error_message = reference_error_message
        )

        next
      }

      max_grad <- tryCatch(
        {
          if (
            !is.null(fit_ref$TMBfn) &&
              is.function(fit_ref$TMBfn$gr)
          ) {
            gradient_values <- as.numeric(
              fit_ref$TMBfn$gr()
            )

            gradient_values <- gradient_values[
              is.finite(gradient_values)
            ]

            if (length(gradient_values) > 0L) {
              max(
                abs(
                  gradient_values
                )
              )
            } else {
              NA_real_
            }
          } else {
            NA_real_
          }
        },
        error = function(e) {
          NA_real_
        }
      )

      converged <-
        isTRUE(fit_ref$convergence) ||
        (is.numeric(fit_ref$convergence) &&
          length(fit_ref$convergence) == 1L &&
          is.finite(fit_ref$convergence) &&
          fit_ref$convergence == 0)

      factor_results[[ref_level]] <- list(
        fit_success = TRUE,
        optimizer_converged = converged,
        max_abs_gradient = max_grad,
        gradient_valid = if (is.finite(max_grad)) {
          max_grad <= gradient_threshold
        } else {
          NA
        },
        logLik = tryCatch(
          as.numeric(
            stats::logLik(
              fit_ref
            )
          ),
          error = function(e) {
            NA_real_
          }
        ),
        warning_messages = reference_warning_messages,
        error_message = reference_error_message
      )
    }

    results[[factor_name]] <- factor_results

    fit_success_values <- vapply(
      factor_results,
      function(x) {
        isTRUE(
          x$fit_success
        )
      },
      logical(1)
    )

    convergence_values <- vapply(
      factor_results,
      function(x) {
        isTRUE(
          x$optimizer_converged
        )
      },
      logical(1)
    )

    gradient_values <- vapply(
      factor_results,
      function(x) {
        x$max_abs_gradient
      },
      numeric(1)
    )

    logLik_values <- vapply(
      factor_results,
      function(x) {
        x$logLik
      },
      numeric(1)
    )

    fit_failure_instability <-
      any(fit_success_values) &&
      any(!fit_success_values)

    convergence_instability <-
      length(
        unique(
          convergence_values
        )
      ) >
        1L

    finite_gradient_values <- gradient_values[
      is.finite(gradient_values)
    ]

    gradient_threshold_instability <-
      length(finite_gradient_values) >= 2L &&
      any(
        finite_gradient_values <= gradient_threshold
      ) &&
      any(
        finite_gradient_values > gradient_threshold
      )

    gradient_floor <- 1e-12

    gradient_ratio <- if (length(finite_gradient_values) >= 2L) {
      max(
        finite_gradient_values
      ) /
        max(
          min(
            finite_gradient_values
          ),
          gradient_floor
        )
    } else {
      NA_real_
    }

    gradient_scale_instability <-
      is.finite(gradient_ratio) &&
      gradient_ratio > 100 &&
      max(
        finite_gradient_values
      ) >
        gradient_threshold / 10

    finite_logLik_values <- logLik_values[
      is.finite(logLik_values)
    ]

    logLik_range <- if (length(finite_logLik_values) >= 2L) {
      diff(
        range(
          finite_logLik_values
        )
      )
    } else {
      NA_real_
    }

    logLik_scale <- if (length(finite_logLik_values) > 0L) {
      max(
        1,
        max(
          abs(
            finite_logLik_values
          )
        )
      )
    } else {
      NA_real_
    }

    effective_logLik_tolerance <-
      reference_logLik_tolerance *
      logLik_scale

    logLik_instability <- if (
      is.finite(logLik_range) &&
        is.finite(
          effective_logLik_tolerance
        )
    ) {
      logLik_range > effective_logLik_tolerance
    } else {
      FALSE
    }

    factor_summaries[[factor_name]] <- list(
      fit_failure_instability = fit_failure_instability,
      convergence_instability = convergence_instability,
      gradient_threshold_instability = gradient_threshold_instability,
      gradient_scale_instability = gradient_scale_instability,
      gradient_ratio = gradient_ratio,
      logLik_instability = logLik_instability,
      logLik_range = logLik_range,
      effective_logLik_tolerance = effective_logLik_tolerance
    )

    sensitivity_flags <- c(
      sensitivity_flags,
      fit_failure_instability,
      convergence_instability,
      gradient_threshold_instability,
      gradient_scale_instability,
      logLik_instability
    )
  }

  all_warning_messages <- unique(
    trimws(
      all_warning_messages
    )
  )

  all_warning_messages <- all_warning_messages[
    !is.na(all_warning_messages) &
      nzchar(all_warning_messages)
  ]

  recycling_warning_pattern <- paste0(
    "longer object length is not a multiple ",
    "of shorter object length"
  )

  zero_row_warning_pattern <-
    "There are rows full of zeros in y"

  has_vector_recycling_warning <- any(
    grepl(
      recycling_warning_pattern,
      all_warning_messages,
      fixed = TRUE
    )
  )

  has_all_zero_row_warning <- any(
    grepl(
      zero_row_warning_pattern,
      all_warning_messages,
      fixed = TRUE
    )
  )

  has_betaH_responses <- any(
    as.character(family) == "betaH"
  )

  classified_warning_indices <-
    grepl(
      recycling_warning_pattern,
      all_warning_messages,
      fixed = TRUE
    ) |
    grepl(
      zero_row_warning_pattern,
      all_warning_messages,
      fixed = TRUE
    )

  unclassified_warning_messages <-
    all_warning_messages[
      !classified_warning_indices
    ]

  structural_warning_messages <-
    all_warning_messages[
      grepl(
        recycling_warning_pattern,
        all_warning_messages,
        fixed = TRUE
      ) |
        all_warning_messages %in%
          unclassified_warning_messages
    ]

  data_warning_messages <-
    all_warning_messages[
      grepl(
        zero_row_warning_pattern,
        all_warning_messages,
        fixed = TRUE
      )
    ]

  all_reference_results <- unlist(
    results,
    recursive = FALSE,
    use.names = TRUE
  )

  n_reference_fits <- length(
    all_reference_results
  )

  n_successful_reference_fits <- sum(
    vapply(
      all_reference_results,
      function(x) {
        isTRUE(
          x$fit_success
        )
      },
      logical(1)
    )
  )

  n_failed_reference_fits <-
    n_reference_fits -
    n_successful_reference_fits

  completed <-
    n_reference_fits > 0L &&
    n_failed_reference_fits == 0L

  list(
    checked = TRUE,
    completed = completed,
    completed_with_warnings = length(all_warning_messages) > 0L,
    completed_with_structural_warnings = length(structural_warning_messages) >
      0L,
    completed_with_data_warnings = length(data_warning_messages) > 0L,
    n_reference_fits = n_reference_fits,
    n_successful_reference_fits = n_successful_reference_fits,
    n_failed_reference_fits = n_failed_reference_fits,
    factor_predictors = factor_names,
    excluded_factor_predictors = excluded_factor_predictors,
    level_limit = reference_stability_max_levels,
    reference_logLik_tolerance = reference_logLik_tolerance,
    reference_stability_seed = reference_stability_seed,
    factor_results = results,
    factor_summaries = factor_summaries,
    warning_messages = all_warning_messages,
    warning_classification = list(
      vector_recycling = has_vector_recycling_warning,
      vector_recycling_in_betaH_model = has_vector_recycling_warning &&
        has_betaH_responses,
      all_zero_response_rows = has_all_zero_row_warning,
      structural = structural_warning_messages,
      data = data_warning_messages,
      unclassified = unclassified_warning_messages
    ),
    reference_sensitive = any(
      sensitivity_flags
    )
  )
}
