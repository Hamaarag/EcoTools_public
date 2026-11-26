# ==============================================================================
# GLLVM cached fitting wrapper
# ==============================================================================

#' Fit and Cache a GLLVM Model
#'
#' Fits a Generalized Linear Latent Variable Model with stable cross-session
#' cache validation, sparse-data screening, convergence checks, and optional
#' diagnostic plotting.
#'
#' @param y Response matrix or data frame, generally species counts.
#' @param X Environmental predictor matrix or data frame.
#' @param studyDesign Study-design data frame used by random effects.
#' @param formula Fixed-effect formula.
#' @param TR Optional trait data.
#' @param model_name Cache/model identifier.
#' @param row.eff Row-effect specification.
#' @param family GLLVM response family, either one family shared by all
#'   responses or a vector with one family per response.
#' @param num.lv Number of latent variables.
#' @param sd.errors Whether to calculate standard errors.
#' @param method Approximation method.
#' @param control Optimizer-control list.
#' @param control.start Starting-value control list.
#' @param cache_dir Cache directory.
#' @param force_refit Whether to ignore a matching cache and refit.
#' @param sparse_action
#' Action when sparse responses are detected. `"proceed"` (the default) records and reports
#' the warning but continues fitting, `"abort"` returns an aborted model object,
#' and `"prompt"` asks for interactive confirmation. The non-interactive
#' `"proceed"` policy is the default.
#' @param abort_on Character vector identifying which sparse-warning categories
#'   activate the configured `sparse_action`. Supported categories are `"all"`,
#'   `"occupancy"`, `"events_per_parameter"`, and `"total_events"`.
#' @param num_cores
#' Optional positive integer specifying the number of TMB/OpenMP threads
#' requested within each GLLVM fit. This argument does not create parallel R
#' workers and does not parallelize the multiple starting values requested
#' through `control.start$n.init`. When `NULL`, the value is selected by
#' [get_auto_cores()]. Actual parallel execution requires an OpenMP-enabled
#' installation of TMB and gllvm, and performance gains depend on model size,
#' approximation method, hardware, and memory bandwidth.
#' @param plot_diag Whether to draw GLLVM diagnostic plots.
#' @param diagnostic_check_status
#' Logical. Passed to [plot_gllvm_diagnostics()].
#' If `TRUE` (default), diagnostic plots are suppressed for model statuses
#' considered unsuitable for residual diagnostics.
#' If `FALSE`, plotting is attempted regardless of model status.
#' @param diagnostic_debug
#' Logical. Passed to [plot_gllvm_diagnostics()].
#' If `TRUE`, diagnostic plotting emits debugging messages describing
#' residual extraction, returned object structure, and related diagnostics.
#' @param gradient_threshold The gradient threshold for trusted model convergence
#' @param coef_threshold Absolute fixed-coefficient warning threshold.
#' @param sparse_screen Whether to run sparse-data screening.
#' @param sparse_min_events_per_parameter Sparse-screen threshold.
#' @param sparse_min_total_events Sparse-screen total-count threshold.
#' @param sparse_min_nonzero_rows Sparse-screen occupied-row threshold.
#' @param sparse_event_basis Sparse-screen basis.
#' @param sparse_count_like_only Restrict screen to count-like responses.
#' @param design_check
#' Logical. If `TRUE`, construct and validate the fixed-effect design matrix
#' before cache evaluation and model fitting.
#'
#' @param design_action
#' Action taken when unused factor levels, all-zero design columns,
#' non-finite values, or rank deficiency are detected. One of:
#'
#' \itemize{
#'   \item `"error"` to stop before fitting;
#'   \item `"warn"` to warn and continue with the unchanged predictor data;
#'   \item `"drop"` to remove unused factor levels from a local copy of `X`,
#'   rebuild the design matrix, and stop if structural problems remain.
#' }
#'
#' The default resolves to `"error"`.
#'
#' @param design_tolerance
#' Non-negative numeric tolerance used by [check_model_matrix()] to identify
#' effectively all-zero design columns and determine numerical matrix rank.
#'
#' @param reference_stability_check
#' Character string controlling automatic reference-level stability
#' assessment for factor predictors. Determines whether the function fits multiple models
#' (one for each eligble factor reference level) and compares their numerical diagnostics. See
#' Details.
#'
#' One of:
#'
#' \itemize{
#'   \item `"none"` Never perform reference-level stability checks.
#'   \item `"candidate"` Perform checks only for models that remain
#'   candidates for interpretation or model selection.
#'   \item `"always"` Perform checks regardless of model status.
#' }
#'
#' The default is `"candidate"`.
#'
#' @param reference_stability_max_levels
#' Integer controlling the maximum number of levels a factor may contain
#' to remain eligible for reference level stability assessment. Defaults to 5.
#'
#' @param reference_logLik_tolerance
#' Non-negative relative tolerance used when comparing maximized
#' log-likelihoods across equivalent factor-reference parameterizations.
#' The effective tolerance is this value multiplied by the larger of one and
#' the largest absolute log-likelihood. Defaults to `1e-6`.
#'
#' @param reference_stability_seed
#' Integer seed used for each alternative-reference fit so that equivalent
#' parameterizations receive reproducible stochastic starting values.
#' Use `NULL` to avoid resetting the seed. The caller's random-number state is
#' restored after the reference-level check. Defaults to `1L`.
#'
#' @param ... Additional arguments passed to gllvm::gllvm().
#'
#' @details
#' Fixed-effect design validation occurs before fingerprint calculation.
#' Therefore, when `design_action = "drop"`, both the cache fingerprint and the
#' fitted model use the cleaned local copy of `X`.
#'
#' Automatic level removal can change factor contrasts, coefficient names, and
#' reference-level interpretation. For reproducible inferential analyses,
#' `design_action = "error"` is recommended so that data cleaning remains
#' explicit in the analysis script.
#'
#' Reference-level stability assessment
#'
#' When fixed-effect predictors contain factors, fit_gllvm_cached()
#' may internally evaluate alternative reference-level parameterizations
#' of those factors as a numerical-stability diagnostic.
#'
#' Ecologically equivalent models should produce identical fitted values,
#' likelihoods, and biological conclusions regardless of factor reference
#' level. Large changes in optimizer convergence, gradient magnitude, or
#' likelihood, or other fitting success across reference levels may
#' indicate parameterization-sensitive optimization and potential
#' identifiability problems.
#'
#' Detected instability is reported through the
#' `reference_stability` component of the returned bundle and may trigger
#' the model status
#' `"suspicious_reference_sensitivity"`.
#'
#' Reference-level stability checks are performed only after a model
#' survives major numerical validation and remains a candidate for
#' inference or model selection, unless
#' `reference_stability_check = "always"`.
#'
#' Sparse-data warnings alone do not prevent reference-level assessment.
#'
#' @seealso
#' [check_model_matrix()],
#' [screen_sparse_gllvm()],
#' [assert_valid_model()],
#' [print_convergence_report()]
#'
#' @return A fitted GLLVM object or a standardized aborted/invalid list.
#' @export
fit_gllvm_cached <- function(
  y,
  X,
  studyDesign,
  formula,
  TR = NULL,
  model_name,
  row.eff = FALSE,
  family = "negative.binomial",
  num.lv = 2,
  sd.errors = TRUE,
  method = "VA",
  control = list(maxit = 10000, optimizer = "nlminb", max.iter = 10000),
  control.start = list(n.init = 20, jitter.var = 0.25, n.init.max = 20),
  cache_dir = "Analysis/",
  force_refit = FALSE,
  sparse_action = c("proceed", "prompt", "abort"),
  abort_on = c("all", "occupancy", "events_per_parameter", "total_events"),
  num_cores = NULL,
  plot_diag = TRUE,
  diagnostic_check_status = TRUE,
  diagnostic_debug = FALSE,
  gradient_threshold = 0.01,
  coef_threshold = 5,
  sparse_screen = TRUE,
  sparse_min_events_per_parameter = 5,
  sparse_min_total_events = 5,
  sparse_min_nonzero_rows = 3,
  sparse_event_basis = c("nonzero", "total", "both"),
  sparse_count_like_only = TRUE,
  design_check = TRUE,
  design_action = c("error", "warn", "drop"),
  design_tolerance = sqrt(.Machine$double.eps),
  reference_stability_check = c("candidate", "none", "always"),
  reference_stability_max_levels = 5L,
  reference_logLik_tolerance = 1e-6,
  reference_stability_seed = 1L,
  ...
) {
  # Check required packages:
  if (!requireNamespace("gllvm", quietly = TRUE)) {
    stop("Package 'gllvm' is required.", call. = FALSE)
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required.", call. = FALSE)
  }

  response_count <- tryCatch(
    ncol(
      as.data.frame(
        y,
        stringsAsFactors = FALSE
      )
    ),
    error = function(e) {
      NA_integer_
    }
  )

  if (
    length(response_count) != 1L ||
      is.na(response_count) ||
      response_count < 1L
  ) {
    stop(
      "'y' must be coercible to a response table with at least one column.",
      call. = FALSE
    )
  }

  response_family <- normalize_gllvm_response_family(
    family = family,
    n_responses = response_count
  )

  # Match arguments:
  sparse_event_basis <- match.arg(sparse_event_basis)
  sparse_action <- match.arg(sparse_action)
  abort_on <- match.arg(abort_on, several.ok = TRUE)
  design_action <- match.arg(design_action)
  reference_stability_check <- match.arg(
    reference_stability_check
  )

  # Verify that the maximum factor levels to check stability for as an integer > 1:
  if (
    !is.numeric(reference_stability_max_levels) ||
      length(reference_stability_max_levels) != 1L ||
      is.na(reference_stability_max_levels) ||
      !is.finite(reference_stability_max_levels) ||
      reference_stability_max_levels < 2L ||
      reference_stability_max_levels !=
        as.integer(reference_stability_max_levels)
  ) {
    stop(
      "'reference_stability_max_levels' must be one integer ",
      "greater than or equal to 2.",
      call. = FALSE
    )
  }

  reference_stability_max_levels <- as.integer(
    reference_stability_max_levels
  )

  if (
    !is.numeric(reference_logLik_tolerance) ||
      length(reference_logLik_tolerance) != 1L ||
      is.na(reference_logLik_tolerance) ||
      !is.finite(reference_logLik_tolerance) ||
      reference_logLik_tolerance < 0
  ) {
    stop(
      paste0(
        "'reference_logLik_tolerance' must be one non-negative ",
        "finite number."
      ),
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

  extra_args <- list(...)

  # Decide the number of TMB/OpenMP threads requested within each GLLVM fit.
  gllvm_tmb_threads <- get_auto_cores(num_cores)

  # Sanitize the model name for caching:
  sanitized_model_name <- stringr::str_replace_all(
    model_name,
    "[^a-zA-Z0-9_.-]",
    "_"
  )

  cache_dir_full <- normalizePath(cache_dir, mustWork = FALSE)
  if (!dir.exists(cache_dir_full)) {
    dir.create(cache_dir_full, recursive = TRUE)
  }
  file_path <- file.path(cache_dir_full, paste0(sanitized_model_name, ".rds"))

  design_check_result <- NULL

  if (isTRUE(design_check)) {
    design_check_result <- check_model_matrix(
      formula = formula,
      data = X,
      action = design_action,
      tolerance = design_tolerance
    )

    X <- design_check_result$data
  }

  response_matrix <- tryCatch(
    as.matrix(
      y
    ),
    error = function(e) {
      stop(
        "'y' could not be converted to a numeric response matrix.",
        call. = FALSE
      )
    }
  )

  if (!is.numeric(response_matrix)) {
    stop(
      "'y' must be coercible to a numeric response matrix.",
      call. = FALSE
    )
  }

  finite_response_counts <- rowSums(
    is.finite(
      response_matrix
    )
  )

  nonzero_response_counts <- rowSums(
    is.finite(
      response_matrix
    ) &
      response_matrix != 0
  )

  all_missing_response_rows <- which(
    finite_response_counts == 0L
  )

  all_zero_response_rows <- which(
    finite_response_counts > 0L &
      nonzero_response_counts == 0L
  )

  if (length(all_missing_response_rows) > 0L) {
    stop(
      paste0(
        length(all_missing_response_rows),
        " response row(s) contain no finite observations: ",
        paste(
          utils::head(
            all_missing_response_rows,
            10L
          ),
          collapse = ", "
        ),
        if (length(all_missing_response_rows) > 10L) {
          paste0(
            ", and ",
            length(all_missing_response_rows) - 10L,
            " additional row(s)"
          )
        } else {
          ""
        },
        ". Model fitting was not started."
      ),
      call. = FALSE
    )
  }

  response_screen_warning <- NULL

  if (length(all_zero_response_rows) > 0L) {
    response_screen_warning <- paste0(
      "GLLVM RESPONSE WARNING: ",
      length(all_zero_response_rows),
      " response row(s) contain only zeros: ",
      paste(
        utils::head(
          all_zero_response_rows,
          10L
        ),
        collapse = ", "
      ),
      if (length(all_zero_response_rows) > 10L) {
        paste0(
          ", and ",
          length(all_zero_response_rows) - 10L,
          " additional row(s)"
        )
      } else {
        ""
      },
      ". Verify that these are genuine sampled communities rather than ",
      "missing surveys, filtering artifacts, or data-join failures."
    )
  }

  response_screen_obj <- list(
    n_rows = nrow(response_matrix),
    n_responses = ncol(response_matrix),
    n_all_zero_rows = length(all_zero_response_rows),
    all_zero_rows = all_zero_response_rows,
    n_all_missing_rows = length(all_missing_response_rows),
    all_missing_rows = all_missing_response_rows
  )

  # ---------------------------------------------------------------------------
  # 1. Stable, cross-session cache fingerprint
  # ---------------------------------------------------------------------------
  dependency_versions <- list(
    gllvm = as.character(utils::packageVersion("gllvm")),
    TMB = as.character(utils::packageVersion("TMB"))
  )

  current_fingerprint <- get_gllvm_fingerprint(
    y = y,
    X = X,
    TR = TR,
    studyDesign = studyDesign,
    formula = formula,
    row.eff = row.eff,
    family = response_family,
    num.lv = num.lv,
    sd.errors = sd.errors,
    method = method,
    control = control,
    control.start = control.start,
    gradient_threshold = gradient_threshold,
    coef_threshold = coef_threshold,
    sparse_screen = sparse_screen,
    sparse_action = sparse_action,
    abort_on = sort(unique(abort_on)),
    sparse_min_events_per_parameter = sparse_min_events_per_parameter,
    sparse_min_total_events = sparse_min_total_events,
    sparse_min_nonzero_rows = sparse_min_nonzero_rows,
    sparse_event_basis = sparse_event_basis,
    sparse_count_like_only = sparse_count_like_only,
    design_check = design_check,
    design_action = design_action,
    design_tolerance = design_tolerance,
    reference_stability_check = reference_stability_check,
    reference_stability_max_levels = as.integer(reference_stability_max_levels),
    reference_logLik_tolerance = reference_logLik_tolerance,
    reference_stability_seed = reference_stability_seed,
    tmb_threads = gllvm_tmb_threads,
    extra_args = extra_args,
    dependency_versions = dependency_versions
  )

  # ---------------------------------------------------------------------------
  # 2. Cache evaluation
  # ---------------------------------------------------------------------------
  if (!force_refit && file.exists(file_path)) {
    cached_model <- tryCatch(
      readRDS(file_path),
      error = function(e) {
        warning(
          paste0(
            "Cached model could not be read and will be refitted: ",
            e$message
          ),
          call. = FALSE
        )
        NULL
      }
    )

    if (!is.null(cached_model)) {
      if (inherits(cached_model, "gllvm")) {
        class(cached_model) <- unique(
          c("ecotools_gllvm", class(cached_model))
        )
      }

      cached_response_count <- tryCatch(
        {
          cached_y <- as.matrix(
            cached_model$y
          )

          if (
            length(dim(cached_y)) == 2L &&
              ncol(cached_y) >= 1L
          ) {
            ncol(cached_y)
          } else {
            NA_integer_
          }
        },
        error = function(e) {
          NA_integer_
        }
      )

      cache_has_response_family <-
        !inherits(cached_model, "gllvm") ||
        (is.finite(cached_response_count) &&
          !is.null(cached_model$response_family) &&
          length(cached_model$response_family) == cached_response_count &&
          !anyNA(cached_model$response_family))

      cached_fingerprint <- cached_model$analysis_fingerprint

      if (
        !is.null(cached_fingerprint) &&
          identical(cached_fingerprint, current_fingerprint) &&
          isTRUE(cache_has_response_family)
      ) {
        message(stringr::str_glue(
          "Found fingerprint-verified cached GLLVM model: {file_path}"
        ))

        # Print the same status information whether the model is fresh or cached.
        print_convergence_report(
          cached_model,
          sanitized_model_name
        )

        print_sparse_gllvm_report(
          sparse_screen_obj = cached_model$sparse_screen,
          prefix = "CACHED SPARSE GLLVM SCREEN"
        )

        if (!is.null(cached_model$warning_msgs)) {
          for (w_msg in unique(cached_model$warning_msgs)) {
            warning(
              paste0("CACHED WARNING: ", w_msg),
              call. = FALSE
            )
          }
        }

        if (plot_diag) {
          plot_gllvm_diagnostics(
            model = cached_model,
            model_name = sanitized_model_name,
            debug = diagnostic_debug,
            diagnostic_check_status = diagnostic_check_status
          )
        }

        return(cached_model)
      }

      if (
        inherits(cached_model, "gllvm") &&
          !isTRUE(cache_has_response_family)
      ) {
        message(
          paste0(
            "\nCACHE INVALIDATED: Cached GLLVM object lacks valid ",
            "'response_family' metadata required by EcoTools 0.1.1. ",
            "Re-fitting the model.\n"
          )
        )
      } else if (is.null(cached_fingerprint)) {
        message(
          "\nCACHE INVALIDATED: Cached object has no fingerprint metadata. ",
          "Evaluating pre-fit checks...\n"
        )
      } else {
        message(
          "\nCACHE INVALIDATED: Model configuration or wrapper version changed. ",
          "Evaluating pre-fit checks...\n"
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 3. Pre-fit sparse-data screen
  # ---------------------------------------------------------------------------
  sparse_screen_obj <- NULL
  captured_warnings <- character()

  if (!is.null(response_screen_warning)) {
    captured_warnings <- c(
      captured_warnings,
      response_screen_warning
    )

    warning(
      response_screen_warning,
      call. = FALSE
    )
  }

  if (isTRUE(sparse_screen)) {
    sparse_screen_obj <- tryCatch(
      screen_sparse_gllvm(
        y = y,
        X = X,
        formula = formula,
        min_events_per_parameter = sparse_min_events_per_parameter,
        min_total_events = sparse_min_total_events,
        min_nonzero_rows = sparse_min_nonzero_rows,
        event_basis = sparse_event_basis,
        count_like_only = sparse_count_like_only
      ),
      error = function(e) {
        list(
          pass = NA,
          warning = paste(
            "SPARSE GLLVM SCREEN ERROR:",
            conditionMessage(e)
          ),
          warning_types = "screen_error",
          details = NULL,
          sparse_details = NULL,
          sparse_species = character(0),
          displayed_sparse_species = character(0),
          n_sparse = NA_integer_
        )
      }
    )

    if (
      length(sparse_screen_obj$pass) == 1L &&
        is.na(sparse_screen_obj$pass)
    ) {
      stop(
        paste0(
          "The GLLVM sparse-data screen failed because of a programming error.\n",
          sparse_screen_obj$warning,
          "\nModel fitting was not started."
        ),
        call. = FALSE
      )
    }

    if (!is.null(sparse_screen_obj$warning)) {
      captured_warnings <- c(
        captured_warnings,
        sparse_screen_obj$warning
      )

      warning(
        paste0(
          "\n=======================================================\n",
          "*** PRE-FIT MODELING RISK DETECTED ***\n",
          sparse_screen_obj$warning,
          "\n=======================================================\n"
        ),
        immediate. = TRUE,
        call. = FALSE
      )

      print_sparse_gllvm_report(
        sparse_screen_obj = sparse_screen_obj,
        prefix = "PRE-FIT SPARSE GLLVM SCREEN"
      )

      warning_types <- sparse_screen_obj$warning_types

      if (is.null(warning_types)) {
        warning_types <- character(0)
      }

      should_abort_type <-
        "all" %in% abort_on || any(warning_types %in% abort_on)

      if (sparse_action == "abort" && should_abort_type) {
        abort_message <- paste0(
          "Aborted automatically due to sparse-data screen. ",
          sparse_screen_obj$warning
        )

        aborted_model <- list(
          model_name = sanitized_model_name,
          response_family = response_family,
          model_status = "aborted",
          valid_for_inference = FALSE,
          valid_for_model_selection = FALSE,
          opt_converged = FALSE,
          se_valid = FALSE,
          diagnostics = list(
            design_matrix = design_check_result,
            response_screen = response_screen_obj
          ),
          warning_msgs = unique(c(
            captured_warnings,
            abort_message
          )),
          sparse_screen = sparse_screen_obj,
          parallel_config = list(
            strategy = "within_fit_tmb_openmp",
            tmb_threads_requested = gllvm_tmb_threads,
            independent_start_workers = 1L
          ),
          analysis_fingerprint = current_fingerprint
        )

        saveRDS(aborted_model, file = file_path)
        message(abort_message)
        return(aborted_model)
      }

      if (sparse_action == "prompt" && should_abort_type) {
        if (!is_rendering_session() && interactive()) {
          ans <- readline(
            prompt = "Do you want to proceed with fitting anyway? (y/n): "
          )

          if (tolower(trimws(ans)) != "y") {
            abort_message <- paste0(
              "Model fitting aborted by user. ",
              sparse_screen_obj$warning
            )

            aborted_model <- list(
              model_name = sanitized_model_name,
              response_family = response_family,
              model_status = "aborted",
              valid_for_inference = FALSE,
              valid_for_model_selection = FALSE,
              opt_converged = FALSE,
              se_valid = FALSE,
              diagnostics = list(
                design_matrix = design_check_result,
                response_screen = response_screen_obj
              ),
              warning_msgs = unique(c(
                captured_warnings,
                abort_message
              )),
              sparse_screen = sparse_screen_obj,
              parallel_config = list(
                strategy = "within_fit_tmb_openmp",
                tmb_threads_requested = gllvm_tmb_threads,
                independent_start_workers = 1L
              ),
              analysis_fingerprint = current_fingerprint
            )

            saveRDS(aborted_model, file = file_path)
            message(abort_message)
            return(aborted_model)
          }
        } else {
          noninteractive_msg <- paste0(
            "Rendering/non-interactive session detected. ",
            "Auto-proceeding because sparse_action = 'prompt' cannot request input."
          )

          captured_warnings <- c(
            captured_warnings,
            noninteractive_msg
          )

          warning(
            noninteractive_msg,
            call. = FALSE
          )
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 4. Fit the full candidate GLLVM
  #
  # This is the only computationally expensive step in the workflow.
  # All validation, sparse screening, and cache checks occur before fitting.
  # ---------------------------------------------------------------------------
  message(stringr::str_glue(
    "Started fitting GLLVM '{sanitized_model_name}' at {Sys.time()} ",
    "requesting {gllvm_tmb_threads} TMB/OpenMP thread(s) within the fit..."
  ))

  engine_warnings <- character()

  model <- with_gllvm_tmb_threads(
    threads = gllvm_tmb_threads,
    code = tryCatch(
      withCallingHandlers(
        {
          gllvm::gllvm(
            y = y,
            X = X,
            TR = TR,
            studyDesign = studyDesign,
            formula = formula,
            row.eff = row.eff,
            family = response_family,
            num.lv = num.lv,
            sd.errors = sd.errors,
            method = method,
            control = control,
            control.start = control.start,
            ...
          )
        },
        warning = function(w) {
          engine_warnings <<- c(
            engine_warnings,
            conditionMessage(w)
          )

          invokeRestart(
            "muffleWarning"
          )
        }
      ),
      error = function(e) {
        error_msg <- paste(
          "Critical optimizer failure:",
          conditionMessage(e)
        )

        captured_warnings <<- c(
          captured_warnings,
          error_msg
        )

        warning(
          error_msg,
          call. = FALSE
        )

        NULL
      }
    )
  )

  engine_warnings <- unique(
    trimws(engine_warnings)
  )

  engine_warnings <- engine_warnings[
    !is.na(engine_warnings) &
      nzchar(engine_warnings)
  ]

  captured_warnings <- unique(
    c(
      captured_warnings,
      engine_warnings
    )
  )

  # Standardized failed-fit object.
  if (is.null(model)) {
    invalid_model <- list(
      model_name = sanitized_model_name,
      response_family = response_family,
      model_status = "invalid",
      valid_for_inference = FALSE,
      valid_for_model_selection = FALSE,
      opt_converged = FALSE,
      se_valid = FALSE,
      diagnostics = list(
        design_matrix = design_check_result,
        response_screen = response_screen_obj
      ),
      warning_msgs = unique(captured_warnings),
      sparse_screen = sparse_screen_obj,
      parallel_config = list(
        strategy = "within_fit_tmb_openmp",
        tmb_threads_requested = gllvm_tmb_threads,
        independent_start_workers = 1L
      ),
      analysis_fingerprint = current_fingerprint
    )

    saveRDS(invalid_model, file = file_path)
    return(invalid_model)
  } else {
    model$call$formula <- formula # Inject the literal formulas back into the call object to save downstream predict() methods
    if (!isFALSE(row.eff)) model$call$row.eff <- row.eff
  }

  # ---------------------------------------------------------------------------
  # 5. Optional Poisson overdispersion warning
  # ---------------------------------------------------------------------------
  if (
    all(response_family == "poisson") &&
      !is.null(model$y) &&
      !is.null(model$fitted)
  ) {
    pearson_res_sq <- tryCatch(
      (model$y - model$fitted)^2 / model$fitted,
      error = function(e) NULL
    )

    if (!is.null(pearson_res_sq)) {
      finite_values <- pearson_res_sq[
        is.finite(pearson_res_sq)
      ]

      if (length(finite_values) > 0) {
        disp_ratio <- mean(
          finite_values,
          na.rm = TRUE
        )

        if (is.finite(disp_ratio) && disp_ratio > 1.5) {
          disp_msg <- paste0(
            "OVERDISPERSION WARNING: Poisson dispersion ratio = ",
            round(disp_ratio, 2),
            ", exceeding the screening threshold of 1.5. ",
            "Consider a negative-binomial sensitivity model."
          )

          captured_warnings <- c(
            captured_warnings,
            disp_msg
          )

          warning(
            disp_msg,
            call. = FALSE
          )
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 6. Numerical validation
  #
  # These checks evaluate:
  #   - optimizer convergence,
  #   - gradient magnitude,
  #   - Hessian and covariance validity,
  #   - standard-error reliability.
  #
  # Together they determine whether the fitted model is suitable for
  # inference and/or model selection.
  # ---------------------------------------------------------------------------
  # 6a. Optimizer convergence
  # ---------------------------------------------------------------------------
  # Accept either the logical value TRUE or the conventional numeric optimizer
  # return code 0 as evidence that the optimizer reported convergence.
  opt_converged <-
    isTRUE(model$convergence) ||
    (is.numeric(model$convergence) &&
      length(model$convergence) == 1L &&
      is.finite(model$convergence) &&
      model$convergence == 0)

  # ---------------------------------------------------------------------------
  # 6b. Gradient assessment
  # ---------------------------------------------------------------------------
  # The optimizer's convergence flag is not sufficient by itself. At a valid
  # optimum, the objective-function gradients should also be close to zero.
  #
  # gradient_threshold is a function argument and must not be overwritten here.

  gradient_values <- tryCatch(
    {
      if (
        !is.null(model$TMBfn) &&
          is.function(model$TMBfn$gr)
      ) {
        as.numeric(model$TMBfn$gr())
      } else {
        NULL
      }
    },
    error = function(e) {
      gradient_error_msg <- paste(
        "GLLVM GRADIENT CHECK ERROR:",
        conditionMessage(e)
      )

      captured_warnings <<- c(
        captured_warnings,
        gradient_error_msg
      )

      NULL
    }
  )

  finite_gradient_values <- if (!is.null(gradient_values)) {
    gradient_values[is.finite(gradient_values)]
  } else {
    numeric(0)
  }

  gradient_available <- length(finite_gradient_values) > 0L

  max_abs_gradient <- if (gradient_available) {
    max(abs(finite_gradient_values))
  } else {
    NA_real_
  }

  gradient_ok <- if (is.finite(max_abs_gradient)) {
    max_abs_gradient <= gradient_threshold
  } else {
    NA
  }

  gradient_severity <- if (!is.finite(max_abs_gradient)) {
    "unassessed"
  } else if (max_abs_gradient <= gradient_threshold) {
    "pass"
  } else if (max_abs_gradient <= 0.05) {
    "mild"
  } else if (max_abs_gradient <= 0.5) {
    "substantial"
  } else {
    "severe"
  }
  # ---------------------------------------------------------------------------
  # 6c. Standard-error availability and finiteness
  # ---------------------------------------------------------------------------
  se_requested <- isTRUE(sd.errors)

  se_object_available <-
    !is.null(model$sd) &&
    is.list(model$sd)

  xcoef_se_available <-
    se_object_available &&
    !is.null(model$sd$Xcoef)

  xcoef_se_total <- if (xcoef_se_available) {
    length(model$sd$Xcoef)
  } else {
    0L
  }

  xcoef_se_nonfinite <- if (xcoef_se_available) {
    sum(!is.finite(model$sd$Xcoef))
  } else {
    NA_integer_
  }

  xcoef_se_finite <-
    xcoef_se_available &&
    identical(xcoef_se_nonfinite, 0L)

  # ---------------------------------------------------------------------------
  # 6d. Inspect every numeric standard-error component
  # ---------------------------------------------------------------------------
  # Checking only model$sd$Xcoef is not sufficient. Some models can have finite
  # Xcoef standard errors, but gllvm could report invalid variance estimates elsewhere
  # in the variance-covariance calculation.
  numeric_se_components <- if (se_object_available) {
    names(model$sd)[vapply(
      model$sd,
      function(x) {
        is.numeric(x) || is.matrix(x) || is.array(x)
      },
      logical(1)
    )]
  } else {
    character(0)
  }

  se_component_summary <- if (length(numeric_se_components) > 0L) {
    as.data.frame(
      data.table::rbindlist(
        lapply(
          numeric_se_components,
          function(component_name) {
            component_values <-
              model$sd[[component_name]]

            component_numeric <-
              as.numeric(component_values)

            finite_values <- component_numeric[
              is.finite(component_numeric)
            ]

            data.table::data.table(
              component = component_name,
              n_values = length(component_numeric),
              n_finite = sum(
                is.finite(component_numeric)
              ),
              n_nonfinite = sum(
                !is.finite(component_numeric)
              ),
              n_zero = sum(
                is.finite(component_numeric) &
                  component_numeric == 0
              ),
              min_value = if (length(finite_values) > 0L) {
                min(finite_values)
              } else {
                NA_real_
              },
              max_value = if (length(finite_values) > 0L) {
                max(finite_values)
              } else {
                NA_real_
              }
            )
          }
        ),
        fill = TRUE,
        use.names = TRUE
      ),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      component = character(0),
      n_values = integer(0),
      n_finite = integer(0),
      n_nonfinite = integer(0),
      n_zero = integer(0),
      min_value = numeric(0),
      max_value = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  required_se_summary_columns <- c(
    "component",
    "n_values",
    "n_finite",
    "n_nonfinite",
    "n_zero",
    "min_value",
    "max_value"
  )

  missing_se_summary_columns <- setdiff(
    required_se_summary_columns,
    names(se_component_summary)
  )

  if (length(missing_se_summary_columns) > 0L) {
    stop(
      paste0(
        "Internal EcoTools error: standard-error component summary is ",
        "missing required column(s): ",
        paste(
          missing_se_summary_columns,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  any_nonfinite_se <- if (nrow(se_component_summary) > 0L) {
    any(se_component_summary$n_nonfinite > 0L)
  } else {
    FALSE
  }

  # ---------------------------------------------------------------------------
  # 6e. Parse variance-covariance and Hessian warnings from gllvm
  # ---------------------------------------------------------------------------
  # Warnings raised during the gllvm fit were captured in captured_warnings.
  # These warnings may identify covariance failures even when model$sd$Xcoef
  # contains only finite numbers.

  negative_variance_warning <- any(
    grepl(
      "negative variance estimate",
      captured_warnings,
      ignore.case = TRUE
    )
  )

  singular_covariance_warning <- any(
    grepl(
      paste0(
        "determinant.*variance[- ]covariance.*(is )?zero|",
        "variance[- ]covariance.*determinant.*zero|",
        "variance[- ]covariance.*singular|",
        "singular.*variance[- ]covariance|",
        "covariance.*not positive definite|",
        "non[- ]positive definite.*covariance|",
        "hessian.*singular|",
        "singular.*hessian"
      ),
      captured_warnings,
      ignore.case = TRUE
    )
  )

  hessian_warning <- any(
    grepl(
      "hessian",
      captured_warnings,
      ignore.case = TRUE
    )
  )

  # Try to extract the number reported in warnings such as:
  # "11 parameter(s) have negative variance estimates ..."
  negative_variance_count <- NA_integer_

  negative_variance_messages <- captured_warnings[
    grepl(
      "negative variance estimate",
      captured_warnings,
      ignore.case = TRUE
    )
  ]

  if (length(negative_variance_messages) > 0L) {
    count_match <- regexec(
      "([0-9]+)\\s+parameter\\(s\\).*negative variance estimate",
      negative_variance_messages[1],
      ignore.case = TRUE
    )

    count_parts <- regmatches(
      negative_variance_messages[1],
      count_match
    )[[1]]

    if (length(count_parts) >= 2L) {
      negative_variance_count <- suppressWarnings(
        as.integer(count_parts[2])
      )
    }
  }

  covariance_warning <-
    negative_variance_warning ||
    singular_covariance_warning ||
    hessian_warning

  # ---------------------------------------------------------------------------
  # 6f. Determine whether standard-error-based inference is valid
  # ---------------------------------------------------------------------------
  # Standard errors are considered valid only when:
  #   1. standard errors were requested;
  #   2. the optimizer reported convergence;
  #   3. the SE object and Xcoef SE matrix exist;
  #   4. Xcoef SEs are finite;
  #   5. every other numeric SE component is finite;
  #   6. gllvm did not report negative variances, singular covariance, or Hessian
  #      problems.
  #
  # Gradient validity is handled separately so the final status can distinguish
  # a covariance failure from a large-gradient failure.
  model$se_valid <- if (!se_requested) {
    NA
  } else {
    opt_converged &&
      se_object_available &&
      xcoef_se_available &&
      xcoef_se_finite &&
      !any_nonfinite_se &&
      !covariance_warning
  }

  # ---------------------------------------------------------------------------
  # 6g. Build detailed failure explanations
  # ---------------------------------------------------------------------------
  inference_failure_reasons <- character(0)

  if (!opt_converged) {
    inference_failure_reasons <- c(
      inference_failure_reasons,
      "optimizer did not report convergence"
    )
  }

  if (se_requested && !se_object_available) {
    inference_failure_reasons <- c(
      inference_failure_reasons,
      "standard-error object is missing"
    )
  }

  if (se_requested && !xcoef_se_available) {
    inference_failure_reasons <- c(
      inference_failure_reasons,
      "fixed-effect standard errors are missing"
    )
  }

  if (
    se_requested &&
      xcoef_se_available &&
      is.finite(xcoef_se_nonfinite) &&
      xcoef_se_nonfinite > 0L
  ) {
    inference_failure_reasons <- c(
      inference_failure_reasons,
      paste0(
        xcoef_se_nonfinite,
        " of ",
        xcoef_se_total,
        " fixed-effect standard errors are non-finite"
      )
    )
  }

  if (any_nonfinite_se) {
    affected_components <- as.character(
      se_component_summary$component[
        se_component_summary$n_nonfinite > 0L
      ]
    )

    affected_components <- unique(
      affected_components[
        !is.na(affected_components) &
          nzchar(affected_components)
      ]
    )

    inference_failure_reasons <- c(
      inference_failure_reasons,
      paste0(
        "non-finite standard errors occur in component(s): ",
        paste(
          affected_components,
          collapse = ", "
        )
      )
    )
  }

  if (negative_variance_warning) {
    inference_failure_reasons <- c(
      inference_failure_reasons,
      if (is.na(negative_variance_count)) {
        "gllvm reported negative variance estimates"
      } else {
        paste0(
          "gllvm reported ",
          negative_variance_count,
          " negative variance estimates"
        )
      }
    )
  }

  if (singular_covariance_warning) {
    inference_failure_reasons <- c(
      inference_failure_reasons,
      paste0(
        "gllvm reported a singular or zero-determinant ",
        "variance-covariance matrix"
      )
    )
  }

  if (
    hessian_warning &&
      !singular_covariance_warning
  ) {
    inference_failure_reasons <- c(
      inference_failure_reasons,
      "gllvm reported a Hessian warning"
    )
  }

  inference_failure_reasons <- unique(
    inference_failure_reasons
  )

  # ---------------------------------------------------------------------------
  # 6h. Store detailed inference-quality diagnostics
  # ---------------------------------------------------------------------------
  # Section 7 will add coef_screen to this same diagnostics list.
  if (
    is.null(model$diagnostics) ||
      !is.list(model$diagnostics)
  ) {
    model$diagnostics <- list()
  }

  model$diagnostics$design_matrix <- design_check_result

  model$diagnostics$inference_quality <- list(
    # optimizer convergence:
    opt_converged = opt_converged,
    # gradient:
    gradient_available = gradient_available,
    max_abs_gradient = max_abs_gradient,
    gradient_threshold = gradient_threshold,
    gradient_ok = gradient_ok,
    # coeffiecient standard errors:
    se_requested = se_requested,
    se_object_available = se_object_available,
    xcoef_se_available = xcoef_se_available,
    xcoef_se_total = xcoef_se_total,
    xcoef_se_nonfinite = xcoef_se_nonfinite,
    xcoef_se_finite = xcoef_se_finite,
    se_component_summary = se_component_summary,
    any_nonfinite_se = any_nonfinite_se,
    se_valid = model$se_valid,
    # variance, covariance and Hessian:
    negative_variance_warning = negative_variance_warning,
    negative_variance_count = negative_variance_count,
    singular_covariance_warning = singular_covariance_warning,
    hessian_warning = hessian_warning,
    covariance_warning = covariance_warning,

    failure_reasons = inference_failure_reasons
  )

  # gradient validity:
  model$gradient_values <- gradient_values
  model$gradient_assessed <- gradient_available
  model$gradient_valid <- gradient_ok
  model$max_abs_gradient <- max_abs_gradient
  model$gradient_threshold <- gradient_threshold
  model$gradient_severity <- gradient_severity

  # ---------------------------------------------------------------------------
  # 6i. Add concise derived warnings to the cached warning log
  # ---------------------------------------------------------------------------
  if (
    se_requested &&
      !isTRUE(model$se_valid)
  ) {
    inference_failure_msg <- if (length(inference_failure_reasons) > 0L) {
      paste0(
        "GLLVM INFERENCE VALIDATION FAILED: ",
        paste(
          inference_failure_reasons,
          collapse = "; "
        ),
        "."
      )
    } else {
      paste0(
        "GLLVM INFERENCE VALIDATION FAILED for an unclassified ",
        "standard-error or variance-covariance reason."
      )
    }

    captured_warnings <- c(
      captured_warnings,
      inference_failure_msg
    )

    warning(
      inference_failure_msg,
      call. = FALSE
    )
  }

  if (identical(gradient_ok, FALSE)) {
    gradient_msg <- paste0(
      "GLLVM GRADIENT WARNING: Maximum absolute gradient = ",
      signif(max_abs_gradient, 6),
      ", exceeding the screening threshold of ",
      gradient_threshold,
      "."
    )

    # Avoid duplicate warnings if another check generated the same message.
    if (!gradient_msg %in% captured_warnings) {
      captured_warnings <- c(
        captured_warnings,
        gradient_msg
      )

      warning(
        gradient_msg,
        call. = FALSE
      )
    }
  }
  # ---------------------------------------------------------------------------
  # 7. Fixed-effect identifiability screening
  #
  # Large-magnitude coefficients often signal:
  #   - quasi-separation,
  #   - sparse responses,
  #   - weak identifiability,
  #   - unstable extrapolation.
  #
  # These checks do not invalidate the model automatically but
  # downgrade the reported model status.
  # ---------------------------------------------------------------------------
  has_extreme_coef <- FALSE
  extreme_coef_values <- numeric(0)
  prop_extreme <- 0

  xcoef <- NULL

  if (!is.null(model$params$Xcoef)) {
    xcoef <- model$params$Xcoef
  } else {
    coef_try <- tryCatch(
      coef(model),
      error = function(e) NULL
    )

    if (is.list(coef_try) && !is.null(coef_try$Xcoef)) {
      xcoef <- coef_try$Xcoef
    }
  }

  if (!is.null(xcoef)) {
    coef_numeric <- as.numeric(as.matrix(xcoef))
    finite_coefs <- coef_numeric[is.finite(coef_numeric)]

    if (length(finite_coefs) > 0) {
      extreme_coef_values <- finite_coefs[
        abs(finite_coefs) > coef_threshold
      ]

      if (length(extreme_coef_values) > 0) {
        has_extreme_coef <- TRUE
        prop_extreme <- length(extreme_coef_values) /
          length(finite_coefs)

        extreme_msg <- stringr::str_glue(
          "EXTREME GLLVM COEFFICIENTS: ",
          "{length(extreme_coef_values)} out of {length(finite_coefs)} ",
          "({round(prop_extreme * 100, 1)}%) fixed-effect coefficients ",
          "exceed |beta| > {coef_threshold}."
        )

        captured_warnings <- c(
          captured_warnings,
          extreme_msg
        )

        warning(
          extreme_msg,
          call. = FALSE
        )
      }
    }
  }

  # Store diagnostic metadata without assuming model$diagnostics already exists.
  if (
    is.null(model$diagnostics) ||
      !is.list(model$diagnostics)
  ) {
    model$diagnostics <- list()
  }

  model$diagnostics$coef_screen <- list(
    n_extreme = length(extreme_coef_values),
    prop_extreme = prop_extreme,
    coef_threshold = coef_threshold,
    values = extreme_coef_values
  )

  # ---------------------------------------------------------------------------
  # 8. Reference-level stability assessment
  # ---------------------------------------------------------------------------
  #
  # Ecologically equivalent parameterizations should yield equivalent
  # optimization results. This section evaluates whether changing the
  # reference level of factor predictors causes materially different
  # convergence diagnostics.
  # ---------------------------------------------------------------------------

  reference_stability_results <- NULL
  reference_sensitive <- FALSE
  reference_check_failed <- FALSE

  factor_predictors <- detect_formula_factors(
    X = X,
    formula = formula
  )

  candidate_for_reference_check <-
    isTRUE(opt_converged) &&
    isTRUE(model$gradient_valid) &&
    (!isTRUE(sd.errors) ||
      isTRUE(model$se_valid))

  run_reference_check <-
    identical(reference_stability_check, "always") ||
    (identical(reference_stability_check, "candidate") &&
      candidate_for_reference_check)

  reference_skip_reason <- if (identical(reference_stability_check, "none")) {
    "disabled"
  } else if (length(factor_predictors) == 0L) {
    "no_factor_predictors"
  } else if (
    identical(reference_stability_check, "candidate") &&
      !isTRUE(opt_converged)
  ) {
    "optimizer_not_converged"
  } else if (
    identical(reference_stability_check, "candidate") &&
      !isTRUE(model$gradient_assessed)
  ) {
    "gradient_not_assessed"
  } else if (
    identical(reference_stability_check, "candidate") &&
      !isTRUE(model$gradient_valid)
  ) {
    "gradient_failed"
  } else if (
    identical(reference_stability_check, "candidate") &&
      isTRUE(sd.errors) &&
      !isTRUE(model$se_valid)
  ) {
    "standard_errors_invalid"
  } else {
    NULL
  }

  if (
    isTRUE(run_reference_check) &&
      length(factor_predictors) > 0L
  ) {
    reference_stability_results <- check_reference_stability_gllvm(
      y = y,
      X = X,
      factor_names = factor_predictors,
      TR = TR,
      studyDesign = studyDesign,
      formula = formula,
      row.eff = row.eff,
      family = response_family,
      num.lv = num.lv,
      method = method,
      control = control,
      control.start = control.start,
      tmb_threads = gllvm_tmb_threads,
      gradient_threshold = gradient_threshold,
      reference_stability_max_levels = reference_stability_max_levels,
      reference_logLik_tolerance = reference_logLik_tolerance,
      reference_stability_seed = reference_stability_seed,
      ...
    )

    reference_check_failed <-
      !is.null(reference_stability_results) &&
      !isTRUE(
        reference_stability_results$completed
      )

    reference_sensitive <- isTRUE(
      reference_stability_results$reference_sensitive
    )

    if (
      isTRUE(
        reference_stability_results$completed_with_warnings
      )
    ) {
      reference_warning_summary <- paste0(
        "GLLVM REFERENCE-STABILITY WARNING: ",
        length(
          reference_stability_results$warning_messages
        ),
        " unique warning message(s) occurred during ",
        "alternative-reference fits."
      )

      if (
        isTRUE(
          reference_stability_results$warning_classification$vector_recycling_in_betaH_model
        )
      ) {
        reference_warning_summary <- paste0(
          reference_warning_summary,
          " gllvm generated a vector-recycling warning while ",
          "fitting a model containing betaH responses."
        )
      }

      if (
        isTRUE(
          reference_stability_results$warning_classification$all_zero_response_rows
        )
      ) {
        reference_warning_summary <- paste0(
          reference_warning_summary,
          " The response matrix contains one or more all-zero rows."
        )
      }

      captured_warnings <- c(
        captured_warnings,
        reference_warning_summary
      )

      warning(
        reference_warning_summary,
        call. = FALSE
      )
    }
  }

  reference_check_warned <-
    !is.null(reference_stability_results) &&
    isTRUE(
      reference_stability_results$completed_with_warnings
    )

  reference_stability_metadata <- list(
    requested_policy = reference_stability_check,
    candidate_eligible = candidate_for_reference_check,
    check_run = !is.null(reference_stability_results),
    check_completed = if (is.null(reference_stability_results)) {
      NA
    } else {
      reference_stability_results$completed
    },
    completed_with_warnings = reference_check_warned,
    completed_with_structural_warnings = if (
      is.null(reference_stability_results)
    ) {
      FALSE
    } else {
      isTRUE(
        reference_stability_results$completed_with_structural_warnings
      )
    },
    completed_with_data_warnings = if (is.null(reference_stability_results)) {
      FALSE
    } else {
      isTRUE(
        reference_stability_results$completed_with_data_warnings
      )
    },
    skip_reason = reference_skip_reason,
    detected_factor_predictors = factor_predictors,
    result = reference_stability_results,
    reference_sensitive = reference_sensitive
  )

  model$reference_sensitive <- reference_sensitive
  model$diagnostics$reference_stability <-
    reference_stability_metadata
  model$reference_stability <-
    reference_stability_metadata

  recycling_warning_pattern <- paste0(
    "longer object length is not a multiple ",
    "of shorter object length"
  )

  primary_vector_recycling_warning <- any(
    grepl(
      recycling_warning_pattern,
      engine_warnings,
      fixed = TRUE
    )
  )

  primary_all_zero_response_warning <- any(
    grepl(
      "There are rows full of zeros in y",
      engine_warnings,
      fixed = TRUE
    )
  )

  model$diagnostics$response_screen <-
    response_screen_obj

  model$diagnostics$engine_warnings <- list(
    vector_recycling = primary_vector_recycling_warning,
    vector_recycling_in_betaH_model = primary_vector_recycling_warning &&
      any(response_family == "betaH"),
    all_zero_response_rows = primary_all_zero_response_warning,
    messages = engine_warnings
  )

  model$diagnostics$wrapper_warnings <- list(
    messages = setdiff(
      unique(
        trimws(
          captured_warnings
        )
      ),
      engine_warnings
    )
  )

  # ---------------------------------------------------------------------------
  # 9. Unified model status classification
  #
  # All diagnostics are combined into a single status code that is used by:
  #   - downstream assertions,
  #   - cache loading messages,
  #   - automated reporting.
  # ---------------------------------------------------------------------------

  model$sd_errors_requested <- isTRUE(sd.errors)
  model$opt_converged <- opt_converged

  model$model_status <- if (!isTRUE(opt_converged)) {
    "invalid"
  } else if (!isTRUE(model$gradient_assessed)) {
    "suspicious_optimization_unassessed"
  } else if (!isTRUE(model$gradient_valid)) {
    "suspicious_optimization_incomplete"
  } else if (
    isTRUE(sd.errors) &&
      !isTRUE(model$se_valid)
  ) {
    "suspicious_inference_unreliable"
  } else if (isTRUE(reference_check_failed)) {
    "suspicious_reference_check_incomplete"
  } else if (isTRUE(reference_sensitive)) {
    "suspicious_reference_sensitivity"
  } else if (isTRUE(has_extreme_coef)) {
    "suspicious_extreme_coefficients"
  } else if (!isTRUE(sd.errors)) {
    "valid_fit_no_standard_errors"
  } else {
    "valid"
  }

  model$valid_for_inference <- identical(
    model$model_status,
    "valid"
  )

  model$valid_for_model_selection <- model$model_status %in%
    c(
      "valid",
      "valid_fit_no_standard_errors"
    )

  # Store model identity and validation metadata before saving.
  model$model_name <- sanitized_model_name
  model$response_family <- response_family

  names(model$response_family) <- colnames(
    as.matrix(y)
  )

  model$sparse_screen <- sparse_screen_obj

  model$sparse_species <- if (
    !is.null(sparse_screen_obj) &&
      !is.null(sparse_screen_obj$sparse_species)
  ) {
    sort(
      unique(
        as.character(sparse_screen_obj$sparse_species)
      )
    )
  } else {
    character(0)
  }

  model$warning_msgs <- if (length(captured_warnings) > 0L) {
    unique(captured_warnings)
  } else {
    NULL
  }

  model$parallel_config <- list(
    strategy = "within_fit_tmb_openmp",
    tmb_threads_requested = gllvm_tmb_threads,
    independent_start_workers = 1L,
    logical_cores_detected = tryCatch(
      parallel::detectCores(
        logical = TRUE
      ),
      error = function(e) {
        NA_integer_
      }
    ),
    physical_cores_detected = tryCatch(
      parallel::detectCores(
        logical = FALSE
      ),
      error = function(e) {
        NA_integer_
      }
    ),
    omp_num_threads_environment = Sys.getenv(
      "OMP_NUM_THREADS",
      unset = NA_character_
    ),
    openblas_num_threads_environment = Sys.getenv(
      "OPENBLAS_NUM_THREADS",
      unset = NA_character_
    ),
    mkl_num_threads_environment = Sys.getenv(
      "MKL_NUM_THREADS",
      unset = NA_character_
    )
  )

  model$analysis_fingerprint <- current_fingerprint
  class(model) <- unique(
    c("ecotools_gllvm", class(model))
  )

  print_convergence_report(
    model,
    sanitized_model_name
  )

  # ---------------------------------------------------------------------------
  # 10. Save newly fitted model
  # ---------------------------------------------------------------------------
  saveRDS(model, file = file_path)

  message(stringr::str_glue(
    "Saved GLLVM model to: {file_path}"
  ))

  # ---------------------------------------------------------------------------
  # 11. Optional diagnostics
  # ---------------------------------------------------------------------------

  if (plot_diag) {
    plot_gllvm_diagnostics(
      model = model,
      model_name = sanitized_model_name,
      debug = diagnostic_debug,
      diagnostic_check_status = diagnostic_check_status
    )
  }

  return(model)
}
