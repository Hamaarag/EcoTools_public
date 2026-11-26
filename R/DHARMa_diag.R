#' Run DHARMa Residual Diagnostics
#'
#' @description
#' Runs simulation-based residual diagnostics for a fitted generalized linear
#' mixed model or generalized additive model using the DHARMa package.
#'
#' The function accepts either a model bundle returned by
#' [fit_and_simulate()] or a directly supplied compatible fitted model object.
#' When a model bundle is supplied, the fitted model, DHARMa residual object,
#' and model-validation status are extracted automatically.
#'
#' The diagnostic workflow includes:
#'
#' \itemize{
#'   \item a general DHARMa residual overview plot;
#'   \item residual plots against fixed-effect predictors;
#'   \item a dispersion test;
#'   \item a zero-inflation test when the observed response contains zeros;
#'   \item optional PNG export of all generated diagnostic plots.
#' }
#'
#' Predictor names and model data are recovered using the insight package when
#' possible. An optional data object may be supplied explicitly when automatic
#' model-data reconstruction is unavailable or unreliable.
#'
#' Diagnostic execution is stopped when the fitted model is missing, the
#' log-likelihood is unavailable, the DHARMa residual object is missing, or the
#' cached model status is `"invalid"`.
#'
#' Models marked `"suspicious_extreme_coefficients"` or
#' `"suspicious_residual_diagnostics"` are not rejected by this diagnostic
#' function. Instead, a warning is emitted and diagnostic plotting continues,
#' because residual output can help identify the reason the model has not been
#' approved for formal inference.
#'
#' @param res_object
#' A model bundle returned by [fit_and_simulate()], containing at least the
#' elements `model`, `residuals`, and `model_status`; or a compatible directly
#' fitted model object. When a directly fitted model is supplied, a previously
#' created DHARMa residual object is not automatically available, so use of a
#' bundle returned by [fit_and_simulate()] is recommended.
#'
#' @param model_title
#' A character scalar used as the title of the main diagnostic plot and as the
#' basis for exported diagnostic file names.
#'
#' @param print_plots
#' Logical. If `TRUE`, diagnostic plots are drawn on the active graphics device.
#' If `FALSE`, plots are not displayed interactively, although they may still be
#' exported when `export_plots = TRUE`.
#'
#' @param export_plots
#' Logical. If `TRUE`, diagnostic plots are exported as PNG files to
#' `output_dir`.
#'
#' @param output_dir
#' A character string specifying the directory in which exported diagnostic
#' plots are saved. Relative paths are resolved using [here::here()]. The
#' directory is created recursively when it does not already exist.
#'
#' @param data
#' An optional data frame, data table, or other object coercible to a data frame
#' containing the observations used to fit the model. When `NULL`, the function
#' attempts to recover the model data using [insight::get_data()]. Supplying
#' `data` explicitly can be useful when the fitted model call cannot be
#' reconstructed in the current environment.
#'
#' @param n_sim
#' Positive integer giving the number of simulated responses used when
#' `res_object` is a directly fitted model rather than an EcoTools bundle.
#'
#' @param tweedie_simulation_max_power
#' Finite numeric threshold strictly between 1 and 2. Simulation-based
#' diagnostics are skipped when a fitted Tweedie model has an estimated power
#' at or above this threshold. The default is `1.99`.
#'
#' @param seed
#' Integer seed used when DHARMa residuals must be generated from a directly
#' fitted model.
#'
#' @param plot_offset_diagnostics
#' Logical. If `TRUE`, residuals are plotted separately against variables used
#' in model offset terms. Offset diagnostics assess whether the assumed
#' exposure correction is adequate; they are not treated as ordinary
#' fixed-effect predictor diagnostics. Defaults to `TRUE`.
#'
#' @param offset_variables
#' Optional character vector naming raw exposure variables to diagnose, such as
#' `"n_traps"` or `"trap_days"`. When `NULL`, the function attempts to recover
#' variables used inside `offset()` terms from the conditional model formula.
#' Supply this argument explicitly when the offset is stored as a transformed
#' expression or passed through a model-specific `offset` argument.
#'
#' @return
#' For Tweedie models blocked by the simulation-safety guard, the returned list
#' has `status = "skipped_tweedie_power_boundary"`,
#' `valid_for_inference = FALSE`, and contains the guard result in
#' `tweedie_simulation_guard`. Statistical test results are unavailable.
#' Invisibly returns a list with the following elements:
#'
#' \describe{
#'   \item{model_title}{
#'     Character string containing the title supplied through `model_title`.
#'     The same value is used to label diagnostic plots and construct exported
#'     file names.
#'   }
#'   \item{dispersion_test}{
#'     The result returned by [DHARMa::testDispersion()], or `NULL` if the test
#'     could not be completed.
#'   }
#'   \item{zero_inflation_test}{
#'     The result returned by [DHARMa::testZeroInflation()], or `NULL` when the
#'     response contains no observed zeros or the test could not be completed.
#'   }
#'   \item{predictors_tested}{
#'     A character vector containing the fixed-effect predictors used in the
#'     predictor-specific residual plots.
#'   }
#'    \item{offset_variables_tested}{
#'     Character vector containing exposure variables used in separate
#'     offset-specific residual plots.
#'   }
#'   \item{offset_expressions}{
#'     Character vector containing offset expressions detected in the model
#'     formula.
#'   }
#' }
#'
#' The function returns `NULL` invisibly when diagnostics are aborted because
#' the fitted model or residual object is unavailable or invalid.
#'
#' @details
#' Exported plots use sanitized versions of `model_title` and predictor names in
#' their file names. The following files may be created:
#'
#' \itemize{
#'   \item `<model-title>_DHARMa_overview.png`;
#'   \item `<model-title>_residuals_vs_<predictor>.png`;
#'   \item `<model-title>_dispersion_test.png`;
#'   \item `<model-title>_zero_inflation_test.png`, when applicable.
#' }
#'
#' The zero-inflation test is skipped when [response_has_zero()] determines that
#' the observed response contains no zero values. This avoids applying a
#' zero-inflation diagnostic to a response for which excess zeros are
#' impossible by construction.
#'
#' Simulation-based diagnostics can vary across DHARMa simulation realizations.
#' For sparse or diagnostically fragile models, use
#' [DHARMa_stability_check()] to assess whether diagnostic conclusions remain
#' stable across repeated simulations.
#'
#' @seealso
#' [fit_and_simulate()],
#' [DHARMa_stability_check()],
#' [assert_valid_model()],
#' [DHARMa::plotResiduals()],
#' [DHARMa::testDispersion()],
#' [DHARMa::testZeroInflation()]
#'
#' @examples
#' \dontrun{
#' diagnostic_results <- DHARMa_diag(
#'   res_object = fitted_model_bundle,
#'   model_title = "Poisson GLMM",
#'   print_plots = TRUE,
#'   export_plots = TRUE,
#'   output_dir = "Model Diagnostics"
#' )
#'
#' # Supply the original data explicitly when automatic data recovery fails.
#' diagnostic_results <- DHARMa_diag(
#'   res_object = fitted_model_bundle,
#'   model_title = "Poisson GLMM",
#'   print_plots = TRUE,
#'   export_plots = FALSE,
#'   data = analysis_data
#' )
#' }
#'
#' @export
DHARMa_diag <- function(
  res_object,
  data = NULL,
  n_sim = 1000L,
  tweedie_simulation_max_power = 1.99,
  seed = 123L,
  model_title = "DHARMa Diagnostics",
  print_plots = TRUE,
  export_plots = FALSE,
  output_dir = "Diagnostics",
  plot_offset_diagnostics = TRUE,
  offset_variables = NULL
) {
  cat(stringr::str_glue("\n--- Running Diagnostics for: {model_title} ---\n"))

  if (
    !is.logical(plot_offset_diagnostics) ||
      length(plot_offset_diagnostics) != 1L ||
      is.na(plot_offset_diagnostics)
  ) {
    stop(
      "'plot_offset_diagnostics' must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.null(offset_variables)) {
    if (
      !is.character(offset_variables) ||
        anyNA(offset_variables) ||
        any(!nzchar(offset_variables))
    ) {
      stop(
        paste0(
          "'offset_variables' must be NULL or a character vector ",
          "of non-empty variable names."
        ),
        call. = FALSE
      )
    }

    offset_variables <- unique(
      offset_variables
    )
  }

  model <- if (is.list(res_object) && "model" %in% names(res_object)) {
    res_object$model
  } else {
    res_object
  }

  res <- if (is.list(res_object) && "residuals" %in% names(res_object)) {
    res_object$residuals
  } else {
    NULL
  }

  tweedie_simulation_guard <- check_tweedie_simulation_guard(
    model = model,
    max_power = tweedie_simulation_max_power
  )

  stored_guard <- if (
    is.list(res_object) &&
      !is.null(res_object$tweedie_simulation_guard)
  ) {
    res_object$tweedie_simulation_guard
  } else {
    NULL
  }

  stored_guard_matches_threshold <- !is.null(stored_guard) &&
    isTRUE(stored_guard$is_tweedie) &&
    !is.null(stored_guard$threshold) &&
    length(stored_guard$threshold) == 1L &&
    is.finite(stored_guard$threshold) &&
    isTRUE(
      all.equal(
        as.numeric(stored_guard$threshold),
        as.numeric(tweedie_simulation_max_power),
        tolerance = sqrt(.Machine$double.eps)
      )
    )

  if (isTRUE(stored_guard_matches_threshold)) {
    tweedie_simulation_guard <- stored_guard
  }

  if (isTRUE(tweedie_simulation_guard$blocked)) {
    warning(
      tweedie_simulation_guard$message,
      call. = FALSE
    )

    return(
      invisible(
        list(
          model_title = model_title,
          status = "skipped_tweedie_power_boundary",
          valid_for_inference = FALSE,
          tweedie_simulation_guard = tweedie_simulation_guard,
          dispersion_test = NULL,
          dispersion_ratio = NA_real_,
          dispersion_p_value = NA_real_,
          zero_inflation_test = NULL,
          zero_inflation_ratio = NA_real_,
          zero_inflation_p_value = NA_real_,
          response_has_zero = NA,
          predictors_tested = character(0),
          offset_variables_tested = character(0),
          offset_expressions = character(0)
        )
      )
    )
  }

  if (is.null(res) && !is.null(model)) {
    if (
      !is.numeric(n_sim) ||
        length(n_sim) != 1L ||
        is.na(n_sim) ||
        !is.finite(n_sim) ||
        n_sim < 1L ||
        n_sim != as.integer(n_sim)
    ) {
      stop(
        "'n_sim' must be one positive integer.",
        call. = FALSE
      )
    }

    if (
      !is.numeric(seed) ||
        length(seed) != 1L ||
        is.na(seed) ||
        !is.finite(seed)
    ) {
      stop(
        "'seed' must be one finite numeric value.",
        call. = FALSE
      )
    }

    res <- tryCatch(
      DHARMa::simulateResiduals(
        fittedModel = model,
        n = as.integer(n_sim),
        seed = as.integer(seed)
      ),
      error = function(e) {
        warning(
          paste0(
            "Could not generate DHARMa residuals from the directly ",
            "supplied model: ",
            conditionMessage(e)
          ),
          call. = FALSE
        )
        NULL
      }
    )
  }

  model_status <- if (
    is.list(res_object) && "model_status" %in% names(res_object)
  ) {
    res_object$model_status
  } else {
    "valid"
  }

  if (is.null(model)) {
    warning(
      paste0(
        "ABORTING DHARMa DIAGNOSTICS: No fitted model was found in the ",
        "supplied object. The bundle itself may exist, but its '$model' ",
        "component is NULL."
      ),
      call. = FALSE
    )

    return(invisible(NULL))
  }

  ll_val <- tryCatch(
    suppressWarnings(stats::logLik(model)),
    error = function(e) NA
  )

  if (is.na(ll_val) || is.null(res) || model_status == "invalid") {
    warning(
      stringr::str_glue(
        "\nABORTING DHARMa DIAGNOSTICS for '{model_title}'\n"
      ),
      call. = FALSE
    )
    cat(
      "  -> Reason: Model status is INVALID, logLik is NA, or residuals are missing.\n"
    )
    cat(
      "  -> Resolve convergence / Hessian / complete-separation issues before residual diagnostics.\n\n"
    )
    return(invisible(NULL))
  }

  if (
    model_status %in%
      c(
        "suspicious_extreme_coefficients",
        "suspicious_residual_diagnostics",
        "valid_fit_no_residual_diagnostics"
      )
  ) {
    diagnostic_status_message <- switch(
      model_status,

      suspicious_extreme_coefficients = paste0(
        "The model has extreme parameter estimates. Residual diagnostics ",
        "may help identify separation or sparse-data problems."
      ),

      suspicious_residual_diagnostics = paste0(
        "The model passed numerical fitting checks, but the stored DHARMa ",
        "residual validation detected lack of fit. Diagnostics will continue ",
        "to show the relevant residual patterns."
      ),

      valid_fit_no_residual_diagnostics = paste0(
        "The model passed numerical fitting checks, but the wrapper did not ",
        "complete residual validation. Diagnostics will continue using the ",
        "available residual object."
      )
    )

    warning(
      paste0(
        "Model '",
        model_title,
        "' has status '",
        model_status,
        "'. ",
        diagnostic_status_message
      ),
      call. = FALSE
    )
  }

  safe_title <- sanitize_filename_part(model_title, max_width = 40)

  if (export_plots) {
    output_dir_full <- here::here(output_dir)
    if (!dir.exists(output_dir_full)) {
      dir.create(output_dir_full, recursive = TRUE)
    }
  }

  save_plot_png <- function(
    filename,
    plot_expr,
    width = 1600,
    height = 1200,
    res_dpi = 150
  ) {
    if (!isTRUE(export_plots)) {
      return(invisible(NULL))
    }

    file_path <- file.path(
      output_dir_full,
      filename
    )

    grDevices::png(
      filename = file_path,
      width = width,
      height = height,
      res = res_dpi
    )

    on.exit(
      grDevices::dev.off(),
      add = TRUE
    )

    eval(
      plot_expr,
      envir = parent.frame()
    )

    invisible(file_path)
  }

  draw_with_model_title <- function(
    plot_expression,
    plot_model_title
  ) {
    old_graphics_parameters <- graphics::par(
      no.readonly = TRUE
    )

    on.exit(
      graphics::par(
        old_graphics_parameters
      ),
      add = TRUE
    )

    current_outer_margin <- graphics::par(
      "oma"
    )

    current_outer_margin[3L] <- max(
      current_outer_margin[3L],
      3
    )

    graphics::par(
      oma = current_outer_margin
    )

    force(
      plot_expression
    )

    graphics::mtext(
      text = stringr::str_wrap(
        plot_model_title,
        width = 60
      ),
      side = 3,
      outer = TRUE,
      line = 1,
      font = 2,
      cex = 1
    )

    invisible(TRUE)
  }

  cat("\nPlotting main DHARMa diagnostic plot...\n")

  if (print_plots) {
    plot(res, title = stringr::str_wrap(model_title, width = 50))
  }

  save_plot_png(
    filename = paste0(safe_title, "_DHARMa_overview.png"),
    plot_expr = bquote(
      plot(
        .(res),
        title = .(stringr::str_wrap(model_title, width = 50))
      )
    )
  )

  # ---------------------------------------------------------------------------
  # Identify offset expressions and their underlying data variables.
  #
  # Offsets are exposure adjustments with coefficients fixed by the model.
  # They are diagnosed separately from estimated fixed-effect predictors.
  # ---------------------------------------------------------------------------

  conditional_formula <- tryCatch(
    stats::formula(
      model
    ),
    error = function(e) {
      NULL
    }
  )

  collect_offset_calls <- function(expression) {
    offset_calls <- list()

    recurse <- function(node) {
      if (!is.call(node)) {
        return(
          invisible(NULL)
        )
      }

      is_offset_call <- if (
        is.symbol(
          node[[1L]]
        )
      ) {
        identical(
          as.character(
            node[[1L]]
          ),
          "offset"
        )
      } else if (
        is.call(
          node[[1L]]
        ) &&
          length(
            node[[1L]]
          ) ==
            3L &&
          as.character(
            node[[1L]][[1L]]
          ) %in%
            c(
              "::",
              ":::"
            )
      ) {
        identical(
          as.character(
            node[[1L]][[3L]]
          ),
          "offset"
        )
      } else {
        FALSE
      }

      if (isTRUE(is_offset_call)) {
        offset_calls[[
          length(offset_calls) + 1L
        ]] <<- node
      }

      if (length(node) >= 2L) {
        for (element_index in seq.int(
          from = 2L,
          to = length(node)
        )) {
          recurse(
            node[[element_index]]
          )
        }
      }

      invisible(NULL)
    }

    recurse(
      expression
    )

    offset_calls
  }

  offset_calls <- if (!is.null(conditional_formula)) {
    collect_offset_calls(
      conditional_formula[[
        length(conditional_formula)
      ]]
    )
  } else {
    list()
  }

  offset_expressions <- if (length(offset_calls) > 0L) {
    vapply(
      offset_calls,
      function(offset_call) {
        paste(
          deparse(
            offset_call,
            width.cutoff = 500L
          ),
          collapse = ""
        )
      },
      character(1)
    )
  } else {
    character(0)
  }

  offset_variables_detected <- if (length(offset_calls) > 0L) {
    unique(
      unlist(
        lapply(
          offset_calls,
          function(offset_call) {
            if (length(offset_call) < 2L) {
              character(0)
            } else {
              all.vars(
                offset_call[[2L]]
              )
            }
          }
        ),
        use.names = FALSE
      )
    )
  } else {
    character(0)
  }

  offset_variables_requested <- if (is.null(offset_variables)) {
    offset_variables_detected
  } else {
    offset_variables
  }

  fixed_vars_raw <- tryCatch(
    {
      unique(insight::find_predictors(
        model,
        effects = "fixed",
        component = "conditional",
        flatten = TRUE,
        verbose = FALSE
      ))
    },
    error = function(e) character(0)
  )

  # NEW: Robust fallback using base R if insight fails to parse the wrapper's environment
  if (length(fixed_vars_raw) == 0L) {
    fixed_vars_raw <- tryCatch(
      {
        if (is.null(conditional_formula)) {
          character(0)
        } else {
          setdiff(
            all.vars(
              conditional_formula[[
                length(conditional_formula)
              ]]
            ),
            offset_variables_detected
          )
        }
      },
      error = function(e) {
        character(0)
      }
    )
  }

  # Prioritize the explicitly provided data object
  model_data <- data

  if (
    is.null(model_data) &&
      is.list(res_object) &&
      !is.null(res_object$analysis_data)
  ) {
    model_data <- res_object$analysis_data
  }

  if (is.null(model_data)) {
    model_data <- tryCatch(
      insight::get_data(model),
      error = function(e) NULL
    )
  }

  if (is.null(model_data) || length(fixed_vars_raw) == 0) {
    warning(
      "Could not extract model data/predictors with insight; skipping predictor residual plots.\n",
      call. = FALSE
    )
    fixed_vars <- character(0)
  } else {
    model_data_df <- as.data.frame(model_data)
    fixed_vars <- setdiff(
      fixed_vars_raw[fixed_vars_raw %in% names(model_data_df)],
      offset_variables_detected
    )

    poly_terms <- grep("poly\\(", fixed_vars_raw, value = TRUE)
    extract_var <- function(term) sub("poly\\(([^,]+),.*\\)", "\\1", term)
    poly_vars <- unique(sapply(poly_terms, extract_var))

    fixed_vars <- unique(c(fixed_vars, poly_vars))

    fixed_vars <- setdiff(
      fixed_vars[fixed_vars %in% names(model_data_df)],
      offset_variables_detected
    )
  }

  offset_variables_tested <- character(0)

  if (
    isTRUE(plot_offset_diagnostics) &&
      !is.null(model_data)
  ) {
    model_data_df <- as.data.frame(model_data)

    missing_offset_variables <- setdiff(
      offset_variables_requested,
      names(model_data_df)
    )

    if (length(missing_offset_variables) > 0L) {
      warning(
        paste0(
          "Offset diagnostic variable(s) were not found in the recovered ",
          "model data: ",
          paste(
            missing_offset_variables,
            collapse = ", "
          ),
          ". Supply the original analysis data and/or specify ",
          "'offset_variables' explicitly."
        ),
        call. = FALSE
      )
    }

    offset_variables_tested <- intersect(
      offset_variables_requested,
      names(model_data_df)
    )
  }

  if (length(fixed_vars) > 0) {
    for (var in fixed_vars) {
      safe_var <- sanitize_filename_part(var)

      cat(stringr::str_glue(
        "\nPlotting residuals against predictor: {var}...\n"
      ))

      pred_vec <- as.data.frame(model_data)[[var]]
      if (is.factor(pred_vec) || is.character(pred_vec)) {
        pred_vec <- droplevels(as.factor(pred_vec))
      }

      if (isTRUE(print_plots)) {
        tryCatch(
          {
            draw_with_model_title(
              plot_expression = {
                DHARMa::plotResiduals(
                  res,
                  form = pred_vec,
                  xlab = var
                )
              },
              plot_model_title = model_title
            )
          },
          error = function(e) {
            warning(
              stringr::str_glue(
                "Could not plot residuals for '{var}': {e$message}"
              ),
              call. = FALSE
            )
          }
        )
      }

      save_plot_png(
        filename = paste0(
          safe_title,
          "_residuals_vs_",
          safe_var,
          ".png"
        ),
        plot_expr = bquote(
          tryCatch(
            {
              draw_with_model_title(
                plot_expression = {
                  DHARMa::plotResiduals(
                    .(res),
                    form = .(pred_vec),
                    xlab = .(var)
                  )
                },
                plot_model_title = .(model_title)
              )
            },
            error = function(e) {
              NULL
            }
          )
        )
      )
    }
  } else {
    cat("No fixed predictors plotted.\n")
  }

  # ---------------------------------------------------------------------------
  # Offset / exposure diagnostics
  # ---------------------------------------------------------------------------

  if (
    isTRUE(plot_offset_diagnostics) &&
      length(offset_variables_tested) > 0L
  ) {
    for (offset_variable in offset_variables_tested) {
      safe_offset_variable <-
        sanitize_filename_part(
          offset_variable
        )

      exposure_vector <-
        as.data.frame(
          model_data
        )[[offset_variable]]

      if (!is.numeric(exposure_vector)) {
        warning(
          paste0(
            "Offset diagnostic variable '",
            offset_variable,
            "' is not numeric; its residual plot was skipped."
          ),
          call. = FALSE
        )

        next
      }

      if (
        anyNA(exposure_vector) ||
          any(!is.finite(exposure_vector))
      ) {
        warning(
          paste0(
            "Offset diagnostic variable '",
            offset_variable,
            "' contains missing or non-finite values; its residual ",
            "plot was skipped."
          ),
          call. = FALSE
        )

        next
      }

      cat(
        stringr::str_glue(
          "\nPlotting residuals against offset/exposure variable: ",
          "{offset_variable}...\n"
        )
      )

      if (isTRUE(print_plots)) {
        tryCatch(
          {
            draw_with_model_title(
              plot_expression = {
                DHARMa::plotResiduals(
                  res,
                  form = exposure_vector,
                  xlab = paste0(
                    offset_variable,
                    " (offset/exposure diagnostic)"
                  )
                )
              },
              plot_model_title = paste0(
                model_title,
                "\nOffset adequacy diagnostic"
              )
            )
          },
          error = function(e) {
            warning(
              paste0(
                "Could not plot residuals against offset/exposure variable '",
                offset_variable,
                "': ",
                conditionMessage(e)
              ),
              call. = FALSE
            )
          }
        )
      }

      save_plot_png(
        filename = paste0(
          safe_title,
          "_residuals_vs_offset_",
          safe_offset_variable,
          ".png"
        ),
        plot_expr = bquote(
          tryCatch(
            {
              draw_with_model_title(
                plot_expression = {
                  DHARMa::plotResiduals(
                    .(res),
                    form = .(exposure_vector),
                    xlab = .(
                      paste0(
                        offset_variable,
                        " (offset/exposure diagnostic)"
                      )
                    )
                  )
                },
                plot_model_title = .(
                  paste0(
                    model_title,
                    "\nOffset adequacy diagnostic"
                  )
                )
              )
            },
            error = function(e) {
              NULL
            }
          )
        )
      )
    }
  } else if (isTRUE(plot_offset_diagnostics)) {
    cat(
      "No offset/exposure variables were available for separate diagnostics.\n"
    )
  }

  # ---------------------------------------------------------------------------
  # Dispersion test
  # ---------------------------------------------------------------------------
  cat("\nRunning Dispersion Test...\n")

  dispersion_ratio <- NA_real_
  dispersion_p_value <- NA_real_

  disp_test <- tryCatch(
    DHARMa::testDispersion(
      res,
      alternative = "two.sided",
      plot = FALSE
    ),
    error = function(e) {
      warning(
        paste0(
          "Could not complete the DHARMa dispersion test: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )

      NULL
    }
  )

  if (!is.null(disp_test)) {
    print(disp_test)

    if (
      !is.null(disp_test$statistic) &&
        length(disp_test$statistic) > 0L
    ) {
      dispersion_ratio <- suppressWarnings(
        as.numeric(disp_test$statistic[1L])
      )
    }

    if (
      !is.null(disp_test$p.value) &&
        length(disp_test$p.value) > 0L
    ) {
      dispersion_p_value <- suppressWarnings(
        as.numeric(disp_test$p.value[1L])
      )
    }

    cat(
      "\nDHARMa dispersion summary:\n",
      "  Observed-to-simulated variance ratio: ",
      if (is.finite(dispersion_ratio)) {
        signif(dispersion_ratio, 5)
      } else {
        "not available"
      },
      "\n",
      "  P-value: ",
      if (is.finite(dispersion_p_value)) {
        format.pval(
          dispersion_p_value,
          digits = 4,
          eps = 0.001
        )
      } else {
        "not available"
      },
      "\n",
      sep = ""
    )
  }

  if (isTRUE(print_plots)) {
    tryCatch(
      {
        draw_with_model_title(
          plot_expression = {
            DHARMa::testDispersion(
              res,
              alternative = "two.sided",
              plot = TRUE
            )
          },
          plot_model_title = model_title
        )
      },
      error = function(e) {
        warning(
          paste0(
            "Could not draw the DHARMa dispersion-test plot: ",
            conditionMessage(e)
          ),
          call. = FALSE
        )

        NULL
      }
    )
  }

  save_plot_png(
    filename = paste0(
      safe_title,
      "_dispersion_test.png"
    ),
    plot_expr = bquote(
      tryCatch(
        {
          draw_with_model_title(
            plot_expression = {
              DHARMa::testDispersion(
                .(res),
                alternative = "two.sided",
                plot = TRUE
              )
            },
            plot_model_title = .(model_title)
          )
        },
        error = function(e) {
          NULL
        }
      )
    )
  )

  # ---------------------------------------------------------------------------
  # Zero-inflation test
  # ---------------------------------------------------------------------------
  zi_test <- NULL
  zero_inflation_ratio <- NA_real_
  zero_inflation_p_value <- NA_real_

  has_zero <- response_has_zero(
    model = model,
    residuals_object = res
  )

  if (isTRUE(has_zero)) {
    cat("\nRunning Zero-Inflation Test...\n")

    zi_test <- tryCatch(
      DHARMa::testZeroInflation(
        res,
        plot = FALSE
      ),
      error = function(e) {
        warning(
          paste0(
            "Could not complete the DHARMa zero-inflation test: ",
            conditionMessage(e)
          ),
          call. = FALSE
        )

        NULL
      }
    )

    if (!is.null(zi_test)) {
      print(zi_test)

      if (
        !is.null(zi_test$statistic) &&
          length(zi_test$statistic) > 0L
      ) {
        zero_inflation_ratio <- suppressWarnings(
          as.numeric(zi_test$statistic[1L])
        )
      }

      if (
        !is.null(zi_test$p.value) &&
          length(zi_test$p.value) > 0L
      ) {
        zero_inflation_p_value <- suppressWarnings(
          as.numeric(zi_test$p.value[1L])
        )
      }
    }

    if (isTRUE(print_plots)) {
      tryCatch(
        {
          draw_with_model_title(
            plot_expression = {
              DHARMa::testZeroInflation(
                res,
                plot = TRUE
              )
            },
            plot_model_title = model_title
          )
        },
        error = function(e) {
          warning(
            paste0(
              "Could not draw the DHARMa zero-inflation-test plot: ",
              conditionMessage(e)
            ),
            call. = FALSE
          )

          NULL
        }
      )
    }

    save_plot_png(
      filename = paste0(
        safe_title,
        "_zero_inflation_test.png"
      ),
      plot_expr = bquote(
        tryCatch(
          {
            draw_with_model_title(
              plot_expression = {
                DHARMa::testZeroInflation(
                  .(res),
                  plot = TRUE
                )
              },
              plot_model_title = .(model_title)
            )
          },
          error = function(e) {
            NULL
          }
        )
      )
    )
  } else {
    cat(
      "\nSkipping Zero-Inflation Test: ",
      "response contains no observed zeros.\n",
      sep = ""
    )
  }

  # ---------------------------------------------------------------------------
  # Return diagnostic results
  # ---------------------------------------------------------------------------
  invisible(
    list(
      model_title = model_title,
      dispersion_test = disp_test,
      dispersion_ratio = dispersion_ratio,
      dispersion_p_value = dispersion_p_value,
      zero_inflation_test = zi_test,
      zero_inflation_ratio = zero_inflation_ratio,
      zero_inflation_p_value = zero_inflation_p_value,
      response_has_zero = isTRUE(has_zero),
      predictors_tested = fixed_vars,
      offset_variables_tested = offset_variables_tested,
      offset_expressions = offset_expressions
    )
  )
}
