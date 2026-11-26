# ==============================================================================
# Sparse-data screening for GLMMs and GAMs
# ==============================================================================

#' Screen a Count Model for Sparse Data
#'
#' @description
#' Screens a fitted GLMM or GAM for sparse response data relative to the number
#' of estimated parameters and, when applicable, within levels of fixed-effect
#' categorical predictors.
#'
#' For count responses, the default event basis is the number of non-zero
#' observations rather than total abundance. A single observation containing
#' many individuals is therefore not treated as equivalent to many independent
#' occupied observations distributed across the model design.
#'
#' @param model
#' A fitted model object.
#'
#' @param data
#' Optional data frame or object coercible to a data frame. When `NULL`, the
#' function attempts to recover the model frame from `model`.
#'
#' @param response
#' Optional character scalar identifying the response variable. When `NULL`,
#' the response is inferred from the model frame or model formula.
#'
#' @param group_vars
#' Optional character vector identifying categorical predictors to screen.
#' When `NULL`, fixed-effect categorical predictors are inferred from the fitted
#' model.
#'
#' @param estimated_parameters
#' Optional positive numeric value giving the approximate number of estimated
#' model parameters. When `NULL`, the value is extracted from the degrees of
#' freedom associated with `logLik(model)`.
#'
#' @param min_events_per_parameter
#' Minimum number of informative events per estimated parameter.
#'
#' @param min_events_per_level
#' Minimum number of informative events within each categorical predictor level.
#'
#' @param event_basis
#' Basis used to define informative events. One of `"nonzero"`, `"total"`, or
#' `"both"`.
#'
#' @param count_like_only
#' Logical. If `TRUE`, screening is restricted to non-negative integer-like
#' responses.
#'
#' @return
#' A list containing:
#'
#' \describe{
#'   \item{pass}{
#'     Logical indicating whether all sparse-data checks passed.
#'   }
#'   \item{warning}{
#'     A warning message, or `NULL` when the screen passed.
#'   }
#'   \item{warning_types}{
#'     Character vector identifying the failed screening categories.
#'   }
#'   \item{details}{
#'     Detailed response-level and categorical-level screening results.
#'   }
#' }
#'
#' @export
screen_sparse_count_model <- function(
  model,
  data = NULL,
  response = NULL,
  group_vars = NULL,
  estimated_parameters = NULL,
  min_events_per_parameter = 5,
  min_events_per_level = 3,
  event_basis = c("nonzero", "total", "both"),
  count_like_only = TRUE
) {
  event_basis <- match.arg(event_basis)

  # ---------------------------------------------------------------------------
  # 1. Validate thresholds
  # ---------------------------------------------------------------------------
  threshold_values <- c(
    min_events_per_parameter = min_events_per_parameter,
    min_events_per_level = min_events_per_level
  )

  if (
    any(!is.finite(threshold_values)) ||
      any(threshold_values < 0)
  ) {
    stop(
      "Sparse-screen thresholds must be finite, non-negative values.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 2. Recover the model frame and response
  # ---------------------------------------------------------------------------
  model_frame <- tryCatch(
    stats::model.frame(model),
    error = function(e) NULL
  )

  supplied_data <- if (!is.null(data)) {
    as.data.frame(
      data,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }

  if (!is.null(response)) {
    response <- as.character(
      unlist(
        response,
        recursive = TRUE,
        use.names = FALSE
      )
    )

    response <- response[
      !is.na(response) &
        nzchar(response)
    ]

    if (length(response) > 1L) {
      response <- response[1L]
    }
  }

  if (is.null(response) || length(response) == 0L) {
    response <- tryCatch(
      {
        model_formula <- stats::formula(model)
        response_name <- all.vars(model_formula[[2L]])

        if (length(response_name) > 0L) {
          response_name[1L]
        } else {
          NULL
        }
      },
      error = function(e) NULL
    )
  }

  if (
    (is.null(response) || length(response) == 0L) &&
      !is.null(model_frame) &&
      ncol(model_frame) > 0L
  ) {
    response <- names(model_frame)[1L]
  }

  response <- if (
    !is.null(response) &&
      length(response) > 0L
  ) {
    as.character(response[1L])
  } else {
    NULL
  }

  # Prefer the supplied data when it contains the inferred response because it
  # has simple user-facing predictor names and avoids engine-specific model-frame
  # columns.
  if (
    !is.null(supplied_data) &&
      !is.null(response) &&
      response %in% names(supplied_data)
  ) {
    screen_data <- supplied_data
  } else if (!is.null(model_frame)) {
    screen_data <- as.data.frame(
      model_frame,
      stringsAsFactors = FALSE
    )
  } else {
    screen_data <- supplied_data
  }

  if (
    is.null(screen_data) ||
      is.null(response) ||
      !response %in% names(screen_data)
  ) {
    return(
      list(
        pass = NA,
        warning = paste0(
          "SPARSE DATA SCREEN ERROR: Could not recover the model response ",
          "and screening data."
        ),
        warning_types = "screen_error",
        details = NULL
      )
    )
  }

  y <- screen_data[[response]]

  # Matrix-valued responses require response-specific handling that this
  # univariate screen does not provide.
  if (
    is.matrix(y) ||
      is.data.frame(y) ||
      is.list(y)
  ) {
    return(
      list(
        pass = NA,
        warning = paste0(
          "SPARSE DATA SCREEN ERROR: Response '",
          response,
          "' is not a simple univariate vector."
        ),
        warning_types = "screen_error",
        details = NULL
      )
    )
  }

  y <- suppressWarnings(
    as.numeric(y)
  )

  # ---------------------------------------------------------------------------
  # 3. Check whether the response is count-like
  # ---------------------------------------------------------------------------
  if (
    isTRUE(count_like_only) &&
      !is_count_like_response(y)
  ) {
    return(
      list(
        pass = TRUE,
        warning = NULL,
        warning_types = character(0),
        details = list(
          reason = "response_not_count_like",
          response = response
        )
      )
    )
  }

  total_events <- sum(
    y,
    na.rm = TRUE
  )

  nonzero_n <- sum(
    y > 0,
    na.rm = TRUE
  )

  n_obs <- sum(
    !is.na(y)
  )

  zero_prop <- if (n_obs > 0L) {
    mean(
      y == 0,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  observed_y <- y[is.finite(y)]

  is_binary_response <- length(observed_y) > 0L &&
    all(observed_y %in% c(0, 1))

  minority_class_n <- if (is_binary_response) {
    min(
      sum(observed_y == 0),
      sum(observed_y == 1)
    )
  } else {
    NA_integer_
  }

  informative_event_n <- if (is_binary_response) {
    minority_class_n
  } else {
    nonzero_n
  }

  # ---------------------------------------------------------------------------
  # 4. Estimate the number of model parameters
  # ---------------------------------------------------------------------------
  if (is.null(estimated_parameters)) {
    estimated_parameters <- tryCatch(
      as.numeric(
        attr(
          stats::logLik(model),
          "df"
        )
      ),
      error = function(e) NA_real_
    )
  }

  estimated_parameters <- suppressWarnings(
    as.numeric(
      unlist(
        estimated_parameters,
        recursive = TRUE,
        use.names = FALSE
      )
    )[1L]
  )

  if (
    length(estimated_parameters) == 0L ||
      !is.finite(estimated_parameters) ||
      estimated_parameters <= 0
  ) {
    estimated_parameters <- NA_real_
  }

  total_events_per_parameter <- if (is.finite(estimated_parameters)) {
    total_events / estimated_parameters
  } else {
    NA_real_
  }

  informative_events_per_parameter <- if (is.finite(estimated_parameters)) {
    informative_event_n / estimated_parameters
  } else {
    NA_real_
  }

  basis_value <- switch(
    event_basis,
    nonzero = informative_events_per_parameter,
    total = total_events_per_parameter,
    both = min(
      total_events_per_parameter,
      informative_events_per_parameter,
      na.rm = TRUE
    )
  )

  if (!is.finite(basis_value)) {
    basis_value <- NA_real_
  }

  # ---------------------------------------------------------------------------
  # 5. Identify fixed-effect categorical predictors
  # ---------------------------------------------------------------------------
  if (is.null(group_vars)) {
    fixed_vars_raw <- tryCatch(
      insight::find_predictors(
        model,
        effects = "fixed",
        component = "conditional"
      ),
      error = function(e) NULL
    )

    fixed_vars <- if (!is.null(fixed_vars_raw)) {
      as.character(
        unlist(
          fixed_vars_raw,
          recursive = TRUE,
          use.names = FALSE
        )
      )
    } else {
      character(0)
    }

    fixed_vars <- unique(
      fixed_vars[
        !is.na(fixed_vars) &
          nzchar(fixed_vars)
      ]
    )

    fixed_vars <- intersect(
      fixed_vars,
      names(screen_data)
    )

    group_vars <- fixed_vars[
      vapply(
        fixed_vars,
        function(variable_name) {
          variable_values <- screen_data[[variable_name]]

          is.factor(variable_values) ||
            is.character(variable_values)
        },
        logical(1)
      )
    ]
  } else {
    group_vars <- as.character(
      unlist(
        group_vars,
        recursive = TRUE,
        use.names = FALSE
      )
    )

    group_vars <- unique(
      group_vars[
        !is.na(group_vars) &
          nzchar(group_vars)
      ]
    )

    missing_group_vars <- setdiff(
      group_vars,
      names(screen_data)
    )

    if (length(missing_group_vars) > 0L) {
      warning(
        paste0(
          "Sparse screen ignored grouping variable(s) not found in data: ",
          paste(missing_group_vars, collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }

    group_vars <- intersect(
      group_vars,
      names(screen_data)
    )
  }

  # ---------------------------------------------------------------------------
  # 6. Screen categorical predictor levels
  # ---------------------------------------------------------------------------
  group_result_list <- list()

  if (length(group_vars) > 0L) {
    for (group_variable in group_vars) {
      group_values <- screen_data[[group_variable]]

      valid_rows <- !is.na(group_values) &
        !is.na(y)

      group_values <- droplevels(
        factor(group_values[valid_rows])
      )

      response_values <- y[valid_rows]

      if (length(group_values) == 0L) {
        next
      }

      split_indices <- split(
        seq_along(group_values),
        group_values,
        drop = TRUE
      )

      variable_results <- lapply(
        names(split_indices),
        function(level_name) {
          level_indices <- split_indices[[level_name]]
          level_response <- response_values[level_indices]

          data.frame(
            grouping_variable = group_variable,
            level = level_name,
            n = length(level_response),
            total_events = sum(
              level_response,
              na.rm = TRUE
            ),
            nonzero_n = sum(
              level_response > 0,
              na.rm = TRUE
            ),
            zero_prop = if (sum(!is.na(level_response)) > 0L) {
              mean(
                level_response == 0,
                na.rm = TRUE
              )
            } else {
              NA_real_
            },
            stringsAsFactors = FALSE
          )
        }
      )

      group_result_list[[group_variable]] <- do.call(
        rbind,
        variable_results
      )
    }
  }

  group_details <- if (length(group_result_list) > 0L) {
    do.call(
      rbind,
      group_result_list
    )
  } else {
    data.frame(
      grouping_variable = character(0),
      level = character(0),
      n = integer(0),
      total_events = numeric(0),
      nonzero_n = integer(0),
      zero_prop = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  row.names(group_details) <- NULL

  if (nrow(group_details) > 0L) {
    group_details[["basis_value"]] <- switch(
      event_basis,
      nonzero = group_details[["nonzero_n"]],
      total = group_details[["total_events"]],
      both = pmin(
        group_details[["nonzero_n"]],
        group_details[["total_events"]]
      )
    )

    group_details[["fails_minimum"]] <-
      group_details[["basis_value"]] < min_events_per_level
  } else {
    group_details[["basis_value"]] <- numeric(0)
    group_details[["fails_minimum"]] <- logical(0)
  }

  group_fail <- nrow(group_details) > 0L &&
    any(
      group_details[["fails_minimum"]],
      na.rm = TRUE
    )

  failed_group_details <- group_details[
    !is.na(group_details[["fails_minimum"]]) &
      group_details[["fails_minimum"]],
    ,
    drop = FALSE
  ]

  # ---------------------------------------------------------------------------
  # 7. Apply the events-per-parameter criterion
  # ---------------------------------------------------------------------------
  epp_fail <- is.finite(basis_value) &&
    basis_value < min_events_per_parameter

  warning_types <- character(0)

  if (epp_fail) {
    warning_types <- c(
      warning_types,
      "events_per_parameter"
    )
  }

  if (group_fail) {
    warning_types <- c(
      warning_types,
      "events_per_level"
    )
  }

  warning_types <- unique(
    warning_types
  )

  # ---------------------------------------------------------------------------
  # 8. Construct the warning text
  # ---------------------------------------------------------------------------
  warning_msg <- NULL

  if (epp_fail || group_fail) {
    failed_level_text <- if (nrow(failed_group_details) > 0L) {
      paste0(
        failed_group_details[["grouping_variable"]],
        " = ",
        failed_group_details[["level"]],
        " (",
        failed_group_details[["basis_value"]],
        " events)"
      )
    } else {
      character(0)
    }

    warning_msg <- paste0(
      "SPARSE DATA WARNING: Count/binary response is sparse for ",
      "inferential modelling. Total abundance/events = ",
      total_events,
      ", non-zero observations = ",
      nonzero_n,
      "/",
      n_obs,
      ", zero proportion = ",
      round(zero_prop, 3),
      if (is.finite(total_events_per_parameter)) {
        paste0(
          ", total events per estimated parameter = ",
          round(total_events_per_parameter, 2)
        )
      } else {
        ", total events per estimated parameter = NA"
      },
      if (is.finite(informative_events_per_parameter)) {
        paste0(
          ", informative events per estimated parameter = ",
          round(informative_events_per_parameter, 2),
          if (is_binary_response) {
            " using the minority outcome class"
          } else {
            " using non-zero observations"
          }
        )
      } else {
        ", informative events per estimated parameter = NA"
      },
      " (warning threshold = ",
      min_events_per_parameter,
      ", basis = '",
      event_basis,
      "')",
      if (length(failed_level_text) > 0L) {
        paste0(
          ". Sparse categorical level(s): ",
          paste(
            failed_level_text,
            collapse = "; "
          ),
          ". Minimum events per level = ",
          min_events_per_level
        )
      } else {
        "."
      },
      " Treat fitted effects and DHARMa diagnostics cautiously; ",
      "consider descriptive summaries, pooling, penalized methods, or ",
      "presence/absence modelling."
    )
  }

  # ---------------------------------------------------------------------------
  # 9. Return results
  # ---------------------------------------------------------------------------
  list(
    pass = is.null(warning_msg),
    warning = warning_msg,
    warning_types = warning_types,
    details = list(
      response = response,
      total_events = total_events,
      nonzero_n = nonzero_n,
      n_obs = n_obs,
      zero_prop = zero_prop,
      estimated_parameters = estimated_parameters,
      total_events_per_parameter = total_events_per_parameter,
      is_binary_response = is_binary_response,
      minority_class_n = minority_class_n,
      informative_event_n = informative_event_n,
      informative_events_per_parameter = informative_events_per_parameter,
      event_basis = event_basis,
      basis_value = basis_value,
      min_events_per_parameter = min_events_per_parameter,
      min_events_per_level = min_events_per_level,
      group_vars = group_vars,
      group_details = group_details,
      failed_group_details = failed_group_details
    )
  )
}
