#' Check Whether Tweedie Simulation Is Computationally Safe
#'
#' Examines fitted Tweedie power parameters before simulation-based
#' diagnostics are started. Tweedie simulation can become prohibitively
#' expensive when the power approaches the Gamma boundary at 2.
#'
#' @param model A fitted model.
#' @param max_power One finite numeric simulation-safety threshold strictly
#'   between 1 and 2.
#' @param fitted_power Optional numeric vector containing already transformed
#'   Tweedie powers. This is used internally for model types such as
#'   `mgcv::twlss()`, for which the power can vary among observations.
#'
#' @return An internal named list containing the model family, fitted power
#'   range, configured threshold, blocking decision, reason code, and
#'   explanatory message.
#'
#' @keywords internal
check_tweedie_simulation_guard <- function(
  model,
  max_power = 1.99,
  fitted_power = NULL
) {
  if (
    !is.numeric(max_power) ||
      length(max_power) != 1L ||
      is.na(max_power) ||
      !is.finite(max_power) ||
      max_power <= 1 ||
      max_power >= 2
  ) {
    stop(
      paste0(
        "'max_power' must be one finite numeric value strictly ",
        "between 1 and 2."
      ),
      call. = FALSE
    )
  }

  result <- list(
    is_tweedie = FALSE,
    model_family = NA_character_,
    power = NA_real_,
    min_power = NA_real_,
    max_power = NA_real_,
    n_power_values = 0L,
    threshold = max_power,
    blocked = FALSE,
    reason = "not_tweedie",
    message = NULL
  )

  if (is.null(model)) {
    result$reason <- "missing_model"
    result$message <- paste0(
      "No fitted model was available for the Tweedie simulation-safety ",
      "check."
    )

    return(result)
  }

  family_object <- tryCatch(
    stats::family(model),
    error = function(e) NULL
  )

  family_name <- if (
    !is.null(family_object) &&
      !is.null(family_object$family) &&
      length(family_object$family) == 1L
  ) {
    as.character(family_object$family)
  } else {
    NA_character_
  }

  result$model_family <- family_name

  is_tweedie_family <- !is.na(family_name) &&
    tolower(family_name) %in%
      c(
        "tweedie",
        "tw",
        "twlss"
      )

  if (!is_tweedie_family && is.null(fitted_power)) {
    return(result)
  }

  result$is_tweedie <- TRUE
  result$reason <- "power_unavailable"

  extracted_power <- fitted_power

  if (is.null(extracted_power) && inherits(model, "glmmTMB")) {
    family_parameters <- tryCatch(
      glmmTMB::family_params(model),
      error = function(e) NULL
    )

    if (
      !is.null(family_parameters) &&
        length(family_parameters) > 0L
    ) {
      parameter_names <- names(family_parameters)

      power_index <- if (!is.null(parameter_names)) {
        grep(
          pattern = "power",
          x = parameter_names,
          ignore.case = TRUE
        )
      } else {
        integer(0)
      }

      if (
        length(power_index) == 0L &&
          length(family_parameters) == 1L
      ) {
        power_index <- 1L
      }

      if (length(power_index) > 0L) {
        extracted_power <- family_parameters[
          power_index[1L]
        ]
      }
    }
  }

  if (
    is.null(extracted_power) &&
      inherits(model, "gam") &&
      !is.null(family_object$getTheta) &&
      is.function(family_object$getTheta) &&
      !identical(tolower(family_name), "twlss")
  ) {
    extracted_power <- tryCatch(
      family_object$getTheta(TRUE),
      error = function(e) NULL
    )
  }

  if (
    is.null(extracted_power) &&
      inherits(model, "gam") &&
      !is.null(family_object$p)
  ) {
    extracted_power <- family_object$p
  }

  if (is.null(extracted_power)) {
    result$message <- paste0(
      "The fitted model uses a Tweedie family, but EcoTools could not ",
      "retrieve its fitted power parameter. Simulation was not blocked ",
      "automatically."
    )

    return(result)
  }

  extracted_power <- suppressWarnings(
    as.numeric(extracted_power)
  )

  extracted_power <- extracted_power[
    is.finite(extracted_power)
  ]

  if (length(extracted_power) == 0L) {
    result$message <- paste0(
      "The fitted model uses a Tweedie family, but all retrieved power ",
      "values were non-finite."
    )

    return(result)
  }

  result$power <- if (length(extracted_power) == 1L) {
    extracted_power
  } else {
    max(extracted_power)
  }

  result$min_power <- min(extracted_power)
  result$max_power <- max(extracted_power)
  result$n_power_values <- length(extracted_power)

  result$blocked <- result$max_power >= max_power

  if (result$blocked) {
    result$reason <- "power_at_or_above_simulation_threshold"

    result$message <- sprintf(
      paste0(
        "DHARMa simulation was skipped because the maximum fitted ",
        "Tweedie power was %.7f, at or above the EcoTools ",
        "simulation-safety threshold of %.4f. Tweedie simulation can ",
        "become prohibitively expensive as the power approaches the ",
        "Gamma boundary at 2. Consider explicitly evaluating a Gamma ",
        "model with a log link when the response is strictly positive. ",
        "The current model has not passed simulation-based validation ",
        "and is not valid for inference under the EcoTools policy."
      ),
      result$max_power,
      max_power
    )
  } else {
    result$reason <- "power_below_simulation_threshold"

    result$message <- sprintf(
      paste0(
        "The maximum fitted Tweedie power was %.7f, below the ",
        "simulation-safety threshold of %.4f."
      ),
      result$max_power,
      max_power
    )
  }

  result
}
