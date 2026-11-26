#' Collapse GLLVM Family Prediction Components
#'
#' Converts family-specific predictions to one expected response value per
#' response. This includes hurdle-beta and ordinal probability outputs and
#' trial-scaled or zero-inflated binomial means.
#'
#' @param predictions Prediction matrix or category-probability array returned
#'   for one family block.
#' @param family A single response family name.
#' @param response_values Numeric response matrix for the requested responses.
#' @param conditional_mean Optional matrix of conditional beta means for
#'   `orderedBeta`.
#' @param n_trials Optional matrix of binomial trial counts.
#' @param inflation_probability Optional zero-inflation probability matrix for
#'   `ZIB` and `ZNIB`.
#' @param upper_inflation_probability Optional upper-inflation probability
#'   matrix for `ZNIB`.
#'
#' @return A numeric matrix of expected response values, one column per
#'   requested response.
#'
#' @keywords internal
collapse_gllvm_family_predictions <- function(
  predictions,
  family,
  response_values,
  conditional_mean = NULL,
  n_trials = NULL,
  inflation_probability = NULL,
  upper_inflation_probability = NULL
) {
  response_values <- as.matrix(response_values)
  n_observations <- nrow(response_values)
  n_responses <- ncol(response_values)

  if (identical(family, "ordinal")) {
    if (
      length(dim(predictions)) != 3L ||
        dim(predictions)[2L] != n_observations ||
        dim(predictions)[3L] != n_responses
    ) {
      stop(
        "Ordinal predictions must be a categories-by-observations-by-responses array.",
        call. = FALSE
      )
    }

    expected_values <- matrix(
      NA_real_,
      nrow = n_observations,
      ncol = n_responses,
      dimnames = dimnames(response_values)
    )

    for (response_index in seq_len(n_responses)) {
      observed_values <- response_values[, response_index]
      observed_values <- observed_values[is.finite(observed_values)]

      if (length(observed_values) == 0L) {
        next
      }

      category_values <- seq.int(
        min(observed_values),
        max(observed_values)
      )

      if (length(category_values) > dim(predictions)[1L]) {
        stop(
          "Ordinal category probabilities do not cover the observed response range.",
          call. = FALSE
        )
      }

      category_probabilities <- matrix(
        predictions[
          seq_along(category_values),
          ,
          response_index,
          drop = FALSE
        ],
        nrow = length(category_values),
        ncol = n_observations
      )

      expected_values[, response_index] <- colSums(
        sweep(
          category_probabilities,
          MARGIN = 1L,
          STATS = category_values,
          FUN = "*"
        ),
        na.rm = TRUE
      )
    }

    return(expected_values)
  }

  if (identical(family, "orderedBeta")) {
    if (
      length(dim(predictions)) != 3L ||
        dim(predictions)[1L] != 3L ||
        dim(predictions)[2L] != n_observations ||
        dim(predictions)[3L] != n_responses
    ) {
      stop(
        "orderedBeta predictions must contain three category probabilities per response.",
        call. = FALSE
      )
    }

    if (
      is.null(conditional_mean) ||
        !identical(
          dim(as.matrix(conditional_mean)),
          c(n_observations, n_responses)
        )
    ) {
      stop(
        "orderedBeta requires one conditional beta mean per response cell.",
        call. = FALSE
      )
    }

    upper_probability <- matrix(
      predictions[3L, , , drop = FALSE],
      nrow = n_observations,
      ncol = n_responses
    )
    interior_probability <- matrix(
      predictions[2L, , , drop = FALSE],
      nrow = n_observations,
      ncol = n_responses
    )

    expected_values <- upper_probability +
      interior_probability * as.matrix(conditional_mean)
    dimnames(expected_values) <- dimnames(response_values)
    return(expected_values)
  }

  predictions <- as.matrix(predictions)
  expected_columns <- if (identical(family, "betaH")) {
    2L * n_responses
  } else {
    n_responses
  }

  expected_dimensions <- c(
    n_observations,
    expected_columns
  )

  if (!identical(dim(predictions), expected_dimensions)) {
    stop(
      paste0(
        "Unexpected prediction dimensions for family '",
        family,
        "': expected ",
        paste(
          expected_dimensions,
          collapse = " x "
        ),
        ", obtained ",
        paste(
          dim(predictions),
          collapse = " x "
        ),
        "."
      ),
      call. = FALSE
    )
  }

  if (!identical(family, "betaH")) {
    response_predictions <- predictions
  } else {
    response_columns <- seq_len(n_responses)
    occurrence_columns <- n_responses + response_columns

    response_predictions <- predictions[, response_columns, drop = FALSE] *
      predictions[, occurrence_columns, drop = FALSE]
  }

  if (family %in% c("binomial", "beta.binomial", "ZIB", "ZNIB")) {
    if (
      is.null(n_trials) ||
        !identical(
          dim(as.matrix(n_trials)),
          c(n_observations, n_responses)
        )
    ) {
      stop(
        paste0(
          "Family '",
          family,
          "' requires one trial count per response cell."
        ),
        call. = FALSE
      )
    }

    n_trials <- as.matrix(n_trials)
  }

  if (family %in% c("binomial", "beta.binomial")) {
    response_predictions <- response_predictions * n_trials
  }

  if (identical(family, "ZIB")) {
    if (
      is.null(inflation_probability) ||
        !identical(
          dim(as.matrix(inflation_probability)),
          c(n_observations, n_responses)
        )
    ) {
      stop(
        "ZIB requires one zero-inflation probability per response cell.",
        call. = FALSE
      )
    }

    response_predictions <- response_predictions *
      n_trials *
      (1 - as.matrix(inflation_probability))
  }

  if (identical(family, "ZNIB")) {
    if (
      is.null(inflation_probability) ||
        is.null(upper_inflation_probability) ||
        !identical(
          dim(as.matrix(inflation_probability)),
          c(n_observations, n_responses)
        ) ||
        !identical(
          dim(as.matrix(upper_inflation_probability)),
          c(n_observations, n_responses)
        )
    ) {
      stop(
        "ZNIB requires zero- and upper-inflation probabilities per response cell.",
        call. = FALSE
      )
    }

    response_predictions <- n_trials *
      (as.matrix(upper_inflation_probability) +
        (1 -
          as.matrix(inflation_probability) -
          as.matrix(upper_inflation_probability)) *
          response_predictions)
  }

  dimnames(response_predictions) <- dimnames(response_values)
  response_predictions
}
