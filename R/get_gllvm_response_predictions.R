#' Extract Response-Scale Predictions from a GLLVM
#'
#' Obtains one response-scale fitted value for every response cell. Species are
#' predicted in homogeneous family blocks to support mixed-family GLLVMs,
#' including hurdle-beta and count responses.
#'
#' @param model A fitted `gllvm` object.
#' @param family A character vector containing one family name per response
#'   column.
#'
#' @return A numeric matrix with the same dimensions and column order as the
#'   model response matrix.
#'
#' @keywords internal
get_gllvm_response_predictions <- function(
  model,
  family
) {
  if (!inherits(model, "gllvm")) {
    stop(
      "'model' must inherit from class 'gllvm'.",
      call. = FALSE
    )
  }

  response_matrix <- tryCatch(
    as.matrix(
      model$y
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    is.null(response_matrix) ||
      length(
        dim(
          response_matrix
        )
      ) !=
        2L
  ) {
    stop(
      "Could not recover the fitted GLLVM response matrix.",
      call. = FALSE
    )
  }

  n_responses <- ncol(
    response_matrix
  )

  if (
    !is.character(family) ||
      length(family) != n_responses ||
      anyNA(family)
  ) {
    stop(
      paste0(
        "'family' must contain one non-missing character value per ",
        "response column."
      ),
      call. = FALSE
    )
  }

  fitted_matrix <- matrix(
    NA_real_,
    nrow = nrow(
      response_matrix
    ),
    ncol = n_responses,
    dimnames = dimnames(
      response_matrix
    )
  )

  family_groups <- split(
    seq_len(
      n_responses
    ),
    family
  )

  expand_response_parameter <- function(
    values,
    species_index,
    parameter_name
  ) {
    if (is.null(values)) {
      stop(
        paste0("Model is missing required parameter '", parameter_name, "'."),
        call. = FALSE
      )
    }

    values_dimensions <- dim(values)
    values_matrix <- as.matrix(values)
    species_count <- length(species_index)

    if (length(values_matrix) == 1L) {
      return(
        matrix(
          values_matrix[1L],
          nrow = nrow(response_matrix),
          ncol = species_count
        )
      )
    }

    if (is.null(values_dimensions)) {
      values_vector <- as.numeric(values)

      if (length(values_vector) == n_responses) {
        return(
          matrix(
            rep(values_vector[species_index], each = nrow(response_matrix)),
            nrow = nrow(response_matrix),
            ncol = species_count
          )
        )
      }

      if (length(values_vector) == nrow(response_matrix)) {
        return(
          matrix(
            rep(values_vector, times = species_count),
            nrow = nrow(response_matrix),
            ncol = species_count
          )
        )
      }

      if (length(values_vector) == nrow(response_matrix) * n_responses) {
        values_matrix <- matrix(
          values_vector,
          nrow = nrow(response_matrix),
          ncol = n_responses
        )
      } else {
        stop(
          paste0(
            "Could not align '",
            parameter_name,
            "' to the response matrix."
          ),
          call. = FALSE
        )
      }
    }

    if (
      identical(
        dim(values_matrix),
        dim(response_matrix)
      )
    ) {
      return(
        values_matrix[,
          species_index,
          drop = FALSE
        ]
      )
    }

    if (
      nrow(values_matrix) == 1L &&
        ncol(values_matrix) == n_responses
    ) {
      return(
        matrix(
          rep(values_matrix[1L, species_index], each = nrow(response_matrix)),
          nrow = nrow(response_matrix),
          ncol = species_count
        )
      )
    }

    if (
      nrow(values_matrix) == nrow(response_matrix) &&
        ncol(values_matrix) == 1L
    ) {
      return(
        matrix(
          rep(values_matrix[, 1L], times = species_count),
          nrow = nrow(response_matrix),
          ncol = species_count
        )
      )
    }

    stop(
      paste0(
        "Could not align '",
        parameter_name,
        "' to the selected response columns."
      ),
      call. = FALSE
    )
  }

  select_family_parameter <- function(values, species_index, parameter_name) {
    if (is.null(values)) {
      stop(
        paste0("Model is missing required parameter '", parameter_name, "'."),
        call. = FALSE
      )
    }

    values <- as.numeric(values)

    if (length(values) == 1L) {
      return(
        rep(values, length(species_index))
      )
    }

    if (length(values) == n_responses) {
      return(values[species_index])
    }

    if (length(values) == length(species_index)) {
      return(values)
    }

    stop(
      paste0(
        "Could not align '",
        parameter_name,
        "' to the selected response columns."
      ),
      call. = FALSE
    )
  }

  for (family_name in names(family_groups)) {
    species_index <- family_groups[[family_name]]

    family_predictions <- tryCatch(
      predict(
        model,
        type = "response",
        spp = species_index
      ),
      error = function(e) {
        stop(
          paste0(
            "Could not obtain response-scale predictions for family '",
            family_name,
            "' and response column(s) ",
            paste(
              species_index,
              collapse = ", "
            ),
            ": ",
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    response_values <- response_matrix[,
      species_index,
      drop = FALSE
    ]

    conditional_mean <- NULL
    if (identical(family_name, "orderedBeta")) {
      link_predictions <- tryCatch(
        predict(
          model,
          type = "link",
          spp = species_index
        ),
        error = function(e) {
          stop(
            paste0(
              "Could not obtain orderedBeta conditional means: ",
              conditionMessage(e)
            ),
            call. = FALSE
          )
        }
      )
      link_predictions <- as.matrix(link_predictions)

      expected_dimensions <- c(
        nrow(response_matrix),
        length(species_index)
      )
      if (!identical(dim(link_predictions), expected_dimensions)) {
        stop(
          "Unexpected orderedBeta link-prediction dimensions.",
          call. = FALSE
        )
      }

      link_names <- model$link
      if (length(link_names) == 1L) {
        link_names <- rep(link_names, n_responses)
      }
      if (length(link_names) != n_responses) {
        stop(
          "Could not align orderedBeta link functions to response columns.",
          call. = FALSE
        )
      }

      conditional_mean <- vapply(
        seq_along(species_index),
        function(index) {
          stats::make.link(
            link_names[species_index[index]]
          )$linkinv(link_predictions[, index])
        },
        numeric(nrow(response_matrix))
      )
      conditional_mean <- matrix(
        conditional_mean,
        nrow = nrow(response_matrix),
        ncol = length(species_index)
      )
    }

    n_trials <- if (
      family_name %in% c("binomial", "beta.binomial", "ZIB", "ZNIB")
    ) {
      expand_response_parameter(
        model$Ntrials,
        species_index,
        "Ntrials"
      )
    } else {
      NULL
    }

    inflation_probability <- NULL
    upper_inflation_probability <- NULL

    if (identical(family_name, "ZIB")) {
      zero_inflation <- select_family_parameter(
        model$params$phi,
        species_index,
        "params$phi"
      )
      inflation_probability <- matrix(
        rep(zero_inflation, each = nrow(response_matrix)),
        nrow = nrow(response_matrix),
        ncol = length(species_index)
      )
    }

    if (identical(family_name, "ZNIB")) {
      zero_parameter <- select_family_parameter(
        model$params$phi,
        species_index,
        "params$phi"
      )
      upper_parameter <- select_family_parameter(
        model$params$ZINB.phi,
        species_index,
        "params$ZINB.phi"
      )
      normalizer <- 1 + zero_parameter + upper_parameter

      inflation_probability <- matrix(
        rep(zero_parameter / normalizer, each = nrow(response_matrix)),
        nrow = nrow(response_matrix),
        ncol = length(species_index)
      )
      upper_inflation_probability <- matrix(
        rep(upper_parameter / normalizer, each = nrow(response_matrix)),
        nrow = nrow(response_matrix),
        ncol = length(species_index)
      )
    }

    family_predictions <- collapse_gllvm_family_predictions(
      predictions = family_predictions,
      family = family_name,
      response_values = response_values,
      conditional_mean = conditional_mean,
      n_trials = n_trials,
      inflation_probability = inflation_probability,
      upper_inflation_probability = upper_inflation_probability
    )

    fitted_matrix[,
      species_index
    ] <- family_predictions
  }

  fitted_matrix
}
