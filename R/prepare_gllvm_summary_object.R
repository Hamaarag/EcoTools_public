#' Prepare a GLLVM Object for Summary
#'
#' Removes internal extra beta-hurdle component rows from species-level
#' parameter arrays on a temporary copy, so `summary.gllvm()` sees one row per
#' observed response. The fitted object supplied by the caller is not modified.
#'
#' @param object A fitted `gllvm` object.
#'
#' @return A copy suitable for passing to `summary.gllvm()`.
#'
#' @keywords internal
prepare_gllvm_summary_object <- function(
  object
) {
  if (!inherits(object, "gllvm")) {
    stop(
      "'object' must inherit from class 'gllvm'.",
      call. = FALSE
    )
  }

  response_matrix <- tryCatch(
    as.matrix(
      object$y
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    is.null(response_matrix) ||
      length(dim(response_matrix)) != 2L ||
      ncol(response_matrix) < 1L
  ) {
    stop(
      "Could not recover the observed GLLVM response matrix.",
      call. = FALSE
    )
  }

  response_count <- ncol(
    response_matrix
  )

  # EcoTools stores exactly one family per observed response. This is the
  # authoritative representation for mixed-family models because gllvm may
  # expand object$family internally for multi-component families such as betaH.
  response_family <- if (!is.null(object$response_family)) {
    normalize_gllvm_response_family(
      family = object$response_family,
      n_responses = response_count
    )
  } else {
    stored_family <- as.character(
      object$family
    )

    if (length(stored_family) == 1L) {
      rep(
        stored_family,
        response_count
      )
    } else if (length(stored_family) == response_count) {
      stored_family
    } else {
      stop(
        paste0(
          "Cannot prepare the GLLVM summary because the model contains ",
          length(stored_family),
          " internal family entries for ",
          response_count,
          " observed responses and has no valid 'response_family' metadata."
        ),
        call. = FALSE
      )
    }
  }

  if (!any(response_family == "betaH")) {
    return(object)
  }

  n_beta_hurdle <- sum(
    response_family == "betaH"
  )

  expected_internal_count <-
    response_count +
    n_beta_hurdle

  stored_family <- as.character(
    object$family
  )

  if (
    length(stored_family) %in%
      c(
        response_count,
        expected_internal_count
      )
  ) {
    # Valid representations:
    #   response_count: already one family per response;
    #   expected_internal_count: gllvm-expanded betaH components.
  } else {
    stop(
      paste0(
        "Cannot prepare the GLLVM summary: the model has ",
        length(stored_family),
        " stored family entries, but expected either ",
        response_count,
        " response-level entries or ",
        expected_internal_count,
        " component-level entries for ",
        n_beta_hurdle,
        " betaH response(s)."
      ),
      call. = FALSE
    )
  }

  trim_species_vector <- function(
    values,
    parameter_name
  ) {
    if (is.null(values)) {
      return(values)
    }

    if (!is.atomic(values)) {
      return(values)
    }

    if (length(values) == expected_internal_count) {
      return(
        values[
          seq_len(response_count)
        ]
      )
    }

    if (length(values) == response_count) {
      return(values)
    }

    stop(
      paste0(
        "Cannot prepare the GLLVM summary: '",
        parameter_name,
        "' has ",
        length(values),
        " entries; expected ",
        response_count,
        " response-level entries or ",
        expected_internal_count,
        " component-level entries."
      ),
      call. = FALSE
    )
  }

  trim_species_rows <- function(
    values,
    parameter_name
  ) {
    if (is.null(values)) {
      return(values)
    }

    if (
      !is.matrix(values) &&
        !is.data.frame(values)
    ) {
      return(values)
    }

    values <- as.matrix(
      values
    )

    if (nrow(values) == expected_internal_count) {
      return(
        values[
          seq_len(response_count),
          ,
          drop = FALSE
        ]
      )
    }

    if (nrow(values) == response_count) {
      return(values)
    }

    stop(
      paste0(
        "Cannot prepare the GLLVM summary: '",
        parameter_name,
        "' has ",
        nrow(values),
        " rows; expected ",
        response_count,
        " response-level rows or ",
        expected_internal_count,
        " component-level rows."
      ),
      call. = FALSE
    )
  }

  if (
    !is.null(object$params) &&
      is.list(object$params)
  ) {
    object$params$beta0 <- trim_species_vector(
      values = object$params$beta0,
      parameter_name = "params$beta0"
    )

    object$params$theta <- trim_species_rows(
      values = object$params$theta,
      parameter_name = "params$theta"
    )
  }

  if (
    !is.null(object$sd) &&
      is.list(object$sd)
  ) {
    object$sd$beta0 <- trim_species_vector(
      values = object$sd$beta0,
      parameter_name = "sd$beta0"
    )

    object$sd$theta <- trim_species_rows(
      values = object$sd$theta,
      parameter_name = "sd$theta"
    )
  }

  # summary.gllvm() should receive one family label per observed response,
  # rather than gllvm's expanded component-level representation.
  object$family <- response_family

  object
}
