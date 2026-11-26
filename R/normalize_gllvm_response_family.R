#' Normalize GLLVM Response Families
#'
#' Expands a shared family name to all response columns, or validates a vector
#' containing one family name per response column.
#'
#' @param family A family name or character vector of family names.
#' @param n_responses Positive integer number of response columns.
#'
#' @return A character vector with one family name per response column.
#'
#' @keywords internal
normalize_gllvm_response_family <- function(
  family,
  n_responses
) {
  if (
    !is.numeric(n_responses) ||
      length(n_responses) != 1L ||
      is.na(n_responses) ||
      !is.finite(n_responses) ||
      n_responses < 1L ||
      n_responses != as.integer(n_responses)
  ) {
    stop(
      "'n_responses' must be one positive integer.",
      call. = FALSE
    )
  }

  family_names <- tryCatch(
    as.character(family),
    error = function(e) {
      stop(
        "'family' must be coercible to character family names.",
        call. = FALSE
      )
    }
  )

  if (
    length(family_names) == 0L ||
      anyNA(family_names) ||
      any(!nzchar(family_names))
  ) {
    stop(
      "'family' must contain non-missing, non-empty family names.",
      call. = FALSE
    )
  }

  if (length(family_names) == 1L) {
    return(
      rep(
        family_names,
        as.integer(n_responses)
      )
    )
  }

  if (length(family_names) != n_responses) {
    stop(
      paste0(
        "'family' must contain either one shared family or one family per ",
        "response column (",
        n_responses,
        ")."
      ),
      call. = FALSE
    )
  }

  family_names
}
