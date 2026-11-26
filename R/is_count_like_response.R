#' Identify a Non-Negative Integer-Like Response
#'
#' Determines whether all non-missing response values are non-negative and
#' sufficiently close to integers.
#'
#' @param x A numeric response vector.
#' @param tolerance Numeric tolerance used when comparing values with rounded
#'   integers.
#'
#' @return A logical scalar.
#'
#' @keywords internal
is_count_like_response <- function(
  x,
  tolerance = 1e-8
) {
  if (
    !is.numeric(tolerance) ||
      length(tolerance) != 1L ||
      is.na(tolerance) ||
      !is.finite(tolerance) ||
      tolerance < 0
  ) {
    stop(
      "'tolerance' must be one finite, non-negative numeric value.",
      call. = FALSE
    )
  }

  if (!is.numeric(x) && !is.integer(x)) {
    return(FALSE)
  }

  # Treat NaN and infinite values as invalid response values. Ordinary NA
  # values are omitted only for the purpose of classifying the observed values.
  if (any(is.nan(x))) {
    return(FALSE)
  }

  observed_values <- x[!is.na(x)]

  if (length(observed_values) == 0L) {
    return(FALSE)
  }

  if (any(!is.finite(observed_values))) {
    return(FALSE)
  }

  if (any(observed_values < 0)) {
    return(FALSE)
  }

  all(
    abs(
      observed_values - round(observed_values)
    ) <=
      tolerance
  )
}
