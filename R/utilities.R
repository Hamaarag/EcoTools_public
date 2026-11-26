# ==============================================================================
# General utility functions
# ==============================================================================

#' Reverse Character Strings
#'
#' Reverses the characters in each supplied string. This utility may be useful
#' for graphics devices that do not correctly render right-to-left text.
#'
#' @param x A character vector or an object coercible to character.
#'
#' @return A character vector containing the reversed strings. Returns `NULL`
#'   when `x` is `NULL`.
#'
#' @details
#' This function reverses Unicode code points in their stored order. It is
#' intended primarily as a compatibility workaround for graphics devices that
#' display Hebrew text from left to right.
#'
#' Character reversal is not a complete implementation of the Unicode
#' bidirectional algorithm. Labels containing a mixture of Hebrew, Latin text,
#' numbers, parentheses, or punctuation may therefore require a graphics device
#' with native bidirectional-text support instead.
#'
#' @references
#' The initial implementation was inspired by a base-R string-reversal idiom
#' posted by Josh O'Brien on Stack Overflow:
#' \url{https://stackoverflow.com/a/13613183}
#'
#' The implementation used here was subsequently rewritten to provide explicit
#' input coercion, `NULL` handling, and type-stable vectorized output.
#'
#' @examples
#' strReverse("abc")
#' strReverse(c("abc", "def"))
#' strReverse(NULL)
#'
#' @export
strReverse <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  x_char <- as.character(x)

  vapply(
    strsplit(x_char, split = "", fixed = TRUE),
    function(z) paste(rev(z), collapse = ""),
    character(1)
  )
}

#' Calculate the Standard Error
#'
#' Calculates the standard error of a numeric vector, with optional removal of
#' missing values.
#'
#' @param x A numeric vector.
#' @param na.rm Logical. If `TRUE`, missing values are removed before
#'   calculation.
#'
#' @return A numeric value. Returns `NA_real_` when fewer than two usable values
#'   are available.
#'
#' @examples
#' std_error(c(1, 2, 3, 4))
#' std_error(c(1, 2, NA, 4), na.rm = TRUE)
#'
#' @export
std_error <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }

  if (length(x) <= 1L) {
    return(NA_real_)
  }

  stats::sd(x) / sqrt(length(x))
}

#' Automatically choose the number of CPU cores for model fitting
#' Uses physical cores where possible, leaves one core free,
#' and caps the absolute maximum to prevent OpenMP stack overflows.
#'
#' @param user_cores Specify the number of cores you want to use, if you want to force this number
#' @param max_cap The maximum number of cores to use, even if your machine has more; this limit is
#'  set to prevent OpenMP stack overflows.
#'
#' @return The number of cores that will be used to fit the model
#'
#' @export
get_auto_cores <- function(user_cores = NULL, max_cap = 6) {
  # If the user explicitly asks for a specific number, respect it
  if (!is.null(user_cores)) {
    return(as.integer(user_cores))
  }

  # Detect physical cores (avoids hyper-threading overhead)
  phys_cores <- parallel::detectCores(logical = FALSE)

  if (is.na(phys_cores) || phys_cores < 1) {
    phys_cores <- parallel::detectCores(logical = TRUE)
  }

  if (is.na(phys_cores) || phys_cores < 1) {
    return(1L)
  }

  # Leave one core free for the OS
  safe_cores <- max(1L, as.integer(phys_cores) - 1L)

  # Apply the OpenMP stability cap
  return(min(safe_cores, as.integer(max_cap)))
}

#' Sanitize Text *or Use in a File Name
#'
#' Replaces runs of non-alphanumeric characters with underscores and optionally
#' truncates the resulting text.
#'
#' @param x A character vector or an object coercible to character.
#' @param max_width Optional positive integer giving the maximum number of
#'   characters retained in each result.
#'
#' @return A character vector suitable for use as file-name components.
#'
#' @examples
#' sanitize_filename_part("Residuals: treatment/site")
#'
#' @export
sanitize_filename_part <- function(x, max_width = NULL) {
  x <- as.character(x)

  if (
    !is.null(max_width) &&
      (!is.numeric(max_width) ||
        length(max_width) != 1L ||
        is.na(max_width) ||
        !is.finite(max_width) ||
        max_width < 1L ||
        max_width != as.integer(max_width))
  ) {
    stop(
      "'max_width' must be NULL or one positive integer.",
      call. = FALSE
    )
  }

  out <- x |>
    stringr::str_replace_all("[^[:alnum:]_]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_replace_all("^_|_$", "")

  if (!is.null(max_width)) {
    out <- stringr::str_trunc(out, width = as.integer(max_width), ellipsis = "")
  }

  out
}

#' Determine Whether a Model Response Contains Zeros
#'
#' Attempts to determine whether the observed response associated with a fitted
#' model or DHARMa residual object contains at least one zero. This is used internally
#' to decide whether zero dispersion should be checked.
#'
#' @param model A fitted model object. The function attempts to obtain the model
#'   response from the model frame.
#' @param residuals_object An optional DHARMa residual object. When the response
#'   cannot be recovered from `model`, the function attempts to use the observed
#'   response stored in this object.
#'
#' @return A logical scalar. Returns `TRUE` when at least one observed response
#'   equals zero, `FALSE` when no zeros are present, and possibly `NA` when the
#'   response cannot be recovered.
#'
#' @keywords internal
response_has_zero <- function(model = NULL, residuals_object = NULL) {
  y <- NULL

  if (!is.null(residuals_object)) {
    y <- tryCatch(
      residuals_object$observedResponse,
      error = function(e) NULL
    )
  }

  if (is.null(y) && !is.null(model)) {
    y <- tryCatch(
      model.response(model.frame(model)),
      error = function(e) NULL
    )
  }

  if (is.null(y)) {
    return(NA)
  }

  y <- suppressWarnings(as.numeric(y))
  y <- y[is.finite(y)]

  if (length(y) == 0) {
    return(NA)
  }

  any(y == 0)
}

#' @keywords internal
is_rendering_session <- function() {
  isTRUE(getOption("knitr.in.progress")) ||
    !is.null(getOption("quarto.version"))
}
