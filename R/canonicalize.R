# ==============================================================================
# Stable fingerprint helpers
# ==============================================================================

#' Canonicalize tabular data for cache hashing
#'
#' Converts a tabular object into a stable representation suitable for
#' cross-session cache fingerprinting.
#' Convert tabular input into a stable plain-data-frame representation before
#' hashing.
#' This removes data.table-specific bookkeeping such as internal self
#' references while preserving values, column names, column classes, and row order.
#'
#' @param x A tabular object.
#'
#' @return A plain data frame or `NULL`.
#'
#' @keywords internal
canonicalize_tabular_for_hash <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  out <- as.data.frame(x, stringsAsFactors = FALSE)

  # Remove volatile or irrelevant data-frame-level attributes.
  attributes(out) <- attributes(out)[
    intersect(
      names(attributes(out)),
      c("names", "row.names", "class")
    )
  ]

  class(out) <- "data.frame"
  row.names(out) <- seq_len(nrow(out))

  out
}

# Recursively canonicalize additional arguments passed through "...".
# This prevents formula environments and data.table internals inside extra
# arguments from destabilizing the fingerprint.
#' Canonicalize Additional Arguments for Hashing
#'
#' Recursively converts objects into stable representations suitable for
#' cross-session cache fingerprints.
#'
#' @param x An R object.
#'
#' @return A canonicalized representation of `x`.
#'
#' @keywords internal
canonicalize_extra_args_for_hash <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.environment(x)) {
    stop(
      "Environments cannot be included in a stable cache fingerprint.",
      call. = FALSE
    )
  }

  if (is.function(x)) {
    stop(
      "Functions passed through additional arguments cannot currently be ",
      "included in a stable cross-session cache fingerprint.",
      call. = FALSE
    )
  }

  if (
    inherits(x, "data.table") ||
      is.data.frame(x) ||
      is.matrix(x)
  ) {
    return(canonicalize_tabular_for_hash(x))
  }

  if (
    inherits(x, "formula") ||
      is.call(x) ||
      is.expression(x)
  ) {
    return(canonicalize_language_for_hash(x))
  }

  if (is.list(x)) {
    if (!is.null(names(x))) {
      x <- x[order(names(x))]
    }

    return(lapply(x, canonicalize_extra_args_for_hash))
  }

  x
}


# Convert formulas and calls to stable one-line text without their environments.
#' Recursively converts objects into stable representations suitable for
#' cross-session cache fingerprints.
#'
#' @param x An R object.
#'
#' @return A canonicalized representation of `x`.
#'
#' @keywords internal
canonicalize_language_for_hash <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (isFALSE(x)) {
    return(FALSE)
  }
  if (isTRUE(x)) {
    return(TRUE)
  }

  if (inherits(x, "formula") || is.call(x) || is.expression(x)) {
    return(paste(deparse(x, width.cutoff = 500L), collapse = ""))
  }

  x
}
