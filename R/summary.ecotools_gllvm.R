#' Summarize an EcoTools GLLVM
#'
#' Delegates to `gllvm::summary.gllvm()`, adapting beta-hurdle parameter arrays
#' on a temporary copy to work around mixed-family row-name handling in some
#' `gllvm` versions.
#'
#' @param object A GLLVM returned by [fit_gllvm_cached()].
#' @param ... Additional arguments passed to `summary.gllvm()`.
#'
#' @return An object of class `summary.gllvm`.
#'
#' @keywords internal
#' @exportS3Method summary ecotools_gllvm
summary.ecotools_gllvm <- function(
  object,
  ...
) {
  summary_method <- utils::getS3method(
    "summary",
    "gllvm"
  )

  summary_method(
    object = prepare_gllvm_summary_object(object),
    ...
  )
}
