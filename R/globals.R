# ------------------------------------------------------------------------------
# Wrapper version
# ------------------------------------------------------------------------------
# Included in cache fingerprints. Increment when wrapper behavior changes in a way
# that should force refitting old cached models.

#' Get Analysis Wrapper Version
#'
#' Retrieves the current version string used for cache fingerprinting.
#'
#' @return A character string.
#' @export
get_wrapper_version <- function() {
  "2026-11-18_v24"
}

#' EcoTools Plotting Defaults
#'
#' Default configuration values used by plotting helpers in EcoTools.
#'
#' @format
#' \describe{
#'   \item{fontname}{A character scalar naming the default plot font.}
#'   \item{fontsize}{A numeric scalar giving the default base font size.}
#'   \item{okabe}{A character vector containing an accessible colour palette.}
#'   \item{pdf_width}{A numeric scalar giving the default PDF width in mm.}
#'   \item{pdf_aspect_ratio}{A numeric scalar giving the default PDF aspect ratio.}
#' }
#'
#' @name plotting_defaults
NULL

#' @rdname plotting_defaults
fontname <- "Almoni ML v5 AAA"

#' @rdname plotting_defaults
fontsize <- 22

#' @rdname plotting_defaults
okabe <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7"
)

#' @rdname plotting_defaults
pdf_width <- 160

#' @rdname plotting_defaults
pdf_aspect_ratio <- 2 / 3

utils::globalVariables(c(
  ".",
  ".N",
  ".SD",
  ":=",
  "ColIndex",
  "Fitted",
  "Residuals",
  "RowIndex",
  "SqrtAbsRes",
  "component",
  "diagnostic",
  "estimate",
  "grouping_variable",
  "heb_name",
  "level",
  "median_p",
  "min_p",
  "n",
  "n_nonfinite",
  "p.value",
  "pred_heb",
  "prop_borderline",
  "prop_significant",
  "q10_p",
  "q90_p",
  "significance_alpha",
  "significant",
  "sparse_flag",
  "std.error",
  "variable",
  "x",
  "ymax",
  "ymin",
  "Family"
))
