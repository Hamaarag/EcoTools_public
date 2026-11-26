# ==============================================================================
# emmeans support for EcoTools model bundles
# ==============================================================================

#' Recover Data from an EcoTools Model Bundle
#'
#' @description
#' Provides `emmeans` with the model data required to construct a reference grid
#' for a model bundle returned by [fit_and_simulate()].
#'
#' The method delegates model-specific reconstruction to the `emmeans` method
#' for the embedded fitted model. When the EcoTools bundle contains
#' `analysis_data`, that stored data frame is supplied explicitly so the result
#' does not depend on objects in the interactive global environment.
#'
#' @param object
#' An object of class `"ecotools_model_bundle"` returned by
#' [fit_and_simulate()].
#'
#' @param ...
#' Additional arguments passed to [emmeans::recover_data()].
#'
#' @return
#' A data frame with the attributes required by `emmeans`.
#'
#' @export
recover_data.ecotools_model_bundle <- function(
  object,
  ...
) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop(
      "Package 'emmeans' is required for estimated marginal means.",
      call. = FALSE
    )
  }

  if (
    !is.list(object) ||
      is.null(object$model)
  ) {
    stop(
      paste0(
        "The EcoTools bundle does not contain a fitted model. ",
        "Estimated marginal means cannot be calculated."
      ),
      call. = FALSE
    )
  }

  fitted_model <- object$model
  dots <- list(...)

  dots$object <- NULL

  if (
    !is.null(object$analysis_data) &&
      !"data" %in% names(dots)
  ) {
    recovery_arguments <- c(
      list(
        object = fitted_model,
        data = as.data.frame(
          object$analysis_data,
          stringsAsFactors = FALSE
        )
      ),
      dots
    )
  } else {
    recovery_arguments <- c(
      list(
        object = fitted_model
      ),
      dots
    )
  }

  do.call(
    emmeans::recover_data,
    recovery_arguments
  )
}


#' Construct an emmeans Basis for an EcoTools Model Bundle
#'
#' @description
#' Supplies the coefficient, covariance, and linear-function basis needed by
#' `emmeans` for a model bundle returned by [fit_and_simulate()].
#'
#' The method delegates to the existing `emmeans` method for the embedded
#' `glmmTMB`, `gam`, or other supported fitted model.
#'
#' @param object
#' An object of class `"ecotools_model_bundle"` returned by
#' [fit_and_simulate()].
#'
#' @param trms
#' Fixed-effect terms supplied by `emmeans`.
#'
#' @param xlev
#' Predictor-level information supplied by `emmeans`.
#'
#' @param grid
#' Reference grid supplied by `emmeans`.
#'
#' @param ...
#' Additional arguments passed to [emmeans::emm_basis()].
#'
#' @return
#' A list containing the linear-function basis and associated inferential
#' information required by `emmeans`.
#'
#' @export
emm_basis.ecotools_model_bundle <- function(
  object,
  trms,
  xlev,
  grid,
  ...
) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop(
      "Package 'emmeans' is required for estimated marginal means.",
      call. = FALSE
    )
  }

  if (
    !is.list(object) ||
      is.null(object$model)
  ) {
    stop(
      paste0(
        "The EcoTools bundle does not contain a fitted model. ",
        "An emmeans basis cannot be constructed."
      ),
      call. = FALSE
    )
  }

  emmeans::emm_basis(
    object$model,
    trms = trms,
    xlev = xlev,
    grid = grid,
    ...
  )
}

# ==============================================================================
# Validated estimated marginal means
# ==============================================================================

#' Calculate Estimated Marginal Means Safely
#'
#' @description
#' Calculates estimated marginal means from an EcoTools model bundle or a
#' directly fitted model object while supplying recoverable model data
#' explicitly.
#'
#' When an EcoTools bundle is supplied, the function validates the bundle for
#' formal inference and uses its stored `analysis_data`. When a directly fitted
#' model is supplied, the function uses the model's stored frame where
#' available.
#'
#' @param object
#' An EcoTools model bundle or a directly fitted model supported by `emmeans`.
#'
#' @param specs
#' Specifications passed to [emmeans::emmeans()].
#'
#' @param data
#' Optional model data. When `NULL`, the function uses bundle
#' `analysis_data` or the fitted model's stored frame.
#'
#' @param validate
#' Logical. If `TRUE`, an EcoTools bundle must pass
#' [assert_valid_model()] for inference. Direct fitted models cannot undergo
#' complete EcoTools bundle validation.
#'
#' @param ...
#' Additional arguments passed to [emmeans::emmeans()].
#'
#' @return
#' An object returned by [emmeans::emmeans()].
#'
#' @export
emmeans_checked <- function(
  object,
  specs,
  data = NULL,
  validate = TRUE,
  ...
) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop(
      "Package 'emmeans' is required.",
      call. = FALSE
    )
  }

  is_bundle <- inherits(
    object,
    "ecotools_model_bundle"
  ) ||
    (is.list(object) &&
      "model" %in% names(object) &&
      "model_status" %in% names(object))

  if (is_bundle) {
    if (isTRUE(validate)) {
      assert_valid_model(
        object,
        purpose = "inference"
      )
    }

    fitted_model <- object$model

    if (
      is.null(data) &&
        !is.null(object$analysis_data)
    ) {
      data <- as.data.frame(
        object$analysis_data,
        stringsAsFactors = FALSE
      )
    }
  } else {
    fitted_model <- object

    if (isTRUE(validate)) {
      warning(
        paste0(
          "A directly fitted model was supplied, so EcoTools bundle ",
          "validity metadata could not be checked."
        ),
        call. = FALSE
      )
    }

    if (
      is.null(data) &&
        !is.null(fitted_model$frame)
    ) {
      data <- as.data.frame(
        fitted_model$frame,
        stringsAsFactors = FALSE
      )
    }
  }

  if (is.null(fitted_model)) {
    stop(
      "No fitted model was found.",
      call. = FALSE
    )
  }

  emmeans_arguments <- list(
    object = fitted_model,
    specs = specs
  )

  if (!is.null(data)) {
    emmeans_arguments$data <- data
  }

  emmeans_arguments <- c(
    emmeans_arguments,
    list(...)
  )

  do.call(
    emmeans::emmeans,
    emmeans_arguments
  )
}
