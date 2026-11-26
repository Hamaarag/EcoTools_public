# ==============================================================================
# Unified Model Validity Gatekeeper
# ==============================================================================

#' Assert That a Model Is Valid for Its Intended Use
#'
#' @description
#' Applies a standardized validity gate to models returned by:
#'
#' 1. `fit_and_simulate()` for GLMMs and GAMs;
#' 2. `fit_gllvm_cached()` for GLLVMs.
#'
#' The `purpose` argument distinguishes a model intended for formal inference
#' from a converged model intended only for exploratory diagnostics or model
#' selection.
#'
#' This function does not repeat engine-specific convergence calculations.
#' Those checks belong inside the fitting wrappers. Instead, this function
#' enforces the standardized `model_status` assigned by those wrappers.
#'
#' @param model_object
#' A GLMM/GAM bundle returned by `fit_and_simulate()`, or a GLLVM model object
#' returned by `fit_gllvm_cached()`. The object must contain `model_status`.
#'
#' @param purpose
#' Intended use of the model:
#'
#' - `"inference"` requires a fully validated model with status `"valid"`.
#' - `"fit"` permits either `"valid"` or
#'   `"valid_fit_no_standard_errors"`.
#'
#' Use `"inference"` before extracting p-values, confidence intervals,
#' significance classifications, or inferential predictions.
#'
#' Use `"fit"` before exploratory diagnostics, ordination, preliminary model
#' assessment, or model selection using a converged fit without standard errors.
#'
#' @param allow_sparse
#' Logical. If `TRUE`, a failed sparse-data screen generates a warning but does
#' not stop execution. If `FALSE`, a failed sparse-data screen stops execution
#' even when the model otherwise passes the requested status gate.
#'
#' @param warn_sparse
#' Logical. If `TRUE`, print the cached sparse-data warning when
#' `allow_sparse = TRUE`. Ignored when no sparse warning is present.
#'
#' @param print_sparse_details
#' Logical. If `TRUE`, print the response/species-level sparse-screen details
#' when available.
#'
#' @return
#' Invisibly returns `TRUE` when the model passes the requested validity gate.
#'
#' @export
assert_valid_model <- function(
  model_object,
  purpose = c("inference", "fit"),
  allow_sparse = TRUE,
  warn_sparse = TRUE,
  print_sparse_details = FALSE
) {
  purpose <- match.arg(purpose)

  # ---------------------------------------------------------------------------
  # 0. Make sure the sparsity policy indicated by the user is clear
  # ---------------------------------------------------------------------------
  if (
    !is.logical(allow_sparse) ||
      length(allow_sparse) != 1L ||
      is.na(allow_sparse)
  ) {
    stop(
      "'allow_sparse' must be a single non-missing logical value.",
      call. = FALSE
    )
  }

  if (
    !is.logical(warn_sparse) ||
      length(warn_sparse) != 1L ||
      is.na(warn_sparse)
  ) {
    stop(
      "'warn_sparse' must be a single non-missing logical value.",
      call. = FALSE
    )
  }

  if (
    !is.logical(print_sparse_details) ||
      length(print_sparse_details) != 1L ||
      is.na(print_sparse_details)
  ) {
    stop(
      "'print_sparse_details' must be a single non-missing logical value.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 1. Validate object structure
  # ---------------------------------------------------------------------------
  if (
    !is.list(model_object) ||
      !"model_status" %in% names(model_object)
  ) {
    stop(
      paste0(
        "Model object lacks a standardized 'model_status'. ",
        "Please refit the model with the current wrapper version."
      ),
      call. = FALSE
    )
  }

  status <- model_object$model_status

  if (
    length(status) != 1L ||
      is.na(status) ||
      !is.character(status) ||
      !nzchar(status)
  ) {
    stop(
      paste0(
        "Model object contains an invalid or ambiguous 'model_status': ",
        paste(status, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 2. Extract cached warnings
  # ---------------------------------------------------------------------------
  stored_warnings <- NULL

  # GLMM/GAM bundles normally use $warnings.
  if ("warnings" %in% names(model_object) && !is.null(model_object$warnings)) {
    stored_warnings <- model_object$warnings
  }

  # GLLVM objects normally use $warning_msgs.
  if (
    is.null(stored_warnings) &&
      "warning_msgs" %in% names(model_object) &&
      !is.null(model_object$warning_msgs)
  ) {
    stored_warnings <- model_object$warning_msgs
  }

  if (!is.null(stored_warnings)) {
    stored_warnings <- unique(
      as.character(stored_warnings)
    )

    stored_warnings <- stored_warnings[
      !is.na(stored_warnings) &
        nzchar(stored_warnings)
    ]
  }

  sparse_screen <- if ("sparse_screen" %in% names(model_object)) {
    model_object$sparse_screen
  } else {
    NULL
  }

  # ---------------------------------------------------------------------------
  # 3. Define statuses allowed for each purpose
  # ---------------------------------------------------------------------------
  permitted_statuses <- switch(
    purpose,

    "inference" = c(
      "valid"
    ),

    "fit" = c(
      "valid",
      "valid_fit_no_standard_errors",
      "valid_fit_no_residual_diagnostics"
    )
  )

  # ---------------------------------------------------------------------------
  # 4. Explain standardized statuses
  # ---------------------------------------------------------------------------
  status_explanation <- switch(
    status,

    "valid" = paste0(
      "The model passed the wrapper's optimizer, gradient, coefficient, ",
      "and uncertainty checks applicable to this model."
    ),

    "valid_fit_no_standard_errors" = paste0(
      "The optimizer and applicable convergence checks passed, but standard ",
      "errors were deliberately not requested. The fit may be used for ",
      "exploratory diagnostics or model selection, but not for formal ",
      "standard-error-based inference."
    ),
    # Legacy status retained for compatibility with older cached model objects
    "suspicious_gradient" = paste0(
      "The optimizer may have reported success, but the maximum absolute ",
      "gradient exceeded the configured threshold. The fitted solution is ",
      "not accepted as a sufficiently stationary optimum."
    ),

    "suspicious_optimization_incomplete" = paste0(
      "The optimizer reported success, but the maximum absolute gradient ",
      "exceeded the configured threshold. The fitted solution is not accepted ",
      "as a sufficiently stationary optimum and must not be used for model ",
      "selection or formal inference."
    ),

    "suspicious_inference_unreliable" = paste0(
      "Point estimation may have converged, but standard-error, Hessian, ",
      "or variance-covariance validation failed. Wald p-values and ",
      "confidence intervals must not be reported."
    ),

    "valid_fit_no_residual_diagnostics" = paste0(
      "The model passed the available numerical fitting checks, but residual ",
      "simulation diagnostics were unavailable or unsuccessful. The fit is not ",
      "accepted for formal inference until residual adequacy has been assessed."
    ),

    "suspicious_residual_diagnostics" = paste0(
      "The model passed numerical optimizer and Hessian checks, but one or ",
      "more DHARMa residual diagnostics detected disagreement between the ",
      "observed data and simulations from the fitted model. Formal inference ",
      "from this model is not accepted without resolving or explicitly ",
      "justifying the residual lack of fit."
    ),

    "suspicious_optimization_unassessed" = paste0(
      "The optimizer reported success, but the final objective-function gradient ",
      "could not be calculated. The stationary-point criterion therefore remains ",
      "unverified."
    ),

    "suspicious_reference_sensitivity" = paste0(
      "Numerically equivalent factor reference-level parameterizations produced ",
      "different optimization outcomes, gradient classifications, or fitting ",
      "success. The model is sensitive to arbitrary coding choices and must not ",
      "be used for formal inference or model selection without resolving that ",
      "instability."
    ),

    "suspicious_extreme_coefficients" = paste0(
      "One or more coefficients exceeded the configured magnitude threshold. ",
      "This may indicate sparsity, separation, weak identification, or a ",
      "boundary estimate and requires manual review."
    ),

    "invalid" = paste0(
      "The optimizer or another fatal fitting check failed. The model must ",
      "not be used for estimation or inference."
    ),

    "aborted" = paste0(
      "The model was not fitted because a pre-fit governance rule aborted ",
      "the analysis."
    ),

    paste0(
      "The model has an unrecognized status code: '",
      status,
      "'. Manual review is required."
    )
  )

  # ---------------------------------------------------------------------------
  # 5. Enforce the requested status gate
  # ---------------------------------------------------------------------------
  if (!status %in% permitted_statuses) {
    purpose_description <- switch(
      purpose,

      "inference" = paste0(
        "formal inference, including coefficient p-values, confidence ",
        "intervals, significance classifications, and inferential predictions"
      ),

      "fit" = paste0(
        "exploratory diagnostics, ordination, preliminary model assessment, ",
        "or model selection"
      )
    )

    warning_text <- if (
      !is.null(stored_warnings) &&
        length(stored_warnings) > 0L
    ) {
      paste(
        stored_warnings,
        collapse = "\n  - "
      )
    } else {
      "No stored warnings."
    }

    stop(
      paste0(
        "Model is not valid for ",
        purpose_description,
        ".\n",
        "Model status: ",
        status,
        "\n",
        "Reason: ",
        status_explanation,
        "\n",
        "Stored warnings:\n  - ",
        warning_text
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 6. Apply sparse-data policy separately
  # ---------------------------------------------------------------------------

  # A sparse-screen result can have three states:
  #
  #   TRUE  = screen completed and passed;
  #   FALSE = screen completed and detected sparse data;
  #   NA    = screen failed or could not be assessed.
  sparse_pass <- if (
    !is.null(sparse_screen) &&
      is.list(sparse_screen) &&
      !is.null(sparse_screen$pass) &&
      length(sparse_screen$pass) == 1L
  ) {
    sparse_screen$pass
  } else {
    NA
  }

  sparse_failed <- identical(
    sparse_pass,
    FALSE
  )

  sparse_unassessed <- length(sparse_pass) == 1L &&
    is.na(sparse_pass) &&
    !is.null(sparse_screen)

  # Also detect an explicitly stored sparse-screen error.
  sparse_messages <- character(0)

  if (
    !is.null(sparse_screen) &&
      !is.null(sparse_screen$warning)
  ) {
    sparse_messages <- c(
      sparse_messages,
      as.character(sparse_screen$warning)
    )
  }

  if (!is.null(stored_warnings)) {
    sparse_messages <- c(
      sparse_messages,
      as.character(stored_warnings)
    )
  }

  sparse_messages <- sparse_messages[
    !is.na(sparse_messages) &
      nzchar(sparse_messages)
  ]

  sparse_error_detected <- any(
    grepl(
      "^SPARSE GLLVM SCREEN ERROR:",
      sparse_messages,
      ignore.case = TRUE
    )
  )

  sparse_unassessed <- sparse_unassessed ||
    sparse_error_detected

  # ---------------------------------------------------------------------------
  # 6a. Handle a failed or unassessed sparse screen
  # ---------------------------------------------------------------------------
  if (sparse_unassessed) {
    sparse_error_text <- if (
      !is.null(sparse_screen) &&
        !is.null(sparse_screen$warning) &&
        length(sparse_screen$warning) > 0L &&
        !is.na(sparse_screen$warning[1])
    ) {
      as.character(sparse_screen$warning[1])
    } else {
      "The sparse-data screen was not completed successfully."
    }

    stop(
      paste0(
        "Model passed the '",
        purpose,
        "' model-status gate, but sparse-data validity could not be assessed.\n",
        sparse_error_text,
        "\nRefit the model with the current wrapper before proceeding."
      ),
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 6b. Handle a completed screen that detected sparse data
  # ---------------------------------------------------------------------------
  if (sparse_failed) {
    sparse_warning <- if (
      !is.null(sparse_screen$warning) &&
        length(sparse_screen$warning) > 0L &&
        !is.na(sparse_screen$warning[1])
    ) {
      as.character(sparse_screen$warning[1])
    } else {
      "The model failed the sparse-data screen."
    }

    if (!isTRUE(allow_sparse)) {
      stop(
        paste0(
          "Model passed the '",
          purpose,
          "' model-status gate but failed the sparse-data screen.\n",
          sparse_warning
        ),
        call. = FALSE
      )
    }

    if (isTRUE(warn_sparse)) {
      warning(
        paste0(
          "Model passed the '",
          purpose,
          "' model-status gate but has a sparse-data warning:\n",
          sparse_warning
        ),
        call. = FALSE
      )
    }

    if (isTRUE(print_sparse_details)) {
      print_sparse_gllvm_report(
        sparse_screen_obj = sparse_screen,
        prefix = paste0(
          "MODEL SPARSE-DATA SCREEN FOR ",
          toupper(purpose)
        )
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 7. Return successfully
  # ---------------------------------------------------------------------------
  invisible(TRUE)
}
