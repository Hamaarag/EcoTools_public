# ==============================================================================
# Convergence Report Helper
# ==============================================================================

#' Print Model Convergence Report
#'
#' @description
#' Prints a detailed convergence and inference-quality report for a GLLVM model.
#'
#' The report distinguishes among:
#'
#' \enumerate{
#'   \item optimizer convergence;
#'   \item gradient assessment;
#'   \item availability and finiteness of fixed-effect standard errors;
#'   \item validity of the broader variance-covariance calculation;
#'   \item extreme fixed-effect coefficients;
#'   \item model suitability for exploratory use or formal Wald inference;
#'   \item sparse-response screening status.
#' }
#'
#' @param model A fitted or standardized GLLVM object returned by
#'   [fit_gllvm_cached()].
#' @param model_name Optional character scalar identifying the model. When
#'   `NULL`, the function uses `model$model_name` when available.
#'
#' @return Invisibly returns the input model object.
#'
#' @export
print_convergence_report <- function(
  model,
  model_name = NULL
) {
  # ---------------------------------------------------------------------------
  # 1. Basic model metadata
  # ---------------------------------------------------------------------------
  if (is.null(model_name)) {
    model_name <- if (
      !is.null(model$model_name) &&
        length(model$model_name) == 1L &&
        !is.na(model$model_name) &&
        nzchar(as.character(model$model_name))
    ) {
      as.character(model$model_name)
    } else {
      "Unnamed GLLVM model"
    }
  }

  model_status <- if (
    !is.null(model$model_status) &&
      length(model$model_status) == 1L
  ) {
    as.character(model$model_status)
  } else {
    "unknown"
  }

  opt_converged <- isTRUE(model$opt_converged)

  inference_quality <- tryCatch(
    model$diagnostics$inference_quality,
    error = function(e) NULL
  )

  coef_screen <- tryCatch(
    model$diagnostics$coef_screen,
    error = function(e) NULL
  )

  # ---------------------------------------------------------------------------
  # 2. Gradient metadata
  # ---------------------------------------------------------------------------
  gradient_assessed <- if (!is.null(model$gradient_assessed)) {
    isTRUE(model$gradient_assessed)
  } else if (
    !is.null(inference_quality) &&
      !is.null(inference_quality$gradient_available)
  ) {
    isTRUE(inference_quality$gradient_available)
  } else {
    FALSE
  }

  gradient_valid <- if (!is.null(model$gradient_valid)) {
    model$gradient_valid
  } else if (
    !is.null(inference_quality) &&
      !is.null(inference_quality$gradient_ok)
  ) {
    inference_quality$gradient_ok
  } else {
    NA
  }

  max_abs_gradient <- if (!is.null(model$max_abs_gradient)) {
    as.numeric(model$max_abs_gradient)[1]
  } else if (
    !is.null(inference_quality) &&
      !is.null(inference_quality$max_abs_gradient)
  ) {
    as.numeric(inference_quality$max_abs_gradient)[1]
  } else {
    NA_real_
  }

  gradient_threshold <- if (!is.null(model$gradient_threshold)) {
    as.numeric(model$gradient_threshold)[1]
  } else if (
    !is.null(inference_quality) &&
      !is.null(inference_quality$gradient_threshold)
  ) {
    as.numeric(inference_quality$gradient_threshold)[1]
  } else {
    NA_real_
  }

  # ---------------------------------------------------------------------------
  # 3. Standard-error metadata
  # ---------------------------------------------------------------------------
  sd_errors_requested <- if (!is.null(model$sd_errors_requested)) {
    isTRUE(model$sd_errors_requested)
  } else if (
    !is.null(inference_quality) &&
      !is.null(inference_quality$se_requested)
  ) {
    isTRUE(inference_quality$se_requested)
  } else {
    !identical(
      model_status,
      "valid_fit_no_standard_errors"
    ) &&
      !is.null(model$sd) &&
      is.list(model$sd)
  }

  # Preserve TRUE, FALSE, and NA rather than collapsing NA into FALSE.
  se_valid <- if (!is.null(model$se_valid)) {
    model$se_valid
  } else if (
    !is.null(inference_quality) &&
      !is.null(inference_quality$se_valid)
  ) {
    inference_quality$se_valid
  } else {
    NA
  }

  se_total <- if (
    !is.null(model$sd) &&
      is.list(model$sd) &&
      !is.null(model$sd$Xcoef)
  ) {
    length(model$sd$Xcoef)
  } else {
    0L
  }

  se_finite <- if (se_total > 0L) {
    sum(is.finite(model$sd$Xcoef))
  } else {
    0L
  }

  # ---------------------------------------------------------------------------
  # 4. Covariance-validation reasons
  # ---------------------------------------------------------------------------
  covariance_failure_reasons <- character(0)

  if (!is.null(inference_quality)) {
    if (isTRUE(inference_quality$negative_variance_warning)) {
      negative_variance_count <-
        inference_quality$negative_variance_count

      covariance_failure_reasons <- c(
        covariance_failure_reasons,
        if (
          is.null(negative_variance_count) ||
            length(negative_variance_count) == 0L ||
            is.na(negative_variance_count)
        ) {
          "negative variance estimates reported"
        } else {
          paste0(
            as.integer(negative_variance_count),
            " negative variance estimates reported"
          )
        }
      )
    }

    if (isTRUE(inference_quality$singular_covariance_warning)) {
      covariance_failure_reasons <- c(
        covariance_failure_reasons,
        paste0(
          "singular or zero-determinant ",
          "variance-covariance matrix"
        )
      )
    }

    if (
      isTRUE(inference_quality$hessian_warning) &&
        !isTRUE(inference_quality$singular_covariance_warning)
    ) {
      covariance_failure_reasons <- c(
        covariance_failure_reasons,
        "Hessian warning reported"
      )
    }

    if (isTRUE(inference_quality$any_nonfinite_se)) {
      covariance_failure_reasons <- c(
        covariance_failure_reasons,
        paste0(
          "at least one numeric standard-error component ",
          "contains non-finite values"
        )
      )
    }
  }

  # Compatibility fallback for model objects created by older wrapper versions.
  if (
    length(covariance_failure_reasons) == 0L &&
      !is.null(model$covariance_failure_reasons)
  ) {
    covariance_failure_reasons <- as.character(
      model$covariance_failure_reasons
    )
  }

  covariance_failure_reasons <- unique(
    covariance_failure_reasons[
      !is.na(covariance_failure_reasons) &
        nzchar(covariance_failure_reasons)
    ]
  )

  # ---------------------------------------------------------------------------
  # 5. Extreme-coefficient status
  # ---------------------------------------------------------------------------
  has_extreme_coef <- if (
    !is.null(coef_screen) &&
      !is.null(coef_screen$n_extreme) &&
      length(coef_screen$n_extreme) == 1L &&
      is.finite(coef_screen$n_extreme)
  ) {
    coef_screen$n_extreme > 0L
  } else {
    identical(
      model_status,
      "suspicious_extreme_coefficients"
    )
  }

  extreme_text <- if (has_extreme_coef) {
    if (
      !is.null(coef_screen$n_extreme) &&
        !is.null(coef_screen$coef_threshold)
    ) {
      paste0(
        "FLAGGED: ",
        coef_screen$n_extreme,
        " coefficient(s) exceeded |beta| > ",
        coef_screen$coef_threshold
      )
    } else {
      "FLAGGED"
    }
  } else {
    "NONE FLAGGED"
  }

  # ---------------------------------------------------------------------------
  # 5b. Reference-level stability diagnostics
  # ---------------------------------------------------------------------------

  reference_stability <- tryCatch(
    model$reference_stability,
    error = function(e) NULL
  )

  reference_sensitive <- isTRUE(
    tryCatch(
      model$reference_sensitive,
      error = function(e) FALSE
    )
  )

  # Current objects store metadata around the actual check result. Older cached
  # objects may store the check result directly, so retain compatibility with both.
  reference_metadata_format <-
    is.list(reference_stability) &&
    any(
      c(
        "requested_policy",
        "candidate_eligible",
        "check_run",
        "skip_reason",
        "detected_factor_predictors",
        "result"
      ) %in%
        names(reference_stability)
    )

  if (reference_metadata_format) {
    reference_check_run <- isTRUE(
      reference_stability$check_run
    )

    reference_skip_reason <- if (
      !is.null(reference_stability$skip_reason) &&
        length(reference_stability$skip_reason) == 1L &&
        !is.na(reference_stability$skip_reason)
    ) {
      as.character(
        reference_stability$skip_reason
      )
    } else {
      NULL
    }

    detected_factor_predictors <- if (
      !is.null(
        reference_stability$detected_factor_predictors
      )
    ) {
      as.character(
        reference_stability$detected_factor_predictors
      )
    } else {
      character(0)
    }

    reference_result <-
      reference_stability$result

    reference_check_warned <- isTRUE(
      reference_stability$completed_with_warnings
    )

    reference_check_completed <- if (
      !is.null(
        reference_stability$check_completed
      ) &&
        length(
          reference_stability$check_completed
        ) ==
          1L
    ) {
      reference_stability$check_completed
    } else if (
      !is.null(reference_result) &&
        !is.null(
          reference_result$completed
        )
    ) {
      reference_result$completed
    } else {
      NA
    }
  } else {
    reference_check_run <-
      !is.null(reference_stability)

    reference_skip_reason <- NULL

    detected_factor_predictors <- if (
      !is.null(reference_stability) &&
        !is.null(
          reference_stability$factor_predictors
        )
    ) {
      as.character(
        reference_stability$factor_predictors
      )
    } else {
      character(0)
    }

    reference_result <-
      reference_stability

    reference_check_warned <- if (
      !is.null(reference_result) &&
        !is.null(
          reference_result$completed_with_warnings
        )
    ) {
      isTRUE(
        reference_result$completed_with_warnings
      )
    } else {
      FALSE
    }

    reference_check_completed <- if (
      !is.null(reference_result) &&
        !is.null(
          reference_result$completed
        )
    ) {
      reference_result$completed
    } else {
      NA
    }
  }

  checked_factor_predictors <- if (
    !is.null(reference_result) &&
      !is.null(reference_result$factor_predictors)
  ) {
    as.character(
      reference_result$factor_predictors
    )
  } else {
    character(0)
  }

  excluded_factor_predictors <- if (
    !is.null(reference_result) &&
      !is.null(reference_result$excluded_factor_predictors)
  ) {
    as.character(
      reference_result$excluded_factor_predictors
    )
  } else {
    character(0)
  }

  reference_level_limit <- if (
    !is.null(reference_result) &&
      !is.null(reference_result$level_limit) &&
      length(reference_result$level_limit) == 1L
  ) {
    reference_result$level_limit
  } else {
    NA_integer_
  }

  reference_text <- if (is.null(reference_stability)) {
    "NOT RUN / NOT AVAILABLE"
  } else if (identical(reference_skip_reason, "disabled")) {
    "NOT RUN: disabled by reference_stability_check = 'none'"
  } else if (identical(reference_skip_reason, "optimizer_not_converged")) {
    "SKIPPED: optimizer did not converge"
  } else if (identical(reference_skip_reason, "gradient_not_assessed")) {
    "SKIPPED: final gradient was not assessed"
  } else if (identical(reference_skip_reason, "gradient_failed")) {
    paste0(
      "SKIPPED: gradient check failed",
      if (length(detected_factor_predictors) > 0L) {
        paste0(
          " [detected factor(s): ",
          paste(
            detected_factor_predictors,
            collapse = ", "
          ),
          "]"
        )
      } else {
        ""
      }
    )
  } else if (identical(reference_skip_reason, "standard_errors_invalid")) {
    paste0(
      "SKIPPED: covariance / standard-error validation failed",
      if (length(detected_factor_predictors) > 0L) {
        paste0(
          " [detected factor(s): ",
          paste(
            detected_factor_predictors,
            collapse = ", "
          ),
          "]"
        )
      } else {
        ""
      }
    )
  } else if (identical(reference_skip_reason, "no_factor_predictors")) {
    "NO FACTOR PREDICTORS DETECTED"
  } else if (!isTRUE(reference_check_run)) {
    "NOT RUN"
  } else if (
    length(checked_factor_predictors) == 0L &&
      length(excluded_factor_predictors) > 0L
  ) {
    paste0(
      "SKIPPED: factor(s) exceeded level limit",
      if (is.finite(reference_level_limit)) {
        paste0(
          " (",
          reference_level_limit,
          ")"
        )
      } else {
        ""
      },
      " [",
      paste(
        excluded_factor_predictors,
        collapse = ", "
      ),
      "]"
    )
  } else if (length(checked_factor_predictors) == 0L) {
    "NO ELIGIBLE FACTOR PREDICTORS CHECKED"
  } else if (
    identical(
      reference_check_completed,
      FALSE
    )
  ) {
    paste0(
      "INCOMPLETE: one or more alternative-reference fits failed",
      if (length(checked_factor_predictors) > 0L) {
        paste0(
          " [",
          paste(
            checked_factor_predictors,
            collapse = ", "
          ),
          "]"
        )
      } else {
        ""
      }
    )
  } else if (reference_sensitive) {
    paste0(
      "FLAGGED [",
      paste(
        checked_factor_predictors,
        collapse = ", "
      ),
      "]"
    )
  } else if (isTRUE(reference_check_warned)) {
    paste0(
      "PASS WITH WARNINGS [",
      paste(
        checked_factor_predictors,
        collapse = ", "
      ),
      "]"
    )
  } else {
    paste0(
      "PASS [",
      paste(
        checked_factor_predictors,
        collapse = ", "
      ),
      "]"
    )
  }

  # ---------------------------------------------------------------------------
  # 6. Sparse-response screening status
  # ---------------------------------------------------------------------------
  sparse_screen_obj <- if (!is.null(model$sparse_screen)) {
    model$sparse_screen
  } else {
    NULL
  }

  sparse_pass <- if (
    !is.null(sparse_screen_obj) &&
      !is.null(sparse_screen_obj$pass) &&
      length(sparse_screen_obj$pass) == 1L
  ) {
    sparse_screen_obj$pass
  } else {
    NA
  }

  sparse_error_messages <- character(0)

  if (
    !is.null(sparse_screen_obj) &&
      !is.null(sparse_screen_obj$warning)
  ) {
    sparse_error_messages <- c(
      sparse_error_messages,
      as.character(sparse_screen_obj$warning)
    )
  }

  if (!is.null(model$warning_msgs)) {
    sparse_error_messages <- c(
      sparse_error_messages,
      as.character(model$warning_msgs)
    )
  }

  sparse_error_messages <- sparse_error_messages[
    !is.na(sparse_error_messages)
  ]

  sparse_screen_error <- any(
    grepl(
      "^SPARSE GLLVM SCREEN ERROR:",
      sparse_error_messages,
      ignore.case = TRUE
    )
  )

  sparse_screen_failed <-
    (!is.null(sparse_screen_obj) &&
      length(sparse_pass) == 1L &&
      is.na(sparse_pass)) ||
    sparse_screen_error

  sparse_species <- if (!is.null(model$sparse_species)) {
    as.character(model$sparse_species)
  } else if (
    !is.null(sparse_screen_obj) &&
      !is.null(sparse_screen_obj$sparse_species)
  ) {
    as.character(sparse_screen_obj$sparse_species)
  } else if (
    !is.null(sparse_screen_obj) &&
      !is.null(sparse_screen_obj$sparse_details)
  ) {
    sparse_details <- as.data.frame(
      sparse_screen_obj$sparse_details,
      stringsAsFactors = FALSE
    )

    if ("response" %in% names(sparse_details)) {
      as.character(sparse_details[["response"]])
    } else {
      character(0)
    }
  } else {
    character(0)
  }

  sparse_species <- sort(
    unique(
      sparse_species[
        !is.na(sparse_species) &
          nzchar(sparse_species)
      ]
    )
  )

  sparse_text <- if (sparse_screen_failed) {
    "SCREEN FAILED / NOT ASSESSED"
  } else if (length(sparse_species) > 0L) {
    shown_sparse_species <- utils::head(
      sparse_species,
      10L
    )

    sparse_suffix <- if (length(sparse_species) > 10L) {
      paste0(
        ", and ",
        length(sparse_species) - 10L,
        " additional species"
      )
    } else {
      ""
    }

    paste0(
      length(sparse_species),
      " FLAGGED [",
      paste(
        shown_sparse_species,
        collapse = ", "
      ),
      sparse_suffix,
      "]"
    )
  } else if (isTRUE(sparse_pass)) {
    "NONE FLAGGED"
  } else if (is.null(sparse_screen_obj)) {
    "NOT RUN / NOT AVAILABLE"
  } else {
    "STATUS UNKNOWN"
  }

  # ---------------------------------------------------------------------------
  # 7. Compose individual report lines
  # ---------------------------------------------------------------------------
  optimizer_text <- if (opt_converged) {
    "SUCCESS"
  } else {
    "FAILED"
  }

  gradient_text <- if (!gradient_assessed) {
    "NOT ASSESSED"
  } else {
    paste0(
      if (isTRUE(gradient_valid)) {
        "PASS"
      } else {
        "FAILED"
      },
      " (max |gradient| = ",
      if (is.finite(max_abs_gradient)) {
        signif(max_abs_gradient, 6)
      } else {
        "NA"
      },
      "; threshold = ",
      if (is.finite(gradient_threshold)) {
        gradient_threshold
      } else {
        "NA"
      },
      ")"
    )
  }

  se_text <- if (!sd_errors_requested) {
    "NOT REQUESTED / NOT COMPUTED"
  } else if (se_total == 0L) {
    "REQUESTED BUT NOT AVAILABLE"
  } else {
    paste0(
      if (se_finite == se_total) {
        "FINITE"
      } else {
        "NON-FINITE VALUES DETECTED"
      },
      " (",
      se_finite,
      "/",
      se_total,
      ")"
    )
  }

  covariance_text <- if (!sd_errors_requested) {
    "NOT REQUESTED / NOT COMPUTED"
  } else if (isTRUE(se_valid)) {
    "PASS"
  } else if (length(covariance_failure_reasons) > 0L) {
    paste0(
      "FAILED: ",
      paste(
        covariance_failure_reasons,
        collapse = "; "
      )
    )
  } else {
    "FAILED FOR AN UNCLASSIFIED REASON"
  }

  # ---------------------------------------------------------------------------
  # 8. Overall status interpretation
  # ---------------------------------------------------------------------------
  overall_text <- switch(
    model_status,

    "valid" = paste0(
      "VALID FOR STANDARD-ERROR-BASED INFERENCE, SUBJECT TO ",
      "RESIDUAL, SPARSITY, AND MODEL-ASSUMPTION DIAGNOSTICS"
    ),

    "valid_fit_no_standard_errors" = paste0(
      "VALID CONVERGED FIT FOR EXPLORATION / MODEL SELECTION;\n",
      "                             STANDARD-ERROR-BASED INFERENCE NOT ASSESSED"
    ),

    "suspicious_optimization_incomplete" = paste0(
      "OPTIMIZER REPORTED SUCCESS, BUT THE GRADIENT\n",
      "                             DID NOT SATISFY THE CONVERGENCE THRESHOLD;\n",
      "                             DO NOT USE FOR MODEL SELECTION OR INFERENCE"
    ),

    "suspicious_optimization_unassessed" = paste0(
      "OPTIMIZER REPORTED SUCCESS, BUT THE FINAL GRADIENT\n",
      "                             COULD NOT BE ASSESSED;\n",
      "                             MANUAL REVIEW IS REQUIRED BEFORE MODEL SELECTION OR INFERENCE"
    ),

    "suspicious_inference_unreliable" = paste0(
      "POINT ESTIMATION CONVERGED, BUT VARIANCE-COVARIANCE /\n",
      "                             STANDARD-ERROR VALIDATION FAILED;\n",
      "                             DO NOT CITE WALD P-VALUES OR CONFIDENCE INTERVALS"
    ),

    "suspicious_reference_sensitivity" = paste0(
      "NUMERICAL DIAGNOSTICS CHANGED UNDER ALTERNATIVE\n",
      "                             FACTOR REFERENCE-LEVEL PARAMETERIZATIONS;\n",
      "                             INTERPRET MODEL STABILITY CAUTIOUSLY"
    ),

    "suspicious_reference_check_incomplete" = paste0(
      "ONE OR MORE ALTERNATIVE REFERENCE-LEVEL FITS FAILED;\n",
      "                             REFERENCE STABILITY COULD NOT BE FULLY ASSESSED"
    ),

    "suspicious_extreme_coefficients" = paste0(
      "EXTREME FIXED-EFFECT COEFFICIENTS DETECTED;\n",
      "                             INSPECT SPARSITY AND IDENTIFIABILITY BEFORE USE"
    ),

    "invalid" = paste0(
      "INVALID OPTIMIZATION; DO NOT USE FOR ESTIMATION,\n",
      "                             MODEL SELECTION, OR INFERENCE"
    ),

    "aborted" = "FIT ABORTED BEFORE MODEL ESTIMATION",

    "UNKNOWN / UNCLASSIFIED MODEL STATUS"
  )

  # ---------------------------------------------------------------------------
  # 9. Print one complete report
  # ---------------------------------------------------------------------------
  report_text <- paste0(
    "\n--------------------------------------------------\n",
    "MODEL CONVERGENCE REPORT: '",
    model_name,
    "'\n",
    "  1. Optimizer Exit Status : ",
    optimizer_text,
    "\n",
    "  2. Gradient Check        : ",
    gradient_text,
    "\n",
    "  3. Fixed-effect SEs      : ",
    se_text,
    "\n",
    "  4. Covariance Validation : ",
    covariance_text,
    "\n",
    "  5. Extreme Coefficients  : ",
    extreme_text,
    "\n",
    "  6. Reference Stability   : ",
    reference_text,
    "\n",
    "  7. Model Status Code     : ",
    model_status,
    "\n",
    "  8. Sparse Responses      : ",
    sparse_text,
    "\n",
    "  -> OVERALL STATUS        : ",
    overall_text,
    "\n",
    "--------------------------------------------------\n"
  )

  if (
    model_status %in%
      c(
        "valid",
        "valid_fit_no_standard_errors"
      )
  ) {
    message(report_text)
  } else {
    warning(
      report_text,
      call. = FALSE,
      immediate. = TRUE
    )
  }

  invisible(model)
}
