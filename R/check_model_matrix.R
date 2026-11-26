# ==============================================================================
# Fixed-effect design-matrix validation
# ==============================================================================

#' Check a Fixed-Effect Model Matrix
#'
#' @description
#' Constructs and validates the fixed-effect design matrix associated with a
#' model formula and predictor data.
#'
#' Mixed-model terms such as `(1 | group)` are permitted. Random-effect
#' terms are validated separately and removed before constructing the
#' fixed-effect design matrix. Every random-effect grouping expression must
#' reference at least one variable present in `data`.
#'
#' The function checks for:
#'
#' \enumerate{
#'   \item factor levels that are defined but absent from the supplied data;
#'   \item design-matrix columns containing only zeroes;
#'   \item exact rank deficiency among fixed-effect columns;
#'   \item non-finite values in the resulting design matrix.
#' }
#'
#' An unused factor level can create a structurally non-estimable coefficient.
#' For example, a retained factor level with no observations may produce a
#' coefficient estimate and standard error of zero, followed by undefined
#' Wald statistics and p-values.
#'
#' @param formula
#' A model formula describing the fixed-effect structure. The response, if
#' present, is ignored when constructing the design matrix.
#'
#' @param data
#' A data frame, data table, matrix, or compatible object containing the
#' predictors referenced by `formula`.
#'
#' @param action
#' Action to take when a design problem is detected:
#'
#' \itemize{
#'   \item `"error"` stops execution and reports the problem;
#'   \item `"warn"` emits a warning and returns the unchanged data and design
#'   matrix;
#'   \item `"drop"` removes unused levels from factor columns in a local copy of
#'   `data`, rebuilds the design matrix, and stops if structural problems remain.
#' }
#'
#' @param tolerance
#' Non-negative numeric tolerance used when identifying design-matrix columns
#' whose entries are effectively all zero. The default is
#' `sqrt(.Machine$double.eps)`.
#'
#' @return
#' A list with the following elements:
#'
#' \describe{
#'   \item{pass}{
#'     Logical indicating whether the final design matrix passed all checks.
#'   }
#'   \item{data}{
#'     The predictor data used to construct the final design matrix. When
#'     `action = "drop"`, this is a copy with unused factor levels removed.
#'   }
#'   \item{model_matrix}{
#'     The final fixed-effect model matrix.
#'   }
#'    \item{fixed_formula}{
#'     Formula remaining after valid random-effect terms have been removed.
#'     This formula is used to construct the fixed-effect design matrix.
#'   }
#'   \item{raw_random_effect_terms}{
#'     Character vector containing bar expressions detected directly in the
#'     original formula before random-effect syntax is normalized.
#'   }
#'   \item{random_effect_terms}{
#'     Character vector containing the random-effect terms detected in the
#'     original formula.
#'   }
#'   \item{random_effect_grouping_variables}{
#'     Character vector containing grouping variables referenced by the
#'     random-effect terms.
#'   }
#'   \item{rank}{
#'     The numerical rank of the final model matrix.
#'   }
#'   \item{n_columns}{
#'     Number of columns in the final model matrix.
#'   }
#'   \item{zero_columns}{
#'     Names of design-matrix columns containing only zeroes within
#'     `tolerance`.
#'   }
#'   \item{dependent_columns}{
#'     Names of columns identified as linearly dependent by the QR
#'     decomposition.
#'   }
#'   \item{unused_levels}{
#'     A named list containing unused levels for each affected factor.
#'   }
#'   \item{nonfinite_columns}{
#'     Names of model-matrix columns containing non-finite values.
#'   }
#'   \item{data_changed}{
#'     Logical indicating whether unused factor levels were dropped.
#'   }
#'   \item{message}{
#'     A human-readable summary of the detected problem, or `NULL` when the
#'     final matrix passed.
#'   }
#' }
#'
#' @details
#' Removing unused factor levels only addresses levels that remain defined after
#' subsetting but have no observations in the supplied data. It does not repair
#' general collinearity, missing interaction support, duplicated predictors, or
#' other causes of rank deficiency.
#'
#' When `action = "drop"`, the function changes only its local copy of `data`.
#' The caller must explicitly use the returned `data` object when fitting the
#' model.
#'
#' For reproducible analyses, `"error"` is the recommended default because it
#' requires factor cleanup and reference-level choices to remain explicit in the
#' analysis script.
#'
#' @examples
#' example_data <- data.frame(
#'   year = 1:6,
#'   group = factor(
#'     c("A", "A", "B", "B", "A", "B"),
#'     levels = c("A", "B", "unused")
#'   )
#' )
#'
#' \dontrun{
#' check_model_matrix(
#'   formula = ~ year + group,
#'   data = example_data,
#'   action = "error"
#' )
#' }
#'
#' cleaned <- check_model_matrix(
#'   formula = ~ year + group,
#'   data = example_data,
#'   action = "drop"
#' )
#'
#' levels(cleaned$data$group)
#' colnames(cleaned$model_matrix)
#'
#' @export
check_model_matrix <- function(
  formula,
  data,
  action = c("error", "warn", "drop"),
  tolerance = sqrt(.Machine$double.eps)
) {
  action <- match.arg(action)

  if (!inherits(formula, "formula")) {
    formula <- tryCatch(
      stats::as.formula(formula),
      error = function(e) {
        stop(
          "'formula' must be a valid formula.",
          call. = FALSE
        )
      }
    )
  }

  if (is.null(data)) {
    stop(
      "'data' must not be NULL.",
      call. = FALSE
    )
  }

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

  data_use <- as.data.frame(
    data,
    stringsAsFactors = FALSE
  )

  if (nrow(data_use) == 0L) {
    stop(
      "'data' contains no rows.",
      call. = FALSE
    )
  }

  if (ncol(data_use) == 0L) {
    stop(
      "'data' contains no columns.",
      call. = FALSE
    )
  }

  factor_columns <- names(data_use)[vapply(
    data_use,
    is.factor,
    logical(1)
  )]

  find_unused_levels <- function(data_object) {
    out <- lapply(
      factor_columns,
      function(variable_name) {
        factor_values <- data_object[[variable_name]]

        observed_values <- unique(
          as.character(
            factor_values[!is.na(factor_values)]
          )
        )

        setdiff(
          levels(factor_values),
          observed_values
        )
      }
    )

    names(out) <- factor_columns

    out[
      lengths(out) > 0L
    ]
  }

  # ---------------------------------------------------------------------------
  # Separate and validate mixed-model terms before fixed-effect validation.
  #
  # stats::model.matrix() does not understand mixed-model bar expressions.
  # reformulas::findbars() identifies "|" and expands "||" expressions.
  # reformulas::nobars() then removes those terms before fixed-effect matrix
  # construction.
  # ---------------------------------------------------------------------------
  if (!requireNamespace("reformulas", quietly = TRUE)) {
    stop(
      paste0(
        "Package 'reformulas' is required to validate formulas ",
        "that may contain random-effect terms."
      ),
      call. = FALSE
    )
  }

  # Inspect the original formula before reformulas transforms "|" or "||".
  #
  # Some formula parsers may simplify or discard constant grouping
  # expressions such as (1 | 1). Those expressions must be rejected before
  # random-effect terms are normalized.
  collect_raw_bar_terms <- function(expression) {
    collected_terms <- list()

    recurse <- function(node) {
      if (!is.call(node)) {
        return(
          invisible(NULL)
        )
      }

      operator_name <- as.character(
        node[[1L]]
      )

      if (
        length(operator_name) == 1L &&
          operator_name %in%
            c(
              "|",
              "||"
            )
      ) {
        collected_terms[[
          length(collected_terms) + 1L
        ]] <<- node
      }

      if (length(node) >= 2L) {
        for (element_index in seq.int(
          from = 2L,
          to = length(node)
        )) {
          recurse(
            node[[element_index]]
          )
        }
      }

      invisible(NULL)
    }

    recurse(
      expression
    )

    collected_terms
  }

  raw_bar_terms <- collect_raw_bar_terms(
    formula[[length(formula)]]
  )

  if (length(raw_bar_terms) > 0L) {
    raw_grouping_expressions <- lapply(
      raw_bar_terms,
      function(bar_term) {
        if (length(bar_term) < 3L) {
          stop(
            paste0(
              "Invalid random-effect term in 'formula': ",
              paste(
                deparse(
                  bar_term,
                  width.cutoff = 500L
                ),
                collapse = ""
              ),
              "."
            ),
            call. = FALSE
          )
        }

        bar_term[[3L]]
      }
    )

    raw_grouping_variables <- lapply(
      raw_grouping_expressions,
      all.vars
    )

    raw_constant_grouping <- vapply(
      seq_along(raw_grouping_expressions),
      function(term_index) {
        grouping_expression <- raw_grouping_expressions[[term_index]]

        grouping_variables <- raw_grouping_variables[[term_index]]

        length(grouping_variables) == 0L ||
          is.numeric(grouping_expression) ||
          is.logical(grouping_expression) ||
          is.character(grouping_expression)
      },
      logical(1)
    )

    if (any(raw_constant_grouping)) {
      invalid_terms <- vapply(
        raw_bar_terms[raw_constant_grouping],
        function(bar_term) {
          paste(
            deparse(
              bar_term,
              width.cutoff = 500L
            ),
            collapse = ""
          )
        },
        character(1)
      )

      stop(
        paste0(
          "Invalid random-effect grouping term(s): ",
          paste(
            invalid_terms,
            collapse = ", "
          ),
          ". A random-effect term must identify a grouping variable ",
          "present in 'data'. Expressions such as '(predictor | 1)' ",
          "or '(predictor || 1)' are invalid because 1 is an ",
          "intercept expression, not a grouping variable."
        ),
        call. = FALSE
      )
    }
  }

  random_effect_terms <- tryCatch(
    reformulas::findbars(
      formula
    ),
    error = function(e) {
      stop(
        paste0(
          "Could not identify random-effect terms in 'formula': ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  random_effect_grouping_variables <- character(0)

  if (length(random_effect_terms) > 0L) {
    grouping_expressions <- lapply(
      random_effect_terms,
      function(random_term) {
        if (
          !is.call(random_term) ||
            length(random_term) < 3L
        ) {
          stop(
            paste0(
              "Internal EcoTools formula-parsing error: an extracted ",
              "random-effect term does not contain both an effect and ",
              "a grouping expression."
            ),
            call. = FALSE
          )
        }
        random_term[[3L]]
      }
    )

    grouping_variables_by_term <- lapply(
      grouping_expressions,
      all.vars
    )

    constant_grouping_terms <- vapply(
      seq_along(grouping_expressions),
      function(term_index) {
        grouping_expression <- grouping_expressions[[term_index]]
        grouping_variables <- grouping_variables_by_term[[term_index]]
        length(grouping_variables) == 0L ||
          is.numeric(grouping_expression) ||
          is.logical(grouping_expression) ||
          is.character(grouping_expression)
      },
      logical(1)
    )

    if (any(constant_grouping_terms)) {
      invalid_terms <- vapply(
        random_effect_terms[
          constant_grouping_terms
        ],
        function(random_term) {
          paste(
            deparse(
              random_term,
              width.cutoff = 500L
            ),
            collapse = ""
          )
        },
        character(1)
      )
      stop(
        paste0(
          "Invalid random-effect grouping term(s): ",
          paste(
            invalid_terms,
            collapse = ", "
          ),
          ". A random-effect term must identify a grouping variable ",
          "present in 'data'. Expressions such as '(predictor | 1)' ",
          "or '(predictor || 1)' are invalid because 1 is an ",
          "intercept expression, not a grouping variable."
        ),
        call. = FALSE
      )
    }

    random_effect_grouping_variables <- unique(
      unlist(
        grouping_variables_by_term,
        use.names = FALSE
      )
    )
    missing_grouping_variables <- setdiff(
      random_effect_grouping_variables,
      names(data_use)
    )

    if (length(missing_grouping_variables) > 0L) {
      stop(
        paste0(
          "Random-effect grouping variable(s) referenced by 'formula' ",
          "are missing from 'data': ",
          paste(
            missing_grouping_variables,
            collapse = ", "
          ),
          "."
        ),
        call. = FALSE
      )
    }
  }

  fixed_formula <- tryCatch(
    reformulas::nobars(
      formula
    ),
    error = function(e) {
      stop(
        paste0(
          "Could not remove random-effect terms before constructing ",
          "the fixed-effect model matrix: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  random_effect_term_labels <- if (length(random_effect_terms) > 0L) {
    vapply(
      random_effect_terms,
      function(random_term) {
        paste(
          deparse(
            random_term,
            width.cutoff = 500L
          ),
          collapse = ""
        )
      },
      character(1)
    )
  } else {
    character(0)
  }

  empty_factorial_cells_before <- detect_empty_factorial_cells(
    formula = fixed_formula,
    data = data_use
  )

  unused_levels_before <- find_unused_levels(data_use)
  data_changed <- FALSE

  if (
    action == "drop" &&
      length(unused_levels_before) > 0L
  ) {
    data_use <- droplevels(data_use)
    factor_columns <- names(data_use)[vapply(
      data_use,
      is.factor,
      logical(1)
    )]
    data_changed <- TRUE
  }

  formula_terms <- tryCatch(
    stats::terms(
      fixed_formula,
      data = data_use
    ),
    error = function(e) {
      stop(
        paste0(
          "Could not construct fixed-effect formula terms: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  predictor_terms <- stats::delete.response(
    formula_terms
  )

  empty_factorial_cells <- detect_empty_factorial_cells(
    formula = fixed_formula,
    data = data_use
  )

  model_matrix <- tryCatch(
    stats::model.matrix(
      predictor_terms,
      data = data_use,
      na.action = stats::na.pass
    ),
    error = function(e) {
      stop(
        paste0(
          "Could not construct the fixed-effect model matrix: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  if (nrow(model_matrix) == 0L) {
    stop(
      "The constructed model matrix contains no rows.",
      call. = FALSE
    )
  }

  if (ncol(model_matrix) == 0L) {
    stop(
      "The constructed model matrix contains no columns.",
      call. = FALSE
    )
  }

  nonfinite_columns <- colnames(model_matrix)[
    vapply(
      seq_len(ncol(model_matrix)),
      function(column_index) {
        any(!is.finite(model_matrix[, column_index]))
      },
      logical(1)
    )
  ]

  zero_columns <- colnames(model_matrix)[
    vapply(
      seq_len(ncol(model_matrix)),
      function(column_index) {
        column_values <- model_matrix[, column_index]

        all(
          is.finite(column_values) &
            abs(column_values) <= tolerance
        )
      },
      logical(1)
    )
  ]

  qr_result <- qr(
    model_matrix,
    tol = tolerance
  )

  design_rank <- qr_result$rank
  design_columns <- ncol(model_matrix)

  dependent_columns <- if (design_rank < design_columns) {
    pivot_order <- qr_result$pivot

    colnames(model_matrix)[
      pivot_order[
        seq.int(
          from = design_rank + 1L,
          to = design_columns
        )
      ]
    ]
  } else {
    character(0)
  }

  unused_levels_after <- {
    refreshed_factor_columns <- names(data_use)[vapply(
      data_use,
      is.factor,
      logical(1)
    )]

    out <- lapply(
      refreshed_factor_columns,
      function(variable_name) {
        factor_values <- data_use[[variable_name]]

        observed_values <- unique(
          as.character(
            factor_values[!is.na(factor_values)]
          )
        )

        setdiff(
          levels(factor_values),
          observed_values
        )
      }
    )

    names(out) <- refreshed_factor_columns

    out[
      lengths(out) > 0L
    ]
  }

  pass <-
    length(zero_columns) == 0L &&
    length(dependent_columns) == 0L &&
    length(nonfinite_columns) == 0L &&
    length(unused_levels_after) == 0L

  problem_parts <- character(0)

  if (length(unused_levels_after) > 0L) {
    unused_text <- vapply(
      names(unused_levels_after),
      function(variable_name) {
        paste0(
          variable_name,
          " = [",
          paste(
            unused_levels_after[[variable_name]],
            collapse = ", "
          ),
          "]"
        )
      },
      character(1)
    )

    problem_parts <- c(
      problem_parts,
      paste0(
        "unused factor levels: ",
        paste(unused_text, collapse = "; ")
      )
    )
  }

  if (length(zero_columns) > 0L) {
    problem_parts <- c(
      problem_parts,
      paste0(
        "all-zero model-matrix columns: ",
        paste(zero_columns, collapse = ", ")
      )
    )
  }

  if (
    length(dependent_columns) > 0L &&
      nrow(empty_factorial_cells) > 0L
  ) {
    empty_cell_text <- paste0(
      empty_factorial_cells$interaction_term,
      " has no observations for [",
      empty_factorial_cells$combination,
      "]"
    )

    problem_parts <- c(
      problem_parts,
      paste0(
        "incomplete factorial support: ",
        paste(
          empty_cell_text,
          collapse = "; "
        )
      )
    )
  }

  if (length(dependent_columns) > 0L) {
    problem_parts <- c(
      problem_parts,
      paste0(
        "QR-selected linearly dependent model-matrix column(s): ",
        paste(
          dependent_columns,
          collapse = ", "
        )
      )
    )
  }

  if (length(nonfinite_columns) > 0L) {
    problem_parts <- c(
      problem_parts,
      paste0(
        "non-finite model-matrix columns: ",
        paste(nonfinite_columns, collapse = ", ")
      )
    )
  }

  result_message <- if (length(problem_parts) > 0L) {
    paste0(
      "FIXED-EFFECT DESIGN-MATRIX CHECK FAILED: ",
      paste(problem_parts, collapse = "; "),
      "."
    )
  } else {
    NULL
  }

  result <- list(
    pass = pass,
    data = data_use,
    model_matrix = model_matrix,
    fixed_formula = fixed_formula,
    raw_random_effect_terms = if (length(raw_bar_terms) > 0L) {
      vapply(
        raw_bar_terms,
        function(bar_term) {
          paste(
            deparse(
              bar_term,
              width.cutoff = 500L
            ),
            collapse = ""
          )
        },
        character(1)
      )
    } else {
      character(0)
    },
    random_effect_terms = random_effect_term_labels,
    random_effect_grouping_variables = random_effect_grouping_variables,
    rank = design_rank,
    n_columns = design_columns,
    zero_columns = zero_columns,
    dependent_columns = dependent_columns,
    unused_levels = unused_levels_after,
    unused_levels_removed = if (data_changed) {
      unused_levels_before
    } else {
      list()
    },
    empty_factorial_cells = empty_factorial_cells,
    empty_factorial_cells_before_cleanup = empty_factorial_cells_before,
    nonfinite_columns = nonfinite_columns,
    data_changed = data_changed,
    action = action,
    tolerance = tolerance,
    message = result_message
  )

  if (!pass) {
    if (action == "error" || action == "drop") {
      design_repair_guidance <- if (
        nrow(empty_factorial_cells) > 0L &&
          length(dependent_columns) > 0L
      ) {
        paste0(
          " The formula requests one or more categorical interaction ",
          "combinations for which no observations are available. The interaction ",
          "is therefore not fully estimable. Consider simplifying the interaction, ",
          "restricting the inferential scope to supported combinations, fitting an ",
          "observed-cell-means parameterization, or defining scientifically ",
          "justified estimable contrasts. The QR-selected dependent column need ",
          "not correspond directly to the empty cell because its identity depends ",
          "on factor contrasts, reference levels, and model-matrix column order."
        )
      } else if (action == "drop") {
        paste0(
          " Unused factor levels were dropped, but the remaining ",
          "design matrix is still invalid."
        )
      } else {
        paste0(
          " Remove unused factor levels or repair the predictor ",
          "structure before fitting."
        )
      }

      stop(
        paste0(
          result_message,
          design_repair_guidance
        ),
        call. = FALSE
      )
    }

    warning(
      paste0(
        result_message,
        " Fitting will continue because action = 'warn'."
      ),
      call. = FALSE
    )
  } else if (data_changed) {
    removed_text <- vapply(
      names(unused_levels_before),
      function(variable_name) {
        paste0(
          variable_name,
          " = [",
          paste(
            unused_levels_before[[variable_name]],
            collapse = ", "
          ),
          "]"
        )
      },
      character(1)
    )

    message(
      paste0(
        "Removed unused factor levels before fitting: ",
        paste(removed_text, collapse = "; "),
        "."
      )
    )
  }

  result
}
