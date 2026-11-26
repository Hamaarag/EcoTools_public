#' Detect Empty Factorial Cells Requested by a Formula
#'
#' Identifies combinations of observed categorical predictor levels that are
#' absent from the data but are requested by interaction terms in a model
#' formula.
#'
#' Only interaction terms composed entirely of direct factor or character
#' columns are assessed. Transformed terms and interactions involving numeric
#' predictors are ignored.
#'
#' @param formula A model formula.
#' @param data A data frame containing the predictors.
#'
#' @return
#' A data frame with one row per empty combination. The `interaction_term`
#' column identifies the formula term, and `combination` gives a readable
#' representation of the missing factor-level combination.
#'
#' @keywords internal
detect_empty_factorial_cells <- function(
  formula,
  data
) {
  if (!inherits(formula, "formula")) {
    formula <- tryCatch(
      stats::as.formula(
        formula
      ),
      error = function(e) {
        stop(
          "'formula' must be a valid formula.",
          call. = FALSE
        )
      }
    )
  }

  data_use <- as.data.frame(
    data,
    stringsAsFactors = FALSE
  )

  predictor_terms <- stats::delete.response(
    stats::terms(
      formula,
      data = data_use
    )
  )

  term_labels <- attr(
    predictor_terms,
    "term.labels"
  )

  if (
    is.null(term_labels) ||
      length(term_labels) == 0L
  ) {
    return(
      data.frame(
        interaction_term = character(0),
        combination = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  result_rows <- list()
  result_index <- 0L

  for (term_label in term_labels) {
    term_formula <- tryCatch(
      stats::as.formula(
        paste0(
          "~ ",
          term_label
        ),
        env = environment(formula)
      ),
      error = function(e) {
        NULL
      }
    )

    if (is.null(term_formula)) {
      next
    }

    term_variables <- all.vars(
      term_formula
    )

    term_variables <- unique(
      term_variables
    )

    if (length(term_variables) < 2L) {
      next
    }

    # Only diagnose direct columns. This avoids misrepresenting transformed
    # terms such as factor(x), cut(x), or interactions involving expressions.
    if (!all(term_variables %in% names(data_use))) {
      next
    }

    categorical_variables <- vapply(
      term_variables,
      function(variable_name) {
        variable_values <- data_use[[variable_name]]

        is.factor(variable_values) ||
          is.character(variable_values) ||
          is.logical(variable_values)
      },
      logical(1)
    )

    if (!all(categorical_variables)) {
      next
    }

    observed_level_lists <- lapply(
      term_variables,
      function(variable_name) {
        variable_values <- data_use[[variable_name]]

        if (is.factor(variable_values)) {
          observed_values <- unique(
            as.character(
              variable_values[
                !is.na(variable_values)
              ]
            )
          )

          factor_levels <- levels(
            variable_values
          )

          factor_levels[
            factor_levels %in% observed_values
          ]
        } else {
          sort(
            unique(
              as.character(
                variable_values[
                  !is.na(variable_values)
                ]
              )
            )
          )
        }
      }
    )

    names(observed_level_lists) <- term_variables

    if (
      any(
        lengths(observed_level_lists) == 0L
      )
    ) {
      next
    }

    possible_combinations <- expand.grid(
      observed_level_lists,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    observed_combinations <- unique(
      data_use[
        stats::complete.cases(
          data_use[
            term_variables
          ]
        ),
        term_variables,
        drop = FALSE
      ]
    )

    for (variable_name in term_variables) {
      observed_combinations[[variable_name]] <- as.character(
        observed_combinations[[variable_name]]
      )

      possible_combinations[[variable_name]] <- as.character(
        possible_combinations[[variable_name]]
      )
    }

    make_combination_key <- function(
      data_object
    ) {
      if (nrow(data_object) == 0L) {
        return(
          character(0)
        )
      }

      do.call(
        paste,
        c(
          data_object,
          sep = "\r"
        )
      )
    }

    possible_keys <- make_combination_key(
      possible_combinations
    )

    observed_keys <- make_combination_key(
      observed_combinations
    )

    empty_rows <- !possible_keys %in% observed_keys

    if (!any(empty_rows)) {
      next
    }

    empty_combinations <- possible_combinations[
      empty_rows,
      ,
      drop = FALSE
    ]

    for (combination_index in seq_len(
      nrow(empty_combinations)
    )) {
      result_index <- result_index + 1L

      combination_values <- vapply(
        term_variables,
        function(variable_name) {
          paste0(
            variable_name,
            " = ",
            empty_combinations[
              combination_index,
              variable_name,
              drop = TRUE
            ]
          )
        },
        character(1)
      )

      result_rows[[result_index]] <- data.frame(
        interaction_term = term_label,
        combination = paste(
          combination_values,
          collapse = ", "
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(result_rows) == 0L) {
    return(
      data.frame(
        interaction_term = character(0),
        combination = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  result <- do.call(
    rbind,
    result_rows
  )

  row.names(result) <- NULL

  unique(
    result
  )
}
