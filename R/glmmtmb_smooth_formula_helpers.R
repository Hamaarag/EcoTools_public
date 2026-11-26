#' Normalize mgcv Smooth Calls for glmmTMB
#'
#' Converts namespaced `mgcv::s()` calls to the unqualified `s()` form expected
#' by glmmTMB's experimental smooth-term interface.
#'
#' @param formula_object A formula or a character representation of a formula.
#'
#' @return A formula with `mgcv::s()` call heads replaced by `s()`.
#'
#' @keywords internal
normalize_glmmtmb_smooth_formula <- function(
  formula_object
) {
  formula_result <- if (inherits(formula_object, "formula")) {
    formula_object
  } else {
    stats::as.formula(
      formula_object,
      env = parent.frame()
    )
  }

  original_environment <- environment(
    formula_result
  )

  rewrite_expression <- function(
    expression_object
  ) {
    if (!is.call(expression_object)) {
      return(expression_object)
    }

    expression_list <- as.list(
      expression_object
    )

    expression_list <- lapply(
      expression_list,
      rewrite_expression
    )

    rewritten_call <- as.call(
      expression_list
    )

    if (
      length(rewritten_call) >= 1L &&
        identical(
          rewritten_call[[1L]],
          quote(mgcv::s)
        )
    ) {
      rewritten_call[[1L]] <- as.name(
        "s"
      )
    }

    rewritten_call
  }

  normalized_formula <- rewrite_expression(
    formula_result
  )

  class(normalized_formula) <- "formula"

  environment(normalized_formula) <-
    original_environment

  normalized_formula
}

#' Detect by Arguments in glmmTMB Smooth Terms
#'
#' Detects `s()` calls containing an explicitly supplied `by` argument.
#' glmmTMB's current experimental smooth implementation does not support these
#' terms.
#'
#' @param formula_object A formula or character representation of a formula.
#'
#' @return A logical scalar.
#'
#' @keywords internal
formula_has_glmmtmb_by_smooth <- function(
  formula_object
) {
  formula_result <- if (inherits(formula_object, "formula")) {
    formula_object
  } else {
    stats::as.formula(
      formula_object,
      env = parent.frame()
    )
  }

  contains_by_smooth <- function(
    expression_object
  ) {
    if (!is.call(expression_object)) {
      return(FALSE)
    }

    call_head <- expression_object[[1L]]

    is_smooth_call <-
      identical(
        call_head,
        as.name("s")
      ) ||
      identical(
        call_head,
        quote(mgcv::s)
      )

    if (is_smooth_call) {
      call_arguments <- as.list(
        expression_object
      )[-1L]

      argument_names <- names(
        call_arguments
      )

      if (
        !is.null(argument_names) &&
          "by" %in% argument_names
      ) {
        return(TRUE)
      }
    }

    nested_expressions <- as.list(
      expression_object
    )[-1L]

    if (length(nested_expressions) == 0L) {
      return(FALSE)
    }

    any(
      vapply(
        nested_expressions,
        contains_by_smooth,
        logical(1)
      )
    )
  }

  contains_by_smooth(
    formula_result
  )
}

#' Detect mgcv Smooth Terms in a Formula
#'
#' @param formula_object A formula or character representation of a formula.
#'
#' @return A logical scalar.
#'
#' @keywords internal
formula_has_mgcv_smooth <- function(
  formula_object
) {
  formula_result <- if (inherits(formula_object, "formula")) {
    formula_object
  } else {
    stats::as.formula(
      formula_object,
      env = parent.frame()
    )
  }

  contains_smooth <- function(
    expression_object
  ) {
    if (!is.call(expression_object)) {
      return(FALSE)
    }

    call_head <- expression_object[[1L]]

    if (
      identical(
        call_head,
        as.name("s")
      ) ||
        identical(
          call_head,
          quote(mgcv::s)
        )
    ) {
      return(TRUE)
    }

    any(
      vapply(
        as.list(expression_object)[-1L],
        contains_smooth,
        logical(1)
      )
    )
  }

  contains_smooth(
    formula_result
  )
}


#' Remove Smooth Terms from a Parametric Design-Check Formula
#'
#' Removes terms containing `s()` from a formula before constructing the
#' ordinary parametric model matrix. glmmTMB performs the smooth-basis and
#' penalty construction later during its own model preprocessing.
#'
#' @param formula_object A formula.
#'
#' @return A formula containing only the ordinary parametric fixed effects.
#'
#' @keywords internal
strip_mgcv_smooth_terms <- function(
  formula_object
) {
  if (!inherits(formula_object, "formula")) {
    formula_object <- stats::as.formula(
      formula_object,
      env = parent.frame()
    )
  }

  formula_terms <- stats::terms(
    formula_object
  )

  term_labels <- attr(
    formula_terms,
    "term.labels"
  )

  smooth_term_pattern <- paste0(
    "(^|:|\\*)",
    "\\s*(mgcv::)?s\\s*\\("
  )

  retained_terms <- term_labels[
    !grepl(
      smooth_term_pattern,
      term_labels,
      perl = TRUE
    )
  ]

  formula_response <- if (length(formula_object) == 3L) {
    paste(
      deparse(
        formula_object[[2L]]
      ),
      collapse = ""
    )
  } else {
    NULL
  }

  parametric_formula <- stats::reformulate(
    termlabels = retained_terms,
    response = formula_response,
    intercept = attr(
      formula_terms,
      "intercept"
    ) ==
      1L,
    env = environment(
      formula_object
    )
  )

  parametric_formula
}
