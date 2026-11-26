detect_formula_factors <- function(
  X,
  formula
) {
  vars_in_formula <-
    all.vars(
      stats::delete.response(
        stats::terms(formula)
      )
    )

  factor_names <- names(X)[
    vapply(
      X,
      function(x) {
        is.factor(x) ||
          is.character(x)
      },
      logical(1)
    )
  ]

  intersect(
    factor_names,
    vars_in_formula
  )
}
