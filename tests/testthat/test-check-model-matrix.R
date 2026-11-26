test_that("empty factorial cells are reported with dependent columns", {
  test_data <- data.frame(
    year_fct = factor(
      c(
        "2013",
        "2013",
        "2015",
        "2015",
        "2015"
      ),
      levels = c(
        "2013",
        "2015"
      )
    ),
    habitat = factor(
      c(
        "Loess",
        "Planting",
        "Agriculture",
        "Loess",
        "Planting"
      ),
      levels = c(
        "Agriculture",
        "Loess",
        "Planting"
      )
    )
  )

  expect_error(
    check_model_matrix(
      formula = ~ year_fct * habitat,
      data = test_data,
      action = "error"
    ),
    regexp = paste0(
      "incomplete factorial support:.*",
      "year_fct = 2013, habitat = Agriculture"
    )
  )
})

test_that("additive models are not rejected solely because a factorial cell is empty", {
  test_data <- data.frame(
    year_fct = factor(
      c(
        "2013",
        "2013",
        "2015",
        "2015",
        "2015"
      ),
      levels = c(
        "2013",
        "2015"
      )
    ),
    habitat = factor(
      c(
        "Loess",
        "Planting",
        "Agriculture",
        "Loess",
        "Planting"
      ),
      levels = c(
        "Agriculture",
        "Loess",
        "Planting"
      )
    )
  )

  result <- check_model_matrix(
    formula = ~ year_fct + habitat,
    data = test_data,
    action = "error"
  )

  expect_true(
    result$pass
  )

  expect_equal(
    nrow(
      result$empty_factorial_cells
    ),
    0L
  )
})

test_that("empty-cell diagnostics distinguish the empty cell from QR aliasing", {
  test_data <- data.frame(
    year_fct = factor(
      c(
        "2013",
        "2013",
        "2015",
        "2015",
        "2015",
        "2023",
        "2023",
        "2023"
      ),
      levels = c(
        "2013",
        "2015",
        "2023"
      )
    ),
    habitat = factor(
      c(
        "Loess",
        "Planting",
        "Agriculture",
        "Loess",
        "Planting",
        "Agriculture",
        "Loess",
        "Planting"
      ),
      levels = c(
        "Agriculture",
        "Loess",
        "Planting"
      )
    )
  )

  captured_warning <- NULL

  result <- withCallingHandlers(
    check_model_matrix(
      formula = ~ year_fct * habitat,
      data = test_data,
      action = "warn"
    ),
    warning = function(w) {
      captured_warning <<- conditionMessage(w)
      invokeRestart(
        "muffleWarning"
      )
    }
  )

  expect_false(
    result$pass
  )

  expect_true(
    any(
      result$empty_factorial_cells$combination ==
        "year_fct = 2013, habitat = Agriculture"
    )
  )

  expect_true(
    length(
      result$dependent_columns
    ) >
      0L
  )

  expect_match(
    captured_warning,
    "incomplete factorial support"
  )

  expect_match(
    captured_warning,
    "QR-selected"
  )
})

test_that("check_model_matrix handles a valid random-intercept term", {
  X <- data.frame(
    yr_scaled = seq(
      -1,
      1,
      length.out = 12L
    ),
    habitat = factor(
      rep(
        c(
          "A",
          "B"
        ),
        each = 6L
      )
    ),
    plot_code = factor(
      rep(
        seq_len(6L),
        each = 2L
      )
    )
  )

  result <- check_model_matrix(
    formula = ~ yr_scaled +
      habitat +
      (1 | plot_code),
    data = X,
    action = "error"
  )

  expect_true(
    result$pass
  )

  expect_identical(
    result$random_effect_grouping_variables,
    "plot_code"
  )

  expect_true(
    any(
      grepl(
        "plot_code",
        result$random_effect_terms,
        fixed = TRUE
      )
    )
  )

  expect_false(
    any(
      grepl(
        "plot_code",
        colnames(
          result$model_matrix
        ),
        fixed = TRUE
      )
    )
  )
})

test_that("check_model_matrix handles valid double-bar random slopes", {
  X <- data.frame(
    yr_scaled = seq(
      -1,
      1,
      length.out = 12L
    ),
    habitat = factor(
      rep(
        c(
          "A",
          "B"
        ),
        each = 6L
      )
    ),
    plot_code = factor(
      rep(
        seq_len(6L),
        each = 2L
      )
    )
  )

  result <- check_model_matrix(
    formula = ~ habitat +
      (0 + yr_scaled || plot_code),
    data = X,
    action = "error"
  )

  expect_true(
    result$pass
  )

  expect_identical(
    result$random_effect_grouping_variables,
    "plot_code"
  )

  expect_true(
    "habitatB" %in%
      colnames(
        result$model_matrix
      )
  )
})

test_that("check_model_matrix rejects a constant random-effect group", {
  X <- data.frame(
    yr_scaled = seq(-1, 1, length.out = 12L),
    habitat = factor(rep(c("A", "B"), each = 6L))
  )
  expect_error(
    check_model_matrix(
      formula = ~ yr_scaled + habitat + (1 | 1),
      data = X,
      action = "error"
    ),
    regexp = "not a grouping variable"
  )

  expect_error(
    check_model_matrix(
      formula = ~ yr_scaled + habitat + (1 | "all_species"),
      data = X,
      action = "error"
    ),
    regexp = "not a grouping variable"
  )

  expect_error(
    check_model_matrix(
      formula = ~ habitat + (0 + yr_scaled || 1),
      data = X,
      action = "error"
    ),
    regexp = "not a grouping variable"
  )
})

test_that("check_model_matrix rejects missing grouping variables", {
  X <- data.frame(
    yr_scaled = seq(
      -1,
      1,
      length.out = 12L
    ),
    habitat = factor(
      rep(
        c(
          "A",
          "B"
        ),
        each = 6L
      )
    )
  )

  expect_error(
    check_model_matrix(
      formula = ~ yr_scaled +
        habitat +
        (1 | plot_code),
      data = X,
      action = "error"
    ),
    regexp = "plot_code"
  )
})

test_that("check_model_matrix accepts nested grouping variables", {
  X <- data.frame(
    yr_scaled = seq(-1, 1, length.out = 12L),
    region = factor(rep(c("north", "south"), each = 6L)),
    plot_code = factor(rep(seq_len(6L), each = 2L))
  )

  result <- check_model_matrix(
    formula = ~ yr_scaled + (1 | region / plot_code),
    data = X,
    action = "error"
  )

  expect_true(
    result$pass
  )

  expect_true(
    all(
      c(
        "region",
        "plot_code"
      ) %in%
        result$random_effect_grouping_variables
    )
  )
})
