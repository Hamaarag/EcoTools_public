test_that("GLLVM response screen warns about all-zero rows", {
  skip_if_not_installed(
    "gllvm"
  )

  X <- data.frame(
    treatment = factor(
      c(
        "A",
        "A",
        "B",
        "B",
        "B",
        "A"
      )
    )
  )

  y <- matrix(
    c(
      1,
      0,
      0,
      0,
      2,
      1,
      0,
      3,
      1,
      1,
      0,
      0
    ),
    nrow = 6L,
    byrow = TRUE
  )

  colnames(
    y
  ) <- c(
    "species_a",
    "species_b"
  )

  fit <- NULL

  captured_warning_messages <- character()
  fit <- withCallingHandlers(
    fit_gllvm_cached(
      model_name = "test_gllvm_all_zero_rows",
      y = y,
      X = X,
      studyDesign = X,
      formula = ~treatment,
      family = rep(
        "poisson",
        2L
      ),
      num.lv = 0L,
      sd.errors = FALSE,
      control.start = list(
        n.init = 1L,
        jitter.var = 0,
        n.init.max = 1L
      ),
      force_refit = TRUE,
      reference_stability_check = "none",
      cache_dir = tempdir(),
      sparse_screen = FALSE,
      plot_diag = FALSE,
      num_cores = 1L
    ),
    warning = function(w) {
      captured_warning_messages <<- c(
        captured_warning_messages,
        conditionMessage(w)
      )
      invokeRestart(
        "muffleWarning"
      )
    }
  )
  expect_true(
    any(
      grepl(
        paste0(
          "GLLVM RESPONSE WARNING: ",
          "2 response row(s) contain only zeros: 2, 6"
        ),
        captured_warning_messages,
        fixed = TRUE
      )
    )
  )

  expect_type(
    fit,
    "list"
  )

  expect_identical(
    fit$diagnostics$response_screen$n_all_zero_rows,
    2L
  )

  expect_identical(
    fit$diagnostics$response_screen$all_zero_rows,
    c(
      2L,
      6L
    )
  )

  expect_true(
    length(captured_warning_messages) >= 1L
  )
})
