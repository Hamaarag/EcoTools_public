test_that("with_gllvm_tmb_threads validates its thread argument", {
  expect_error(
    EcoTools:::with_gllvm_tmb_threads(
      threads = 0L,
      code = TRUE
    ),
    "positive integer"
  )

  expect_error(
    EcoTools:::with_gllvm_tmb_threads(
      threads = NA_integer_,
      code = TRUE
    ),
    "positive integer"
  )

  expect_error(
    EcoTools:::with_gllvm_tmb_threads(
      threads = 1.5,
      code = TRUE
    ),
    "positive integer"
  )
})

test_that("with_gllvm_tmb_threads configures and restores the gllvm DLL", {
  skip_if_not_installed("gllvm")
  skip_if_not_installed("TMB")

  loadNamespace("gllvm")

  original_threads <- as.integer(
    TMB::openmp(
      DLL = "gllvm"
    )
  )

  on.exit(
    TMB::openmp(
      original_threads,
      DLL = "gllvm"
    ),
    add = TRUE
  )

  target_threads <- if (original_threads == 1L) {
    2L
  } else {
    1L
  }

  observed_inside <- EcoTools:::with_gllvm_tmb_threads(
    threads = target_threads,
    code = {
      as.integer(
        TMB::openmp(
          DLL = "gllvm"
        )
      )
    }
  )

  observed_after <- as.integer(
    TMB::openmp(
      DLL = "gllvm"
    )
  )

  expect_identical(
    observed_inside,
    target_threads
  )

  expect_identical(
    observed_after,
    original_threads
  )
})

test_that("with_gllvm_tmb_threads restores state after an error", {
  skip_if_not_installed("gllvm")
  skip_if_not_installed("TMB")

  loadNamespace("gllvm")

  original_threads <- as.integer(
    TMB::openmp(
      DLL = "gllvm"
    )
  )

  on.exit(
    TMB::openmp(
      original_threads,
      DLL = "gllvm"
    ),
    add = TRUE
  )

  target_threads <- if (original_threads == 1L) {
    2L
  } else {
    1L
  }

  expect_error(
    EcoTools:::with_gllvm_tmb_threads(
      threads = target_threads,
      code = {
        stop(
          "deliberate test error"
        )
      }
    ),
    "deliberate test error"
  )

  observed_after <- as.integer(
    TMB::openmp(
      DLL = "gllvm"
    )
  )

  expect_identical(
    observed_after,
    original_threads
  )
})

test_that("GLLVM fingerprint changes when requested TMB threads change", {
  skip_if_not_installed("digest")

  y <- matrix(
    c(
      1,
      2,
      3,
      4,
      5,
      6
    ),
    nrow = 3L,
    byrow = TRUE
  )

  X <- data.frame(
    treatment = factor(
      c(
        "A",
        "B",
        "A"
      )
    )
  )

  common_arguments <- list(
    y = y,
    X = X,
    studyDesign = X,
    formula = ~treatment,
    TR = NULL,
    row.eff = FALSE,
    family = "poisson",
    num.lv = 1L,
    sd.errors = FALSE,
    method = "VA",
    control = list(),
    control.start = list(),
    gradient_threshold = 0.01,
    coef_threshold = 5,
    sparse_screen = FALSE,
    sparse_action = "proceed",
    abort_on = "all",
    sparse_min_events_per_parameter = 5,
    sparse_min_total_events = 5,
    sparse_min_nonzero_rows = 3,
    sparse_event_basis = "nonzero",
    sparse_count_like_only = TRUE,
    design_check = TRUE,
    design_action = "error",
    design_tolerance = sqrt(.Machine$double.eps),
    reference_stability_check = "none",
    reference_stability_max_levels = 5L,
    dependency_versions = list(
      gllvm = "test",
      TMB = "test"
    ),
    extra_args = list()
  )

  one_thread_arguments <- common_arguments
  one_thread_arguments$tmb_threads <- 1L

  two_thread_arguments <- common_arguments
  two_thread_arguments$tmb_threads <- 2L

  one_thread_fingerprint <- do.call(
    EcoTools:::get_gllvm_fingerprint,
    one_thread_arguments
  )

  two_thread_fingerprint <- do.call(
    EcoTools:::get_gllvm_fingerprint,
    two_thread_arguments
  )

  expect_false(
    identical(
      one_thread_fingerprint,
      two_thread_fingerprint
    )
  )
})
