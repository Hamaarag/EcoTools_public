# E11 (2026-10-22): the family component of the analysis fingerprint must not
# depend on process state. Before the fix, `family$a` partial-matched the
# glmmTMB `aic` stub and its serialised closure environment changed after the
# first glmmTMB fit and after summary() in the same process.

fingerprint_spec <- function() {
  set.seed(1)
  d <- data.frame(
    y = stats::rnbinom(120, mu = 4, size = 2),
    g = factor(rep(1:3, 40)),
    s = factor(rep(1:12, each = 10))
  )
  list(
    data = d,
    model_formula = "y ~ g + (1 | s)",
    zi_formula = ~0,
    disp_formula = ~1,
    family = glmmTMB::nbinom2(),
    weights = NULL,
    offset = NULL,
    package = "glmmTMB",
    dependency_versions = list(
      glmmTMB = "x",
      mgcv = "x",
      mgcViz = "x",
      DHARMa = "x"
    )
  )
}

test_that("nbinom2 fingerprint is identical before and after a fit and summary() in one process", {
  skip_if_not_installed("glmmTMB")
  before <- do.call(get_analysis_fingerprint, fingerprint_spec())
  spec <- fingerprint_spec()
  m <- suppressWarnings(
    glmmTMB::glmmTMB(
      y ~ g + (1 | s),
      family = glmmTMB::nbinom2(),
      data = spec$data
    )
  )
  invisible(summary(m))
  after <- do.call(get_analysis_fingerprint, fingerprint_spec())
  expect_identical(after, before)
})

test_that("nbinom2 fingerprint is identical in a second clean process", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("callr")

  here_value <- do.call(
    get_analysis_fingerprint,
    fingerprint_spec()
  )

  candidate_source_root <- normalizePath(
    testthat::test_path(
      "..",
      ".."
    ),
    winslash = "/",
    mustWork = FALSE
  )

  use_source_tree <- file.exists(
    file.path(
      candidate_source_root,
      "DESCRIPTION"
    )
  )

  if (use_source_tree) {
    skip_if_not_installed("pkgload")
  }

  other <- callr::r(
    function(
      use_source_tree,
      source_root
    ) {
      if (isTRUE(use_source_tree)) {
        pkgload::load_all(
          source_root,
          quiet = TRUE
        )
      } else {
        loadNamespace(
          "EcoTools"
        )
      }

      fingerprint_function <- getFromNamespace(
        "get_analysis_fingerprint",
        "EcoTools"
      )

      set.seed(1)

      d <- data.frame(
        y = stats::rnbinom(
          120,
          mu = 4,
          size = 2
        ),
        g = factor(
          rep(
            1:3,
            40
          )
        ),
        s = factor(
          rep(
            1:12,
            each = 10
          )
        )
      )

      specification <- list(
        data = d,
        model_formula = "y ~ g + (1 | s)",
        zi_formula = ~0,
        disp_formula = ~1,
        family = glmmTMB::nbinom2(),
        weights = NULL,
        offset = NULL,
        package = "glmmTMB",
        dependency_versions = list(
          glmmTMB = "x",
          mgcv = "x",
          mgcViz = "x",
          DHARMa = "x"
        )
      )

      do.call(
        fingerprint_function,
        specification
      )
    },
    args = list(
      use_source_tree = use_source_tree,
      source_root = candidate_source_root
    ),
    libpath = .libPaths()
  )

  expect_identical(
    other,
    here_value
  )
})

test_that("a function-valued family element is not hashed", {
  fam <- glmmTMB::nbinom2()
  expect_null(
    fam[["a"]]
  )

  expect_true(
    "aic" %in% names(fam)
  )

  expect_true(
    is.function(
      fam[["aic"]]
    )
  )
  spec <- fingerprint_spec()
  fam2 <- fam
  fam2$a <- function(...) NULL
  spec2 <- spec
  spec2$family <- fam2
  expect_identical(
    do.call(get_analysis_fingerprint, spec2),
    do.call(get_analysis_fingerprint, spec)
  )
})
