test_that("glmmTMB character family names are normalized", {
  normalized_family <- normalize_model_family(
    family = "nbinom2",
    package = "glmmTMB"
  )

  expect_true(
    is.list(normalized_family)
  )

  expect_identical(
    normalized_family$family,
    "nbinom2"
  )

  expect_identical(
    normalized_family$link,
    "log"
  )
})


test_that("base family names are normalized for glmmTMB", {
  normalized_family <- normalize_model_family(
    family = "poisson",
    package = "glmmTMB"
  )

  expect_true(
    is.list(normalized_family)
  )

  expect_identical(
    normalized_family$family,
    "poisson"
  )

  expect_identical(
    normalized_family$link,
    "log"
  )
})


test_that("base family names are normalized for mgcv", {
  normalized_family <- normalize_model_family(
    family = "gaussian",
    package = "mgcv"
  )

  expect_true(
    is.list(normalized_family)
  )

  expect_identical(
    normalized_family$family,
    "gaussian"
  )

  expect_identical(
    normalized_family$link,
    "identity"
  )
})


test_that("existing family objects are returned unchanged", {
  original_family <- stats::poisson()

  normalized_family <- normalize_model_family(
    family = original_family,
    package = "glmmTMB"
  )

  expect_identical(
    normalized_family,
    original_family
  )
})


test_that("unsupported family names fail informatively", {
  expect_error(
    normalize_model_family(
      family = "not_a_real_family",
      package = "glmmTMB"
    ),
    regexp = "Unsupported glmmTMB family name"
  )

  expect_error(
    normalize_model_family(
      family = "not_a_real_family",
      package = "mgcv"
    ),
    regexp = "Unsupported mgcv family name"
  )
})


test_that("invalid character family specifications are rejected", {
  expect_error(
    normalize_model_family(
      family = character(0),
      package = "glmmTMB"
    ),
    regexp = "one non-empty character value"
  )

  expect_error(
    normalize_model_family(
      family = NA_character_,
      package = "glmmTMB"
    ),
    regexp = "one non-empty character value"
  )

  expect_error(
    normalize_model_family(
      family = "",
      package = "glmmTMB"
    ),
    regexp = "one non-empty character value"
  )
})


test_that("unknown engines reject character family specifications", {
  expect_error(
    normalize_model_family(
      family = "poisson",
      package = "unknown_engine"
    ),
    regexp = "not implemented for engine"
  )
})
