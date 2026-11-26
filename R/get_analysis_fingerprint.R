# ==============================================================================
# Fingerprinting for GLMM/GAM cache validation
# ==============================================================================

#' Generate an Analysis Fingerprint
#'
#' @description
#' Creates a stable cache fingerprint from the model data, formulae, family,
#' fitting engine, additional arguments, sparse-data settings, fixed-effect
#' design-matrix settings, DHARMa simulation settings, residual-validation
#' settings, and the current EcoTools wrapper version.
#'
#' The fingerprint is used by [fit_and_simulate()] to determine whether a
#' cached GLMM or GAM was fitted and validated using the same analysis
#' specification as the current request.
#'
#' @param data
#' A data frame, data table, matrix, or compatible object containing the model
#' response and predictors.
#'
#' @param model_formula
#' The main model formula, or an engine-specific list of formulae.
#'
#' @param zi_formula
#' The zero-inflation formula used by a compatible model-fitting engine.
#'
#' @param disp_formula
#' The dispersion formula used by a compatible model-fitting engine.
#'
#' @param family
#' The model family object.
#'
#' @param weights
#' Optional model weights.
#'
#' @param offset
#' Optional model offset.
#'
#' @param package
#' Character string identifying the model-fitting engine.
#'
#' @param extra_args
#' A named list of additional arguments passed to the model-fitting engine.
#'
#' @param separation_threshold
#' Numeric threshold used to identify unusually large absolute parametric
#' coefficients.
#'
#' @param sparse_screen
#' Logical indicating whether sparse-data screening is enabled.
#'
#' @param sparse_min_events_per_parameter
#' Minimum number of informative events per estimated parameter used by the
#' sparse-data screen.
#'
#' @param sparse_min_events_per_level
#' Minimum number of informative events within each focal categorical level
#' used by the sparse-data screen.
#'
#' @param sparse_event_basis
#' Character string specifying whether the sparse-data screen uses nonzero
#' observations, total abundance, or both.
#'
#' @param sparse_count_like_only
#' Logical indicating whether sparse-data screening is restricted to
#' non-negative integer-like responses.
#'
#' @param design_check
#' Logical indicating whether the fixed-effect design matrix is validated
#' before model fitting.
#'
#' @param design_action
#' Character string specifying the action used by [check_model_matrix()] when a
#' design problem is detected. One of `"error"`, `"warn"`, or `"drop"`.
#'
#' @param design_tolerance
#' Non-negative numeric tolerance used to identify effectively zero columns and
#' determine the numerical rank of the fixed-effect design matrix.
#'
#' @param dharma_n_sim
#' Positive integer giving the number of simulated responses used to construct
#' the cached DHARMa residual object.
#'
#' @param dharma_seed
#' Integer random seed used for reproducible DHARMa residual simulation.
#'
#' @param tweedie_simulation_max_power
#' Finite numeric threshold strictly between 1 and 2 used to prevent
#' computationally unsafe Tweedie residual simulation near the Gamma boundary.
#' The threshold is included in the cache fingerprint.
#'
#' @param residual_alpha
#' Numeric significance threshold used by the formal residual-validation gate.
#'
#' @param dependency_versions
#' A named list containing the versions of the underlying modeling packages.
#'
#' @return
#' A character scalar containing the serialized digest of the analysis and
#' validation specification.
#'
#' @details
#' The function attempts to restrict the fingerprinted data to variables
#' referenced by the main, zero-inflation, and dispersion formulae. If those
#' variables cannot be identified reliably, the complete data object is used.
#'
#' Family metadata are converted to a stable list containing the family name,
#' link, class, and available family-specific parameters.
#'
#' DHARMa settings are included because changing the number of simulations,
#' random seed, or residual significance threshold may change the cached
#' residual diagnostics and final model-validity status.
#'
#' @seealso
#' [fit_and_simulate()],
#' [check_model_matrix()]
#'
#' @keywords internal
get_analysis_fingerprint <- function(
  data,
  model_formula,
  zi_formula,
  disp_formula,
  family,
  weights,
  offset,
  package,
  mgcv_fitter = "gam",
  extra_args = list(),
  separation_threshold = 5,
  sparse_screen = TRUE,
  sparse_min_events_per_parameter = 5,
  sparse_min_events_per_level = 3,
  sparse_event_basis = "nonzero",
  sparse_count_like_only = TRUE,
  design_check = TRUE,
  design_action = c("error", "warn", "drop"),
  design_tolerance = sqrt(.Machine$double.eps),
  dharma_n_sim = 1000L,
  dharma_seed = 123L,
  tweedie_simulation_max_power = 1.99,
  residual_alpha = 0.05,
  dependency_versions
) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop(
      "Package 'digest' is required for strict cache validation.",
      call. = FALSE
    )
  }

  design_action <- match.arg(design_action)

  # Validate the input of the mgcv_fitter argument:
  if (
    !is.character(mgcv_fitter) ||
      length(mgcv_fitter) != 1L ||
      is.na(mgcv_fitter) ||
      !mgcv_fitter %in%
        c(
          "gam",
          "bam"
        )
  ) {
    stop(
      "'mgcv_fitter' must be either 'gam' or 'bam'.",
      call. = FALSE
    )
  }

  data_df <- canonicalize_tabular_for_hash(data)

  formula_vars <- tryCatch(
    {
      if (is.list(model_formula)) {
        unique(
          unlist(
            lapply(
              model_formula,
              function(f) {
                all.vars(
                  stats::as.formula(f)
                )
              }
            )
          )
        )
      } else {
        all.vars(
          stats::as.formula(model_formula)
        )
      }
    },
    error = function(e) {
      character(0)
    }
  )

  all_formula_vars <- tryCatch(
    {
      unique(
        c(
          formula_vars,
          all.vars(
            stats::as.formula(zi_formula)
          ),
          all.vars(
            stats::as.formula(disp_formula)
          )
        )
      )
    },
    error = function(e) {
      character(0)
    }
  )

  data_subset <- if (
    length(all_formula_vars) > 0L &&
      all(all_formula_vars %in% names(data_df))
  ) {
    data_df[,
      all_formula_vars,
      drop = FALSE
    ]
  } else {
    data_df
  }

  data_subset <- canonicalize_tabular_for_hash(data_subset)

  exact_family_element <- function(
    family_object,
    element_name,
    default = NULL
  ) {
    if (
      is.null(family_object) ||
        !is.list(family_object) ||
        !element_name %in% names(family_object)
    ) {
      return(default)
    }

    family_object[[element_name]]
  }

  stable_non_function <- function(
    value
  ) {
    if (
      is.null(value) ||
        is.function(value) ||
        is.environment(value)
    ) {
      return(NULL)
    }

    value
  }

  family_name <- tryCatch(
    {
      value <- exact_family_element(
        family,
        "family",
        default = NA_character_
      )

      if (length(value) == 0L) {
        NA_character_
      } else {
        as.character(value[[1L]])
      }
    },
    error = function(e) {
      NA_character_
    }
  )

  family_link <- tryCatch(
    {
      value <- exact_family_element(
        family,
        "link",
        default = NA_character_
      )

      if (length(value) == 0L) {
        NA_character_
      } else {
        as.character(value[[1L]])
      }
    },
    error = function(e) {
      NA_character_
    }
  )

  family_class <- tryCatch(
    class(family),
    error = function(e) NA_character_
  )

  # Use exact family-element matching and exclude functions and environments
  # from the serialized fingerprint. This prevents partial matches such as
  # family$a matching glmmTMB's family$aic closure.

  family_stable <- list(
    name = family_name,
    link = family_link,
    class = family_class,
    p = tryCatch(
      stable_non_function(
        exact_family_element(
          family,
          "p"
        )
      ),
      error = function(e) {
        NULL
      }
    ),
    theta = tryCatch(
      {
        get_theta <- exact_family_element(
          family,
          "getTheta"
        )

        if (is.function(get_theta)) {
          stable_non_function(
            get_theta()
          )
        } else {
          stable_non_function(
            exact_family_element(
              family,
              "theta"
            )
          )
        }
      },
      error = function(e) {
        NULL
      }
    ),
    a = tryCatch(
      stable_non_function(
        exact_family_element(
          family,
          "a"
        )
      ),
      error = function(e) {
        NULL
      }
    ),
    b = tryCatch(
      stable_non_function(
        exact_family_element(
          family,
          "b"
        )
      ),
      error = function(e) {
        NULL
      }
    )
  )

  fingerprint_list <- list(
    wrapper_version = get_wrapper_version(),
    dependency_versions = dependency_versions,
    data_subset = data_subset,
    data_dimensions = dim(data_df),
    model_formula = canonicalize_language_for_hash(
      model_formula
    ),
    zero_inflation_formula = canonicalize_language_for_hash(
      zi_formula
    ),
    dispersion_formula = canonicalize_language_for_hash(
      disp_formula
    ),
    family = family_stable,
    weights = canonicalize_extra_args_for_hash(
      weights
    ),
    offset = canonicalize_extra_args_for_hash(
      offset
    ),
    package = as.character(package),
    mgcv_fitter = if (identical(package, "mgcv")) {
      as.character(
        mgcv_fitter
      )
    } else {
      NA_character_
    },
    extra_args = canonicalize_extra_args_for_hash(
      extra_args
    ),
    separation_threshold = as.numeric(
      separation_threshold
    ),
    sparse_screen = isTRUE(sparse_screen),
    sparse_min_events_per_parameter = as.numeric(
      sparse_min_events_per_parameter
    ),
    sparse_min_events_per_level = as.numeric(
      sparse_min_events_per_level
    ),
    sparse_event_basis = as.character(
      sparse_event_basis
    ),
    sparse_count_like_only = isTRUE(
      sparse_count_like_only
    ),
    design_check = isTRUE(design_check),
    design_action = as.character(design_action),
    design_tolerance = as.numeric(
      design_tolerance
    ),
    dharma_n_sim = as.integer(
      dharma_n_sim
    ),
    dharma_seed = as.integer(
      dharma_seed
    ),
    tweedie_simulation_max_power = tweedie_simulation_max_power,
    residual_alpha = as.numeric(
      residual_alpha
    )
  )

  digest::digest(
    fingerprint_list,
    algo = "xxhash64",
    serialize = TRUE
  )
}
