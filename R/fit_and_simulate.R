# ==============================================================================
# GLMM/GAM fitting wrapper with DHARMa simulation
# ==============================================================================
#' Fit a GLMM or GAM and Generate DHARMa Residuals
#'
#' @description
#' Fits a generalized linear mixed model using [glmmTMB::glmmTMB()] or a
#' generalized additive model using [mgcv::gam()], performs convergence and
#' coefficient checks, optionally screens count-like responses for sparse-data
#' risk, generates DHARMa simulation residuals when appropriate, and caches the
#' resulting model bundle to disk.
#'
#' Cached models are reused only when their stored analysis fingerprint matches
#' the current data, formula, family, model arguments, screening thresholds, and
#' wrapper version. Otherwise, the cache is invalidated and the model is refit.
#'
#' @param model_name
#' A character scalar used to identify the model and construct the cached RDS
#' file name.
#'
#' @param package
#' Character string specifying the model-fitting engine. Supported values are
#' `"glmmTMB"` and `"mgcv"`.
#'
#' @param mgcv_fitter
#' Character string selecting the mgcv fitting function when
#' `package = "mgcv"`. One of:
#'
#' \itemize{
#'   \item `"gam"` uses [mgcv::gam()] and is the default;
#'   \item `"bam"` uses [mgcv::bam()] for large-data GAM fitting.
#' }
#'
#' The argument is not applicable when `package = "glmmTMB"`.
#' EcoTools does not select `bam()` automatically because `gam()` and `bam()`
#' can use different numerical algorithms and smoothing-parameter defaults.
#'
#' For `mgcv_fitter = "bam"`, the wrapper passes `num_cores` through the
#' direct `nthreads` argument of [mgcv::bam()]. Additional options such as
#' `discrete`, `method`, `chunk.size`, `rho`, and `AR.start` may be supplied
#' through `...`.
#'
#' @param model_formula
#' A model formula, character representation of a formula, or an engine-specific
#' list of formulas accepted by the selected model-fitting package.
#'
#' @param zi_formula
#' A formula specifying the zero-inflation component for a
#' [glmmTMB::glmmTMB()] model. The default `~0` disables zero inflation.
#'
#' @param disp_formula
#' A formula specifying the dispersion component for a
#' [glmmTMB::glmmTMB()] model. The default is `~1`.
#'
#' @param family
#' A family object accepted by the selected fitting engine, or a supported
#' character family name. Examples include `stats::poisson()`,
#' `glmmTMB::nbinom2()`, `"poisson"`, and `"nbinom2"`. Character family names
#' are resolved inside the appropriate package namespace and do not require the
#' modeling package to be attached to the search path.
#'
#' @param data
#' A data frame, data table, or compatible object containing the response and
#' predictors used in the model.
#'
#' @param weights
#' Optional model weights. When `NULL`, no weights argument is included in the
#' underlying model call.
#'
#' @param offset
#' Optional model offset. When `NULL`, no offset argument is included in the
#' underlying model call.
#'
#' @param output_dir
#' Character string specifying the cache directory. Relative paths are resolved
#' using [here::here()]. The directory is created recursively when necessary.
#'
#' @param force_refit
#' Logical. If `TRUE`, ignore any existing cached model and fit the requested
#' model again. If `FALSE`, reuse a cache only when its fingerprint matches the
#' current specification.
#'
#' @param num_cores
#' Optional positive integer giving the number of processor cores or threads
#' supplied to the fitting engine. When `NULL`, the number is selected by
#' [get_auto_cores()].
#'
#' @param separation_threshold
#' Positive numeric threshold used to flag unusually large absolute parametric
#' coefficients. Exceeding this threshold produces a suspicious model status
#' rather than automatically treating the optimizer as failed.
#'
#' @param sparse_screen
#' Logical. If `TRUE`, run [screen_sparse_count_model()] for count-like or binary
#' responses.
#'
#' @param sparse_min_events_per_parameter
#' Positive numeric threshold used by the sparse-data screen for the minimum
#' number of informative observations per approximate estimated model parameter.
#'
#' @param sparse_min_events_per_level
#' Positive numeric threshold used by the sparse-data screen for the minimum
#' number of informative events in each focal categorical-predictor level.
#'
#' @param sparse_event_basis
#' Character string controlling how the sparse-data screen defines informative
#' events. One of `"nonzero"`, `"total"`, or `"both"`. The default selection is
#' `"nonzero"`.
#'
#' @param sparse_count_like_only
#' Logical. If `TRUE`, sparse-data screening is limited to non-negative,
#' integer-like responses.
#'
#' @param design_check
#' Logical. If `TRUE`, validate the conditional fixed-effect design matrix before
#' cache evaluation and model fitting. The current pre-fit implementation is
#' applied automatically to `glmmTMB` models.
#'
#' @param design_action
#' Action taken when unused factor levels, all-zero design columns,
#' non-finite values, or rank deficiency are detected. One of:
#'
#' \itemize{
#'   \item `"error"` to stop before fitting;
#'   \item `"warn"` to warn and continue with the unchanged data;
#'   \item `"drop"` to remove unused factor levels from a local copy of `data`,
#'   rebuild the design matrix, and stop if structural problems remain.
#' }
#'
#' The default resolves to `"error"`.
#'
#' @param design_tolerance
#' Non-negative numeric tolerance used by [check_model_matrix()] to identify
#' effectively all-zero columns and determine numerical matrix rank.
#'
#' @param dharma_n_sim
#' Positive integer giving the number of simulated responses generated for the
#' cached DHARMa residual object. Larger values reduce Monte Carlo variability
#' but increase computation time. The default is `1000`.
#'
#' @param dharma_seed
#' Integer random seed passed to [DHARMa::simulateResiduals()] and used by the
#' custom `twlss` simulation branch. Including the seed in the cache fingerprint
#' makes the stored residual-validation status reproducible.
#'
#' @param residual_alpha
#' Numeric significance threshold used by the formal DHARMa residual-validation
#' gate. The value must be strictly between zero and one. The default is `0.05`.
#'
#' @param tweedie_simulation_max_power
#' Finite numeric threshold strictly between 1 and 2. For fitted Tweedie
#' models, DHARMa simulation and repeated diagnostic stability checks are
#' skipped when the estimated power is at or above this threshold. The default is
#' `1.99`. Models blocked by this guard are marked invalid for inference
#' because simulation-based validation was not completed.
#'
#' @param ...
#' Additional named arguments passed to the selected model-fitting engine.
#' Core arguments managed directly by this wrapper, such as `formula`, `data`,
#' `family`, `weights`, `offset`, and `control`, must not be supplied through
#' `...`.
#'
#' @details
#' For `glmmTMB` models, fixed-effect design validation should be performed on
#' the conditional formula after random-effect terms have been removed.
#' Zero-inflation and dispersion formulae are fingerprinted but are not currently
#' subjected to the same pre-fit design-matrix check.
#'
#' Experimental glmmTMB models containing [mgcv::s()] terms receive special
#' call-storage handling. Their evaluated analysis data are embedded in the
#' stored glmmTMB call because population-level prediction and simulation may
#' reconstruct the smooth-generated fixed and random-effect structures after
#' the caller's original data object no longer exists. This can increase the
#' serialized model size. Smooth support in glmmTMB remains experimental, and
#' EcoTools currently rejects `s(..., by = ...)` terms because the underlying
#' glmmTMB implementation does not support by-variable smooths.
#'
#' For `mgcv` models, the wrapper checks available optimizer-convergence
#' information and applies the coefficient threshold only to parametric
#' coefficients. Smooth-basis coefficients are excluded from the
#' extreme-coefficient screen.
#'
#' Sparse-data warnings are cached but do not automatically invalidate a model.
#' Use [assert_valid_model()] to enforce the validity policy required by a
#' downstream analysis.
#'
#' When `design_action = "drop"`, the wrapper fits the model using a cleaned
#' local copy of `data` and stores that evaluated data frame in
#' `analysis_data`. The fitted model's engine-generated `call$data` expression
#' is preserved rather than replaced with the caller's original data symbol,
#' because caller-local symbols may no longer exist after model fitting or
#' serialization.
#'
#' For ordinary `mgcv` GAM objects, DHARMa response simulation uses the
#' `simulate.gam()` support provided by the mgcViz package. EcoTools loads the
#' mgcViz namespace when GAM residual simulation is requested.
#'
#' @return
#' A list containing model results and cached validation metadata. Depending on
#' fitting success and model type, the returned object may include:
#'
#' \describe{
#'   \item{model}{
#'     The fitted `glmmTMB` or `mgcv` model, or `NULL` when fitting failed.
#'   }
#'   \item{residuals}{
#'     A DHARMa residual object, or `NULL` when simulation was unavailable,
#'     inappropriate, or skipped because the model was invalid.
#'   }
#'   \item{converged}{
#'     Logical indicating whether the numerical optimizer, likelihood, and
#'     Hessian checks passed. This field does not imply that residual
#'     diagnostics passed.
#'   }
#'   \item{model_status}{
#'     Standardized model-validity status. Possible values include:
#'
#'     \itemize{
#'       \item `"valid"` when numerical and residual checks pass;
#'       \item `"valid_fit_no_residual_diagnostics"` when numerical checks pass
#'       but residual validity could not be assessed;
#'       \item `"suspicious_residual_diagnostics"` when numerical checks pass
#'       but one or more primary DHARMa residual tests fail;
#'       \item `"suspicious_extreme_coefficients"` when extreme coefficients
#'       indicate possible separation or weak identification;
#'       \item `"invalid"` when optimizer, likelihood, or Hessian validation
#'       fails.
#'     }
#'   }
#'   \item{valid_for_inference}{
#'     Logical indicating whether the wrapper considers the model eligible for
#'     downstream formal inference.
#'   }
#'   \item{diagnose}{
#'     Optional detailed diagnostic output produced for a problematic
#'     `glmmTMB` model.
#'   }
#'   \item{warnings}{
#'     Character vector of fitting, convergence, coefficient, sparse-data, or
#'     simulation warnings captured by the wrapper.
#'   }
#'   \item{sparse_screen}{
#'     Detailed output returned by [screen_sparse_count_model()], or `NULL` when
#'     sparse screening was disabled or unavailable.
#'   }
#'   \item{data_dim}{
#'     Dimensions of the supplied model data.
#'   }
#'   \item{simulation_error}{
#'     Detailed metadata for an error encountered while generating DHARMa
#'     residual simulations, including the condition message, class, and
#'     failing call. This element is `NULL` when simulation succeeds or is not
#'     attempted.
#'   }
#'   \item{residual_diagnostics}{
#'     A list containing the primary DHARMa test results, extracted p-values,
#'     failed-test names, and the observed-to-simulated dispersion ratio.
#'   }
#'   \item{residual_diagnostics_valid}{
#'     Logical indicating whether all required residual tests were available and
#'     passed the configured significance threshold.
#'   }
#'   \item{design_check}{
#'     Detailed result returned by [check_model_matrix()], or `NULL` when fixed-
#'     effect design validation was disabled or not applicable.
#'   }
#'  \item{valid_for_model_selection}{
#'     Logical indicating whether the wrapper permits the model to remain under
#'     consideration for exploratory comparison or model selection. This field
#'     is separate from `valid_for_inference`.
#'   }
#'  \item{fit_error}{
#'     Detailed metadata for a fitting error, including the condition message,
#'     condition classes, failing call, and captured call stack. This element is
#'     `NULL` for successful fits.
#'   }
#'  \item{analysis_data}{
#'     The data frame used by the fitted model after any local preprocessing or
#'     unused-level removal performed by the wrapper.
#'   }
#'  \item{data_expression}{
#'     Character representation of the original expression supplied through
#'     the `data` argument. This is retained for provenance and is not injected
#'     into the fitted model call.
#'   }
#'   \item{analysis_fingerprint}{
#'     Character hash identifying the model specification and relevant inputs.
#'   }
#'   \item{session_info}{
#'     Session information recorded when the model bundle was fitted.
#'   }
#' }
#'
#' @seealso
#' [DHARMa_diag()],
#' [DHARMa_stability_check()],
#' [screen_sparse_count_model()],
#' [assert_valid_model()],
#' [glmmTMB::glmmTMB()],
#' [mgcv::gam()],
#' [DHARMa::simulateResiduals()]
#'
#' @examples
#' \dontrun{
#' model_bundle <- fit_and_simulate(
#'   model_name = "example_poisson_glmm",
#'   package = "glmmTMB",
#'   model_formula = abundance ~ year + (1 | plot_id),
#'   family = poisson(),
#'   data = analysis_data,
#'   output_dir = "Analysis/Models",
#'   force_refit = FALSE
#' )
#'
#' assert_valid_model(
#'   model_bundle,
#'   purpose = "inference"
#' )
#' }
#'
#' @export
fit_and_simulate <- function(
  model_name,
  package = c("glmmTMB", "mgcv"),
  mgcv_fitter = c("gam", "bam"),
  model_formula,
  zi_formula = ~0,
  disp_formula = ~1,
  family,
  data,
  weights = NULL,
  offset = NULL,
  output_dir = "Models",
  force_refit = FALSE,
  num_cores = NULL,
  separation_threshold = 5,
  sparse_screen = TRUE,
  sparse_min_events_per_parameter = 5,
  sparse_min_events_per_level = 3,
  sparse_event_basis = c("nonzero", "total", "both"),
  sparse_count_like_only = TRUE,
  design_check = TRUE,
  design_action = c("error", "warn", "drop"),
  design_tolerance = sqrt(.Machine$double.eps),
  dharma_n_sim = 1000L,
  dharma_seed = 123L,
  tweedie_simulation_max_power = 1.99,
  residual_alpha = 0.05,
  ...
) {
  # CAPTURE THE TRUE DATASET NAME
  # This grabs the literal text the user typed into the 'data' argument
  # e.g., it captures the symbol `biodiv_coast` instead of the evaluated data frame.
  true_data_name <- substitute(data)

  # Capture the family expression as a text
  family_expression <- substitute(family)

  # Match arguments:
  sparse_event_basis <- match.arg(sparse_event_basis)
  design_action <- match.arg(design_action)
  package <- match.arg(package)
  mgcv_fitter <- match.arg(mgcv_fitter)
  # Error message when a user attempts to change mgcv_fitter argument when the choesn package is glmmTMB:
  if (
    identical(package, "glmmTMB") &&
      !identical(mgcv_fitter, "gam")
  ) {
    stop(
      "'mgcv_fitter' is applicable only when package = 'mgcv'.",
      call. = FALSE
    )
  }

  bare_family_name <- if (
    is.call(family_expression) &&
      length(family_expression) >= 1L &&
      is.symbol(family_expression[[1L]])
  ) {
    as.character(family_expression[[1L]])
  } else {
    NULL
  }

  known_engine_family_names <- if (identical(package, "glmmTMB")) {
    c(
      "nbinom1",
      "nbinom2",
      "nbinom12",
      "compois",
      "genpois",
      "truncated_poisson",
      "truncated_nbinom1",
      "truncated_nbinom2",
      "truncated_compois",
      "truncated_genpois",
      "beta_family",
      "betabinomial",
      "tweedie",
      "lognormal",
      "ziGamma",
      "t_family",
      "ordbeta"
    )
  } else {
    c(
      "betar",
      "cnorm",
      "cox.ph",
      "gevlss",
      "gumbls",
      "multinom",
      "nb",
      "ocat",
      "scat",
      "shash",
      "tw",
      "twlss",
      "ziP"
    )
  }

  engine_namespace <- if (identical(package, "glmmTMB")) {
    "glmmTMB"
  } else {
    "mgcv"
  }

  if (
    !is.null(bare_family_name) &&
      bare_family_name %in% known_engine_family_names &&
      !exists(
        bare_family_name,
        mode = "function",
        envir = parent.frame(),
        inherits = TRUE
      )
  ) {
    stop(
      paste0(
        "Family function '",
        bare_family_name,
        "()' is not available on the current search path. ",
        "Use family = \"",
        bare_family_name,
        "\" or family = ",
        engine_namespace,
        "::",
        bare_family_name,
        "()."
      ),
      call. = FALSE
    )
  }

  # Force the evaluation of the arguments early to detect any problems
  force(model_name)
  force(model_formula)
  force(output_dir)

  glmmTMB_has_smooth <- FALSE

  if (identical(package, "glmmTMB")) {
    if (is.list(model_formula)) {
      stop(
        "glmmTMB smooth support currently requires one model formula, ",
        "not a list of formulas.",
        call. = FALSE
      )
    }

    normalized_glmmTMB_formula <-
      normalize_glmmtmb_smooth_formula(
        model_formula
      )

    glmmTMB_has_smooth <-
      formula_has_mgcv_smooth(
        normalized_glmmTMB_formula
      )

    if (isTRUE(glmmTMB_has_smooth)) {
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop(
          "Package 'mgcv' is required when a glmmTMB formula contains s().",
          call. = FALSE
        )
      }

      if (
        formula_has_glmmtmb_by_smooth(
          normalized_glmmTMB_formula
        )
      ) {
        stop(
          paste0(
            "UNSUPPORTED GLMMTMB SMOOTH: the formula contains an s() term ",
            "with a 'by' argument. The current experimental glmmTMB smooth ",
            "implementation cannot handle by-variable smooths. EcoTools can ",
            "currently pass ordinary s(x, ...) terms to glmmTMB, but cannot fit ",
            "separate smooth trajectories using s(x, by = factor, ...). ",
            "Consider a common smooth plus the factor main effect, a parametric ",
            "time-by-factor interaction, or an explicitly constructed spline ",
            "basis interaction."
          ),
          call. = FALSE
        )
      }

      model_formula <- normalized_glmmTMB_formula
    }
  }

  if (
    identical(package, "glmmTMB") &&
      (formula_has_mgcv_smooth(zi_formula) ||
        formula_has_mgcv_smooth(disp_formula))
  ) {
    stop(
      paste0(
        "EcoTools currently supports glmmTMB s() terms only in the ",
        "conditional model formula. Smooth terms in zi_formula or ",
        "disp_formula are not currently supported."
      ),
      call. = FALSE
    )
  }

  family <- normalize_model_family(
    family = family,
    package = package
  )

  # DEFENSIVELY COERCE TO BASE R DATA.FRAME
  # This permanently protects glmmTMB, mgcv, and all downstream screeners
  # from data.table and tibble dimension-dropping bugs.
  data <- as.data.frame(data)

  if (
    length(dharma_n_sim) != 1L ||
      is.na(dharma_n_sim) ||
      !is.numeric(dharma_n_sim) ||
      !is.finite(dharma_n_sim) ||
      dharma_n_sim < 1L ||
      dharma_n_sim != as.integer(dharma_n_sim)
  ) {
    stop(
      "'dharma_n_sim' must be one positive integer.",
      call. = FALSE
    )
  }

  dharma_n_sim <- as.integer(dharma_n_sim)

  if (
    length(dharma_seed) != 1L ||
      is.na(dharma_seed) ||
      !is.numeric(dharma_seed) ||
      !is.finite(dharma_seed)
  ) {
    stop(
      "'dharma_seed' must be one finite numeric value.",
      call. = FALSE
    )
  }

  dharma_seed <- as.integer(dharma_seed)

  if (
    !is.numeric(tweedie_simulation_max_power) ||
      length(tweedie_simulation_max_power) != 1L ||
      is.na(tweedie_simulation_max_power) ||
      !is.finite(tweedie_simulation_max_power) ||
      tweedie_simulation_max_power <= 1 ||
      tweedie_simulation_max_power >= 2
  ) {
    stop(
      paste0(
        "'tweedie_simulation_max_power' must be one finite numeric ",
        "value strictly between 1 and 2."
      ),
      call. = FALSE
    )
  }

  if (
    length(residual_alpha) != 1L ||
      is.na(residual_alpha) ||
      !is.numeric(residual_alpha) ||
      !is.finite(residual_alpha) ||
      residual_alpha <= 0 ||
      residual_alpha >= 1
  ) {
    stop(
      "'residual_alpha' must be one finite numeric value strictly ",
      "between 0 and 1.",
      call. = FALSE
    )
  }

  design_check_result <- NULL

  if (
    isTRUE(design_check) &&
      identical(package, "glmmTMB")
  ) {
    design_formula <- if (is.list(model_formula)) {
      model_formula[[1L]]
    } else {
      model_formula
    }

    design_formula <- stats::as.formula(
      design_formula,
      env = parent.frame()
    )

    fixed_design_formula <- reformulas::nobars(
      design_formula
    )

    design_formula_used <- if (isTRUE(glmmTMB_has_smooth)) {
      strip_mgcv_smooth_terms(
        fixed_design_formula
      )
    } else {
      fixed_design_formula
    }

    design_check_result <- check_model_matrix(
      formula = design_formula_used,
      data = data,
      action = design_action,
      tolerance = design_tolerance
    )

    if (isTRUE(glmmTMB_has_smooth)) {
      smooth_term_labels <- attr(
        stats::terms(
          fixed_design_formula
        ),
        "term.labels"
      )

      smooth_term_labels <- smooth_term_labels[
        grepl(
          "(mgcv::)?s\\s*\\(",
          smooth_term_labels,
          perl = TRUE
        )
      ]

      design_check_result$smooth_terms_detected <-
        smooth_term_labels

      design_check_result$smooth_terms_deferred_to_glmmTMB <-
        TRUE

      design_check_result$parametric_formula_checked <-
        design_formula_used

      design_check_result$smooth_design_note <- paste0(
        "Ordinary parametric fixed effects were checked before fitting. ",
        "The smooth basis and smoothing penalty were not passed to ",
        "stats::model.matrix(); they were deferred to glmmTMB's experimental ",
        "mgcv smooth-term preprocessing."
      )
    }

    data <- design_check_result$data
  }

  wrapper_call <- match.call()
  data_was_inline_expression <- !is.null(wrapper_call$data) &&
    !is.symbol(wrapper_call$data)
  extra_args <- list(...)

  formula_contains_offset <- tryCatch(
    {
      formulas_to_check <- if (is.list(model_formula)) {
        model_formula
      } else {
        list(
          model_formula
        )
      }

      any(
        vapply(
          formulas_to_check,
          function(formula_object) {
            formula_text <- paste(
              deparse(
                formula_object,
                width.cutoff = 500L
              ),
              collapse = ""
            )

            grepl(
              "(^|[^[:alnum:]_.])offset\\s*\\(",
              formula_text,
              perl = TRUE
            )
          },
          logical(1)
        )
      )
    },
    error = function(e) {
      FALSE
    }
  )

  if (
    !is.null(offset) &&
      isTRUE(formula_contains_offset)
  ) {
    stop(
      paste0(
        "Specify the model offset either inside 'model_formula' using ",
        "offset(...), or through the separate 'offset' argument, but not ",
        "both. Providing both would apply the offset twice."
      ),
      call. = FALSE
    )
  }

  if (
    identical(package, "glmmTMB") &&
      isTRUE(glmmTMB_has_smooth)
  ) {
    if (!"REML" %in% names(extra_args)) {
      extra_args$REML <- TRUE

      message(
        paste0(
          "glmmTMB smooth term detected. EcoTools is using REML = TRUE, ",
          "as recommended by the experimental glmmTMB smooth-term interface."
        )
      )
    } else if (!isTRUE(extra_args$REML)) {
      warning(
        paste0(
          "A glmmTMB smooth term was detected with REML = FALSE. ",
          "glmmTMB currently recommends REML = TRUE for models containing ",
          "mgcv smooths. Proceeding with the explicitly requested setting."
        ),
        call. = FALSE
      )
    }
  }

  forbidden_extra_args <- c(
    "formula",
    "data",
    "family",
    "weights",
    "offset",
    "control",
    "ziformula",
    "dispformula",
    "nthreads",
    "cluster"
  )

  bad_extra_args <- intersect(names(extra_args), forbidden_extra_args)

  if (length(bad_extra_args) > 0) {
    stop(
      "Do not pass the following arguments through '...': ",
      paste(bad_extra_args, collapse = ", "),
      ". These are managed directly by fit_and_simulate().",
      call. = FALSE
    )
  }

  auto_cores <- get_auto_cores(num_cores)

  output_dir_full <- here::here(output_dir)
  file_path <- file.path(output_dir_full, paste0(model_name, ".rds"))

  if (!dir.exists(output_dir_full)) {
    dir.create(output_dir_full, recursive = TRUE)
  }

  # Record package versions:
  dependency_versions <- list(
    glmmTMB = tryCatch(
      as.character(utils::packageVersion("glmmTMB")),
      error = function(e) NA_character_
    ),
    mgcv = tryCatch(
      as.character(utils::packageVersion("mgcv")),
      error = function(e) NA_character_
    ),
    mgcViz = tryCatch(
      as.character(
        utils::packageVersion("mgcViz")
      ),
      error = function(e) NA_character_
    ),
    DHARMa = tryCatch(
      as.character(utils::packageVersion("DHARMa")),
      error = function(e) NA_character_
    )
  )

  current_fingerprint <- get_analysis_fingerprint(
    data = data,
    model_formula = model_formula,
    zi_formula = zi_formula,
    disp_formula = disp_formula,
    family = family,
    weights = weights,
    offset = offset,
    package = package,
    mgcv_fitter = mgcv_fitter,
    extra_args = extra_args,
    separation_threshold = separation_threshold,
    sparse_screen = sparse_screen,
    sparse_min_events_per_parameter = sparse_min_events_per_parameter,
    sparse_min_events_per_level = sparse_min_events_per_level,
    sparse_event_basis = sparse_event_basis,
    sparse_count_like_only = sparse_count_like_only,
    design_check = design_check,
    design_action = design_action,
    design_tolerance = design_tolerance,
    dharma_n_sim = dharma_n_sim,
    dharma_seed = dharma_seed,
    tweedie_simulation_max_power = tweedie_simulation_max_power,
    residual_alpha = residual_alpha,
    dependency_versions = dependency_versions
  )

  # ---------------------------------------------------------------------------
  # 1. Cache inspection and validation
  # ---------------------------------------------------------------------------
  if (!force_refit && file.exists(file_path)) {
    output_bundle <- readRDS(file_path)
    if (
      is.list(output_bundle) &&
        "model" %in% names(output_bundle) &&
        !inherits(output_bundle, "ecotools_model_bundle")
    ) {
      class(output_bundle) <- c(
        "ecotools_model_bundle",
        "list"
      )
    }
    cached_fingerprint <- if (!is.null(output_bundle$analysis_fingerprint)) {
      output_bundle$analysis_fingerprint
    } else {
      NULL
    }

    if (is.null(cached_fingerprint)) {
      warning(
        "\nCached model lacks fingerprint metadata. Re-fitting to ensure reproducibility.\n\n",
        call. = FALSE
      )
    } else if (!identical(cached_fingerprint, current_fingerprint)) {
      warning(stringr::str_glue(
        "\nCACHE INVALIDATED: Data, formula, family, wrapper version, thresholds, or arguments changed for '{model_name}'.\n",
        "  -> Re-fitting model...\n\n"
      ))
    } else {
      cat(stringr::str_glue(
        "\nLoading existing model '{model_name}' ({package}) from disk [Fingerprint Verified]. Skipping fit.\n"
      ))

      status_value <- if (
        !is.null(output_bundle$model_status) &&
          length(output_bundle$model_status) == 1L &&
          !is.na(output_bundle$model_status)
      ) {
        as.character(output_bundle$model_status)
      } else {
        "unknown"
      }

      status_message <- switch(
        status_value,

        valid = paste0(
          "SUCCESS: Model status is VALID FOR INFERENCE, subject to the ",
          "stored residual, sparsity, and model-assumption diagnostics."
        ),

        valid_fit_no_standard_errors = paste0(
          "CAUTION: The model fit is available, but standard errors were ",
          "not assessed."
        ),

        valid_fit_no_residual_diagnostics = paste0(
          "CAUTION: The model passed numerical checks, but residual ",
          "diagnostics were unavailable."
        ),

        suspicious_residual_diagnostics = paste0(
          "WARNING: The model passed optimizer and Hessian checks, but ",
          "DHARMa residual validation detected lack of fit."
        ),

        suspicious_extreme_coefficients = paste0(
          "WARNING: The model has extreme coefficients or possible separation."
        ),

        suspicious_optimization_incomplete = paste0(
          "WARNING: Optimization did not satisfy the configured gradient criterion."
        ),

        suspicious_optimization_unassessed = paste0(
          "WARNING: Optimization could not be fully assessed."
        ),

        suspicious_inference_unreliable = paste0(
          "WARNING: The fitted model has unreliable inferential uncertainty."
        ),

        invalid = if (
          !is.null(output_bundle$tweedie_simulation_guard) &&
            isTRUE(output_bundle$tweedie_simulation_guard$blocked)
        ) {
          paste0(
            "INVALID FOR INFERENCE: Simulation-based validation was skipped ",
            "because the fitted Tweedie power was at or above the configured ",
            "simulation-safety threshold."
          )
        } else {
          paste0(
            "INVALID: Optimizer, likelihood, Hessian, or another required ",
            "validation check failed."
          )
        },

        aborted = "ABORTED: Model fitting was not completed.",

        paste0(
          "UNKNOWN MODEL STATUS: ",
          status_value
        )
      )

      if (identical(status_value, "valid")) {
        message(
          paste0(
            "  -> ",
            status_message,
            "\n"
          )
        )
      } else {
        warning(
          paste0(
            "  -> ",
            status_message
          ),
          call. = FALSE,
          immediate. = TRUE
        )
      }

      if (!is.null(output_bundle$warnings)) {
        cat("Associated Warnings:\n")
        lapply(output_bundle$warnings, function(w) cat(paste0("  - ", w, "\n")))
        cat("\n")
      }

      if (
        identical(status_value, "invalid") &&
          !is.null(output_bundle$fit_error) &&
          !is.null(output_bundle$fit_error$message)
      ) {
        warning(
          paste0(
            "  -> Stored fitting error: ",
            output_bundle$fit_error$message
          ),
          call. = FALSE,
          immediate. = TRUE
        )
      }

      return(output_bundle)
    }
  }

  # ---------------------------------------------------------------------------
  # 2. Pre-fit checks
  # ---------------------------------------------------------------------------
  captured_warnings <- character()

  if (data_was_inline_expression) {
    inline_msg <- paste0(
      "DATA PROVENANCE NOTICE: data was supplied as an inline expression. ",
      "The evaluated analysis data are stored in the model bundle, but a named ",
      "data object generally produces a more readable model call and audit trail."
    )

    captured_warnings <- c(
      captured_warnings,
      inline_msg
    )

    message(inline_msg)
  }

  family_name <- tryCatch(family$family, error = function(e) NA_character_)
  is_tweedie_like <- tolower(family_name) %in%
    c(
      "tweedie",
      "tw",
      "twlss"
    )

  if (is_tweedie_like) {
    if (!requireNamespace("tweedie", quietly = TRUE)) {
      stop("Package 'tweedie' is required.")
    }
    if (!requireNamespace("statmod", quietly = TRUE)) {
      stop("Package 'statmod' is required.")
    }

    resp_var <- tryCatch(
      {
        all.vars(as.formula(
          if (is.list(model_formula)) model_formula[[1]] else model_formula
        ))[1]
      },
      error = function(e) NULL
    )

    if (
      !is.null(resp_var) &&
        resp_var %in% names(data) &&
        any(data[[resp_var]] < 0, na.rm = TRUE)
    ) {
      stop(
        stringr::str_glue(
          "FATAL: Tweedie response '{resp_var}' contains negative values, ",
          "which violates Tweedie distributional support (y >= 0)."
        ),
        call. = FALSE
      )
    }
  }

  # Report the start of the model fitting process:
  fitter_label <- if (identical(package, "mgcv")) {
    paste0(
      "mgcv::",
      mgcv_fitter
    )
  } else {
    "glmmTMB::glmmTMB"
  }

  cat(
    stringr::str_glue(
      "\nStarted fitting '{model_name}' via {fitter_label} at {Sys.time()}, ",
      "requesting up to {auto_cores} fitting thread(s)...\n"
    )
  )

  caller_env <- parent.frame()

  # Retain the exact evaluated analysis data in an environment carried by the
  # fitted formula. This allows the stored model call to use a compact data symbol
  # and remain usable after the wrapper returns or the model is serialized.
  model_env <- new.env(
    parent = caller_env
  )

  if (
    identical(package, "glmmTMB") &&
      isTRUE(glmmTMB_has_smooth)
  ) {
    assign(
      x = "s",
      value = mgcv::s,
      envir = model_env
    )
  }

  data_call_name <- if (
    is.symbol(true_data_name) &&
      length(true_data_name) == 1L &&
      nzchar(as.character(true_data_name))
  ) {
    as.character(true_data_name)
  } else {
    ".ecotools_analysis_data"
  }

  assign(
    x = data_call_name,
    value = data,
    envir = model_env
  )

  # Preserve separately supplied offsets in the persistent model environment.
  #
  # glmmTMB may internally refer to an evaluated offset vector through the
  # temporary symbol '..offset'. Prediction, simulation, model-frame
  # reconstruction, and DHARMa diagnostics can fail after the wrapper returns
  # unless that symbol remains available in the fitted formula environment.
  evaluated_offset <- NULL

  if (!is.null(offset)) {
    evaluated_offset <- tryCatch(
      {
        offset_values <- as.numeric(offset)

        if (length(offset_values) == 1L) {
          offset_values <- rep(
            offset_values,
            nrow(data)
          )
        }

        if (length(offset_values) != nrow(data)) {
          stop(
            paste0(
              "The evaluated 'offset' must have length 1 or match the ",
              "number of rows in 'data' (",
              nrow(data),
              "). Received length ",
              length(offset_values),
              "."
            ),
            call. = FALSE
          )
        }

        if (
          anyNA(offset_values) ||
            any(!is.finite(offset_values))
        ) {
          stop(
            "The evaluated 'offset' contains missing or non-finite values.",
            call. = FALSE
          )
        }

        offset_values
      },
      error = function(e) {
        stop(
          paste0(
            "Could not evaluate and validate the model offset: ",
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    assign(
      x = "..offset",
      value = evaluated_offset,
      envir = model_env
    )

    assign(
      x = ".ecotools_offset",
      value = evaluated_offset,
      envir = model_env
    )
  }

  set_formula_environment <- function(
    formula_object,
    target_environment
  ) {
    formula_result <- if (inherits(formula_object, "formula")) {
      formula_object
    } else {
      stats::as.formula(
        formula_object,
        env = target_environment
      )
    }

    environment(formula_result) <- target_environment

    formula_result
  }

  eval_formula <- if (is.list(model_formula)) {
    lapply(
      model_formula,
      set_formula_environment,
      target_environment = model_env
    )
  } else {
    set_formula_environment(
      formula_object = model_formula,
      target_environment = model_env
    )
  }

  # Zero inflation formula
  eval_zi_formula <- set_formula_environment(
    formula_object = zi_formula,
    target_environment = model_env
  )

  eval_disp_formula <- set_formula_environment(
    formula_object = disp_formula,
    target_environment = model_env
  )

  # ---------------------------------------------------------------------------
  # 3. Fit model
  # ---------------------------------------------------------------------------
  fit_error <- NULL
  simulation_error <- NULL

  tweedie_simulation_guard <- list(
    is_tweedie = FALSE,
    model_family = NA_character_,
    power = NA_real_,
    min_power = NA_real_,
    max_power = NA_real_,
    n_power_values = 0L,
    threshold = tweedie_simulation_max_power,
    blocked = FALSE,
    reason = "not_evaluated",
    message = NULL
  )

  dharma_status <- "not_attempted"

  fitted_model <- tryCatch(
    withCallingHandlers(
      {
        if (identical(package, "glmmTMB")) {
          fit_args <- list(
            formula = eval_formula,
            ziformula = eval_zi_formula,
            dispformula = eval_disp_formula,
            family = family,
            data = as.name(
              data_call_name
            ),
            control = glmmTMB::glmmTMBControl(
              parallel = list(
                n = auto_cores,
                autopar = TRUE
              )
            )
          )

          if (!is.null(weights)) {
            fit_args$weights <- weights
          }

          if (!is.null(evaluated_offset)) {
            fit_args$offset <- evaluated_offset
          }

          fit_args <- c(
            fit_args,
            extra_args
          )

          do.call(
            what = glmmTMB::glmmTMB,
            args = fit_args,
            quote = FALSE,
            envir = model_env
          )
        } else if (identical(package, "mgcv")) {
          # mgcv fitting branch
          if (identical(mgcv_fitter, "gam")) {
            # mgcv::gam() fitting branch
            fit_args <- list(
              formula = eval_formula,
              family = family,
              data = as.name(data_call_name),
              control = mgcv::gam.control(
                nthreads = auto_cores
              )
            )

            if (!is.null(weights)) {
              fit_args$weights <- weights
            }

            if (!is.null(evaluated_offset)) {
              fit_args$offset <- evaluated_offset
            }

            fit_args <- c(
              fit_args,
              extra_args
            )

            do.call(
              what = mgcv::gam,
              args = fit_args,
              quote = FALSE,
              envir = model_env
            )
          } else if (identical(mgcv_fitter, "bam")) {
            # mgcv::bam() fitting branch
            fit_args <- list(
              formula = eval_formula,
              family = family,
              data = as.name(data_call_name),
              control = mgcv::gam.control(),
              nthreads = auto_cores
            )

            if (!is.null(weights)) {
              fit_args$weights <- weights
            }

            if (!is.null(evaluated_offset)) {
              fit_args$offset <- evaluated_offset
            }

            fit_args <- c(
              fit_args,
              extra_args
            )

            do.call(
              what = mgcv::bam,
              args = fit_args,
              quote = FALSE,
              envir = model_env
            )
          } else {
            stop(
              paste0(
                "Unsupported mgcv fitter: '",
                mgcv_fitter,
                "'."
              ),
              call. = FALSE
            )
          }
        }
      },
      warning = function(w) {
        warning_message <- conditionMessage(w)

        captured_warnings <<- c(
          captured_warnings,
          warning_message
        )

        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_call <- conditionCall(e)

      fit_error <<- list(
        message = conditionMessage(e),
        class = class(e),
        call = if (!is.null(error_call)) {
          paste(
            deparse(error_call),
            collapse = " "
          )
        } else {
          NULL
        },
        calls = vapply(
          sys.calls(),
          function(call_object) {
            paste(
              deparse(call_object),
              collapse = " "
            )
          },
          character(1)
        )
      )

      NULL
    }
  )

  if (is.null(fitted_model)) {
    if (
      identical(package, "glmmTMB") &&
        isTRUE(glmmTMB_has_smooth) &&
        !is.null(fit_error$message)
    ) {
      fit_error$message <- paste0(
        fit_error$message,
        " EcoTools detected an experimental glmmTMB smooth-term model. ",
        "Use an unqualified s() term, verify that current versions of glmmTMB ",
        "and mgcv are installed, use REML = TRUE unless a documented analysis ",
        "requires otherwise, and confirm the smooth basis dimension is supported ",
        "by the available unique covariate values."
      )
    }
    error_message <- if (
      !is.null(fit_error$message) &&
        nzchar(fit_error$message)
    ) {
      fit_error$message
    } else {
      "No specific fitting error was captured."
    }

    error_call_text <- if (
      !is.null(fit_error$call) &&
        nzchar(fit_error$call)
    ) {
      paste0(
        "\n  Failing call: ",
        fit_error$call
      )
    } else {
      ""
    }

    fitting_failure_message <- paste0(
      "\nFATAL ERROR fitting '",
      model_name,
      "' via ",
      package,
      ".\n",
      "  Underlying error: ",
      error_message,
      error_call_text,
      "\n",
      "  The failed fit has been cached for inspection, but it is not ",
      "valid for inference.\n"
    )

    captured_warnings <- unique(
      c(
        captured_warnings,
        paste0(
          "FITTING ERROR: ",
          error_message
        )
      )
    )

    warning(
      fitting_failure_message,
      call. = FALSE,
      immediate. = TRUE
    )

    tweedie_simulation_guard$reason <- "fit_failed"
    tweedie_simulation_guard$message <- paste0(
      "The Tweedie simulation-safety check was not performed because model ",
      "fitting failed."
    )

    dharma_status <- "not_attempted_fit_failed"

    output_bundle <- list(
      model = NULL,
      residuals = NULL,
      tweedie_simulation_guard = tweedie_simulation_guard,
      dharma_status = dharma_status,
      converged = FALSE,
      numerically_converged = FALSE,
      model_status = "invalid",
      valid_for_inference = FALSE,
      valid_for_model_selection = FALSE,
      diagnose = NULL,
      fit_error = fit_error,
      warnings = captured_warnings,
      sparse_screen = NULL,
      design_check = design_check_result,
      analysis_data = data,
      offset_metadata = list(
        supplied_separately = !is.null(evaluated_offset),
        expression = if (!is.null(wrapper_call$offset)) {
          paste(
            deparse(
              wrapper_call$offset,
              width.cutoff = 500L
            ),
            collapse = ""
          )
        } else {
          NULL
        },
        length = if (is.null(evaluated_offset)) {
          0L
        } else {
          length(evaluated_offset)
        },
        finite = if (is.null(evaluated_offset)) {
          NA
        } else {
          all(
            is.finite(
              evaluated_offset
            )
          )
        },
        storage_symbol = if (is.null(evaluated_offset)) {
          NULL
        } else {
          ".ecotools_offset"
        }
      ),
      data_expression = paste(
        deparse(true_data_name),
        collapse = " "
      ),
      data_call_name = data_call_name,
      model_call_data_storage = if (
        identical(package, "glmmTMB") &&
          isTRUE(glmmTMB_has_smooth)
      ) {
        "embedded_evaluated_data"
      } else {
        "persistent_environment_symbol"
      },
      glmmTMB_smooth_terms = if (identical(package, "glmmTMB")) {
        isTRUE(
          glmmTMB_has_smooth
        )
      } else {
        FALSE
      },
      fitting_engine = package,
      mgcv_fitter = if (identical(package, "mgcv")) {
        mgcv_fitter
      } else {
        NULL
      },
      data_dim = dim(data),
      simulation_error = simulation_error,
      analysis_fingerprint = current_fingerprint,
      session_info = utils::sessionInfo()
    )

    class(output_bundle) <- c(
      "ecotools_model_bundle",
      "list"
    )

    saveRDS(
      output_bundle,
      file = file_path
    )

    return(output_bundle)
  }

  # -------------------------------------------------------------------------
  # glmmTMB serialization repair
  #
  # glmmTMB resets formula environments during fitting.
  # Reattach the analysis-data environment so that:
  #
  #   * predict()
  #   * DHARMa::simulateResiduals()
  #   * saveRDS()/readRDS()
  #
  # continue to work when the original data object was local to a function.
  # -------------------------------------------------------------------------
  if (identical(package, "glmmTMB")) {
    if (!is.null(evaluated_offset)) {
      assign(
        x = "..offset",
        value = evaluated_offset,
        envir = model_env
      )

      assign(
        x = ".ecotools_offset",
        value = evaluated_offset,
        envir = model_env
      )
    }

    if (
      !is.null(fitted_model$call$formula) &&
        inherits(
          fitted_model$call$formula,
          "formula"
        )
    ) {
      environment(
        fitted_model$call$formula
      ) <- model_env
    }

    if (
      !is.null(fitted_model$call$ziformula) &&
        inherits(
          fitted_model$call$ziformula,
          "formula"
        )
    ) {
      environment(
        fitted_model$call$ziformula
      ) <- model_env
    }

    if (
      !is.null(fitted_model$call$dispformula) &&
        inherits(
          fitted_model$call$dispformula,
          "formula"
        )
    ) {
      environment(
        fitted_model$call$dispformula
      ) <- model_env
    }
  }

  if (package == "mgcv" && identical(fitted_model$family$family, "twlss")) {
    class(fitted_model) <- c("gam_twlss", class(fitted_model))
  }

  # ---------------------------------------------------------------------------
  # 4. Sparse-data screen
  # ---------------------------------------------------------------------------
  sparse_screen_obj <- NULL

  if (isTRUE(sparse_screen)) {
    sparse_screen_obj <- tryCatch(
      {
        screen_sparse_count_model(
          model = fitted_model,
          data = as.data.frame(data),
          min_events_per_parameter = sparse_min_events_per_parameter,
          min_events_per_level = sparse_min_events_per_level,
          event_basis = sparse_event_basis,
          count_like_only = sparse_count_like_only
        )
      },
      error = function(e) {
        list(
          pass = NA,
          warning = paste(
            "SPARSE DATA SCREEN ERROR:",
            conditionMessage(e)
          ),
          details = NULL
        )
      }
    )

    if (
      !is.null(sparse_screen_obj$pass) &&
        length(sparse_screen_obj$pass) == 1L &&
        is.na(sparse_screen_obj$pass)
    ) {
      stop(
        paste0(
          "The sparse-data screen failed and model validity could not be assessed.\n",
          sparse_screen_obj$warning
        ),
        call. = FALSE
      )
    }

    if (!is.null(sparse_screen_obj$warning)) {
      captured_warnings <- c(
        captured_warnings,
        sparse_screen_obj$warning
      )

      warning(
        sparse_screen_obj$warning,
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 5. Convergence, Hessian, and extreme-coefficient checks
  # ---------------------------------------------------------------------------
  is_converged <- TRUE
  has_extreme <- FALSE

  if (package == "glmmTMB") {
    if (
      !is.null(fitted_model$fit$convergence) &&
        fitted_model$fit$convergence != 0
    ) {
      is_converged <- FALSE
      opt_msg <- stringr::str_glue(
        "OPTIMIZER FAILURE: Code {fitted_model$fit$convergence} - {fitted_model$fit$message}"
      )
      captured_warnings <- c(captured_warnings, opt_msg)
      warning(opt_msg, call. = FALSE)
    }

    pd_hess <- if (!is.null(fitted_model$sdr)) {
      isTRUE(fitted_model$sdr$pdHess)
    } else {
      FALSE
    }

    ll_val <- suppressWarnings(tryCatch(
      logLik(fitted_model),
      error = function(e) NA
    ))

    if (!pd_hess || is.na(ll_val)) {
      is_converged <- FALSE
      hess_msg <- "HESSIAN FAILURE: Non-positive definite Hessian or NA log-likelihood."
      captured_warnings <- c(captured_warnings, hess_msg)
      warning(hess_msg, call. = FALSE)
    }

    fix_ef <- tryCatch(glmmTMB::fixef(fitted_model)$cond, error = function(e) {
      NULL
    })

    if (!is.null(fix_ef)) {
      extreme_ef <- fix_ef[abs(fix_ef) > separation_threshold]

      if (length(extreme_ef) > 0) {
        has_extreme <- TRUE
        sep_msg <- stringr::str_glue(
          "SEPARATION / SPARSITY WARNING: Extreme parameters (|beta| > {separation_threshold}):"
        )
        captured_warnings <- c(captured_warnings, sep_msg)
        warning(sep_msg, call. = FALSE)
        print(extreme_ef)
      }
    }

    re_var <- tryCatch(
      {
        vc <- as.data.frame(glmmTMB::VarCorr(fitted_model)$cond)
        if (!is.null(vc$vcov)) vc$vcov else numeric(0)
      },
      error = function(e) numeric(0)
    )

    if (length(re_var) > 0 && any(re_var < 1e-8, na.rm = TRUE)) {
      singular_msg <- "SINGULAR FIT: One or more random-effect variances are near zero."
      captured_warnings <- c(captured_warnings, singular_msg)
      warning(paste0("\n", singular_msg, "\n"), call. = FALSE)
    }
  } else if (package == "mgcv") {
    if (!is.null(fitted_model$converged) && !fitted_model$converged) {
      is_converged <- FALSE
      opt_msg <- "OPTIMIZER FAILURE: GAM optimizer failed to converge."
      captured_warnings <- c(captured_warnings, opt_msg)
      warning(opt_msg, call. = FALSE)
    }

    if (!is.null(fitted_model$outer.info)) {
      outer_conv <- tryCatch(fitted_model$outer.info$conv, error = function(e) {
        NULL
      })

      if (
        !is.null(outer_conv) &&
          !grepl("full convergence", outer_conv, ignore.case = TRUE)
      ) {
        mgcv_msg <- stringr::str_glue(
          "MGCV WARNING: Outer optimizer convergence message: {outer_conv}"
        )
        captured_warnings <- c(captured_warnings, mgcv_msg)
        warning(mgcv_msg, call. = FALSE)
      }
    }

    gam_coefs <- tryCatch(coef(fitted_model), error = function(e) NULL)

    if (!is.null(gam_coefs)) {
      coef_names <- names(gam_coefs)
      is_smooth_coef <- grepl("^(s|te|ti|t2)\\(", coef_names)
      param_coefs <- gam_coefs[!is_smooth_coef]
      extreme_gam <- param_coefs[
        is.finite(param_coefs) & abs(param_coefs) > separation_threshold
      ]

      if (length(extreme_gam) > 0) {
        has_extreme <- TRUE
        sep_msg <- stringr::str_glue(
          "SEPARATION / SPARSITY WARNING in GAM parametric terms (|beta| > {separation_threshold}):"
        )
        captured_warnings <- c(captured_warnings, sep_msg)
        warning(paste0("\n", sep_msg, "\n"), call. = FALSE)
        print(extreme_gam)
      }
    }
  }

  model_status <- if (!is_converged) {
    "invalid"
  } else if (has_extreme) {
    "suspicious_extreme_coefficients"
  } else {
    "valid"
  }

  diag_obj <- NULL

  if (
    identical(package, "glmmTMB") &&
      (!isTRUE(is_converged) ||
        isTRUE(has_extreme))
  ) {
    diag_obj <- tryCatch(
      glmmTMB::diagnose(fitted_model),
      error = function(e) NULL
    )
  }

  # ---------------------------------------------------------------------------
  # 6. Repair model call for downstream tools
  # ---------------------------------------------------------------------------
  if (
    package %in%
      c(
        "glmmTMB",
        "mgcv"
      )
  ) {
    # do.call() may store the fitting function and evaluated data directly in the
    # model call. Restore a readable function head, formula, and compact data symbol.
    # The formula environment contains the exact evaluated analysis data under this
    # symbol, so call reconstruction remains possible after serialization.
    if (!is.null(fitted_model$call)) {
      if (identical(package, "glmmTMB")) {
        fitted_model$call[[1L]] <-
          quote(glmmTMB::glmmTMB)
      } else if (
        identical(package, "mgcv") &&
          identical(mgcv_fitter, "gam")
      ) {
        fitted_model$call[[1L]] <-
          quote(mgcv::gam)
      } else if (
        identical(package, "mgcv") &&
          identical(mgcv_fitter, "bam")
      ) {
        fitted_model$call[[1L]] <-
          quote(mgcv::bam)
      }
    }

    if (!is.null(fitted_model$call$formula)) {
      fitted_model$call$formula <- eval_formula
    }

    if (
      identical(package, "glmmTMB") &&
        !is.null(fitted_model$call$ziformula)
    ) {
      fitted_model$call$ziformula <- eval_zi_formula
    }

    if (
      identical(package, "glmmTMB") &&
        !is.null(fitted_model$call$dispformula)
    ) {
      fitted_model$call$dispformula <- eval_disp_formula
    }

    if (!is.null(fitted_model$call)) {
      if (
        identical(package, "glmmTMB") &&
          isTRUE(glmmTMB_has_smooth)
      ) {
        # glmmTMB population-level prediction and some simulation paths rebuild
        # their model frames from the stored call. Embedding the evaluated data
        # prevents reconstruction from depending on a caller-local data symbol
        # that may no longer exist after the wrapper returns or after RDS
        # serialization.
        fitted_model$call$data <- data
      } else {
        fitted_model$call$data <- as.name(
          data_call_name
        )
      }
    }

    if (is.null(weights)) {
      fitted_model$call$weights <- NULL
    } else if (!is.null(wrapper_call$weights)) {
      fitted_model$call$weights <-
        wrapper_call$weights
    }

    if (is.null(evaluated_offset)) {
      fitted_model$call$offset <- NULL
    } else {
      fitted_model$call$offset <- as.name(".ecotools_offset")
    }
  }

  smooth_reconstruction_check <- NULL

  if (
    identical(package, "glmmTMB") &&
      isTRUE(glmmTMB_has_smooth)
  ) {
    smooth_reconstruction_check <- tryCatch(
      {
        reconstruction_predictions <- stats::predict(
          fitted_model,
          type = "response",
          re.form = ~0
        )

        list(
          pass = length(reconstruction_predictions) == nrow(data) &&
            all(
              is.finite(
                reconstruction_predictions
              )
            ),
          n_predictions = length(
            reconstruction_predictions
          ),
          error = NULL
        )
      },
      error = function(e) {
        list(
          pass = FALSE,
          n_predictions = NA_integer_,
          error = conditionMessage(e)
        )
      }
    )

    if (!isTRUE(smooth_reconstruction_check$pass)) {
      reconstruction_message <- paste0(
        "GLMMTMB SMOOTH RECONSTRUCTION FAILED: population-level prediction ",
        "could not reconstruct the experimental smooth model",
        if (!is.null(smooth_reconstruction_check$error)) {
          paste0(
            ": ",
            smooth_reconstruction_check$error
          )
        } else {
          "."
        }
      )

      captured_warnings <- unique(
        c(
          captured_warnings,
          reconstruction_message
        )
      )

      warning(
        reconstruction_message,
        call. = FALSE
      )
    }
  }

  if (!identical(fitted_model$family$family, "twlss")) {
    tweedie_simulation_guard <- check_tweedie_simulation_guard(
      model = fitted_model,
      max_power = tweedie_simulation_max_power
    )
  }

  # ---------------------------------------------------------------------------
  # 7. DHARMa residual simulation
  # ---------------------------------------------------------------------------
  is_twlss_model <-
    identical(package, "mgcv") &&
    identical(
      fitted_model$family$family,
      "twlss"
    )

  mgcv_simulation_available <- TRUE
  mgcv_simulation_error <- NULL

  if (
    identical(package, "mgcv") &&
      !isTRUE(is_twlss_model)
  ) {
    mgcv_simulation_available <- requireNamespace(
      "mgcViz",
      quietly = TRUE
    )

    if (!isTRUE(mgcv_simulation_available)) {
      mgcv_simulation_error <- paste0(
        "Package 'mgcViz' is required to simulate responses from ",
        "ordinary mgcv GAM objects for DHARMa diagnostics."
      )

      captured_warnings <- c(
        captured_warnings,
        mgcv_simulation_error
      )

      warning(
        mgcv_simulation_error,
        call. = FALSE
      )
    }
  }

  ll_check <- suppressWarnings(tryCatch(
    logLik(fitted_model),
    error = function(e) NA
  ))

  if (
    is.na(ll_check) ||
      identical(model_status, "invalid")
  ) {
    model_res <- NULL
    dharma_status <- "not_attempted_invalid_fit"
  } else if (isTRUE(is_twlss_model)) {
    preds <- stats::predict(
      fitted_model,
      type = "response"
    )

    mu <- preds[, 1L]
    phi <- exp(preds[, 2L])
    theta <- preds[, 3L]

    fam_a <- if (!is.null(fitted_model$family$a)) {
      fitted_model$family$a
    } else {
      1.01
    }

    fam_b <- if (!is.null(fitted_model$family$b)) {
      fitted_model$family$b
    } else {
      1.99
    }

    fitted_power <- (fam_a +
      fam_b * exp(theta)) /
      (1 +
        exp(theta))

    tweedie_simulation_guard <- check_tweedie_simulation_guard(
      model = fitted_model,
      max_power = tweedie_simulation_max_power,
      fitted_power = fitted_power
    )

    if (isTRUE(tweedie_simulation_guard$blocked)) {
      captured_warnings <- unique(
        c(
          captured_warnings,
          tweedie_simulation_guard$message
        )
      )

      warning(
        tweedie_simulation_guard$message,
        call. = FALSE
      )

      model_res <- NULL
      dharma_status <- "skipped_tweedie_power_boundary"
    } else {
      message(
        "-> Running direct twlss simulation engine via tweedie::rtweedie..."
      )

      set.seed(dharma_seed)

      n_observations <- length(mu)

      sim_matrix <- matrix(
        NA_real_,
        nrow = n_observations,
        ncol = dharma_n_sim
      )

      if (length(unique(fitted_power)) == 1L) {
        power_scalar <- fitted_power[1L]

        for (simulation_index in seq_len(dharma_n_sim)) {
          sim_matrix[, simulation_index] <- tweedie::rtweedie(
            n_observations,
            mu = mu,
            phi = phi,
            power = power_scalar
          )
        }
      } else {
        for (simulation_index in seq_len(dharma_n_sim)) {
          sim_matrix[, simulation_index] <- mapply(
            FUN = function(mean_value, phi_value, power_value) {
              tweedie::rtweedie(
                1L,
                mu = mean_value,
                phi = phi_value,
                power = power_value
              )
            },
            mean_value = mu,
            phi_value = phi,
            power_value = fitted_power
          )
        }
      }

      model_res <- DHARMa::createDHARMa(
        simulatedResponse = sim_matrix,
        observedResponse = fitted_model$y,
        fittedPredictedResponse = mu,
        integerResponse = FALSE
      )

      dharma_status <- "completed"
    }
  } else if (
    identical(package, "mgcv") &&
      !isTRUE(mgcv_simulation_available)
  ) {
    model_res <- NULL
    dharma_status <- "unavailable_mgcViz_missing"
  } else if (isTRUE(tweedie_simulation_guard$blocked)) {
    captured_warnings <- unique(
      c(
        captured_warnings,
        tweedie_simulation_guard$message
      )
    )

    warning(
      tweedie_simulation_guard$message,
      call. = FALSE
    )

    model_res <- NULL
    dharma_status <- "skipped_tweedie_power_boundary"
  } else if (
    identical(package, "glmmTMB") &&
      isTRUE(glmmTMB_has_smooth)
  ) {
    if (
      tolower(
        fitted_model$modelInfo$family$family
      ) %in%
        c(
          "binomial",
          "betabinomial"
        )
    ) {
      stop(
        paste0(
          "Direct DHARMa construction for experimental glmmTMB smooth models ",
          "is not yet implemented for binomial-family responses."
        ),
        call. = FALSE
      )
    }
    model_res <- tryCatch(
      {
        set.seed(
          dharma_seed
        )

        simulated_responses <- stats::simulate(
          fitted_model,
          nsim = dharma_n_sim,
          seed = dharma_seed
        )

        simulated_matrix <- do.call(
          cbind,
          lapply(
            simulated_responses,
            function(simulated_response) {
              if (is.matrix(simulated_response)) {
                simulated_response[, 1L]
              } else {
                as.numeric(
                  simulated_response
                )
              }
            }
          )
        )

        observed_response <- stats::model.response(
          stats::model.frame(
            fitted_model
          )
        )

        fitted_response <- stats::predict(
          fitted_model,
          type = "response"
        )

        integer_response <- all(
          is.finite(observed_response) &
            observed_response >= 0 &
            abs(
              observed_response -
                round(observed_response)
            ) <=
              sqrt(.Machine$double.eps)
        )

        residual_object <- DHARMa::createDHARMa(
          simulatedResponse = simulated_matrix,
          observedResponse = observed_response,
          fittedPredictedResponse = fitted_response,
          integerResponse = integer_response
        )

        dharma_status <-
          "completed_unconditional_glmmTMB_smooth"

        residual_object
      },
      error = function(e) {
        simulation_error <<- list(
          message = conditionMessage(e),
          class = class(e),
          call = {
            error_call <- conditionCall(e)

            if (is.null(error_call)) {
              NULL
            } else {
              paste(
                deparse(error_call),
                collapse = " "
              )
            }
          }
        )

        simulation_error_message <- paste0(
          "DHARMa simulation error for experimental glmmTMB smooth model '",
          model_name,
          "': ",
          simulation_error$message
        )

        captured_warnings <<- c(
          captured_warnings,
          simulation_error_message
        )

        dharma_status <<-
          "simulation_failed_glmmTMB_smooth"

        warning(
          simulation_error_message,
          call. = FALSE
        )

        NULL
      }
    )
  } else {
    model_res <- tryCatch(
      {
        residual_object <- DHARMa::simulateResiduals(
          fittedModel = fitted_model,
          n = dharma_n_sim,
          seed = dharma_seed
        )

        dharma_status <- "completed"

        residual_object
      },
      error = function(e) {
        simulation_error <<- list(
          message = conditionMessage(e),
          class = class(e),
          call = {
            error_call <- conditionCall(e)

            if (is.null(error_call)) {
              NULL
            } else {
              paste(
                deparse(error_call),
                collapse = " "
              )
            }
          }
        )

        simulation_error_message <- paste0(
          "DHARMa simulation error for '",
          model_name,
          "': ",
          simulation_error$message
        )

        captured_warnings <<- c(
          captured_warnings,
          simulation_error_message
        )

        dharma_status <<- "simulation_failed"

        warning(
          simulation_error_message,
          call. = FALSE
        )

        NULL
      }
    )
  }

  # ---------------------------------------------------------------------------
  # 7b. Formal DHARMa residual-validation gate
  # ---------------------------------------------------------------------------

  residual_diagnostics <- NULL
  residual_diagnostics_valid <- NA
  required_tests_available <- FALSE

  if (
    !is.null(model_res) &&
      !identical(model_status, "invalid")
  ) {
    dharma_uniformity <- tryCatch(
      DHARMa::testUniformity(
        model_res,
        plot = FALSE
      ),
      error = function(e) NULL
    )

    dharma_dispersion <- tryCatch(
      DHARMa::testDispersion(
        model_res,
        alternative = "two.sided",
        plot = FALSE
      ),
      error = function(e) NULL
    )

    dharma_outliers <- tryCatch(
      DHARMa::testOutliers(
        model_res,
        nBoot = min(
          100L,
          dharma_n_sim
        ),
        plot = FALSE
      ),
      error = function(e) NULL
    )

    fitted_predictions <- tryCatch(
      model_res$fittedPredictedResponse,
      error = function(e) NULL
    )

    if (
      is.null(fitted_predictions) ||
        length(fitted_predictions) != length(model_res$scaledResiduals)
    ) {
      fitted_predictions <- tryCatch(
        stats::predict(
          fitted_model,
          type = "response"
        ),
        error = function(e) NULL
      )
    }

    dharma_quantile_warning <- NULL

    dharma_quantiles <- if (
      !is.null(fitted_predictions) &&
        length(fitted_predictions) == length(model_res$scaledResiduals)
    ) {
      quantile_result <- tryCatch(
        withCallingHandlers(
          DHARMa::testQuantiles(
            model_res,
            predictor = fitted_predictions,
            plot = FALSE
          ),
          warning = function(w) {
            dharma_quantile_warning <<-
              conditionMessage(w)

            invokeRestart(
              "muffleWarning"
            )
          }
        ),
        error = function(e) {
          dharma_quantile_warning <<-
            conditionMessage(e)

          NULL
        }
      )

      if (!is.null(dharma_quantile_warning)) {
        quantile_warning_message <- paste0(
          "DHARMa predicted-quantile test was unavailable or unreliable: ",
          dharma_quantile_warning
        )

        captured_warnings <- unique(
          c(
            captured_warnings,
            quantile_warning_message
          )
        )

        NULL
      } else {
        quantile_result
      }
    } else {
      NULL
    }

    extract_test_p <- function(test_object) {
      if (
        is.null(test_object) ||
          is.null(test_object$p.value) ||
          length(test_object$p.value) == 0L
      ) {
        return(NA_real_)
      }

      p_value <- suppressWarnings(
        as.numeric(test_object$p.value[1])
      )

      if (
        length(p_value) != 1L ||
          !is.finite(p_value)
      ) {
        return(NA_real_)
      }

      p_value
    }

    residual_p_values <- c(
      uniformity = extract_test_p(
        dharma_uniformity
      ),
      dispersion = extract_test_p(
        dharma_dispersion
      ),
      outliers = extract_test_p(
        dharma_outliers
      ),
      predicted_quantiles = extract_test_p(
        dharma_quantiles
      )
    )

    residual_tests_available <- is.finite(
      residual_p_values
    )

    residual_p_values_adjusted <- stats::p.adjust(
      residual_p_values,
      method = "holm"
    )

    names(residual_p_values_adjusted) <- names(
      residual_p_values
    )

    residual_tests_failed <- residual_tests_available &
      residual_p_values_adjusted < residual_alpha
    dispersion_ratio <- if (
      !is.null(dharma_dispersion) &&
        !is.null(dharma_dispersion$statistic) &&
        length(dharma_dispersion$statistic) > 0L
    ) {
      suppressWarnings(
        as.numeric(
          dharma_dispersion$statistic[1]
        )
      )
    } else {
      NA_real_
    }

    required_residual_tests <- c(
      "uniformity",
      "dispersion",
      "outliers",
      "predicted_quantiles"
    )

    required_tests_available <- all(
      residual_tests_available[
        required_residual_tests
      ]
    )

    residual_diagnostics_valid <-
      required_tests_available &&
      !any(
        residual_tests_failed[
          required_residual_tests
        ]
      )

    residual_diagnostics <- list(
      valid = residual_diagnostics_valid,
      alpha = residual_alpha,
      n_sim = dharma_n_sim,
      seed = dharma_seed,
      p_values = residual_p_values,
      p_values_adjusted = residual_p_values_adjusted,
      adjustment_method = "holm",
      tests_available = residual_tests_available,
      required_tests = required_residual_tests,
      required_tests_available = required_tests_available,
      tests_failed = names(
        residual_p_values
      )[residual_tests_failed],
      dispersion_ratio = dispersion_ratio,
      uniformity_test = dharma_uniformity,
      dispersion_test = dharma_dispersion,
      outlier_test = dharma_outliers,
      predicted_quantile_test = dharma_quantiles,
      predicted_quantile_warning = dharma_quantile_warning
    )

    if (any(residual_tests_failed)) {
      residual_failure_parts <- paste0(
        names(residual_p_values)[
          residual_tests_failed
        ],
        " raw p = ",
        format.pval(
          residual_p_values[
            residual_tests_failed
          ],
          digits = 4,
          eps = 0.001
        ),
        ", Holm-adjusted p = ",
        format.pval(
          residual_p_values_adjusted[
            residual_tests_failed
          ],
          digits = 4,
          eps = 0.001
        )
      )

      residual_failure_msg <- paste0(
        "DHARMa RESIDUAL VALIDATION FAILED: ",
        paste(
          residual_failure_parts,
          collapse = "; "
        ),
        if (is.finite(dispersion_ratio)) {
          paste0(
            ". Observed-to-simulated dispersion ratio = ",
            signif(dispersion_ratio, 5)
          )
        } else {
          ""
        },
        "."
      )

      captured_warnings <- c(
        captured_warnings,
        residual_failure_msg
      )

      warning(
        residual_failure_msg,
        call. = FALSE
      )
    }
  }

  # ---------------------------------------------------------------------------
  # 7c. Final status after residual validation
  # ---------------------------------------------------------------------------

  model_status <- if (!is_converged) {
    "invalid"
  } else if (has_extreme) {
    "suspicious_extreme_coefficients"
  } else if (
    is.null(model_res) ||
      !isTRUE(required_tests_available)
  ) {
    "valid_fit_no_residual_diagnostics"
  } else if (!isTRUE(residual_diagnostics_valid)) {
    "suspicious_residual_diagnostics"
  } else {
    "valid"
  }

  if (isTRUE(tweedie_simulation_guard$blocked)) {
    model_status <- "invalid"

    captured_warnings <- unique(
      c(
        captured_warnings,
        tweedie_simulation_guard$message
      )
    )
  }

  # ---------------------------------------------------------------------------
  # 8. Save model bundle
  # ---------------------------------------------------------------------------

  # This is the successful post-fit path. A suspicious diagnostic status does not
  # mean that the fitted model should be removed from the returned bundle.
  if (is.null(fitted_model)) {
    stop(
      paste0(
        "INTERNAL WRAPPER ERROR: Bundle construction was reached with ",
        "'fitted_model = NULL'."
      ),
      call. = FALSE
    )
  }

  output_bundle <- list(
    model = fitted_model,
    residuals = model_res,
    tweedie_simulation_guard = tweedie_simulation_guard,
    dharma_status = dharma_status,
    converged = is_converged,
    numerically_converged = is_converged,
    model_status = model_status,
    valid_for_inference = identical(
      model_status,
      "valid"
    ),
    valid_for_model_selection = identical(
      model_status,
      "valid"
    ),
    diagnose = diag_obj,
    fit_error = NULL,
    design_check = design_check_result,
    residual_diagnostics = residual_diagnostics,
    residual_diagnostics_valid = residual_diagnostics_valid,
    alpha = residual_alpha,
    n_sim = dharma_n_sim,
    seed = dharma_seed,
    analysis_data = data,
    offset_metadata = list(
      supplied_separately = !is.null(evaluated_offset),
      expression = if (!is.null(wrapper_call$offset)) {
        paste(
          deparse(
            wrapper_call$offset,
            width.cutoff = 500L
          ),
          collapse = ""
        )
      } else {
        NULL
      },
      length = if (is.null(evaluated_offset)) {
        0L
      } else {
        length(evaluated_offset)
      },
      finite = if (is.null(evaluated_offset)) {
        NA
      } else {
        all(
          is.finite(
            evaluated_offset
          )
        )
      },
      storage_symbol = if (is.null(evaluated_offset)) {
        NULL
      } else {
        ".ecotools_offset"
      }
    ),
    data_expression = paste(
      deparse(true_data_name),
      collapse = " "
    ),
    data_call_name = data_call_name,
    model_call_data_storage = if (
      identical(package, "glmmTMB") &&
        isTRUE(glmmTMB_has_smooth)
    ) {
      "embedded_evaluated_data"
    } else {
      "persistent_environment_symbol"
    },
    model_class = class(
      fitted_model
    ),
    glmmTMB_smooth_reconstruction = smooth_reconstruction_check,
    glmmTMB_smooth_terms = if (identical(package, "glmmTMB")) {
      isTRUE(
        glmmTMB_has_smooth
      )
    } else {
      FALSE
    },
    fitting_engine = package,
    mgcv_fitter = if (identical(package, "mgcv")) {
      mgcv_fitter
    } else {
      NULL
    },
    warnings = if (length(captured_warnings) > 0L) {
      unique(captured_warnings)
    } else {
      NULL
    },
    sparse_screen = sparse_screen_obj,
    data_dim = dim(data),
    simulation_error = simulation_error,
    analysis_fingerprint = current_fingerprint,
    session_info = utils::sessionInfo()
  )

  class(output_bundle) <- c(
    "ecotools_model_bundle",
    "list"
  )

  # Defensive checks prevent a successful fit from being cached without its model.
  stopifnot(
    !is.null(output_bundle$model),
    inherits(output_bundle$model, "glmmTMB") ||
      inherits(output_bundle$model, "gam")
  )

  # Make sure the fitted model has retained its class:
  if (
    identical(package, "glmmTMB") &&
      !inherits(output_bundle$model, "glmmTMB")
  ) {
    stop(
      "Internal error: fitted model lost its 'glmmTMB' class.",
      call. = FALSE
    )
  }

  if (
    identical(package, "mgcv") &&
      !inherits(output_bundle$model, "gam")
  ) {
    stop(
      "Internal error: fitted model lost its 'gam' class.",
      call. = FALSE
    )
  }

  saveRDS(
    output_bundle,
    file = file_path
  )

  return(output_bundle)
}
