#' Create a Stable GLLVM Cache Fingerprint
#'
#' @keywords internal
get_gllvm_fingerprint <- function(
  y,
  X,
  studyDesign,
  formula,
  TR = NULL,
  row.eff = FALSE,
  family,
  num.lv,
  sd.errors,
  method,
  control,
  control.start,
  gradient_threshold,
  coef_threshold,
  sparse_screen,
  sparse_action,
  abort_on,
  sparse_min_events_per_parameter,
  sparse_min_total_events,
  sparse_min_nonzero_rows,
  sparse_event_basis,
  sparse_count_like_only,
  design_check = TRUE,
  design_action = c("error", "warn", "drop"),
  design_tolerance = sqrt(.Machine$double.eps),
  reference_stability_check = "candidate",
  reference_stability_max_levels = 5L,
  reference_logLik_tolerance = 1e-6,
  reference_stability_seed = 1L,
  tmb_threads = 1L,
  dependency_versions,
  extra_args = list()
) {
  # Validate the number of TMB threads:
  if (
    !is.numeric(tmb_threads) ||
      length(tmb_threads) != 1L ||
      is.na(tmb_threads) ||
      !is.finite(tmb_threads) ||
      tmb_threads < 1L ||
      tmb_threads != as.integer(tmb_threads)
  ) {
    stop(
      "'tmb_threads' must be one positive integer.",
      call. = FALSE
    )
  }

  tmb_threads <- as.integer(tmb_threads)

  if (
    !is.numeric(reference_logLik_tolerance) ||
      length(reference_logLik_tolerance) != 1L ||
      is.na(reference_logLik_tolerance) ||
      !is.finite(reference_logLik_tolerance) ||
      reference_logLik_tolerance < 0
  ) {
    stop(
      paste0(
        "'reference_logLik_tolerance' must be one non-negative ",
        "finite number."
      ),
      call. = FALSE
    )
  }

  reference_logLik_tolerance <- as.numeric(
    reference_logLik_tolerance
  )

  if (!is.null(reference_stability_seed)) {
    if (
      !is.numeric(reference_stability_seed) ||
        length(reference_stability_seed) != 1L ||
        is.na(reference_stability_seed) ||
        !is.finite(reference_stability_seed) ||
        reference_stability_seed != as.integer(reference_stability_seed)
    ) {
      stop(
        paste0(
          "'reference_stability_seed' must be NULL or one finite ",
          "integer."
        ),
        call. = FALSE
      )
    }

    reference_stability_seed <- as.integer(
      reference_stability_seed
    )
  }

  if (!requireNamespace("digest", quietly = TRUE)) {
    stop(
      "Package 'digest' is required for cache fingerprinting.",
      call. = FALSE
    )
  }

  canonical_tabular <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }

    out <- as.data.frame(
      x,
      stringsAsFactors = FALSE
    )

    attr(out, ".internal.selfref") <- NULL
    row.names(out) <- seq_len(nrow(out))

    out
  }

  canonical_language <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }

    if (isFALSE(x)) {
      return(FALSE)
    }

    if (isTRUE(x)) {
      return(TRUE)
    }

    if (
      inherits(x, "formula") ||
        is.call(x) ||
        is.expression(x)
    ) {
      return(
        paste(
          deparse(x, width.cutoff = 500L),
          collapse = ""
        )
      )
    }

    x
  }

  canonical_extra <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }

    if (
      inherits(x, "data.table") ||
        is.data.frame(x) ||
        is.matrix(x)
    ) {
      return(canonical_tabular(x))
    }

    if (
      inherits(x, "formula") ||
        is.call(x) ||
        is.expression(x)
    ) {
      return(canonical_language(x))
    }

    if (is.environment(x)) {
      stop(
        "An environment was passed through '...'. Environments cannot be ",
        "included in a stable cross-session cache fingerprint.",
        call. = FALSE
      )
    }

    if (is.list(x)) {
      return(lapply(x, canonical_extra))
    }

    x
  }

  fingerprint_input <- list(
    wrapper_version = get_wrapper_version(),
    dependency_versions = dependency_versions,
    y = canonical_tabular(y),
    X = canonical_tabular(X),
    TR = canonical_tabular(TR),
    studyDesign = canonical_tabular(studyDesign),
    formula = canonical_language(formula),
    row.eff = canonical_language(row.eff),
    family = as.character(family),
    num.lv = as.integer(num.lv),
    sd.errors = isTRUE(sd.errors),
    method = as.character(method),
    control = canonical_extra(control),
    control.start = canonical_extra(control.start),
    gradient_threshold = as.numeric(gradient_threshold),
    coef_threshold = as.numeric(coef_threshold),
    sparse_screen = isTRUE(sparse_screen),
    sparse_action = as.character(sparse_action),
    abort_on = sort(unique(as.character(abort_on))),
    sparse_min_events_per_parameter = as.numeric(
      sparse_min_events_per_parameter
    ),
    sparse_min_total_events = as.numeric(
      sparse_min_total_events
    ),
    sparse_min_nonzero_rows = as.numeric(
      sparse_min_nonzero_rows
    ),
    sparse_event_basis = as.character(
      sparse_event_basis
    ),
    sparse_count_like_only = isTRUE(
      sparse_count_like_only
    ),
    design_check = isTRUE(design_check),
    design_action = as.character(design_action),
    design_tolerance = as.numeric(design_tolerance),
    reference_stability_check = as.character(reference_stability_check),
    reference_stability_max_levels = as.integer(reference_stability_max_levels),
    reference_logLik_tolerance = reference_logLik_tolerance,
    reference_stability_seed = reference_stability_seed,
    tmb_threads = tmb_threads,
    extra_args = canonical_extra(extra_args)
  )

  digest::digest(
    fingerprint_input,
    algo = "xxhash64",
    serialize = TRUE
  )
}
