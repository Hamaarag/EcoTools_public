#' Screen GLLVM Responses for Sparse Data
#'
#' Screens each response column in a GLLVM response matrix for low total
#' abundance, low occupancy, and too few informative observations relative to
#' the approximate number of fixed-effect parameters.
#'
#' @param y A response matrix or data frame. Columns represent responses or
#'   species and rows represent observations.
#' @param X Optional predictor matrix or data frame.
#' @param formula Optional fixed-effect formula used to estimate the number of
#'   fixed-effect parameters.
#' @param min_events_per_parameter Minimum informative observations per
#'   approximate fixed-effect parameter.
#' @param min_total_events Minimum total abundance per response.
#' @param min_nonzero_rows Minimum number of occupied or non-zero rows.
#' @param event_basis Basis used for events-per-parameter screening. One of
#'   "nonzero", "total", or "both".
#' @param count_like_only If TRUE, screen only non-negative integer-like
#'   responses.
#'
#' @return A list containing the full screening table, sparse-only rows, and
#'   alphabetically sorted sparse response names.
#'
#' @export
screen_sparse_gllvm <- function(
  y,
  X = NULL,
  formula = NULL,
  min_events_per_parameter = 5,
  min_total_events = 5,
  min_nonzero_rows = 3,
  event_basis = c("nonzero", "total", "both"),
  count_like_only = TRUE
) {
  event_basis <- match.arg(event_basis)

  if (is.null(y)) {
    stop("'y' must not be NULL.", call. = FALSE)
  }

  # Convert the response object to an ordinary matrix.
  # This avoids data.table nonstandard-evaluation and internal-attribute issues.
  y_mat <- as.matrix(
    as.data.frame(
      y,
      stringsAsFactors = FALSE
    )
  )

  if (length(dim(y_mat)) != 2L) {
    stop("'y' must be coercible to a two-dimensional matrix.", call. = FALSE)
  }

  if (ncol(y_mat) == 0L) {
    stop("'y' contains no response columns.", call. = FALSE)
  }

  if (is.null(colnames(y_mat))) {
    colnames(y_mat) <- paste0(
      "response_",
      seq_len(ncol(y_mat))
    )
  }

  # Estimate the approximate number of fixed-effect coefficients per response.
  # model.matrix() accounts for intercepts, factor contrasts, and interactions.
  parameter_count <- tryCatch(
    {
      if (!is.null(formula) && !is.null(X)) {
        X_df <- as.data.frame(
          X,
          stringsAsFactors = FALSE
        )

        formula_terms <- stats::terms(
          formula,
          data = X_df
        )

        predictor_terms <- stats::delete.response(
          formula_terms
        )

        ncol(
          stats::model.matrix(
            predictor_terms,
            data = X_df,
            na.action = stats::na.pass
          )
        )
      } else if (!is.null(X)) {
        ncol(
          as.data.frame(X)
        ) +
          1L
      } else {
        1L
      }
    },
    error = function(e) {
      warning(
        paste0(
          "Could not calculate the fixed-effect parameter count from the ",
          "formula; using the number of X columns plus an intercept. Reason: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )

      if (!is.null(X)) {
        ncol(as.data.frame(X)) + 1L
      } else {
        1L
      }
    }
  )

  parameter_count <- max(
    1L,
    as.integer(parameter_count)
  )

  result_rows <- vector(
    mode = "list",
    length = ncol(y_mat)
  )

  for (j in seq_len(ncol(y_mat))) {
    response_name <- colnames(y_mat)[j]

    yy <- suppressWarnings(
      as.numeric(y_mat[, j])
    )

    yy_valid <- yy[
      !is.na(yy)
    ]

    # A response with no observed values cannot be screened or fitted
    # meaningfully. Retain it in the report and flag it as structurally invalid.
    if (length(yy_valid) == 0L) {
      result_rows[[j]] <- data.frame(
        response = response_name,
        count_like = NA,
        total_events = NA_real_,
        nonzero_n = NA_integer_,
        n_obs = 0L,
        zero_prop = NA_real_,
        estimated_parameters = parameter_count,
        total_events_per_parameter = NA_real_,
        nonzero_per_parameter = NA_real_,
        event_basis = event_basis,
        basis_value = NA_real_,
        fails_no_observations = TRUE,
        fails_nonfinite_values = FALSE,
        fails_total_events = TRUE,
        fails_nonzero_rows = TRUE,
        fails_events_per_parameter = TRUE,
        sparse_flag = TRUE,
        stringsAsFactors = FALSE
      )

      next
    }

    contains_nonfinite_values <- any(
      !is.finite(yy_valid)
    )

    count_like <- !contains_nonfinite_values &&
      all(yy_valid >= 0) &&
      all(
        abs(
          yy_valid - round(yy_valid)
        ) <=
          1e-8
      )

    # Non-count responses are retained in the full report but are not flagged
    # when count_like_only is TRUE. Non-finite response values remain structural
    # failures rather than being classified as ordinary non-count responses.
    if (
      isTRUE(count_like_only) &&
        !isTRUE(count_like) &&
        !contains_nonfinite_values
    ) {
      result_rows[[j]] <- data.frame(
        response = response_name,
        count_like = FALSE,
        total_events = NA_real_,
        nonzero_n = NA_integer_,
        n_obs = length(yy_valid),
        zero_prop = NA_real_,
        estimated_parameters = parameter_count,
        total_events_per_parameter = NA_real_,
        nonzero_per_parameter = NA_real_,
        event_basis = event_basis,
        basis_value = NA_real_,
        fails_no_observations = FALSE,
        fails_nonfinite_values = FALSE,
        fails_total_events = FALSE,
        fails_nonzero_rows = FALSE,
        fails_events_per_parameter = FALSE,
        sparse_flag = FALSE,
        stringsAsFactors = FALSE
      )

      next
    }

    if (contains_nonfinite_values) {
      result_rows[[j]] <- data.frame(
        response = response_name,
        count_like = FALSE,
        total_events = NA_real_,
        nonzero_n = NA_integer_,
        n_obs = length(yy_valid),
        zero_prop = NA_real_,
        estimated_parameters = parameter_count,
        total_events_per_parameter = NA_real_,
        nonzero_per_parameter = NA_real_,
        event_basis = event_basis,
        basis_value = NA_real_,
        fails_no_observations = FALSE,
        fails_nonfinite_values = TRUE,
        fails_total_events = TRUE,
        fails_nonzero_rows = TRUE,
        fails_events_per_parameter = TRUE,
        sparse_flag = TRUE,
        stringsAsFactors = FALSE
      )

      next
    }

    total_events <- sum(
      yy,
      na.rm = TRUE
    )

    nonzero_n <- sum(
      yy > 0,
      na.rm = TRUE
    )

    n_obs <- sum(
      !is.na(yy)
    )

    zero_prop <- if (n_obs > 0L) {
      mean(
        yy == 0,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }

    total_events_per_parameter <- total_events /
      parameter_count

    nonzero_per_parameter <- nonzero_n /
      parameter_count

    basis_value <- switch(
      event_basis,
      nonzero = nonzero_per_parameter,
      total = total_events_per_parameter,
      both = min(
        total_events_per_parameter,
        nonzero_per_parameter,
        na.rm = TRUE
      )
    )

    if (!is.finite(basis_value)) {
      basis_value <- NA_real_
    }

    fails_total_events <- is.finite(total_events) &&
      total_events < min_total_events

    fails_nonzero_rows <- nonzero_n < min_nonzero_rows

    fails_events_per_parameter <- is.finite(basis_value) &&
      basis_value < min_events_per_parameter

    sparse_flag <- isTRUE(fails_total_events) ||
      isTRUE(fails_nonzero_rows) ||
      isTRUE(fails_events_per_parameter)

    result_rows[[j]] <- data.frame(
      response = response_name,
      count_like = count_like,
      total_events = total_events,
      nonzero_n = nonzero_n,
      n_obs = n_obs,
      zero_prop = zero_prop,
      estimated_parameters = parameter_count,
      total_events_per_parameter = total_events_per_parameter,
      nonzero_per_parameter = nonzero_per_parameter,
      event_basis = event_basis,
      basis_value = basis_value,
      fails_no_observations = FALSE,
      fails_nonfinite_values = FALSE,
      fails_total_events = fails_total_events,
      fails_nonzero_rows = fails_nonzero_rows,
      fails_events_per_parameter = fails_events_per_parameter,
      sparse_flag = sparse_flag,
      stringsAsFactors = FALSE
    )
  }

  # Combine all response-level results without using data.table NSE.
  details <- do.call(
    rbind,
    result_rows
  )

  row.names(details) <- NULL

  required_columns <- c(
    "response",
    "sparse_flag"
  )

  missing_columns <- setdiff(
    required_columns,
    names(details)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Sparse-screen result is missing required column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  # Select flagged responses with ordinary data-frame indexing.
  sparse_details <- details[
    !is.na(details$sparse_flag) &
      details$sparse_flag,
    ,
    drop = FALSE
  ]

  sparse_species <- sort(
    unique(
      as.character(
        sparse_details$response
      )
    )
  )

  sparse_species <- sparse_species[
    !is.na(sparse_species) &
      nzchar(sparse_species)
  ]

  n_sparse <- length(sparse_species)

  displayed_sparse_species <- utils::head(
    sparse_species,
    10L
  )

  warning_types <- character(0)
  warning_msg <- NULL

  if (n_sparse > 0L) {
    if (
      "fails_no_observations" %in%
        names(sparse_details) &&
        any(
          sparse_details$fails_no_observations,
          na.rm = TRUE
        )
    ) {
      warning_types <- c(
        warning_types,
        "no_observations"
      )
    }

    if (
      "fails_nonfinite_values" %in%
        names(sparse_details) &&
        any(
          sparse_details$fails_nonfinite_values,
          na.rm = TRUE
        )
    ) {
      warning_types <- c(
        warning_types,
        "nonfinite_response"
      )
    }

    displayed_text <- paste(
      displayed_sparse_species,
      collapse = ", "
    )

    if (n_sparse > 10L) {
      displayed_text <- paste0(
        displayed_text,
        ", and ",
        n_sparse - 10L,
        " additional species"
      )
    }

    if (
      any(
        sparse_details$fails_nonzero_rows,
        na.rm = TRUE
      )
    ) {
      warning_types <- c(
        warning_types,
        "occupancy"
      )
    }

    if (
      any(
        sparse_details$fails_events_per_parameter,
        na.rm = TRUE
      )
    ) {
      warning_types <- c(
        warning_types,
        "events_per_parameter"
      )
    }

    if (
      any(
        sparse_details$fails_total_events,
        na.rm = TRUE
      )
    ) {
      warning_types <- c(
        warning_types,
        "total_events"
      )
    }

    warning_types <- unique(
      warning_types
    )

    warning_msg <- paste0(
      "SPARSE GLLVM WARNING: ",
      n_sparse,
      " response column(s) have sparse count data. ",
      "Sparse species/responses: ",
      displayed_text,
      ". Thresholds: minimum total abundance/events = ",
      min_total_events,
      ", minimum occupied/non-zero rows = ",
      min_nonzero_rows,
      ", minimum informative observations per approximate fixed-effect ",
      "parameter = ",
      min_events_per_parameter,
      " using basis '",
      event_basis,
      "'. Treat affected species coefficients and species-specific inference ",
      "cautiously."
    )
  }

  return(
    list(
      pass = n_sparse == 0L,
      warning = warning_msg,
      warning_types = warning_types,
      details = details,
      sparse_details = sparse_details,
      sparse_species = sparse_species,
      displayed_sparse_species = displayed_sparse_species,
      n_sparse = n_sparse,
      estimated_parameters = parameter_count,
      event_basis = event_basis,
      thresholds = list(
        min_events_per_parameter = min_events_per_parameter,
        min_total_events = min_total_events,
        min_nonzero_rows = min_nonzero_rows
      )
    )
  )
}
