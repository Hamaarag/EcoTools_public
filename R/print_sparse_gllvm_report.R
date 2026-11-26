# ==============================================================================
# Sparse GLLVM report helper
# ==============================================================================

#' Print Sparse-Species Screening Results
#'
#' @description
#' Prints the species or response columns flagged by
#' [screen_sparse_gllvm()], together with the numerical quantities used by the
#' sparse-data screen.
#'
#' The function first uses the precomputed `sparse_details` table returned by
#' [screen_sparse_gllvm()]. For compatibility with older cached objects, it can
#' reconstruct the flagged subset from the complete `details` table when
#' `sparse_details` is unavailable.
#'
#' @param sparse_screen_obj
#' An object returned by [screen_sparse_gllvm()].
#'
#' @param prefix
#' Optional character scalar used as the report heading.
#'
#' @return
#' Invisibly returns a data frame containing the flagged responses. Returns
#' `NULL` when no sparse-screen object or detailed results are available.
#'
#' @export
print_sparse_gllvm_report <- function(
  sparse_screen_obj,
  prefix = "SPARSE GLLVM SCREEN"
) {
  # ---------------------------------------------------------------------------
  # 1. Validate the supplied screening object
  # ---------------------------------------------------------------------------
  if (
    is.null(sparse_screen_obj) ||
      !is.list(sparse_screen_obj)
  ) {
    return(invisible(NULL))
  }

  # A pass value of NA indicates that the screen itself could not be completed.
  if (
    !is.null(sparse_screen_obj$pass) &&
      length(sparse_screen_obj$pass) == 1L &&
      is.na(sparse_screen_obj$pass)
  ) {
    screen_error_text <- if (
      !is.null(sparse_screen_obj$warning) &&
        length(sparse_screen_obj$warning) > 0L
    ) {
      as.character(sparse_screen_obj$warning[1])
    } else {
      "The sparse-data screen could not be completed."
    }

    warning(
      paste0(
        prefix,
        ": SCREEN FAILED / NOT ASSESSED\n",
        screen_error_text
      ),
      call. = FALSE
    )

    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # 2. Obtain the sparse-only details
  # ---------------------------------------------------------------------------
  # Prefer the sparse-only table created by the current screener.
  flagged_details <- if (!is.null(sparse_screen_obj$sparse_details)) {
    as.data.frame(
      sparse_screen_obj$sparse_details,
      stringsAsFactors = FALSE
    )
  } else if (!is.null(sparse_screen_obj$details)) {
    # Compatibility route for older cached objects that contain only the full
    # details table.
    all_details <- as.data.frame(
      sparse_screen_obj$details,
      stringsAsFactors = FALSE
    )

    if (
      nrow(all_details) == 0L ||
        !"sparse_flag" %in% names(all_details)
    ) {
      return(invisible(NULL))
    }

    sparse_rows <- !is.na(all_details[["sparse_flag"]]) &
      all_details[["sparse_flag"]]

    all_details[
      sparse_rows,
      ,
      drop = FALSE
    ]
  } else {
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # 3. Handle a successful screen with no flagged responses
  # ---------------------------------------------------------------------------
  if (nrow(flagged_details) == 0L) {
    message(
      prefix,
      ": no response columns were flagged as sparse."
    )

    return(invisible(flagged_details))
  }

  # ---------------------------------------------------------------------------
  # 4. Print the report heading
  # ---------------------------------------------------------------------------
  cat(
    "\n",
    prefix,
    ": FLAGGED RESPONSE COLUMNS\n",
    sep = ""
  )

  cat(
    "The following ",
    nrow(flagged_details),
    " species/response column(s) triggered at least one ",
    "sparse-data criterion:\n\n",
    sep = ""
  )

  # ---------------------------------------------------------------------------
  # 5. Arrange useful reporting columns
  # ---------------------------------------------------------------------------
  preferred_columns <- c(
    "response",
    "total_events",
    "nonzero_n",
    "n_obs",
    "zero_prop",
    "estimated_parameters",
    "total_events_per_parameter",
    "nonzero_per_parameter",
    "event_basis",
    "basis_value",
    "fails_total_events",
    "fails_nonzero_rows",
    "fails_events_per_parameter",
    "sparse_reason",
    "sparse_flag"
  )

  columns_to_print <- preferred_columns[
    preferred_columns %in% names(flagged_details)
  ]

  # Retain additional columns after the preferred columns so future diagnostic
  # fields returned by the screener are not silently omitted.
  additional_columns <- setdiff(
    names(flagged_details),
    columns_to_print
  )

  columns_to_print <- c(
    columns_to_print,
    additional_columns
  )

  print(
    flagged_details[,
      columns_to_print,
      drop = FALSE
    ],
    row.names = FALSE
  )

  # ---------------------------------------------------------------------------
  # 6. Print a compact species-name list
  # ---------------------------------------------------------------------------
  if ("response" %in% names(flagged_details)) {
    flagged_names <- sort(
      unique(
        as.character(
          flagged_details[["response"]]
        )
      )
    )

    flagged_names <- flagged_names[
      !is.na(flagged_names) &
        nzchar(flagged_names)
    ]

    if (length(flagged_names) > 0L) {
      cat("\nFlagged species/response names:\n")
      cat(
        paste0("  - ", flagged_names),
        sep = "\n"
      )
      cat("\n")
    }
  }

  invisible(flagged_details)
}
