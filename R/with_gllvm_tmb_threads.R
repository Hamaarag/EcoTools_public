#' Evaluate Code with a Temporary GLLVM TMB Thread Setting
#'
#' Configures the OpenMP thread count for the TMB dynamic library used by
#' gllvm, evaluates an expression, and restores the caller's previous thread
#' setting afterward.
#'
#' This helper controls threads within one gllvm fit. It does not create
#' multiple R worker processes and does not parallelize independent starting
#' values.
#'
#' @param threads Positive integer number of TMB/OpenMP threads.
#' @param code Expression to evaluate.
#'
#' @return The value returned by `code`.
#'
#' @keywords internal
with_gllvm_tmb_threads <- function(
  threads,
  code
) {
  if (
    !is.numeric(threads) ||
      length(threads) != 1L ||
      is.na(threads) ||
      !is.finite(threads) ||
      threads < 1L ||
      threads != as.integer(threads)
  ) {
    stop(
      "'threads' must be one positive integer.",
      call. = FALSE
    )
  }

  threads <- as.integer(threads)

  if (!requireNamespace("gllvm", quietly = TRUE)) {
    stop(
      "Package 'gllvm' is required.",
      call. = FALSE
    )
  }

  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop(
      "Package 'TMB' is required to configure GLLVM fitting threads.",
      call. = FALSE
    )
  }

  # Loading the gllvm namespace ensures that its compiled TMB DLL is available
  # before TMB::openmp() attempts to inspect or change its thread setting.
  loadNamespace("gllvm")

  previous_threads <- tryCatch(
    {
      as.integer(
        TMB::openmp(
          DLL = "gllvm"
        )
      )
    },
    error = function(e) {
      NA_integer_
    }
  )

  if (
    length(previous_threads) != 1L ||
      !is.finite(previous_threads) ||
      previous_threads < 1L
  ) {
    previous_threads <- 1L
  }

  on.exit(
    {
      try(
        TMB::openmp(
          previous_threads,
          DLL = "gllvm"
        ),
        silent = TRUE
      )
    },
    add = TRUE
  )

  configured_threads <- tryCatch(
    {
      TMB::openmp(
        threads,
        DLL = "gllvm"
      )

      as.integer(
        TMB::openmp(
          DLL = "gllvm"
        )
      )
    },
    error = function(e) {
      stop(
        paste0(
          "Could not configure TMB/OpenMP threads for the gllvm DLL: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  if (
    length(configured_threads) != 1L ||
      !is.finite(configured_threads) ||
      configured_threads != threads
  ) {
    warning(
      paste0(
        "Requested ",
        threads,
        " TMB/OpenMP thread(s) for gllvm, but the DLL reported ",
        if (
          length(configured_threads) == 1L &&
            is.finite(configured_threads)
        ) {
          configured_threads
        } else {
          "an unavailable thread count"
        },
        ". The installed compiler or package build may not support the ",
        "requested OpenMP configuration."
      ),
      call. = FALSE
    )
  }

  eval.parent(
    substitute(code)
  )
}
