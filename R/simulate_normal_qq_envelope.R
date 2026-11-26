#' Calculate a Simulated Normal Q-Q Envelope
#'
#' Generates a pointwise simulation envelope for a normal Q-Q plot by repeatedly
#' simulating standard-normal samples of the same size as the observed residual
#' vector.
#'
#' @param residuals A numeric vector of residuals.
#' @param n_sim Number of simulated normal samples.
#' @param level Confidence level for the pointwise envelope.
#' @param seed Optional integer seed used to make the displayed diagnostic
#'   envelope reproducible.
#'
#' @return A data frame containing theoretical quantiles, observed ordered
#'   residuals, and lower and upper pointwise simulation bounds.
#'
#' @keywords internal
simulate_normal_qq_envelope <- function(
  residuals,
  n_sim = 1000L,
  level = 0.95,
  seed = 123L
) {
  residuals <- as.numeric(residuals)
  residuals <- residuals[is.finite(residuals)]

  n_obs <- length(residuals)

  if (n_obs < 3L) {
    stop(
      "At least three finite residuals are required for a Q-Q envelope.",
      call. = FALSE
    )
  }

  if (
    length(n_sim) != 1L ||
      !is.finite(n_sim) ||
      n_sim < 100L
  ) {
    stop(
      "'n_sim' must be a single finite number of at least 100.",
      call. = FALSE
    )
  }

  if (
    length(level) != 1L ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1
  ) {
    stop(
      "'level' must be a single number strictly between 0 and 1.",
      call. = FALSE
    )
  }

  n_sim <- as.integer(n_sim)

  if (!is.null(seed)) {
    old_seed_exists <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )

    if (old_seed_exists) {
      old_seed <- get(
        ".Random.seed",
        envir = .GlobalEnv,
        inherits = FALSE
      )
    }

    on.exit(
      {
        if (old_seed_exists) {
          assign(
            ".Random.seed",
            old_seed,
            envir = .GlobalEnv
          )
        } else if (
          exists(
            ".Random.seed",
            envir = .GlobalEnv,
            inherits = FALSE
          )
        ) {
          rm(
            ".Random.seed",
            envir = .GlobalEnv
          )
        }
      },
      add = TRUE
    )

    set.seed(seed)
  }

  theoretical_quantiles <- stats::qnorm(
    stats::ppoints(n_obs)
  )

  observed_ordered <- sort(residuals)

  simulated_ordered <- replicate(
    n = n_sim,
    expr = sort(
      stats::rnorm(n_obs)
    )
  )

  alpha <- 1 - level

  lower_bound <- apply(
    simulated_ordered,
    MARGIN = 1L,
    FUN = stats::quantile,
    probs = alpha / 2,
    names = FALSE,
    na.rm = TRUE,
    type = 8
  )

  upper_bound <- apply(
    simulated_ordered,
    MARGIN = 1L,
    FUN = stats::quantile,
    probs = 1 - alpha / 2,
    names = FALSE,
    na.rm = TRUE,
    type = 8
  )

  data.frame(
    theoretical = theoretical_quantiles,
    observed = observed_ordered,
    lower = lower_bound,
    upper = upper_bound,
    stringsAsFactors = FALSE
  )
}
