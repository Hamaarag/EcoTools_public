#' Normalize a Model Family Specification
#'
#' @description
#' Converts a supported character family name into an engine-specific family
#' object. Existing family objects are returned unchanged.
#'
#' @param family
#' A family object or a character scalar naming a supported family.
#'
#' @param package
#' Character scalar identifying the fitting engine.
#'
#' @return
#' A family object accepted by the selected fitting engine.
#'
#' @keywords internal
normalize_model_family <- function(
  family,
  package
) {
  if (!is.character(family)) {
    return(family)
  }

  if (
    length(family) != 1L ||
      is.na(family) ||
      !nzchar(family)
  ) {
    stop(
      paste0(
        "When supplied as text, 'family' must be one non-empty ",
        "character value."
      ),
      call. = FALSE
    )
  }

  family_name <- family

  resolve_exported_family <- function(
    namespace_name,
    function_name
  ) {
    exported_names <- getNamespaceExports(
      namespace_name
    )

    if (!function_name %in% exported_names) {
      stop(
        paste0(
          "Family function '",
          function_name,
          "' is not exported by package '",
          namespace_name,
          "' in the installed package version."
        ),
        call. = FALSE
      )
    }

    family_function <- getExportedValue(
      namespace_name,
      function_name
    )

    if (!is.function(family_function)) {
      stop(
        paste0(
          "Exported object '",
          namespace_name,
          "::",
          function_name,
          "' is not a function."
        ),
        call. = FALSE
      )
    }

    family_function()
  }

  if (identical(package, "glmmTMB")) {
    glmmtmb_families <- c(
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

    base_families <- c(
      "binomial",
      "Gamma",
      "gaussian",
      "inverse.gaussian",
      "poisson"
    )

    if (family_name %in% glmmtmb_families) {
      return(
        resolve_exported_family(
          namespace_name = "glmmTMB",
          function_name = family_name
        )
      )
    }

    if (family_name %in% base_families) {
      return(
        resolve_exported_family(
          namespace_name = "stats",
          function_name = family_name
        )
      )
    }

    stop(
      paste0(
        "Unsupported glmmTMB family name: '",
        family_name,
        "'. Supported glmmTMB-specific names include: ",
        paste(
          glmmtmb_families,
          collapse = ", "
        ),
        ". Supported base-family names include: ",
        paste(
          base_families,
          collapse = ", "
        ),
        ". Alternatively, supply an explicit family object such as ",
        "glmmTMB::nbinom2()."
      ),
      call. = FALSE
    )
  }

  if (identical(package, "mgcv")) {
    base_families <- c(
      "binomial",
      "Gamma",
      "gaussian",
      "inverse.gaussian",
      "poisson",
      "quasibinomial",
      "quasipoisson"
    )

    mgcv_families <- c(
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

    if (family_name %in% base_families) {
      return(
        resolve_exported_family(
          namespace_name = "stats",
          function_name = family_name
        )
      )
    }

    if (family_name %in% mgcv_families) {
      return(
        resolve_exported_family(
          namespace_name = "mgcv",
          function_name = family_name
        )
      )
    }

    stop(
      paste0(
        "Unsupported mgcv family name: '",
        family_name,
        "'. Supply a supported character name or an explicit family ",
        "object such as mgcv::nb()."
      ),
      call. = FALSE
    )
  }

  stop(
    paste0(
      "Character family normalization is not implemented for engine '",
      package,
      "'."
    ),
    call. = FALSE
  )
}
