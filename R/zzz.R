.onLoad <- function(
  libname,
  pkgname
) {
  if (
    requireNamespace(
      "emmeans",
      quietly = TRUE
    )
  ) {
    emmeans::.emm_register(
      classes = "ecotools_model_bundle",
      pkgname = pkgname
    )
  }

  invisible(NULL)
}
