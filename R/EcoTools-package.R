#' EcoTools: Utilities for Ecological Model Workflows
#'
#' EcoTools provides reproducible validation workflows for ecological
#' generalized linear mixed models, generalized additive models, and
#' generalized linear latent variable models. The package supports
#' fixed-effect design checks, sparse-response screening, cache validation,
#' convergence and uncertainty assessment, simulation-based residual
#' diagnostics, repeated diagnostic-stability assessment, and validated
#' downstream inference.
#'
#' @section Tweedie simulation safety:
#' Tweedie simulation may become prohibitively expensive when the fitted power
#' approaches the Gamma boundary at 2. EcoTools checks the fitted power before
#' starting DHARMa simulation or repeated residual-diagnostic stability
#' assessment. By default, these procedures are skipped when the fitted power
#' is at or above 1.99. The affected model is not considered valid for inference
#' because simulation-based validation was not completed. For strictly
#' positive responses, users should consider and explicitly evaluate a Gamma
#' model with a log link.
#'
#' @keywords internal
#' @importFrom grDevices dev.off png
#' @importFrom stats as.formula coef dnorm logLik model.frame model.response
#'   ppoints predict qnorm residuals sd
#' @importFrom utils sessionInfo
#' @importFrom rlang .data
"_PACKAGE"
