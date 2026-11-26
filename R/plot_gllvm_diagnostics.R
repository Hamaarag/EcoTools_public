# ==============================================================================
# GLLVM Diagnostic Plotting Helper
# ==============================================================================

#' Plot GLLVM Residual Diagnostics
#'
#' @description
#' Produces a set of diagnostic plots from the Dunn-Smyth residuals and linear
#' predictors of a fitted GLLVM model.
#'
#' The diagnostic output includes:
#'
#' \itemize{
#'   \item residuals against fitted linear predictors;
#'   \item a normal Q-Q plot with a simulated pointwise envelope;
#'   \item residuals against observation-row index;
#'   \item residuals against response-column index;
#'   \item a scale-location plot.
#' }
#'
#' The normal Q-Q envelope is generated from repeated standard-normal samples
#' having the same number of observations as the finite residual vector. The
#' envelope is pointwise rather than simultaneous and should therefore be used
#' as a visual diagnostic rather than a formal global hypothesis test.
#'
#' @param model
#' A fitted GLLVM model object for which Dunn-Smyth residuals and linear
#' predictors can be extracted.
#'
#' @param model_name
#' A character scalar identifying the model in diagnostic plot titles.
#'
#' @param diagnostic_check_status
#' Logical. Defaults to `TRUE` (meaning model status is checked before plotting).
#' Set to `FALSE` to force model plotting regardless of model status. This allows
#' plotting diagnostics even to invalid models.
#'
#' @param debug
#' Logical. Defaults to `FALSE` (meaning no debugging messages would appear).
#' Set to `TRUE` if you want to see these debugging messages.
#'
#' @return
#' Invisibly returns `NULL`. The function is called for its plotting side
#' effects.
#'
#' @details
#' If the simulated Q-Q envelope cannot be generated, the function falls back
#' to a conventional normal Q-Q plot and displays the envelope error in the
#' plot caption.
#'
#' Any fatal error encountered while extracting residuals or constructing the
#' diagnostic plots is converted into a warning so that diagnostic plotting
#' does not interrupt the surrounding analysis workflow.
#'
#' @seealso
#' [fit_gllvm_cached()],
#' [simulate_normal_qq_envelope()]
#'
#' @export
plot_gllvm_diagnostics <- function(
  model,
  model_name = NULL,
  diagnostic_check_status = TRUE,
  debug = FALSE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for GLLVM diagnostic plots.",
      call. = FALSE
    )
  }

  plotted_statuses <- c(
    "valid",
    "valid_fit_no_standard_errors",
    "suspicious_extreme_coefficients",
    "suspicious_inference_unreliable",
    "suspicious_optimization_incomplete",
    "suspicious_optimization_unassessed"
  )

  # Support for native GLLVM models:
  ## Extract model name from the model object if model name is null
  if (is.null(model_name)) {
    model_name <- if (!is.null(model$model_name)) {
      as.character(model$model_name)
    } else {
      "gllvm_model"
    }
  }

  ## Generate model status if there is none (as in native gllvm models):
  wrapper_model <- !is.null(model$model_status)

  model_status <- if (wrapper_model) {
    model$model_status
  } else {
    NA_character_
  }

  # Generate diagnostic captions for the graphs, according to the model status
  diagnostic_caption <- if (
    identical(
      model_status,
      "suspicious_optimization_incomplete"
    )
  ) {
    paste0(
      "Exploratory diagnostic only: the fitted model failed the configured ",
      "gradient-convergence criterion."
    )
  } else if (
    identical(
      model_status,
      "suspicious_optimization_unassessed"
    )
  ) {
    paste0(
      "Exploratory diagnostic only: gradient convergence was not assessed."
    )
  } else {
    NULL
  }

  ## Generate gradient validity if there is none (as in native gllvm models):
  gradient_valid <- if (!is.null(model$gradient_valid)) {
    model$gradient_valid
  } else {
    NA
  }

  ## Generate sd errors if there are none (as in native gllvm models):
  sd_errors_requested <- if (!is.null(model$sd_errors_requested)) {
    model$sd_errors_requested
  } else {
    NA
  }

  # Inform the user why the diagnostic plots were omitted:
  if (
    isTRUE(diagnostic_check_status) &&
      wrapper_model
  ) {
    if (!model_status %in% plotted_statuses) {
      warning(
        stringr::str_glue(
          "Diagnostic plots were intentionally omitted because ",
          "model status = '{model_status}'. ",
          "Use diagnostic_check_status = FALSE to force plotting."
        ),
        call. = FALSE
      )

      return(
        invisible(NULL)
      )
    }
  }

  diagnostic_result <- NULL

  tryCatch(
    {
      # -----------------------------------------------------------------------
      # 1. Extract residual and fitted-value matrices
      # -----------------------------------------------------------------------

      # Debugging messages
      if (isTRUE(debug)) {
        message("Attempting residual extraction...")
      }

      res_obj <- tryCatch(
        stats::residuals(model),
        error = function(e) {
          offset_info <- tryCatch(
            {
              if (is.null(model$offset)) {
                "NULL"
              } else {
                paste0(
                  "class = ",
                  paste(class(model$offset), collapse = ", "),
                  "; length = ",
                  length(model$offset),
                  "; dim = ",
                  if (is.null(dim(model$offset))) {
                    "NULL"
                  } else {
                    paste(dim(model$offset), collapse = " x ")
                  }
                )
              }
            },
            error = function(x) {
              "unable to inspect offset"
            }
          )

          stop(
            paste0(
              "Failed to extract GLLVM residuals.\n",
              "Underlying error: ",
              conditionMessage(e),
              "\n",
              "Offset information: ",
              offset_info
            ),
            call. = FALSE
          )
        }
      )

      if (isTRUE(debug)) {
        message("Diagnostic residual object structure:")
        utils::str(res_obj, max.level = 1)

        message(
          "Residual object class: ",
          paste(class(res_obj), collapse = ", ")
        )

        if (is.list(res_obj)) {
          message(
            "Residual object names: ",
            paste(names(res_obj), collapse = ", ")
          )
        }
      }

      if (!is.list(res_obj)) {
        stop(
          paste0(
            "residuals(model) returned class ",
            paste(
              class(res_obj),
              collapse = ", "
            ),
            " instead of the expected list containing a residual matrix."
          ),
          call. = FALSE
        )
      }

      if (is.null(res_obj$residuals)) {
        stop(
          "The GLLVM residual method did not return a 'residuals' component.",
          call. = FALSE
        )
      }

      res_mat <- as.matrix(
        res_obj$residuals
      )

      if (
        length(
          dim(res_mat)
        ) !=
          2L
      ) {
        stop(
          paste0(
            "The extracted GLLVM residual object is not a two-dimensional ",
            "matrix."
          ),
          call. = FALSE
        )
      }

      response_family_all <- if (!is.null(model$response_family)) {
        as.character(
          model$response_family
        )
      } else {
        NULL
      }

      if (
        !is.null(response_family_all) &&
          ncol(res_mat) < length(response_family_all)
      ) {
        warning(
          paste0(
            "GLLVM returned residuals for ",
            ncol(res_mat),
            " of ",
            length(response_family_all),
            " response columns. Diagnostic plots will include only the ",
            "responses present in the residual matrix."
          ),
          call. = FALSE
        )
      }

      response_family <- NULL
      response_names_available <-
        !is.null(names(model$response_family)) &&
        !is.null(colnames(res_mat))

      if (
        !is.null(response_family_all) &&
          response_names_available
      ) {
        response_index <- match(
          colnames(res_mat),
          names(model$response_family)
        )

        if (all(!is.na(response_index))) {
          response_family <- response_family_all[response_index]
        }
      }

      if (
        is.null(response_family) &&
          !is.null(response_family_all) &&
          !response_names_available &&
          length(response_family_all) == ncol(res_mat)
      ) {
        response_family <- response_family_all
      }

      # Use internal linear predictors when their response columns can be
      # aligned to the residual columns. This is commonly true for
      # single-component families, but not always for hurdle families such as
      # betaH.
      internal_linpred <- if (!is.null(res_obj$linpred)) {
        tryCatch(
          as.matrix(
            res_obj$linpred
          ),
          error = function(e) {
            NULL
          }
        )
      } else {
        NULL
      }

      internal_linpred_aligned <- FALSE

      if (
        !is.null(internal_linpred) &&
          nrow(internal_linpred) == nrow(res_mat)
      ) {
        if (
          !is.null(colnames(internal_linpred)) &&
            !is.null(colnames(res_mat))
        ) {
          response_index <- match(
            colnames(res_mat),
            colnames(internal_linpred)
          )

          if (all(!is.na(response_index))) {
            internal_linpred <- internal_linpred[,
              response_index,
              drop = FALSE
            ]
            internal_linpred_aligned <- TRUE
          }
        }

        if (
          !internal_linpred_aligned &&
            identical(dim(internal_linpred), dim(res_mat)) &&
            (is.null(colnames(internal_linpred)) ||
              is.null(colnames(res_mat)))
        ) {
          internal_linpred_aligned <- TRUE
        }
      }

      response_predictions <- NULL
      response_prediction_error <- NULL

      if (!isTRUE(internal_linpred_aligned)) {
        response_predictions <- tryCatch(
          {
            prediction_result <- if (!is.null(response_family_all)) {
              get_gllvm_response_predictions(
                model = model,
                family = response_family_all
              )
            } else {
              stats::predict(
                model,
                type = "response"
              )
            }

            prediction_matrix <- as.matrix(
              prediction_result
            )

            if (
              nrow(prediction_matrix) == nrow(res_mat) &&
                !is.null(colnames(prediction_matrix)) &&
                !is.null(colnames(res_mat))
            ) {
              response_index <- match(
                colnames(res_mat),
                colnames(prediction_matrix)
              )

              if (all(!is.na(response_index))) {
                prediction_matrix <- prediction_matrix[,
                  response_index,
                  drop = FALSE
                ]
              } else {
                stop(
                  "Response-scale prediction columns could not be matched to residual columns by name.",
                  call. = FALSE
                )
              }
            }

            if (
              !identical(
                dim(prediction_matrix),
                dim(res_mat)
              )
            ) {
              stop(
                paste0(
                  "response-scale predictions have dimensions ",
                  paste(
                    dim(prediction_matrix),
                    collapse = " x "
                  ),
                  ", but residuals have dimensions ",
                  paste(
                    dim(res_mat),
                    collapse = " x "
                  ),
                  "."
                ),
                call. = FALSE
              )
            }

            prediction_matrix
          },
          error = function(e) {
            response_prediction_error <<-
              conditionMessage(e)

            NULL
          }
        )
      }

      if (isTRUE(internal_linpred_aligned)) {
        fit_mat <- internal_linpred
        fitted_scale <- "link"
        fitted_axis_label <- "Linear predictor"
      } else if (!is.null(response_predictions)) {
        fit_mat <- response_predictions
        fitted_scale <- "response"
        fitted_axis_label <- "Response-scale fitted value"
      } else {
        stop(
          paste0(
            "GLLVM residuals could not be paired with fitted values. ",
            "Internal linear-predictor dimensions were ",
            if (is.null(internal_linpred)) {
              "unavailable"
            } else {
              paste(
                dim(internal_linpred),
                collapse = " x "
              )
            },
            ", residual dimensions were ",
            paste(
              dim(res_mat),
              collapse = " x "
            ),
            ", and response-scale prediction failed",
            if (
              !is.null(response_prediction_error) &&
                nzchar(response_prediction_error)
            ) {
              paste0(
                ": ",
                response_prediction_error
              )
            } else {
              "."
            }
          ),
          call. = FALSE
        )
      }

      if (isTRUE(debug)) {
        message(
          "Residual dimensions: ",
          paste(
            dim(res_mat),
            collapse = " x "
          )
        )

        message(
          "Internal linear-predictor dimensions: ",
          if (is.null(internal_linpred)) {
            "unavailable"
          } else {
            paste(
              dim(internal_linpred),
              collapse = " x "
            )
          }
        )

        message(
          "Diagnostic fitted-value source: ",
          fitted_scale
        )

        message(
          "Diagnostic fitted-value dimensions: ",
          paste(
            dim(fit_mat),
            collapse = " x "
          )
        )
      }

      # -----------------------------------------------------------------------
      # 2. Flatten matrices into a plotting data frame
      # -----------------------------------------------------------------------
      response_names <- colnames(
        res_mat
      )

      if (
        is.null(response_names) ||
          length(response_names) != ncol(res_mat)
      ) {
        response_names <- paste0(
          "Response_",
          seq_len(
            ncol(res_mat)
          )
        )
      }

      response_family_for_plot <- if (!is.null(response_family)) {
        response_family
      } else {
        rep(
          "unspecified",
          ncol(res_mat)
        )
      }

      plot_data <- data.frame(
        Fitted = as.vector(
          fit_mat
        ),
        Residuals = as.vector(
          res_mat
        ),
        RowIndex = as.vector(
          row(res_mat)
        ),
        ColIndex = as.vector(
          col(res_mat)
        ),
        Response = rep(
          response_names,
          each = nrow(res_mat)
        ),
        Family = rep(
          response_family_for_plot,
          each = nrow(res_mat)
        ),
        FittedScale = fitted_scale,
        stringsAsFactors = FALSE
      )

      plot_data <- plot_data[
        is.finite(plot_data[["Fitted"]]) &
          is.finite(plot_data[["Residuals"]]),
        ,
        drop = FALSE
      ]

      if (nrow(plot_data) < 3L) {
        stop(
          paste0(
            "Fewer than three finite residual and fitted-value pairs ",
            "were available for diagnostic plotting."
          ),
          call. = FALSE
        )
      }

      plot_data[["SqrtAbsRes"]] <- sqrt(
        abs(plot_data[["Residuals"]])
      )

      # -----------------------------------------------------------------------
      # 3. Calculate a simulated pointwise normal Q-Q envelope
      # -----------------------------------------------------------------------
      envelope_error_message <- NULL

      qq_envelope_data <- tryCatch(
        simulate_normal_qq_envelope(
          residuals = plot_data[["Residuals"]],
          n_sim = 1000L,
          level = 0.95,
          seed = 123L
        ),
        error = function(e) {
          envelope_error_message <<- conditionMessage(e)
          NULL
        }
      )

      envelope_available <-
        !is.null(qq_envelope_data) &&
        is.data.frame(qq_envelope_data) &&
        nrow(qq_envelope_data) > 0L

      # -----------------------------------------------------------------------
      # 4. Residuals against fitted linear predictors
      # -----------------------------------------------------------------------
      p1 <- ggplot2::ggplot(
        data = plot_data,
        mapping = ggplot2::aes(
          x = .data$Fitted,
          y = .data$Residuals
        )
      ) +
        ggplot2::geom_point(
          alpha = 0.35,
          color = "#0072B2"
        ) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dashed",
          color = "#D55E00",
          linewidth = 0.8
        ) +
        ggplot2::geom_smooth(
          method = "loess",
          formula = y ~ x,
          se = FALSE,
          color = "black"
        ) +
        ggplot2::facet_wrap(
          ggplot2::vars(
            Family
          ),
          scales = "free_x"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste(
            model_name,
            "Residuals vs Fitted"
          ),
          subtitle = if (identical(fitted_scale, "response")) {
            paste0(
              "Response-scale fitted values are shown separately by family."
            )
          } else {
            paste0(
              "Internal linear predictors are shown separately by family."
            )
          },
          x = fitted_axis_label,
          y = "Dunn-Smyth residual",
          caption = diagnostic_caption
        )

      # -----------------------------------------------------------------------
      # 5. Normal Q-Q plot
      # -----------------------------------------------------------------------
      if (envelope_available) {
        p2 <- ggplot2::ggplot(
          data = qq_envelope_data,
          mapping = ggplot2::aes(
            x = .data$theoretical,
            y = .data$observed
          )
        ) +
          ggplot2::geom_ribbon(
            mapping = ggplot2::aes(
              x = .data$theoretical,
              ymin = .data$lower,
              ymax = .data$upper
            ),
            inherit.aes = FALSE,
            alpha = 0.20,
            fill = "#0072B2",
            color = NA
          ) +
          ggplot2::geom_abline(
            intercept = 0,
            slope = 1,
            color = "black",
            linewidth = 0.8
          ) +
          ggplot2::geom_point(
            alpha = 0.50,
            color = "#0072B2"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::labs(
            title = paste(
              model_name,
              "Normal Q-Q"
            ),
            x = "Theoretical Normal Quantiles",
            y = "Ordered Dunn-Smyth Residuals",
            caption = paste(
              paste0(
                "Shaded region: 95% pointwise envelope from 1,000 ",
                "standard-normal simulations."
              ),
              diagnostic_caption,
              sep = "\n"
            )
          )
      } else {
        p2 <- ggplot2::ggplot(
          data = plot_data,
          mapping = ggplot2::aes(
            sample = .data$Residuals
          )
        ) +
          ggplot2::stat_qq(
            alpha = 0.50,
            color = "#0072B2"
          ) +
          ggplot2::stat_qq_line(
            color = "black"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::labs(
            title = paste(
              model_name,
              "Normal Q-Q"
            ),
            x = "Theoretical Normal Quantiles",
            y = "Ordered Dunn-Smyth Residuals",
            caption = paste0(
              "WARNING: Simulated Q-Q envelope omitted. ",
              if (
                !is.null(envelope_error_message) &&
                  nzchar(envelope_error_message)
              ) {
                envelope_error_message
              } else {
                "The envelope could not be calculated."
              },
              diagnostic_caption
            )
          )
      }

      # -----------------------------------------------------------------------
      # 6. Residuals against observation-row index
      # -----------------------------------------------------------------------
      p3 <- ggplot2::ggplot(
        data = plot_data,
        mapping = ggplot2::aes(
          x = .data$RowIndex,
          y = .data$Residuals
        )
      ) +
        ggplot2::geom_point(
          alpha = 0.5,
          color = "#009E73"
        ) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dashed",
          color = "#D55E00",
          linewidth = 1
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste(
            model_name,
            "Residuals vs Row"
          ),
          x = "Site Index (Row)",
          y = "Dunn-Smyth Residuals",
          caption = diagnostic_caption
        )

      # -----------------------------------------------------------------------
      # 7. Residuals against response-column index
      # -----------------------------------------------------------------------
      p4 <- ggplot2::ggplot(
        data = plot_data,
        mapping = ggplot2::aes(
          x = .data$ColIndex,
          y = .data$Residuals,
          color = .data$Family
        )
      ) +
        ggplot2::geom_jitter(
          width = 0.12,
          height = 0,
          alpha = 0.40
        ) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dashed",
          color = "#D55E00",
          linewidth = 0.8
        ) +
        ggplot2::scale_x_continuous(
          breaks = sort(
            unique(
              plot_data[["ColIndex"]]
            )
          )
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "bottom"
        ) +
        ggplot2::labs(
          title = paste(
            model_name,
            "Residuals vs Response Column"
          ),
          x = "Species/response index",
          y = "Dunn-Smyth residual",
          color = "Response family",
          caption = diagnostic_caption
        )

      # -----------------------------------------------------------------------
      # 8. Scale-location plot
      # -----------------------------------------------------------------------
      p5 <- ggplot2::ggplot(
        data = plot_data,
        mapping = ggplot2::aes(
          x = .data$Fitted,
          y = .data$SqrtAbsRes
        )
      ) +
        ggplot2::geom_point(
          alpha = 0.35,
          color = "#0072B2"
        ) +
        ggplot2::geom_smooth(
          method = "loess",
          formula = y ~ x,
          se = FALSE,
          color = "black"
        ) +
        ggplot2::facet_wrap(
          ggplot2::vars(
            Family
          ),
          scales = "free_x"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste(
            model_name,
            "Scale-Location"
          ),
          subtitle = paste0(
            "Panels show response families separately."
          ),
          x = fitted_axis_label,
          y = "\u221A|Dunn-Smyth residual|",
          caption = diagnostic_caption
        )
      # -----------------------------------------------------------------------
      # 9. Print diagnostic plots
      # -----------------------------------------------------------------------
      print(p1)
      print(p2)
      print(p3)
      print(p4)
      print(p5)

      diagnostic_result <<- list(
        completed = TRUE,
        residual_dimensions = dim(
          res_mat
        ),
        fitted_dimensions = dim(
          fit_mat
        ),
        fitted_scale = fitted_scale,
        fitted_axis_label = fitted_axis_label,
        response_families = sort(
          unique(
            response_family
          )
        ),
        n_finite_pairs = nrow(
          plot_data
        )
      )
    },
    error = function(e) {
      diagnostic_result <<- list(
        completed = FALSE,
        error = conditionMessage(e),
        model_status = model_status,
        gradient_valid = gradient_valid
      )
      warning(
        paste0(
          "Could not generate diagnostic plots for '",
          model_name,
          "'.\n\n",
          "Model status: ",
          if (!is.null(model_status)) {
            model_status
          } else {
            "unknown"
          },
          "\n",
          "Gradient valid: ",
          if (!is.null(gradient_valid)) {
            gradient_valid
          } else {
            NA
          },
          "\n",
          "Standard errors requested: ",
          if (!is.null(sd_errors_requested)) {
            sd_errors_requested
          } else {
            NA
          },
          "\n",
          "Random effects: ",
          if (!is.null(model$row.eff)) {
            deparse(model$row.eff)
          } else {
            "unknown"
          },
          "\n",
          "Offset present: ",
          !is.null(model$offset),
          "\n\n",
          "Underlying error: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  invisible(
    diagnostic_result
  )
}
