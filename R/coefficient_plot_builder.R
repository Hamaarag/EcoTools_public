# ==============================================================================
# Coefficient plot builder
# ==============================================================================

#' Build a Coefficient Plot
#'
#' Creates a coefficient plot from a table containing model estimates and
#' standard errors. The function uses an ordinary data.frame internally to avoid
#' data.table nonstandard-evaluation problems when called from a package
#' namespace.
#'
#' @param data_table A data.frame, data.table, or compatible tabular object.
#' @param reverse_heb Logical. Reverse Hebrew character strings for graphics
#'   devices that do not handle right-to-left text correctly.
#' @param title_txt Plot title.
#' @param x_txt X-axis title.
#' @param y_txt Y-axis title.
#' @param color_txt Color-legend title.
#' @param font_family Font family.
#' @param base_size Base text size.
#' @param color_palette Optional named or unnamed color vector.
#' @param significant_only Logical. If TRUE, plot only rows for which
#'   is_significant is TRUE. If FALSE, plot all rows and distinguish significance
#'   using alpha.
#' @return A ggplot object.
#' @export
build_coef_plot <- function(
  data_table,
  reverse_heb = FALSE,
  title_txt = "Don't forget to insert a title!",
  x_txt = "Forgot to name the X axis?",
  y_txt = "Remember to name the y axis too!",
  color_txt = "Color had to be named as well, you know",
  font_family = "sans",
  base_size = 12,
  color_palette = NULL,
  significant_only = TRUE
) {
  # ---------------------------------------------------------------------------
  # 1. Convert input to a plain data frame
  # ---------------------------------------------------------------------------
  # ggplot2 does not require a data.table. Using a plain data.frame here avoids
  # ':=' and data.table NSE problems inside package namespaces.
  plot_data <- as.data.frame(
    data.table::copy(data_table),
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # 2. Validate required columns
  # ---------------------------------------------------------------------------
  required_cols <- c(
    "is_significant",
    "heb_name",
    "pred_heb",
    "estimate",
    "std.error"
  )

  missing_cols <- setdiff(required_cols, names(plot_data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in coefficient plot data: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Ensure estimate and standard-error columns are numeric.
  if (!is.numeric(plot_data$estimate)) {
    stop("'estimate' must be numeric.", call. = FALSE)
  }

  if (!is.numeric(plot_data$std.error)) {
    stop("'std.error' must be numeric.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 3. Create a stable significance indicator
  # ---------------------------------------------------------------------------
  # Using $ assignment rather than data.table ':=' is deliberate.
  plot_data$significance_alpha <- factor(
    !is.na(plot_data$is_significant) &
      plot_data$is_significant %in% TRUE,
    levels = c(FALSE, TRUE)
  )

  # Optionally retain only statistically significant coefficients.
  if (isTRUE(significant_only)) {
    plot_data <- plot_data[
      !is.na(plot_data$is_significant) &
        plot_data$is_significant %in% TRUE,
      ,
      drop = FALSE
    ]
  }

  if (nrow(plot_data) == 0) {
    stop(
      "No coefficients remain to plot. ",
      if (isTRUE(significant_only)) {
        "No rows had 'is_significant == TRUE'."
      } else {
        "The input table contains no usable rows."
      },
      call. = FALSE
    )
  }

  # Exclude rows with unusable coefficient values.
  unusable_rows <- !is.finite(plot_data$estimate) |
    !is.finite(plot_data$std.error)

  if (any(unusable_rows)) {
    warning(
      sum(unusable_rows),
      " row(s) with non-finite estimates or standard errors were removed.",
      call. = FALSE
    )

    plot_data <- plot_data[!unusable_rows, , drop = FALSE]
  }

  if (nrow(plot_data) == 0) {
    stop(
      "No rows with finite estimates and standard errors remain.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # 4. Optionally reverse Hebrew strings
  # ---------------------------------------------------------------------------
  if (isTRUE(reverse_heb)) {
    plot_data$heb_name <- strReverse(plot_data$heb_name)
    plot_data$pred_heb <- strReverse(plot_data$pred_heb)

    title_txt <- strReverse(title_txt)
    x_txt <- strReverse(x_txt)
    y_txt <- strReverse(y_txt)
    color_txt <- strReverse(color_txt)
  }

  # ---------------------------------------------------------------------------
  # 5. Resolve optional global plotting defaults
  # ---------------------------------------------------------------------------
  palette_use <- if (!is.null(color_palette)) {
    color_palette
  } else if (exists("okabe", inherits = TRUE)) {
    get("okabe", inherits = TRUE)
  } else {
    NULL
  }

  font_use <- if (!identical(font_family, "sans")) {
    font_family
  } else if (exists("fontname", inherits = TRUE)) {
    get("fontname", inherits = TRUE)
  } else {
    "sans"
  }

  size_use <- if (!identical(base_size, 12)) {
    base_size
  } else if (exists("fontsize", inherits = TRUE)) {
    get("fontsize", inherits = TRUE)
  } else {
    12
  }

  # Guard against unusable inherited plotting values.
  if (length(font_use) != 1L || is.na(font_use)) {
    font_use <- "sans"
  }

  if (length(size_use) != 1L || !is.finite(size_use)) {
    size_use <- 12
  }

  # ---------------------------------------------------------------------------
  # 6. Construct coefficient plot
  # ---------------------------------------------------------------------------
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      y = forcats::fct_reorder(heb_name, estimate),
      x = estimate,
      color = pred_heb
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = estimate - 1.96 * std.error,
        xmax = estimate + 1.96 * std.error,
        alpha = significance_alpha
      ),
      orientation = "y",
      height = 0.25
    ) +
    ggplot2::geom_point(
      ggplot2::aes(alpha = significance_alpha)
    ) +
    ggplot2::scale_alpha_manual(
      values = c(
        "FALSE" = 0.30,
        "TRUE" = 0.85
      ),
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::theme_minimal()

  # Apply the optional color palette.
  if (!is.null(palette_use)) {
    p <- p +
      ggplot2::scale_color_manual(
        values = palette_use
      )
  }

  p <- p +
    ggplot2::labs(
      title = title_txt,
      x = x_txt,
      y = y_txt,
      color = color_txt
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      text = ggplot2::element_text(
        family = font_use,
        size = size_use
      ),
      axis.text = ggplot2::element_text(
        size = size_use - 1
      ),
      panel.grid.major.y = ggplot2::element_line(
        color = "grey90",
        linetype = 5,
        linewidth = 1
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )

  return(p)
}
