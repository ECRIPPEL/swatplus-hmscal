# ==============================================================================
# plots.R — Publication-ready scientific plots
# Morris screening, hydrographs, scatter, boxplot
#
# All plots use Times New Roman (BASE_FAMILY). On Windows the font is registered
# via windowsFonts(); on Linux/macOS ggplot resolves "TimesNR" against the
# system "Times New Roman" face if installed (otherwise falls back gracefully).
#
# KNOWN NAMING INCONSISTENCY (flagged 2026-05-19) in plot_morris():
# The cutoff variables `mu_med` and `sigma_med` are computed with mean(), not
# median(). The "_med" suffix is misleading but the behaviour is correct: the
# mean is the intended cutoff for the dashed reference lines on the Morris
# mu* vs sigma plot. Rename to `mu_mean`/`sigma_mean` in a future cleanup
# pass, coordinated with R/morris_classify.R, which uses the same names.
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

# Register Times New Roman as "TimesNR" on Windows.
if (.Platform$OS.type == "windows") {
  windowsFonts(TimesNR = windowsFont("Times New Roman"))
} else if (requireNamespace("sysfonts", quietly = TRUE) &&
           requireNamespace("showtext", quietly = TRUE)) {
  # Optional non-Windows path: register Times New Roman if available.
  try(sysfonts::font_add("TimesNR", regular = "Times New Roman.ttf"), silent = TRUE)
  try(showtext::showtext_auto(), silent = TRUE)
}

BASE_FAMILY <- "TimesNR"

# Centralised theme — every plot inherits this so typography is consistent.
.theme_pub <- function(base_size = 12) {
  theme_bw(base_family = BASE_FAMILY, base_size = base_size) +
    theme(
      text             = element_text(family = BASE_FAMILY, colour = "black"),
      plot.title       = element_text(family = BASE_FAMILY, face = "bold",
                                       size = base_size + 2, hjust = 0.5,
                                       margin = margin(b = 10)),
      plot.subtitle    = element_text(family = BASE_FAMILY, size = base_size,
                                       hjust = 0.5, colour = "grey25",
                                       margin = margin(b = 8)),
      axis.title       = element_text(family = BASE_FAMILY, size = base_size + 1,
                                       colour = "black"),
      axis.text        = element_text(family = BASE_FAMILY, size = base_size,
                                       colour = "black"),
      legend.title     = element_text(family = BASE_FAMILY, size = base_size,
                                       colour = "black"),
      legend.text      = element_text(family = BASE_FAMILY, size = base_size - 1,
                                       colour = "black"),
      legend.position  = "top",
      legend.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.30),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.55),
      axis.ticks       = element_line(colour = "black", linewidth = 0.40),
      axis.ticks.length = unit(2.5, "pt"),
      plot.margin      = margin(8, 12, 8, 8)
    )
}

# ==============================================================================
# MORRIS SCREENING PLOT
# ==============================================================================

#' Morris screening plot (mu* vs sigma)
#'
#' @param ranking Classified tibble (with grupo column, from classify_morris)
#' @param iter_label Screening label (e.g., "iteration1")
#' @return ggplot object
plot_morris <- function(ranking, iter_label = "") {
  mu_med    <- mean(ranking$mu_star, na.rm = TRUE)
  sigma_med <- mean(ranking$sigma,   na.rm = TRUE)

  ranking <- ranking %>%
    mutate(plot_label = if_else(mu_star > 0 | sigma > 0, parameter, ""))

  ggplot(ranking, aes(x = mu_star, y = sigma)) +
    geom_hline(yintercept = sigma_med, linetype = "dashed",
               colour = "grey55", linewidth = 0.35) +
    geom_vline(xintercept = mu_med, linetype = "dashed",
               colour = "grey55", linewidth = 0.35) +
    geom_point(
      aes(fill = grupo, shape = grupo),
      colour = "black", size = 3.6, stroke = 0.45, alpha = 0.85
    ) +
    geom_text_repel(
      aes(label = plot_label),
      family = BASE_FAMILY, size = 3.6, colour = "black",
      box.padding = 0.45, point.padding = 0.30,
      segment.colour = "grey60", segment.size = 0.35,
      min.segment.length = 0.20, max.overlaps = Inf
    ) +
    scale_fill_manual(values = c(
      "Low importance"                = "grey85",
      "High importance"               = "grey55",
      "High importance + interaction" = "black"
    )) +
    scale_shape_manual(values = c(
      "Low importance"                = 21,
      "High importance"               = 22,
      "High importance + interaction" = 24
    )) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(
      title = paste("Morris sensitivity screening —", iter_label),
      x = expression(mu^"*"),
      y = expression(sigma),
      fill = NULL, shape = NULL
    ) +
    .theme_pub(base_size = 12)
}

# ==============================================================================
# HYDROGRAPH - observed vs simulated (all runs + optional accepted highlight)
# ==============================================================================

#' Hydrograph plot: observed vs simulated streamflow
#'
#' @param comparativo data.frame with Date, Flow, and run_* columns
#' @param behavioral_runs Character vector of satisfactory run IDs (NULL = no highlight)
#' @param title_label Plot title
#' @return ggplot object
plot_hydrograph <- function(comparativo, behavioral_runs = NULL,
                            title_label = "Observed vs simulated streamflow") {

  has_behav <- !is.null(behavioral_runs) && length(behavioral_runs) > 0
  sim_runs  <- setdiff(names(comparativo), c("Date", "Flow"))

  sim_long <- comparativo %>%
    select(-Flow) %>%
    pivot_longer(cols = -Date, names_to = "run_id", values_to = "sim") %>%
    mutate(category = if (has_behav) {
      ifelse(run_id %in% behavioral_runs, "Satisfactory", "Non-satisfactory")
    } else {
      "Non-satisfactory"
    })

  obs_df <- comparativo %>%
    select(Date, Flow) %>%
    mutate(category = "Observed")

  # Three categories total; we plot in z-order (back to front):
  #   1. Non-satisfactory  : light grey, thin
  #   2. Satisfactory      : muted blue, slightly thicker
  #   3. Observed          : black DASHED, on top
  pal <- c("Non-satisfactory" = "grey80",
           "Satisfactory"     = "#3B7DB1",
           "Observed"         = "black")

  g <- ggplot()

  # Layer 1 — non-satisfactory sims (always present)
  g <- g +
    geom_line(
      data = sim_long %>% filter(category == "Non-satisfactory"),
      aes(x = Date, y = sim, group = run_id, colour = category),
      alpha = 0.55, linewidth = 0.22
    )

  # Layer 2 — satisfactory sims (only if any)
  if (has_behav) {
    g <- g +
      geom_line(
        data = sim_long %>% filter(category == "Satisfactory"),
        aes(x = Date, y = sim, group = run_id, colour = category),
        alpha = 0.85, linewidth = 0.32
      )
  }

  # Layer 3 — observed (dashed, on top)
  g <- g +
    geom_line(
      data = obs_df,
      aes(x = Date, y = Flow, colour = category),
      linewidth = 0.55, linetype = "dashed"
    )

  # Legend keys: order Observed, Satisfactory, Non-satisfactory; show
  # solid lines for sims and dashed for observed.
  legend_breaks <- if (has_behav) {
    c("Observed", "Satisfactory", "Non-satisfactory")
  } else {
    c("Observed", "Non-satisfactory")
  }
  legend_lty <- if (has_behav) c("dashed", "solid", "solid")
                else            c("dashed", "solid")
  legend_lwd <- if (has_behav) c(0.55, 0.40, 0.35)
                else            c(0.55, 0.35)

  g +
    scale_colour_manual(values = pal, breaks = legend_breaks) +
    guides(colour = guide_legend(
      override.aes = list(linetype = legend_lty, linewidth = legend_lwd)
    )) +
    labs(title  = title_label, x = NULL, colour = NULL,
         y      = expression(Streamflow ~ (m^3 ~ s^{-1}))) +
    .theme_pub(base_size = 12)
}

# ==============================================================================
# SCATTER — NSE vs KGE performance space
# ==============================================================================

#' Performance space scatter plot (NSE vs KGE)
#'
#' @param metrics Tibble with NSE, KGE, satisfactory_run
#' @param nse_col Name of the NSE column (e.g., "NSE_mon")
#' @param kge_col Name of the KGE column (e.g., "KGE_mon")
#' @param threshold_nse NSE reference line
#' @param threshold_kge KGE reference line
#' @param title_label Plot title
#' @return ggplot object
plot_performance_scatter <- function(metrics, nse_col, kge_col,
                                     threshold_nse, threshold_kge,
                                     title_label = "Performance space") {
  ggplot() +
    annotate("rect",
             xmin = threshold_nse, xmax = Inf,
             ymin = threshold_kge, ymax = Inf,
             fill = "grey94", alpha = 0.70) +
    geom_point(
      data = metrics %>% filter(!satisfactory_run),
      aes(x = .data[[nse_col]], y = .data[[kge_col]]),
      shape = 21, size = 2.4, stroke = 0.20,
      colour = "grey25", fill = "grey85", alpha = 0.75
    ) +
    geom_point(
      data = metrics %>% filter(satisfactory_run),
      aes(x = .data[[nse_col]], y = .data[[kge_col]]),
      shape = 21, size = 2.9, stroke = 0.30,
      colour = "black", fill = "#3B7DB1", alpha = 0.92
    ) +
    geom_vline(xintercept = threshold_nse, linetype = "dashed",
               colour = "grey35", linewidth = 0.40) +
    geom_hline(yintercept = threshold_kge, linetype = "dashed",
               colour = "grey35", linewidth = 0.40) +
    labs(title = title_label, x = "NSE", y = "KGE") +
    .theme_pub(base_size = 12)
}

# ==============================================================================
# SCATTER - acceptance domain (RC × obj_bfi_mdae)
# ==============================================================================

#' Acceptance domain scatter plot (BFI × Q_total/P)
#'
#' Layout adopted from the Cap 3 "Band composition diagnostic" figure
#' (CAP3/scripts/composition_diagnostic.R), with axes:
#'   x: median BFI across the N target subbasins (per run, simulated)
#'   y: Q_total / P (basin) = (wateryld + baseflow) / P
#'
#' The box width on each axis encodes the acceptance tolerance:
#'   x-tolerance: median(BFI_obs) ± threshold_obj_bfi
#'                 (proxy for the median absolute error filter Φ_BFI ≤ T_BFI)
#'   y-tolerance: RC_target ± threshold_obj_rc
#'
#' Both axes are physical quantities (not error metrics) — analogous to the
#' Cap 3 plot where X = ET/P and Y = WYLD/P. Run colour reflects the rigorous
#' Φ_BFI MdAE filter applied to the satisfactory_run column.
#'
#' @param ranking            Tibble with rc_sim (or rto_q), bfi_sim_ch* per
#'                           target channel, satisfactory_run
#' @param rc_target          Numeric — observed Q_total / P (basin)
#' @param bfi_obs_median     Numeric — median(BFI_obs) across target subbasins
#' @param threshold_obj_rc   Numeric — RC tolerance half-width
#' @param threshold_obj_bfi  Numeric — BFI tolerance half-width (= filter T_BFI)
#' @param title_label        Plot title
#' @return ggplot object
plot_balance_scatter <- function(ranking,
                                 rc_target,
                                 bfi_obs_median,
                                 threshold_obj_rc,
                                 threshold_obj_bfi,
                                 title_label = "Acceptance domain") {
  y_var <- if ("rc_sim" %in% names(ranking)) "rc_sim" else "rto_q"
  if (!y_var %in% names(ranking)) {
    stop("plot_balance_scatter: expected 'rc_sim' or 'rto_q' in ranking.")
  }

  bfi_cols <- grep("^bfi_sim_ch[0-9]+$", names(ranking), value = TRUE)
  if (length(bfi_cols) == 0) {
    stop("plot_balance_scatter: no 'bfi_sim_chXX' columns found in ranking.")
  }
  ranking$bfi_sim_median <- apply(ranking[, bfi_cols, drop = FALSE],
                                   1, median, na.rm = TRUE)

  # Axis limits — anchor on target ± 1.4× tolerance, expand to fit data
  x_data_min <- min(ranking$bfi_sim_median, na.rm = TRUE)
  x_data_max <- max(ranking$bfi_sim_median, na.rm = TRUE)
  x_lo <- min(bfi_obs_median - 1.4 * threshold_obj_bfi, x_data_min)
  x_hi <- max(bfi_obs_median + 1.4 * threshold_obj_bfi, x_data_max)

  y_data_min <- min(ranking[[y_var]], na.rm = TRUE)
  y_data_max <- max(ranking[[y_var]], na.rm = TRUE)
  y_lo <- min(rc_target - 1.4 * threshold_obj_rc, y_data_min)
  y_hi <- max(rc_target + 1.4 * threshold_obj_rc, y_data_max)

  n_sat <- sum(ranking$satisfactory_run, na.rm = TRUE)
  n_tot <- nrow(ranking)
  n_gauges <- length(bfi_cols)

  ggplot(ranking, aes(x = bfi_sim_median, y = .data[[y_var]],
                      fill = satisfactory_run)) +
    # Tolerance box: acceptance region in the (BFI, RC) plane
    annotate("rect",
             xmin = bfi_obs_median - threshold_obj_bfi,
             xmax = bfi_obs_median + threshold_obj_bfi,
             ymin = rc_target - threshold_obj_rc,
             ymax = rc_target + threshold_obj_rc,
             fill = "grey85", alpha = 0.40,
             colour = "grey50", linewidth = 0.30) +
    # Target × at (BFI_obs_median, RC_obs)
    annotate("point", x = bfi_obs_median, y = rc_target,
             shape = 4, size = 3, stroke = 1, colour = "black") +
    geom_point(shape = 21, size = 2.4, alpha = 0.85,
               stroke = 0.4, colour = "black") +
    scale_fill_manual(
      values = c("FALSE" = "#F5C089", "TRUE" = "#A8C5E8"),
      labels = c("FALSE" = "Non-satisfactory", "TRUE" = "Satisfactory"),
      drop   = FALSE
    ) +
    labs(
      title    = title_label,
      subtitle = sprintf(
        paste("Each point = 1 parameter set. Box = acceptance region",
              "(BFI ± %.2f, RC ± %.2f around ×).  %d / %d satisfactory  (%d gauges)."),
        threshold_obj_bfi, threshold_obj_rc, n_sat, n_tot, n_gauges),
      x = expression(median ~ BFI[sim] ~ "across target subbasins"),
      y = expression(Q[total] / P ~ "(basin)")
    ) +
    coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(y_lo, y_hi), expand = FALSE) +
    .theme_pub(base_size = 12) +
    theme(legend.position = c(0.85, 0.20))
}

# ==============================================================================
# BOXPLOT - normalized accepted parameters
# ==============================================================================

#' Normalized boxplot of accepted parameter distributions
#'
#' @param satisfactory_results Tibble of satisfactory runs
#' @param param_info Tibble with short, par_name, min, max
#' @param title_label Plot title
#' @return ggplot object
plot_param_boxplot <- function(satisfactory_results, param_info,
                               title_label = "Accepted parameter distributions") {
  prior_df <- param_info %>%
    transmute(parameter = short, par_name, prior_min = min, prior_max = max)

  param_long <- satisfactory_results %>%
    select(all_of(param_info$par_name)) %>%
    pivot_longer(cols = everything(), names_to = "par_name", values_to = "value") %>%
    left_join(param_info %>% select(short, par_name), by = "par_name") %>%
    rename(parameter = short) %>%
    filter(!is.na(parameter), !is.na(value))

  param_norm <- param_long %>%
    left_join(prior_df, by = c("parameter", "par_name")) %>%
    mutate(
      value_norm = case_when(
        is.na(prior_min) | is.na(prior_max) ~ NA_real_,
        prior_max <= prior_min ~ NA_real_,
        TRUE ~ (value - prior_min) / (prior_max - prior_min)
      )
    ) %>%
    filter(!is.na(value_norm))

  param_norm$parameter <- factor(param_norm$parameter, levels = param_info$short)

  ggplot(param_norm, aes(x = parameter, y = value_norm)) +
    geom_jitter(width = 0.22, alpha = 0.30, colour = "black", size = 1.1) +
    geom_boxplot(
      width = 0.28, fill = "white", colour = "black",
      linewidth = 0.40, outlier.shape = NA, alpha = 0.85
    ) +
    geom_hline(yintercept = c(0, 1), linetype = "dashed",
               colour = "grey45", linewidth = 0.40) +
    geom_hline(yintercept = 0.5, linetype = "dotted",
               colour = "grey65", linewidth = 0.35) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title_label, x = NULL, y = "Normalised parameter value") +
    .theme_pub(base_size = 12) +
    theme(
      axis.text.x = element_text(family = BASE_FAMILY, angle = 45,
                                  hjust = 1, face = "italic",
                                  size = 12, colour = "black")
    )
}

# ==============================================================================
# HELPER — Save plot as TIFF (publication-ready)
# ==============================================================================

#' Save plot as TIFF for publication
#'
#' Default 300 DPI (publication quality). LZW compression. Times New Roman
#' is inherited from the centralised .theme_pub() applied to every plot
#' function in this file.
#'
#' @param plot     ggplot object
#' @param filename File name (without path)
#' @param save_path Destination folder
#' @param width    Width in cm (default 18)
#' @param height   Height in cm (default 15)
#' @param dpi      Resolution in dots per inch (default 300 — publication)
save_tiff <- function(plot, filename, save_path,
                      width = 18, height = 15, dpi = 300) {
  ggsave(
    filename    = file.path(save_path, filename),
    plot        = plot,
    width       = width,
    height      = height,
    units       = "cm",
    dpi         = dpi,
    compression = "lzw",
    device      = "tiff"
  )
}
