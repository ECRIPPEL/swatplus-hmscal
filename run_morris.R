# ==============================================================================
# run_morris.R - Main Morris screening script (generic for any iteration)
#
# USAGE:
#   Rscript run_morris.R config/iteration1.yaml       # from terminal
#
#   ITERATION_CONFIG <- "config/iteration1.yaml"        # from RStudio/VSCode
#   source("run_morris.R")
# ==============================================================================

.bak_config <- if (exists("ITERATION_CONFIG")) ITERATION_CONFIG else NULL
.bak_master <- if (exists(".MASTER_STATE")) .MASTER_STATE else NULL
rm(list = ls())
gc()
if (!is.null(.bak_config)) ITERATION_CONFIG <- .bak_config
if (!is.null(.bak_master)) .MASTER_STATE <- .bak_master
rm(.bak_config, .bak_master)

library(SWATrunR)
library(yaml)
library(dplyr)
library(tibble)
library(lhs)
library(remotes)
library(purrr)


# ------------------------------------------------------------------------------
# 1. READ CONFIG
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (exists("ITERATION_CONFIG")) {
  config_path <- ITERATION_CONFIG
} else if (length(args) >= 1) {
  config_path <- args[1]
} else {
  stop("Usage: Rscript run_morris.R config/iterationN.yaml\n  Or set ITERATION_CONFIG before source().")
}

cfg <- yaml::read_yaml(config_path)
cat("=== MORRIS SCREENING - before Iteration", cfg$iteration, "===\n")
cat("Config:", config_path, "\n")

# ------------------------------------------------------------------------------
# 2. LOAD MODULES
# ------------------------------------------------------------------------------
source("R/versioning.R")
source("R/lhoat_engine.R")
source("R/param_inherit.R")
source("R/run_filter.R")
source("R/metrics.R")
source("R/morris_classify.R")
source("R/plots.R")
source("R/write_run_info.R")
source("R/build_outputs.R")
source("R/bfi_eckhardt.R")
source("R/baseflow_comparison.R")
if (isTRUE(cfg$check_range$enabled)) source("R/check_range.R")

set.seed(cfg$seed)

# ------------------------------------------------------------------------------
# 3. VERSIONING
# ------------------------------------------------------------------------------
morris_tag_prefix <- cfg$morris$tag_prefix
if (is.null(morris_tag_prefix) || !nzchar(morris_tag_prefix)) {
  stop("config error: morris.tag_prefix is required (e.g., \"MORRIS1\", \"MORRIS2\", \"MORRIS3\")")
}
ver <- create_versioned_folder(cfg$project_path, paste0(morris_tag_prefix, "_v"))
cat("TAG:", ver$tag, "\nFOLDER:", ver$path, "\n")

# ------------------------------------------------------------------------------
# 4. PARAMETERS
# ------------------------------------------------------------------------------
if (cfg$iteration == 1) {
  # Iteration 1: parameters defined in YAML
  param_info <- tibble(
    short    = sapply(cfg$parameters, `[[`, "short"),
    par_name = sapply(cfg$parameters, `[[`, "par_name"),
    min      = sapply(cfg$parameters, `[[`, "min"),
    max      = sapply(cfg$parameters, `[[`, "max")
  )
} else {
  # Iterations 2-4: inherit ranges from the previous iteration
  if (is.null(cfg$inherit$from) || !nzchar(cfg$inherit$from)) {
    stop("config error: inherit.from is required for iteration ", cfg$iteration,
         " (e.g., \"ITER1\", \"ITER2\", \"ITER3\")")
  }
  inherited <- find_inherited_ranges(cfg$project_path, cfg$inherit$from)

  param_map <- unlist(cfg$param_map)

  extra_params <- NULL
  if (!is.null(cfg$extra_params)) {
    extra_params <- tibble(
      short    = sapply(cfg$extra_params, `[[`, "short"),
      par_name = sapply(cfg$extra_params, `[[`, "par_name"),
      min      = sapply(cfg$extra_params, `[[`, "min"),
      max      = sapply(cfg$extra_params, `[[`, "max"),
      grupo    = NA_character_,
      q_low    = NA_real_,
      q_high   = NA_real_,
      post_min = NA_real_,
      post_max = NA_real_,
      mediana  = NA_real_
    )
  }

  param_info <- build_param_info(
    ranges_data     = inherited$data,
    param_map       = param_map,
    extra_params    = extra_params,
    drop            = cfg$inherit$drop,
    override_ranges = cfg$inherit$override_ranges
  )
}

p <- nrow(param_info)
cat("Parameters:", p, "\n")
print(param_info)

# ------------------------------------------------------------------------------
# 5. CHECK PARAMETER RANGE (Iteration 1 only)
# ------------------------------------------------------------------------------
if (isTRUE(cfg$check_range$enabled)) {
  check_result <- check_parameter_range(param_info, cfg$project_path, cfg$check_range$map)
}

# ------------------------------------------------------------------------------
# 6. GENERATE LH-OAT TRAJECTORIES
# ------------------------------------------------------------------------------
morris_start  <- if (!is.null(cfg$morris$start_date)) cfg$morris$start_date else cfg$start_date
morris_end    <- if (!is.null(cfg$morris$end_date))   cfg$morris$end_date   else cfg$end_date
morris_warmup <- if (!is.null(cfg$morris$warmup))     cfg$morris$warmup     else cfg$warmup

design_norm <- generate_lhoat_trajectories(
  m_traj      = cfg$morris$m_traj,
  p           = p,
  delta       = cfg$morris$delta,
  param_names = param_info$short
)

scaled <- scale_parameters(design_norm, param_info)

cat("Expected runs:", nrow(scaled$swat_params), "\n")

# ------------------------------------------------------------------------------
# 7. CLEAN CACHE
# ------------------------------------------------------------------------------
cache_folders <- list.files(cfg$project_path, pattern = "run_|lhoat_",
                            full.names = TRUE, include.dirs = TRUE)
cache_folders <- cache_folders[file.info(cache_folders)$isdir]
if (length(cache_folders) > 0) unlink(cache_folders, recursive = TRUE, force = TRUE)

# ------------------------------------------------------------------------------
# 8. DEFINE SWAT+ OUTPUTS
#    Base block: cfg$morris_outputs (always required).
#    Optional bulk block: cfg$extra_outputs (mirrors swatplus-calibrator's
#    h3_outputs pattern — channels, LSUs, HRUs, aquifers).
# ------------------------------------------------------------------------------
outputs <- build_outputs(cfg$morris_outputs, cfg$extra_outputs)

# ------------------------------------------------------------------------------
# 9. RUN SWAT+
# ------------------------------------------------------------------------------
sim <- run_swatplus(
  project_path = cfg$project_path,
  output       = outputs,
  parameter    = scaled$swat_params,
  start_date   = morris_start,
  end_date     = morris_end,
  years_skip   = morris_warmup,
  n_thread     = cfg$threads,
  save_file    = ver$tag,
  keep_folder  = TRUE
)

# ------------------------------------------------------------------------------
# 10. VALID RUNS
# ------------------------------------------------------------------------------
valid_runs <- get_valid_runs(sim$simulation, names(outputs))
cat("Valid runs:", length(valid_runs), "\n")

# ------------------------------------------------------------------------------
# 11. METRICS
# ------------------------------------------------------------------------------

# BFI multisite (loaded once for all iterations)
bfi_metrics <- NULL
if (!is.null(cfg$bfi_targets_csv)) {
  bfi_csv_path <- file.path(cfg$project_path, cfg$bfi_targets_csv)
  bfi_targets  <- read_bfi_targets(bfi_csv_path)
  cat("[bfi] targets loaded :", nrow(bfi_targets), "channels from",
      cfg$bfi_targets_csv, "\n")

  bfi_output_name  <- cfg$bfi$multi_day_output %||% "vazao_multi_day"
  bfi_period_start <- cfg$bfi$period_start     %||% "1980-01-01"
  bfi_period_end   <- cfg$bfi$period_end       %||% "1989-12-31"
  bfi_eps          <- cfg$bfi$eps              %||% 0.05
  bfi_min_days     <- cfg$bfi$min_days         %||% 365
  bfi_sim_source   <- cfg$bfi$sim_source       %||% "swat_internal_bf"

  bfi_metrics <- compute_bfi_multisite(
    sim          = sim,
    output_name  = bfi_output_name,
    bfi_targets  = bfi_targets,
    valid_runs   = valid_runs,
    sim_source   = bfi_sim_source,
    project_path = cfg$project_path,
    period_start = bfi_period_start,
    period_end   = bfi_period_end,
    eps          = bfi_eps,
    min_days     = bfi_min_days
  )
  cat(sprintf("[bfi] sim_source=%s — computed for %d runs x %d channels\n",
              bfi_sim_source, nrow(bfi_metrics), nrow(bfi_targets)))
}

if (cfg$iteration == 1) {
  # Iteration 1: annual RC + multisite BFI proxy
  if (is.null(bfi_metrics)) {
    stop("The Iteration 1 screening requires bfi_targets_csv in YAML.")
  }
  if (is.null(cfg$morris_targets$rc)) {
    stop("The Iteration 1 screening requires morris_targets.rc (= observed Q_total/P).")
  }
  metrics <- calc_iteration1_metrics(sim, valid_runs,
                                      rc_target      = cfg$morris_targets$rc,
                                      bfi_metrics    = bfi_metrics,
                                      project_path   = cfg$project_path,
                                      outlet_channel = cfg$channel_exutorio)
} else {
  # Iterations 2 to 4: streamflow metrics (NSE/KGE/PBIAS) + RC + BFI append
  obs_file_key <- if (cfg$temporal_scale == "daily") "daily" else "monthly"
  obs_path <- file.path(cfg$project_path, cfg$obs_files[[obs_file_key]])
  obs_data <- read.csv(obs_path, stringsAsFactors = FALSE) %>%
    mutate(Date = as.Date(Date))

  if (cfg$temporal_scale == "monthly") {
    obs_data <- obs_data %>%
      mutate(Date = lubridate::floor_date(Date, "month"))
  }

  obs_data <- obs_data %>%
    filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(morris_end)) %>%
    arrange(Date)

  floor_mon <- cfg$temporal_scale == "monthly"
  comparativo <- build_comparativo(sim$simulation$vazao_exutorio, obs_data, valid_runs, floor_mon)

  filt <- filter_invalid_runs(comparativo, valid_runs)
  valid_runs <- filt$runs_ok
  comparativo <- comparativo %>% select(Date, Flow, all_of(valid_runs))

  metrics <- calc_hard_metrics(comparativo, valid_runs, suffix = cfg$metric_suffix)
  metrics <- metrics %>% mutate(run_index = as.numeric(sub("^run_0*([0-9]+).*$", "\\1", run_id)))

  # Append RC + BFI metrics so Morris EE also ranks on those constraints.
  # obj_total is extended so Morris ranks parameters across all 5 metrics
  # (NSE, KGE, PBIAS, obj_rc, obj_bfi_mdae).
  if (!is.null(bfi_metrics) && !is.null(cfg$morris_targets$rc)) {
    rc_df <- calc_obj_rc(sim, valid_runs, cfg$morris_targets$rc,
                         project_path   = cfg$project_path,
                         outlet_channel = cfg$channel_exutorio) %>%
      select(run_id, rc_sim, obj_rc, rto_et, rto_wyld, rto_baseflow, rto_q,
             rc_hru_level, channel_losses)
    metrics <- metrics %>%
      left_join(rc_df,                                        by = "run_id") %>%
      left_join(bfi_metrics %>% select(run_id, obj_bfi_mdae,
                                        obj_bfi_mae, n_passed_eps05),
                by = "run_id") %>%
      mutate(obj_total = obj_total + obj_rc + obj_bfi_mdae)
  }
}

# ------------------------------------------------------------------------------
# 12. ATTACH TRAJECTORY + COMPUTE EE
# ------------------------------------------------------------------------------
metrics <- design_norm %>%
  select(run_index, traj_id, step_id, all_of(param_info$short)) %>%
  inner_join(metrics, by = "run_index") %>%
  arrange(traj_id, step_id)

ee_all <- calc_all_ee(metrics, cfg$morris_metrics, param_info$short)
sens   <- compute_sensitivity(ee_all)

cat("\n=== GLOBAL RANKING ===\n")
print(sens$ranking)

# ------------------------------------------------------------------------------
# 13. MORRIS SCREENING PLOT
# ------------------------------------------------------------------------------
ranking_classified <- classify_morris(sens$ranking)
g_morris <- plot_morris(ranking_classified, ver$tag)
print(g_morris)
save_tiff(g_morris, paste0("morris_screening_", ver$tag, ".tif"), ver$path)

# ------------------------------------------------------------------------------
# 14. SAVE RESULTS
# ------------------------------------------------------------------------------
iter_label <- cfg$iteration

saveRDS(
  list(
    script_version    = cfg$script_version,
    config            = cfg,
    tag               = ver$tag,
    sim               = sim,
    ranking           = sens$ranking,
    sensitivity       = sens$sensitivity,
    metrics           = metrics,
    ee_all            = ee_all,
    param_info        = param_info,
    lhoat_design_norm = design_norm,
    swat_params       = scaled$swat_params
  ),
  file.path(ver$path, paste0("morris_results_", ver$tag, ".rds"))
)

write.csv(sens$ranking,
          file.path(ver$path, paste0("morris_ranking_", ver$tag, ".csv")),
          row.names = FALSE)
write.csv(sens$sensitivity,
          file.path(ver$path, paste0("morris_sensitivity_", ver$tag, ".csv")),
          row.names = FALSE)
write.csv(metrics,
          file.path(ver$path, paste0("morris_metrics_", ver$tag, ".csv")),
          row.names = FALSE)
write.csv(ee_all,
          file.path(ver$path, paste0("morris_ee_", ver$tag, ".csv")),
          row.names = FALSE)
write.csv(param_info,
          file.path(ver$path, paste0("param_info_", ver$tag, ".csv")),
          row.names = FALSE)

# ---- Run metadata (run_info.txt) ----
top_param <- if (nrow(sens$ranking) > 0) sens$ranking$parameter[1] else NA
top_mu    <- if (nrow(sens$ranking) > 0) sens$ranking$mu_star[1]   else NA
write_run_info(
  path       = ver$path,
  cfg        = cfg,
  run_type   = "Morris screening",
  param_info = param_info,
  summary    = list(
    n_simulations  = nrow(design_norm),
    n_valid        = length(valid_runs),
    top_param      = top_param,
    top_param_mu   = top_mu,
    n_metrics      = length(cfg$morris_metrics %||% character())
  )
)

cat("\nFiles saved to:", ver$path, "\n")
cat("Morris screening", iter_label, "COMPLETED.\n")
