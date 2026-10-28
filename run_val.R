# ==============================================================================
# run_val.R — Split-sample validation of an HMSCal calibrated ensemble
#
# Loads parameter sets from a previous calibration (typically the satisfactory
# subset of an iteration) and propagates them through SWAT+ over an independent
# validation window. Computes the SAME metric set used at calibration:
#   - Daily NSE / KGE / |PBIAS| (Moriasi)
#   - obj_rc (channel-routed Q_outlet/P deviation, Eq. 1)
#   - obj_bfi_mdae (multi-site SWAT internal aquifer-to-channel BFI, Eq. 4)
#
# Classifies each propagated run as "still_satisfactory" if it meets ALL FIVE
# the same thresholds at validation, in addition to having met them at calibration.
#
# IMPORTANT: the BFI targets and RC target MUST be recomputed for the
# validation window (BFI_obs, Q_obs/P_obs are climate-period-dependent).
#
# USAGE:
#   Rscript run_val.R config/validation.yaml
#
#   ITERATION_CONFIG <- "config/validation.yaml"
#   source("run_val.R")
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
library(tidyr)
library(purrr)
library(hydroGOF)
library(lubridate)

# ------------------------------------------------------------------------------
# 1. READ CONFIG
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (exists("ITERATION_CONFIG")) {
  config_path <- ITERATION_CONFIG
} else if (length(args) >= 1) {
  config_path <- args[1]
} else {
  stop("Usage: Rscript run_val.R config/validation.yaml\n",
       "  Or set ITERATION_CONFIG before source().")
}

cfg <- yaml::read_yaml(config_path)
cat("=== VAL — split-sample validation ===\n")
cat("Config:", config_path, "\n")

if (is.null(cfg$validation$param_sets_csv)) {
  stop("config error: validation.param_sets_csv is required (path to ",
       "calibration satisfactory CSV with SWATrunR-style :: parameter columns).")
}

# ------------------------------------------------------------------------------
# 2. LOAD MODULES
# ------------------------------------------------------------------------------
source("R/versioning.R")
source("R/run_filter.R")
source("R/metrics.R")
source("R/plots.R")
source("R/write_run_info.R")
source("R/build_outputs.R")
source("R/bfi_eckhardt.R")
source("R/baseflow_comparison.R")

set.seed(cfg$seed)

# ------------------------------------------------------------------------------
# 3. VERSIONING
# ------------------------------------------------------------------------------
val_tag_prefix <- cfg$validation$tag_prefix %||% "VAL"
ver <- create_versioned_folder(cfg$project_path, paste0(val_tag_prefix, "_v"))
cat("TAG   :", ver$tag, "\nFOLDER:", ver$path, "\n")

# ------------------------------------------------------------------------------
# 4. LOAD PARAMETER SETS FROM CALIBRATION CSV
# ------------------------------------------------------------------------------
param_sets_path <- cfg$validation$param_sets_csv
if (!file.exists(param_sets_path)) {
  # Try resolving relative to project_path
  alt <- file.path(cfg$project_path, param_sets_path)
  if (file.exists(alt)) {
    param_sets_path <- alt
  } else {
    stop("validation.param_sets_csv not found at either:\n",
         "  literal     : ", cfg$validation$param_sets_csv, "\n",
         "  project_rel : ", alt, "\n",
         "Use an absolute path or copy the file under <project_path>.")
  }
}

param_sets <- read.csv(param_sets_path, stringsAsFactors = FALSE,
                        check.names = FALSE) %>% as_tibble()

par_cols <- grep("::", names(param_sets), value = TRUE, fixed = TRUE)
if (length(par_cols) == 0) {
  stop("No SWATrunR-style parameter columns (containing '::') in: ",
       param_sets_path)
}

val_params <- param_sets[, par_cols, drop = FALSE]
n_sim      <- nrow(val_params)

# Preserve original calibration run ids for traceability
orig_run_ids <- if ("run_id" %in% names(param_sets)) {
  as.character(param_sets$run_id)
} else {
  sprintf("run_%03d", seq_len(n_sim))
}

cat("Param sets       :", basename(param_sets_path), "\n")
cat("Ensemble size    :", n_sim, "parameter sets\n")
cat("Parameters/set   :", length(par_cols), "\n")

# ------------------------------------------------------------------------------
# 5. DEFINE SWAT+ OUTPUTS
# ------------------------------------------------------------------------------
outputs <- build_outputs(cfg$cal_outputs, cfg$extra_outputs)

# ------------------------------------------------------------------------------
# 6. CLEAN CACHE
# ------------------------------------------------------------------------------
cache_folders <- list.files(cfg$project_path,
                             pattern = "run_|lhoat_|soft_|hard_|val_",
                             full.names = TRUE, include.dirs = TRUE)
cache_folders <- cache_folders[file.info(cache_folders)$isdir]
if (length(cache_folders) > 0) unlink(cache_folders, recursive = TRUE,
                                       force = TRUE)

# ------------------------------------------------------------------------------
# 7. RUN SWAT+ ON VALIDATION WINDOW
# ------------------------------------------------------------------------------
val_start  <- cfg$cal$start_date %||% cfg$start_date
val_end    <- cfg$cal$end_date   %||% cfg$end_date
val_warmup <- cfg$cal$warmup     %||% cfg$warmup

cat("Validation window:", val_start, "->", val_end,
    "  (warmup", val_warmup, "yr)\n")

sim <- run_swatplus(
  project_path = cfg$project_path,
  output       = outputs,
  parameter    = val_params,
  start_date   = val_start,
  end_date     = val_end,
  years_skip   = val_warmup,
  n_thread     = cfg$threads,
  save_file    = ver$tag,
  keep_folder  = TRUE
)

# ------------------------------------------------------------------------------
# 8. VALID RUNS
# ------------------------------------------------------------------------------
valid_runs <- get_valid_runs(sim$simulation, names(outputs))
cat("Valid runs       :", length(valid_runs), "/", n_sim, "\n")

# ------------------------------------------------------------------------------
# 9. BFI MULTI-SITE (same logic as run_cal.R)
# ------------------------------------------------------------------------------
bfi_metrics <- NULL
bfi_targets <- NULL
if (!is.null(cfg$bfi_targets_csv)) {
  bfi_csv_path <- file.path(cfg$project_path, cfg$bfi_targets_csv)
  bfi_targets  <- read_bfi_targets(bfi_csv_path)
  cat("[bfi] targets    :", nrow(bfi_targets), "channels from",
      cfg$bfi_targets_csv, "\n")

  bfi_output_name  <- cfg$bfi$multi_day_output %||% "vazao_multi_day"
  bfi_period_start <- cfg$bfi$period_start     %||% val_start
  bfi_period_end   <- cfg$bfi$period_end       %||% val_end
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

# Diagnostic: SWAT+ internal aqu_flo_cha vs Eckhardt (same as run_cal.R)
baseflow_internal <- NULL
if (!is.null(bfi_metrics) &&
    isTRUE(cfg$extra_outputs$enabled) &&
    isTRUE(cfg$extra_outputs$aquifer$enabled)) {
  baseflow_internal <- tryCatch({
    compute_baseflow_comparison(
      sim              = sim,
      project_path     = cfg$project_path,
      bfi_targets      = bfi_targets,
      valid_runs       = valid_runs,
      multi_day_output = bfi_output_name
    )
  }, error = function(e) {
    warning("[baseflow] comparison skipped: ", conditionMessage(e))
    NULL
  })
}

# ------------------------------------------------------------------------------
# 10. STREAMFLOW METRICS — daily + monthly diagnostic
# ------------------------------------------------------------------------------
suffix <- if (cfg$temporal_scale == "daily") "_day" else "_mon"

# Primary scale (daily for split-sample validation by convention)
primary_key    <- if (cfg$temporal_scale == "daily") "daily" else "monthly"
primary_output <- if (cfg$temporal_scale == "daily") "vazao_exutorio_day" else "vazao_exutorio_mon"

obs_primary <- read.csv(file.path(cfg$project_path, cfg$obs_files[[primary_key]]),
                         stringsAsFactors = FALSE) %>%
  mutate(Date = as.Date(Date))

if (cfg$temporal_scale == "monthly") {
  obs_primary <- obs_primary %>% mutate(Date = floor_date(Date, "month"))
}

obs_primary <- obs_primary %>%
  filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(val_end)) %>%
  arrange(Date)

floor_mon <- cfg$temporal_scale == "monthly"
comp_primary <- build_comparativo(sim$simulation[[primary_output]], obs_primary,
                                   valid_runs, floor_mon)
filt <- filter_invalid_runs(comp_primary, valid_runs)
valid_runs <- filt$runs_ok
comp_primary <- comp_primary %>% select(Date, Flow, all_of(valid_runs))

metrics <- calc_hard_metrics(comp_primary, valid_runs, suffix = suffix)

# ------------------------------------------------------------------------------
# 11. JOIN RC + BFI INTO METRIC TABLE (same join order as run_cal.R hard branch)
# ------------------------------------------------------------------------------
v2_mode <- !is.null(bfi_metrics) && !is.null(cfg$cal$rc_target)
if (v2_mode) {
  rc_df <- calc_obj_rc(sim, valid_runs, cfg$cal$rc_target,
                       project_path   = cfg$project_path,
                       outlet_channel = cfg$channel_exutorio) %>%
    select(run_id, rto_q, rc_sim, obj_rc, rto_et, rto_wyld, rto_baseflow,
           rc_target, P_sim, ET_sim, WYLD_sim, BASEFLOW_sim, Q_TOTAL_sim,
           rc_hru_level, channel_losses)
  metrics <- metrics %>%
    left_join(rc_df,        by = "run_id") %>%
    left_join(bfi_metrics,  by = "run_id") %>%
    mutate(obj_total = obj_total + obj_rc + obj_bfi_mdae)
}

# ------------------------------------------------------------------------------
# 12. CLASSIFY still_satisfactory (joint criterion at the validation window)
# ------------------------------------------------------------------------------
metrics <- classify_hard_runs(
  metrics,
  threshold_nse        = cfg$cal$threshold_nse,
  threshold_kge        = cfg$cal$threshold_kge,
  threshold_abs_pbias  = cfg$cal$threshold_abs_pbias,
  threshold_obj_rc     = cfg$cal$threshold_obj_rc,
  threshold_obj_bfi    = cfg$cal$threshold_obj_bfi,
  suffix               = suffix
) %>%
  # Keep `satisfactory_run` (used by plot helpers that read this column name)
  # AND expose a semantically clearer alias for the validation context.
  mutate(still_satisfactory = satisfactory_run)

# Map back to original calibration run ids for traceability
val_params_idx <- val_params %>% mutate(run_id = sprintf("run_%d", row_number()))
run_map <- tibble(
  run_id      = valid_runs,
  orig_run_id = orig_run_ids[as.integer(sub("^run_0*", "", valid_runs))]
)

results <- metrics %>%
  left_join(run_map,         by = "run_id") %>%
  left_join(val_params_idx,  by = "run_id")

still_results <- results %>% filter(still_satisfactory)
n_still <- nrow(still_results)

cat("\n=== TEMPORAL STABILITY ===\n")
cat("Still satisfactory at validation:", n_still, "/", length(valid_runs),
    sprintf(" (%.1f%%)\n", 100 * n_still / max(length(valid_runs), 1)))
if (v2_mode) {
  moriasi_only <- sum(metrics[[paste0("NSE", suffix)]] >= cfg$cal$threshold_nse &
                       metrics[[paste0("KGE", suffix)]] >= cfg$cal$threshold_kge &
                       metrics$absPBIAS <= cfg$cal$threshold_abs_pbias,
                       na.rm = TRUE)
  cat("Moriasi-only at validation       :", moriasi_only, "/",
      length(valid_runs),
      sprintf(" (%.1f%%)\n", 100 * moriasi_only / max(length(valid_runs), 1)))
}

# ------------------------------------------------------------------------------
# 13. MONTHLY DIAGNOSTIC (secondary scale)
# ------------------------------------------------------------------------------
secondary_key    <- if (cfg$temporal_scale == "daily") "monthly" else "daily"
secondary_output <- if (cfg$temporal_scale == "daily") "vazao_exutorio_mon" else "vazao_exutorio_day"

if (secondary_output %in% names(sim$simulation) &&
    !is.null(cfg$obs_files[[secondary_key]])) {
  obs_secondary <- read.csv(file.path(cfg$project_path,
                                       cfg$obs_files[[secondary_key]]),
                             stringsAsFactors = FALSE) %>%
    mutate(Date = as.Date(Date))
  if (secondary_key == "monthly") {
    obs_secondary <- obs_secondary %>% mutate(Date = floor_date(Date, "month"))
  }
  obs_secondary <- obs_secondary %>%
    filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(val_end)) %>%
    arrange(Date)
  sec_suffix <- if (secondary_key == "monthly") "_mon" else "_day"
  comp_sec <- build_comparativo(sim$simulation[[secondary_output]],
                                 obs_secondary, valid_runs,
                                 floor_month = (secondary_key == "monthly"))
  metrics_sec <- calc_hard_metrics(comp_sec, valid_runs, suffix = sec_suffix)
  sec_metric_cols <- grep(paste0("(NSE|KGE|PBIAS)", sec_suffix),
                           names(metrics_sec), value = TRUE)
  results <- results %>% left_join(
    metrics_sec %>% select(run_id, all_of(sec_metric_cols)), by = "run_id")
}

# ------------------------------------------------------------------------------
# 14. PLOTS
# ------------------------------------------------------------------------------
nse_col <- paste0("NSE", suffix)
kge_col <- paste0("KGE", suffix)

g_scatter <- plot_performance_scatter(
  metrics, nse_col, kge_col,
  cfg$cal$threshold_nse, cfg$cal$threshold_kge,
  paste(ver$tag, "-", cfg$temporal_scale, "validation"))
print(g_scatter)
save_tiff(g_scatter, paste0("performance_scatter_", ver$tag, ".tif"), ver$path)

primary_suffix_tag <- if (cfg$temporal_scale == "daily") "day" else "mon"
g_hydro <- plot_hydrograph(
  comparativo     = comp_primary,
  behavioral_runs = if (n_still > 0) still_results$run_id else NULL,
  title_label     = paste(ver$tag, "—", cfg$temporal_scale, "streamflow (val)")
)
print(g_hydro)
save_tiff(g_hydro,
           paste0("hydrograph_", ver$tag, "_", primary_suffix_tag, ".tif"),
           ver$path, width = 20, height = 10)

# Acceptance-domain scatter (same axes as cal: median BFI vs RC_sim)
if (v2_mode && "obj_bfi_mdae" %in% names(results) &&
    "rc_sim" %in% names(results)) {
  bfi_obs_median <- if (!is.null(bfi_targets)) median(bfi_targets$bfi_obs)
                    else NA_real_
  g_balance <- plot_balance_scatter(
    ranking            = results,
    rc_target          = cfg$cal$rc_target,
    bfi_obs_median     = bfi_obs_median,
    threshold_obj_rc   = cfg$cal$threshold_obj_rc,
    threshold_obj_bfi  = cfg$cal$threshold_obj_bfi,
    title_label        = paste(ver$tag, "— validation acceptance domain")
  )
  print(g_balance)
  save_tiff(g_balance, paste0("balance_scatter_", ver$tag, ".tif"), ver$path)
}

# ==============================================================================
# 15. SAVE
# ==============================================================================

# Full RDS
saveRDS(
  list(
    script_version    = cfg$script_version %||% "validation",
    config            = cfg,
    tag               = ver$tag,
    sim               = sim,
    param_sets_source = param_sets_path,
    bfi_metrics       = bfi_metrics,
    baseflow_internal = baseflow_internal,
    results           = results,
    still_results     = still_results
  ),
  file.path(ver$path, paste0("simulation_results_", ver$tag, ".rds"))
)

# Run-level CSVs
write.csv(results,
          file.path(ver$path, paste0("simulation_results_", ver$tag, ".csv")),
          row.names = FALSE)
if (n_still > 0) {
  write.csv(still_results,
            file.path(ver$path,
                       paste0("simulation_results_still_satisfactory_",
                              ver$tag, ".csv")),
            row.names = FALSE)
}

# Baseflow comparison CSV
if (!is.null(baseflow_internal) && !is.null(bfi_metrics)) {
  bf_compare <- join_baseflow_with_eckhardt(
    bfi_metrics       = bfi_metrics,
    baseflow_internal = baseflow_internal,
    bfi_targets       = bfi_targets
  )
  write.csv(bf_compare,
            file.path(ver$path,
                       paste0("baseflow_comparison_", ver$tag, ".csv")),
            row.names = FALSE)
}

# Run info
write_run_info(
  path       = ver$path,
  cfg        = cfg,
  run_type   = "VAL",
  param_info = NULL,
  summary    = list(
    tag                  = ver$tag,
    param_sets_source    = basename(param_sets_path),
    n_propagated         = n_sim,
    n_valid              = length(valid_runs),
    n_still_satisfactory = n_still,
    pct_still_satisfactory = sprintf("%.1f%%",
                                      100 * n_still / max(length(valid_runs), 1)),
    rc_target            = cfg$cal$rc_target %||% "(not set)",
    threshold_obj_rc     = cfg$cal$threshold_obj_rc %||% "(not set)",
    threshold_obj_bfi    = cfg$cal$threshold_obj_bfi %||% "(not set)"
  )
)

cat("\nFiles saved to:", ver$path, "\n")
cat("VAL", ver$tag, "COMPLETED.\n")
