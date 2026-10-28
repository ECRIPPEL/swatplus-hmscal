# ==============================================================================
# run_cal.R - Main calibration script (generic for any iteration)
#
# USAGE:
#   Rscript run_cal.R config/iteration1.yaml       # from terminal
#
#   ITERATION_CONFIG <- "config/iteration2.yaml"        # from RStudio/VSCode
#   source("run_cal.R")
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
  stop("Usage: Rscript run_cal.R config/iterationN.yaml\n  Or set ITERATION_CONFIG before source().")
}

cfg <- yaml::read_yaml(config_path)
cat("=== CALIBRATION - Iteration", cfg$iteration, "- type:", cfg$cal$type, "===\n")

# ------------------------------------------------------------------------------
# 2. LOAD MODULES
# ------------------------------------------------------------------------------
source("R/versioning.R")
source("R/param_inherit.R")
source("R/run_filter.R")
source("R/metrics.R")
source("R/morris_classify.R")
source("R/plots.R")
source("R/write_run_info.R")
source("R/build_outputs.R")
source("R/bfi_eckhardt.R")
source("R/baseflow_comparison.R")

set.seed(cfg$seed)

# ------------------------------------------------------------------------------
# 3. VERSIONING
# ------------------------------------------------------------------------------
cal_tag_prefix <- cfg$cal$tag_prefix
if (is.null(cal_tag_prefix) || !nzchar(cal_tag_prefix)) {
  stop("config error: cal.tag_prefix is required ",
       "(e.g., \"ITER1\", \"ITER2\", \"ITER3\", \"ITER4\")")
}
ver <- create_versioned_folder(cfg$project_path, paste0(cal_tag_prefix, "_v"))
cat("TAG:", ver$tag, "\nFOLDER:", ver$path, "\n")

# ------------------------------------------------------------------------------
# 4. LOAD MORRIS SCREENING RESULTS (ranking)
# ------------------------------------------------------------------------------
morris_tag_prefix <- cfg$morris$tag_prefix
if (is.null(morris_tag_prefix) || !nzchar(morris_tag_prefix)) {
  stop("config error: morris.tag_prefix is required to locate the Morris ranking")
}
morris_results <- find_morris_results(cfg$project_path, morris_tag_prefix)

cat("Morris screening results from:", morris_results$folder, "\n")

# ------------------------------------------------------------------------------
# 5. PARAMETERS
# ------------------------------------------------------------------------------
# Init auto_drop info so the SAVE block can reference it for either branch
auto_drop_info <- list(drop = character(0),
                       detail = tibble::tibble(parameter = character(),
                                                mu_star = numeric(),
                                                zero_flag = logical()))

if (cfg$iteration == 1) {
  # Iteration 1: parameters from YAML + extras for calibration
  morris_params <- tibble(
    short    = sapply(cfg$parameters, `[[`, "short"),
    par_name = sapply(cfg$parameters, `[[`, "par_name"),
    min      = sapply(cfg$parameters, `[[`, "min"),
    max      = sapply(cfg$parameters, `[[`, "max")
  )

  if (!is.null(cfg$cal_extra_params)) {
    extras <- tibble(
      short    = sapply(cfg$cal_extra_params, `[[`, "short"),
      par_name = sapply(cfg$cal_extra_params, `[[`, "par_name"),
      min      = sapply(cfg$cal_extra_params, `[[`, "min"),
      max      = sapply(cfg$cal_extra_params, `[[`, "max")
    )
    param_info <- bind_rows(morris_params, extras)
  } else {
    param_info <- morris_params
  }
} else {
  # Iterations 2-4: inherit param_info DIRECTLY from the previous iteration via
  # cfg$inherit$from. (Earlier versions read morris_results$param_info, but that
  # silently broke ITER4, which has no new screening between ITER3 and itself,
  # so it would re-use MORRIS3's prior instead of ITER3's narrowed ranges.)
  if (is.null(cfg$inherit$from) || !nzchar(cfg$inherit$from)) {
    stop("config error: inherit.from is required for iteration ", cfg$iteration,
         " calibration (e.g., \"ITER1\", \"ITER2\", \"ITER3\")")
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

  # ---- Validate manual drop strictly (catches typos in inherit.drop) ---------
  manual_drop <- cfg$inherit$drop %||% character(0)
  if (length(manual_drop) > 0) {
    .all_shorts <- c(inherited$data$parameter,
                     if (!is.null(extra_params)) extra_params$short else character())
    .missing <- setdiff(manual_drop, .all_shorts)
    if (length(.missing) > 0) {
      stop("inherit.drop references unknown parameters: ",
           paste(.missing, collapse = ", "),
           "\nAvailable shorts in inherited + extras: ",
           paste(.all_shorts, collapse = ", "))
    }
  }

  # ---- AUTO-DROP by Morris sensitivity (μ*(obj_total) ≤ threshold) ----------
  # Filters at the calibration stage (not the screening) so morris_screening_<TAG>.tif retains the
  # full ranking — that figure is the evidence justifying the drops.
  auto_drop_info <- compute_auto_drop(
    ranking_morris = morris_results$ranking,
    threshold      = cfg$inherit$auto_drop_threshold,
    protect        = cfg$inherit$auto_drop_protect %||% character(0)
  )

  if (length(auto_drop_info$drop) > 0) {
    cat(sprintf("\n[auto_drop] threshold = %g  |  %d candidate(s):\n",
                cfg$inherit$auto_drop_threshold, length(auto_drop_info$drop)))
    for (i in seq_len(nrow(auto_drop_info$detail))) {
      d <- auto_drop_info$detail[i, ]
      flag <- if (isTRUE(d$zero_flag)) "  [WARN: mu* literal zero]" else ""
      cat(sprintf("  - %-12s (mu* = %g)%s\n", d$parameter, d$mu_star, flag))
    }
    if (any(auto_drop_info$detail$zero_flag)) {
      zeros <- auto_drop_info$detail$parameter[auto_drop_info$detail$zero_flag]
      warning("[auto_drop] Literal-zero mu* for: ", paste(zeros, collapse = ", "),
              ". Possible (a) physical insensitivity, (b) SWATrunR silent fail, ",
              "or (c) SWAT+ initialization-order issue (e.g., channel rating ",
              "curves precomputed in proc_cha before proc_cal applies ",
              "calibration.cal in SWAT+ rev 61.0.2). Verify before publishing.")
    }
  }

  combined_drop <- unique(c(manual_drop, auto_drop_info$drop))

  param_info <- build_param_info(
    ranges_data     = inherited$data,
    param_map       = param_map,
    extra_params    = extra_params,
    drop            = combined_drop,
    override_ranges = cfg$inherit$override_ranges,
    strict_drop     = FALSE  # auto-drop may include params already absent (lenient)
  ) %>% as_tibble() %>% select(short, par_name, min, max)
}

cat("Calibration parameters:", nrow(param_info), "\n")
print(param_info[, c("short", "min", "max")])

# ------------------------------------------------------------------------------
# 6. GENERATE PARAMETER SAMPLES (runif)
# ------------------------------------------------------------------------------
n_sim <- cfg$cal$n_simulations

cal_params <- as_tibble(setNames(
  lapply(seq_len(nrow(param_info)), function(i) {
    runif(n_sim, min = param_info$min[i], max = param_info$max[i])
  }),
  param_info$par_name
))

cat("Simulations:", n_sim, "\n")

# ------------------------------------------------------------------------------
# 7. DEFINE SWAT+ OUTPUTS
#    Base block: cfg$cal_outputs (always required).
#    Optional bulk block: cfg$extra_outputs (mirrors swatplus-calibrator's
#    h3_outputs pattern — channels, LSUs, HRUs, aquifers).
# ------------------------------------------------------------------------------
outputs <- build_outputs(cfg$cal_outputs, cfg$extra_outputs)

# ------------------------------------------------------------------------------
# 8. CLEAN CACHE
# ------------------------------------------------------------------------------
cache_folders <- list.files(cfg$project_path, pattern = "run_|lhoat_|soft_|hard_",
                            full.names = TRUE, include.dirs = TRUE)
cache_folders <- cache_folders[file.info(cache_folders)$isdir]
if (length(cache_folders) > 0) unlink(cache_folders, recursive = TRUE, force = TRUE)

# ------------------------------------------------------------------------------
# 9. RUN SWAT+
# ------------------------------------------------------------------------------
cal_start  <- cfg$cal$start_date
cal_end    <- cfg$cal$end_date
cal_warmup <- cfg$cal$warmup

sim <- run_swatplus(
  project_path = cfg$project_path,
  output       = outputs,
  parameter    = cal_params,
  start_date   = cal_start,
  end_date     = cal_end,
  years_skip   = cal_warmup,
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
# 10b. BFI MULTISITE - load targets CSV and apply Eckhardt to all valid runs
#      The targets CSV is looked up under `cfg$project_path` (same convention as
#      obs_files). Required columns: channel_id, bfi (= BFI_obs), a, bfimax_q9050.
#      Output name in cal_outputs that holds daily Q at the multi-site channels
#      is configurable via cfg$bfi$multi_day_output (default "vazao_multi_day").
# ------------------------------------------------------------------------------
bfi_metrics <- NULL
bfi_targets <- NULL
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

# ------------------------------------------------------------------------------
# 10c. SWAT+ INTERNAL BASEFLOW vs ECKHARDT — diagnostic comparison
#      Gated on extra_outputs.aquifer.enabled (need aqu_flo_cha_<i> outputs).
#      Computes BFI_swat_internal per gauge (aqu_flo_cha aggregated upstream)
#      and joins it next to BFI_eckhardt_sim. Saved separately as
#      baseflow_comparison_<TAG>.csv. Useful for:
#        - detecting equifinality (aqu-active vs latq-dominated regimes)
#        - validating internal partition realism (Section 2.5.X of paper)
#      Channels without a directly-mapped floodplain aquifer in aquifer.con
#      will produce bfi_swat_internal ~ 0; this is documented in §2.3.5.
# ------------------------------------------------------------------------------
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
  if (!is.null(baseflow_internal)) {
    cat(sprintf("[baseflow] computed SWAT+ internal partition for %d runs x %d gauges\n",
                nrow(baseflow_internal), nrow(bfi_targets)))
  }
}

# ==============================================================================
# 11. METRICS + FILTERING (depends on calibration type)
# ==============================================================================

if (cfg$cal$type == "soft") {
  # ---------- ITERATION 1: annual water balance, RC + multisite BFI proxy ----------
  if (is.null(bfi_metrics)) {
    stop("Iteration 1 requires bfi_targets_csv ",
         "in YAML and a multi-channel daily output (default: vazao_multi_day).")
  }
  if (is.null(cfg$cal$rc_target)) {
    stop("Iteration 1 requires cal.rc_target (= observed Q_total/P).")
  }
  if (is.null(cfg$cal$threshold_obj_rc) || is.null(cfg$cal$threshold_obj_bfi)) {
    stop("Iteration 1 requires cal.threshold_obj_rc and ",
         "cal.threshold_obj_bfi in YAML.")
  }

  ranking_balance <- calc_soft_metrics(
    sim                 = sim,
    valid_runs          = valid_runs,
    rc_target           = cfg$cal$rc_target,
    threshold_obj_rc    = cfg$cal$threshold_obj_rc,
    threshold_obj_bfi   = cfg$cal$threshold_obj_bfi,
    bfi_metrics         = bfi_metrics,
    project_path        = cfg$project_path,
    outlet_channel      = cfg$channel_exutorio
  )
  satisfactory <- ranking_balance %>% filter(satisfactory_run)

  cat("Satisfactory runs:", nrow(satisfactory), "\n")

  # Monthly + daily streamflow diagnostics — built for ALL valid runs (not
  # just satisfactory) so the hydrographs can show satisfactory vs. non-
  # satisfactory ensembles side by side.
  comp_mon <- NULL
  comp_day <- NULL
  monthly_metrics <- tibble(run_id = character(), NSE_mon = numeric(),
                            KGE_mon = numeric())

  if ("vazao_exutorio" %in% names(sim$simulation)) {
    obs_mon <- read.csv(file.path(cfg$project_path, cfg$obs_files$monthly),
                         stringsAsFactors = FALSE) %>%
      mutate(Date = as.Date(Date)) %>%
      filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(cal_end)) %>%
      arrange(Date)

    comp_mon <- build_comparativo(sim$simulation$vazao_exutorio, obs_mon,
                                   valid_runs, floor_month = TRUE)

    monthly_metrics <- map_dfr(valid_runs, function(run) {
      tibble(
        run_id  = run,
        NSE_mon = round(NSE(sim = comp_mon[[run]], obs = comp_mon$Flow), 3),
        KGE_mon = round(KGE(sim = comp_mon[[run]], obs = comp_mon$Flow), 3)
      )
    })
  }

  # Daily diagnostic at outlet (ch33). Uses `vazao_exutorio_day` if declared
  # in cal_outputs; falls back to `vazao_multi_day_33` only if 33 happens to
  # be in the BFI multi-channel list (typically NOT — outlet is regulated
  # and excluded from the BFI selection).
  daily_outlet_key <- if ("vazao_exutorio_day" %in% names(sim$simulation)) {
    "vazao_exutorio_day"
  } else if ("vazao_multi_day_33" %in% names(sim$simulation)) {
    "vazao_multi_day_33"
  } else {
    NA_character_
  }
  if (!is.na(daily_outlet_key) &&
      file.exists(file.path(cfg$project_path, cfg$obs_files$daily))) {
    obs_day <- read.csv(file.path(cfg$project_path, cfg$obs_files$daily),
                         stringsAsFactors = FALSE) %>%
      mutate(Date = as.Date(Date)) %>%
      filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(cal_end)) %>%
      arrange(Date)
    comp_day <- build_comparativo(sim$simulation[[daily_outlet_key]],
                                   obs_day, valid_runs, floor_month = FALSE)
  }

  # Consolidated results
  params_out <- cal_params %>% mutate(run_index = row_number())
  results <- ranking_balance %>%
    left_join(params_out, by = "run_index") %>%
    left_join(monthly_metrics, by = "run_id")

  satisfactory_results <- results %>% filter(satisfactory_run)

  # Acceptance-domain plot, physical quantities on both axes:
  # X = median BFI across target subbasins, Y = Q_total/P (basin).
  # Box width on X = 2 × threshold_obj_bfi (Φ_BFI MdAE filter visualised
  # as a tolerance band around median(BFI_obs)).
  bfi_obs_median <- if (!is.null(bfi_targets)) median(bfi_targets$bfi_obs)
                    else NA_real_
  g_scatter <- plot_balance_scatter(
    ranking            = ranking_balance,
    rc_target          = cfg$cal$rc_target,
    bfi_obs_median     = bfi_obs_median,
    threshold_obj_rc   = cfg$cal$threshold_obj_rc,
    threshold_obj_bfi  = cfg$cal$threshold_obj_bfi,
    title_label        = paste(ver$tag, "- acceptance domain")
  )
  print(g_scatter)
  save_tiff(g_scatter, paste0("balance_scatter_", ver$tag, ".tif"), ver$path)

  # Hydrographs at BOTH scales (monthly + daily). Always plotted (regardless
  # of how many runs are satisfactory); accepted runs are highlighted in
  # blue when present, observed is dashed black on top.
  if (!is.null(comp_mon)) {
    g_hydro_mon <- plot_hydrograph(
      comparativo      = comp_mon,
      behavioral_runs  = if (nrow(satisfactory) > 0) satisfactory$run_id else NULL,
      title_label      = paste(ver$tag, "— monthly streamflow")
    )
    print(g_hydro_mon)
    save_tiff(g_hydro_mon, paste0("hydrograph_", ver$tag, "_mon.tif"),
              ver$path, width = 20, height = 10)
  }
  if (!is.null(comp_day)) {
    g_hydro_day <- plot_hydrograph(
      comparativo      = comp_day,
      behavioral_runs  = if (nrow(satisfactory) > 0) satisfactory$run_id else NULL,
      title_label      = paste(ver$tag, "— daily streamflow")
    )
    print(g_hydro_day)
    save_tiff(g_hydro_day, paste0("hydrograph_", ver$tag, "_day.tif"),
              ver$path, width = 20, height = 10)
  }

  # Morris ranking for narrowing (used by unified inheritance block at the end)
  ranking_morris <- morris_results$ranking %>%
    filter(if ("metric" %in% names(.)) metric == "obj_total" else TRUE)

} else {
  # ---------- ITERATIONS 2-4: streamflow criteria ----------
  suffix <- cfg$metric_suffix  # "_mon" or "_day"

  # Primary scale (calibration target)
  primary_key    <- if (cfg$temporal_scale == "daily") "daily" else "monthly"
  primary_output <- if (cfg$temporal_scale == "daily") "vazao_exutorio_day" else "vazao_exutorio_mon"

  obs_primary <- read.csv(file.path(cfg$project_path, cfg$obs_files[[primary_key]]),
                           stringsAsFactors = FALSE) %>%
    mutate(Date = as.Date(Date))

  if (cfg$temporal_scale == "monthly") {
    obs_primary <- obs_primary %>% mutate(Date = floor_date(Date, "month"))
  }

  obs_primary <- obs_primary %>%
    filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(cal_end)) %>%
    arrange(Date)

  floor_mon <- cfg$temporal_scale == "monthly"
  comp_primary <- build_comparativo(sim$simulation[[primary_output]], obs_primary, valid_runs, floor_mon)

  filt <- filter_invalid_runs(comp_primary, valid_runs)
  valid_runs <- filt$runs_ok
  comp_primary <- comp_primary %>% select(Date, Flow, all_of(valid_runs))

  # Primary metrics
  metrics <- calc_hard_metrics(comp_primary, valid_runs, suffix = suffix)

  # Join RC + BFI multisite into the metric table BEFORE classification,
  # so satisfactory_run is the AND of (NSE,KGE,PBIAS) ∧ (obj_rc) ∧ (obj_bfi_mdae).
  # obj_total is also extended to include the RC + BFI components so Morris
  # ranking captures sensitivity across all 5 metrics.
  v2_mode <- !is.null(bfi_metrics) && !is.null(cfg$cal$rc_target)
  if (v2_mode) {
    rc_df <- calc_obj_rc(sim, valid_runs, cfg$cal$rc_target,
                         project_path   = cfg$project_path,
                         outlet_channel = cfg$channel_exutorio) %>%
      select(run_id, rto_q, rc_sim, obj_rc, rto_et, rto_wyld, rto_baseflow,
             rc_target, P_sim, ET_sim, WYLD_sim, BASEFLOW_sim, Q_TOTAL_sim,
             rc_hru_level, channel_losses)
    metrics <- metrics %>%
      left_join(rc_df,                                        by = "run_id") %>%
      left_join(bfi_metrics,                                  by = "run_id") %>%
      mutate(obj_total = obj_total + obj_rc + obj_bfi_mdae)
  }

  metrics <- classify_hard_runs(
    metrics,
    threshold_nse        = cfg$cal$threshold_nse,
    threshold_kge        = cfg$cal$threshold_kge,
    threshold_abs_pbias  = cfg$cal$threshold_abs_pbias,
    threshold_obj_rc     = cfg$cal$threshold_obj_rc,    # NULL = skip the RC predicate
    threshold_obj_bfi    = cfg$cal$threshold_obj_bfi,   # NULL = skip the BFI predicate
    suffix               = suffix
  )

  # Consolidated results
  params_out <- cal_params %>% mutate(run_index = seq_len(n())) %>% relocate(run_index)
  results <- metrics %>% left_join(params_out, by = "run_index")
  satisfactory_results <- results %>% filter(satisfactory_run)

  cat("Satisfactory runs:", nrow(satisfactory_results), "\n")

  # Secondary scale diagnostics
  secondary_key    <- if (cfg$temporal_scale == "daily") "monthly" else "daily"
  secondary_output <- if (cfg$temporal_scale == "daily") "vazao_exutorio_mon" else "vazao_exutorio_day"

  if (secondary_output %in% names(sim$simulation)) {
    obs_secondary <- read.csv(file.path(cfg$project_path, cfg$obs_files[[secondary_key]]),
                               stringsAsFactors = FALSE) %>%
      mutate(Date = as.Date(Date))

    if (secondary_key == "monthly") {
      obs_secondary <- obs_secondary %>% mutate(Date = floor_date(Date, "month"))
    }

    obs_secondary <- obs_secondary %>%
      filter(Date >= as.Date(cfg$obs_start_date) & Date <= as.Date(cal_end)) %>%
      arrange(Date)

    sec_suffix <- if (secondary_key == "monthly") "_mon" else "_day"
    comp_sec <- build_comparativo(sim$simulation[[secondary_output]], obs_secondary,
                                   valid_runs, floor_month = (secondary_key == "monthly"))
    metrics_sec <- calc_hard_metrics(comp_sec, valid_runs, suffix = sec_suffix)
    sec_metric_cols <- grep(paste0("(NSE|KGE|PBIAS)", sec_suffix), names(metrics_sec), value = TRUE)
    results <- results %>% left_join(metrics_sec %>% select(run_id, all_of(sec_metric_cols)), by = "run_id")
  }

  # Water balance diagnostic (daily iterations)
  # rto_* columns were already joined via calc_obj_rc above; skip the
  # legacy diagnostic to avoid duplicate columns.
  if (!v2_mode && !is.null(cfg$diagnostic_balance) && "precip" %in% names(sim$simulation)) {
    diag_balance <- calc_balance_diagnostic(sim, valid_runs, cfg$diagnostic_balance)
    results <- results %>% left_join(diag_balance, by = "run_id")
  }

  satisfactory_results <- results %>% filter(satisfactory_run)

  # Plots
  nse_col <- paste0("NSE", suffix)
  kge_col <- paste0("KGE", suffix)

  g_scatter <- plot_performance_scatter(
    metrics, nse_col, kge_col,
    cfg$cal$threshold_nse, cfg$cal$threshold_kge,
    paste(ver$tag, "-", cfg$temporal_scale, "performance space"))
  print(g_scatter)
  save_tiff(g_scatter, paste0("performance_scatter_", ver$tag, ".tif"), ver$path)

  # Hydrographs — primary (calibration target) + secondary (diagnostic).
  # Each saved with scale suffix (_day / _mon) so both files coexist.
  primary_suffix <- if (cfg$temporal_scale == "daily") "day" else "mon"
  g_hydro <- plot_hydrograph(
    comparativo     = comp_primary,
    behavioral_runs = if (nrow(satisfactory_results) > 0) satisfactory_results$run_id else NULL,
    title_label     = paste(ver$tag, "—", cfg$temporal_scale, "streamflow")
  )
  print(g_hydro)
  save_tiff(g_hydro, paste0("hydrograph_", ver$tag, "_", primary_suffix, ".tif"),
            ver$path, width = 20, height = 10)

  # Secondary scale hydrograph (only if comp_sec was built)
  if (exists("comp_sec") && !is.null(comp_sec)) {
    sec_suffix_tag <- if (secondary_key == "monthly") "mon" else "day"
    sec_label      <- if (secondary_key == "monthly") "monthly" else "daily"
    g_hydro_sec <- plot_hydrograph(
      comparativo     = comp_sec,
      behavioral_runs = if (nrow(satisfactory_results) > 0) satisfactory_results$run_id else NULL,
      title_label     = paste(ver$tag, "—", sec_label, "streamflow (diagnostic)")
    )
    print(g_hydro_sec)
    save_tiff(g_hydro_sec,
              paste0("hydrograph_", ver$tag, "_", sec_suffix_tag, ".tif"),
              ver$path, width = 20, height = 10)
  }

  # Parameter boxplot
  if (nrow(satisfactory_results) > 0) {
    g_box <- plot_param_boxplot(satisfactory_results, param_info,
                                 paste(ver$tag, "- accepted parameters"))
    print(g_box)
    save_tiff(g_box, paste0("params_boxplot_", ver$tag, ".tif"),
              ver$path, width = 25, height = 10)
  }

  # Acceptance-domain plot for the streamflow iterations (ITER2, the daily iterations).
  # Same axes as Iteration 1, physical quantities on both: median BFI across
  # target subbasins (X) and Q_total/P (Y). Gated on bfi_targets being present.
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
      title_label        = paste(ver$tag, "- acceptance domain")
    )
    print(g_balance)
    save_tiff(g_balance, paste0("balance_scatter_", ver$tag, ".tif"), ver$path)
  }

  # Morris ranking for narrowing (used by unified inheritance block at the end)
  ranking_morris <- morris_results$ranking %>%
    filter(if ("metric" %in% names(.)) metric == "obj_total" else TRUE)
}

# ==============================================================================
# 12. SAVE RESULTS — standardised filenames `<descriptor>_<TAG>.<ext>`
# ==============================================================================
iter_label <- cfg$iteration

# ---- (a) Full RDS (run-level state for downstream tooling) ----
saveRDS(
  list(
    script_version       = cfg$script_version,
    config               = cfg,
    tag                  = ver$tag,
    morris_source           = morris_results$folder,
    sim                  = sim,
    param_info           = param_info,
    ranking_morris       = ranking_morris,
    results              = results,
    satisfactory_results = satisfactory_results
  ),
  file.path(ver$path, paste0("simulation_results_", ver$tag, ".rds"))
)

# ---- (b) Run-level CSVs: ALL valid + SATISFACTORIES (always when non-empty) -
write.csv(results,
          file.path(ver$path, paste0("simulation_results_", ver$tag, ".csv")),
          row.names = FALSE)

if (nrow(satisfactory_results) > 0) {
  write.csv(satisfactory_results,
            file.path(ver$path,
                      paste0("simulation_results_satisfactories_", ver$tag, ".csv")),
            row.names = FALSE)
}

# ---- (c) Baseflow diagnostic CSV — Eckhardt sim vs SWAT+ internal partition --
# Joined with BFI eckhardt sim for direct per-gauge comparison. Includes a
# `gap_internal_eckhardt_chXX = bfi_swat_internal − bfi_eckhardt_sim` column
# per channel. Only written when extra_outputs.aquifer.enabled = true.
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
  cat("[baseflow] saved baseflow_comparison_", ver$tag, ".csv (",
      nrow(bf_compare), " runs x ", nrow(bfi_targets), " gauges)\n", sep = "")
}

# ==============================================================================
# 13. RUN INFO — written BEFORE the inheritance step so it always lands on disk
#     even when the inheritance step stops loud (0 satisfactory + "satisfactory"
#     mode). The CAL itself completed; only the post-CAL narrowing is missing.
# ==============================================================================
inh_source <- cfg$cal$inheritance_source %||% "satisfactory"

# Render auto_drop info for the run_info summary (concise, one line)
auto_drop_str <- if (length(auto_drop_info$drop) > 0) {
  paste(sprintf("%s(mu*=%g)", auto_drop_info$detail$parameter,
                                auto_drop_info$detail$mu_star), collapse = ", ")
} else "none"

write_run_info(
  path       = ver$path,
  cfg        = cfg,
  run_type   = "CAL",
  param_info = param_info,
  summary    = list(
    tag                  = ver$tag,
    cal_type             = cfg$cal$type %||% "(unset)",
    inheritance_source   = inh_source,
    auto_drop_threshold  = cfg$inherit$auto_drop_threshold %||% "(disabled)",
    auto_dropped         = auto_drop_str,
    n_simulations        = cfg$cal$n_simulations %||% nrow(results),
    n_valid              = length(valid_runs),
    n_satisfactory       = nrow(satisfactory_results),
    morris_source           = morris_results$folder
  )
)

# ==============================================================================
# 14. INHERITED RANGES — driven by cal.inheritance_source ("satisfactory" | "all_valid")
# ==============================================================================
if (inh_source == "satisfactory") {
  if (nrow(satisfactory_results) > 0) {
    inheritance_input <- satisfactory_results
    cat("[inheritance] source = satisfactory (",
        nrow(satisfactory_results), " runs)\n", sep = "")
  } else {
    stop(sprintf(
      paste(
        "%s produced 0 satisfactory runs (n_valid=%d, n_sim=%d).",
        "Cannot generate inherited_ranges_%s.csv.",
        "Action: widen prior, increase n_simulations, relax thresholds, OR set",
        "  cal.inheritance_source: \"all_valid\"  in the YAML to use all valid runs",
        "  (exploratory mode -- ranges will NOT reflect acceptance filtering).",
        "(run_info.txt was already written; only the inheritance file is missing.)",
        sep = "\n"
      ),
      ver$tag, length(valid_runs), n_sim, ver$tag
    ))
  }
} else if (inh_source == "all_valid") {
  if (nrow(results) == 0) {
    stop(ver$tag, " produced 0 valid runs -- cannot proceed even in all_valid mode.")
  }
  inheritance_input <- results
  warning("[inheritance] source = all_valid (", nrow(results), " runs). ",
          "Inherited ranges are NOT filtered on acceptance ",
          "(only ", nrow(satisfactory_results), " runs would have passed).")
} else {
  stop("Invalid cal.inheritance_source: '", inh_source,
       "'. Use 'satisfactory' or 'all_valid'.")
}

new_ranges <- calc_new_ranges(
  inheritance_input,
  param_info$par_name,
  param_info,
  ranking_morris
)
print(new_ranges)

write.csv(new_ranges,
          file.path(ver$path, paste0("inherited_ranges_", ver$tag, ".csv")),
          row.names = FALSE)
cat("[OK] Wrote inherited_ranges_", ver$tag, ".csv (",
    nrow(new_ranges), " parameters).\n", sep = "")

# ==============================================================================
# 15. UPDATE RUN INFO with final inheritance result
# ==============================================================================
write_run_info(
  path       = ver$path,
  cfg        = cfg,
  run_type   = "CAL",
  param_info = param_info,
  summary    = list(
    tag                  = ver$tag,
    cal_type             = cfg$cal$type %||% "(unset)",
    inheritance_source   = inh_source,
    auto_drop_threshold  = cfg$inherit$auto_drop_threshold %||% "(disabled)",
    auto_dropped         = auto_drop_str,
    n_simulations        = cfg$cal$n_simulations %||% nrow(results),
    n_valid              = length(valid_runs),
    n_satisfactory       = nrow(satisfactory_results),
    n_used_inheritance   = nrow(inheritance_input),
    morris_source           = morris_results$folder
  )
)

cat("\nFiles saved to:", ver$path, "\n")
cat("CAL", ver$tag, "COMPLETED.\n")
