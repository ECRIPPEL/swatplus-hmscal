# ==============================================================================
# metrics.R — Performance metrics and water balance computation
# ==============================================================================

library(dplyr)
library(tibble)
library(purrr)
library(hydroGOF)

#' Compute NSE, KGE, PBIAS, and obj_total for each run
#'
#' @param comparativo data.frame with Date, Flow, and run_* columns
#' @param valid_runs Character vector of run names
#' @param suffix Metric column suffix (e.g., "_mon", "_day")
#' @return Tibble with metrics per run
calc_hard_metrics <- function(comparativo, valid_runs, suffix = "_mon") {
  map_dfr(valid_runs, function(run_nm) {
    sim_vals <- comparativo[[run_nm]]
    obs_vals <- comparativo$Flow

    if (all(is.na(sim_vals)) || length(na.omit(sim_vals)) < 2 ||
        sd(sim_vals, na.rm = TRUE) == 0) {
      out <- tibble(
        run_id    = run_nm,
        NSE       = NA_real_,
        KGE       = NA_real_,
        PBIAS     = NA_real_,
        absPBIAS  = NA_real_,
        err_nse   = NA_real_,
        err_kge   = NA_real_,
        err_pbias = NA_real_,
        obj_total = NA_real_,
        run_index = as.numeric(sub("^run_0*([0-9]+).*$", "\\1", run_nm))
      )
      names(out) <- gsub("^(NSE|KGE|PBIAS)$", paste0("\\1", suffix), names(out))
      return(out)
    }

    nse_val   <- NSE(sim = sim_vals, obs = obs_vals)
    kge_val   <- KGE(sim = sim_vals, obs = obs_vals)
    pbias_val <- pbias(sim = sim_vals, obs = obs_vals)

    nse_n   <- as.numeric(nse_val)
    kge_n   <- as.numeric(kge_val)
    pbias_n <- as.numeric(pbias_val)

    out <- tibble(
      run_id    = run_nm,
      NSE       = nse_n,
      KGE       = kge_n,
      PBIAS     = pbias_n,
      absPBIAS  = abs(pbias_n),
      err_nse   = pmax(0, 1 - pmin(nse_n, 1)),
      err_kge   = pmax(0, 1 - pmin(kge_n, 1)),
      err_pbias = abs(pbias_n) / 100,
      obj_total = pmax(0, 1 - pmin(nse_n, 1)) + pmax(0, 1 - pmin(kge_n, 1)) + abs(pbias_n) / 100,
      run_index = as.numeric(sub("^run_0*([0-9]+).*$", "\\1", run_nm))
    )

    names(out) <- gsub("^(NSE|KGE|PBIAS)$", paste0("\\1", suffix), names(out))
    out
  })
}

#' Classify runs as satisfactory (streamflow thresholds)
#'
#' Optionally ANDs in obj_rc + obj_bfi_mdae thresholds when those columns
#' are present. If `threshold_obj_rc` / `threshold_obj_bfi` are NULL the v1
#' predicate (NSE/KGE/PBIAS only) is used.
#'
#' @param metrics Tibble of metrics (output of calc_hard_metrics, optionally
#'                joined with obj_rc and obj_bfi_mdae columns)
#' @param threshold_nse Minimum NSE threshold
#' @param threshold_kge Minimum KGE threshold
#' @param threshold_abs_pbias Maximum |PBIAS| threshold
#' @param threshold_obj_rc Optional max for obj_rc (= |RC_obs - RC_sim|)
#' @param threshold_obj_bfi Optional max for obj_bfi_mdae (median |BFI_obs - BFI_sim|)
#' @param suffix Column suffix (e.g., "_mon", "_day")
#' @return Tibble with added column satisfactory_run
classify_hard_runs <- function(metrics, threshold_nse, threshold_kge,
                               threshold_abs_pbias,
                               threshold_obj_rc  = NULL,
                               threshold_obj_bfi = NULL,
                               suffix = "_mon") {
  nse_col <- paste0("NSE", suffix)
  kge_col <- paste0("KGE", suffix)

  base_pass <- metrics[[nse_col]] >= threshold_nse &
               metrics[[kge_col]] >= threshold_kge &
               metrics$absPBIAS  <= threshold_abs_pbias

  rc_pass <- if (!is.null(threshold_obj_rc) && "obj_rc" %in% names(metrics)) {
    metrics$obj_rc <= threshold_obj_rc
  } else {
    rep(TRUE, nrow(metrics))
  }

  bfi_pass <- if (!is.null(threshold_obj_bfi) && "obj_bfi_mdae" %in% names(metrics)) {
    metrics$obj_bfi_mdae <= threshold_obj_bfi
  } else {
    rep(TRUE, nrow(metrics))
  }

  metrics %>%
    mutate(satisfactory_run = base_pass & rc_pass & bfi_pass) %>%
    arrange(desc(satisfactory_run), desc(.data[[nse_col]]),
            desc(.data[[kge_col]]), absPBIAS)
}

#' Compute RC (runoff coefficient) and obj_rc per run
#'
#' RC = (channel-routed Q at outlet) / P_basin, dimensionless.
#'
#' Reads `flo_out` (m³/s, annual avg) at the outlet directly from the
#' SWAT+ `channel_sd_aa.txt` output, captured by SWATrunR as
#' `sim$simulation$cha_aa_flo_out_<outlet_channel>` when
#' `extra_outputs.channels.enabled = TRUE` and `outlet_channel` is in
#' the channels list. The conversion to mm/yr uses the upstream
#' contributing area from `aquifer.con` topology:
#'
#'   Q_mm_yr = flo_out_m3s × 86400 × 365 / (A_up_ha × 10)
#'   RC_sim  = Q_mm_yr / P_mm
#'
#' This matches what the ANA gauge measures — accounts for channel +
#' reservoir transmission losses (evap, seep, ΔS) that are subtracted
#' between HRU output and channel outlet (~20 % of HRU output in the
#' Uruguay project). Calibrating against HRU-level (wateryld + baseflow)
#' would silently accept partitions that fail the actual gauged Q.
#'
#' WYLD_sim and BASEFLOW_sim are reported as diagnostic columns
#' (HRU-level partition reference) but do NOT enter rc_sim.
#'
#' Falls back to HRU-level sum with a loud warning if the outlet
#' channel_sd_aa output is missing (e.g., extra_outputs.channels disabled
#' or the outlet not in the channels list).
#'
#' @param sim            SWATrunR simulation object
#' @param valid_runs     Character vector of run names
#' @param rc_target      Numeric — observed Q_outlet / P_basin
#' @param project_path   SWAT+ project path (for topology / area lookup)
#' @param outlet_channel Outlet channel id (e.g., 33 for Iraí)
#' @return Tibble with run_id, P_sim, ET_sim, WYLD_sim, BASEFLOW_sim,
#'         Q_TOTAL_sim, rc_sim, rto_et, rto_wyld, rto_baseflow,
#'         rc_hru_level (diagnostic), channel_losses (diagnostic), obj_rc
calc_obj_rc <- function(sim, valid_runs, rc_target,
                        project_path   = NULL,
                        outlet_channel = NULL) {

  has_wateryld <- !is.null(sim$simulation$wateryld)
  has_baseflow <- !is.null(sim$simulation$baseflow)
  if (!has_wateryld) {
    warning("calc_obj_rc: 'wateryld' output missing; HRU-level diagnostic ",
            "(rto_wyld, rc_hru_level, channel_losses) will read NA.")
  }
  if (!has_baseflow) {
    warning("calc_obj_rc: 'baseflow' output missing; HRU-level diagnostic ",
            "(rto_baseflow) will read zero.")
  }

  # ---- Read flo_out (m³/s, annual avg) at outlet from channel_sd_aa ----
  outlet_key <- if (!is.null(outlet_channel)) {
    paste0("cha_aa_flo_out_", outlet_channel)
  } else NA_character_

  has_outlet <- !is.na(outlet_key) && !is.null(sim$simulation[[outlet_key]])

  basin_area_ha <- NA_real_
  if (has_outlet && !is.null(project_path)) {
    # build_aqu_upstream_map() lives in R/baseflow_comparison.R
    topo <- tryCatch(
      build_aqu_upstream_map(project_path,
                              gauge_channels = outlet_channel,
                              cache = TRUE),
      error = function(e) { warning(conditionMessage(e)); NULL })
    if (!is.null(topo)) {
      ups <- topo$upstream_aqus_all[[as.character(outlet_channel)]]
      if (length(ups) > 0) {
        basin_area_ha <- sum(topo$aqu_area_ha[as.character(ups)])
      }
    }
  }

  if (!has_outlet || is.na(basin_area_ha)) {
    warning("calc_obj_rc: cannot derive channel-routed Q at outlet (need ",
            "extra_outputs.channels with the outlet in units AND project_path ",
            "with aquifer.con topology). Falling back to HRU-level RC ",
            "(wateryld + baseflow), which OVERESTIMATES RC by ~20% due to ",
            "channel/reservoir losses.")
  }

  base <- tibble(
    run_id   = valid_runs,
    P_sim    = as.numeric(sim$simulation$precip[nrow(sim$simulation$precip), valid_runs]),
    ET_sim   = as.numeric(sim$simulation$et[nrow(sim$simulation$et), valid_runs]),
    WYLD_sim = if (has_wateryld) {
      as.numeric(sim$simulation$wateryld[nrow(sim$simulation$wateryld), valid_runs])
    } else {
      rep(NA_real_, length(valid_runs))
    },
    BASEFLOW_sim = if (has_baseflow) {
      as.numeric(sim$simulation$baseflow[nrow(sim$simulation$baseflow), valid_runs])
    } else {
      rep(0, length(valid_runs))
    }
  )

  if (has_outlet && !is.na(basin_area_ha)) {
    # channel_sd_aa flo_out is m³/s annual avg (single row per run)
    df_q <- sim$simulation[[outlet_key]]
    flo_out_m3s <- as.numeric(df_q[nrow(df_q), valid_runs])
    # m³/s × 86400 × 365 = m³/yr; ÷ (ha × 10) = mm/yr
    base$Q_TOTAL_sim <- flo_out_m3s * 86400 * 365 / (basin_area_ha * 10)
  } else {
    # Fallback: HRU-level sum (overestimates by channel-loss amount)
    base$Q_TOTAL_sim <- base$WYLD_sim + base$BASEFLOW_sim
  }

  base %>%
    mutate(
      rto_q        = Q_TOTAL_sim  / P_sim,         # = RC_sim (channel outlet)
      rc_sim       = rto_q,
      rto_et       = ET_sim       / P_sim,         # diagnostic — ET is escape valve
      rto_wyld     = WYLD_sim     / P_sim,         # diagnostic: HRU surq+latq only
      rto_baseflow = BASEFLOW_sim / P_sim,         # diagnostic: HRU gw_q only
      # HRU-level RC for closure-comparison diagnostic (will be > rc_sim due
      # to channel losses; the gap quantifies network transmission losses).
      rc_hru_level   = (WYLD_sim + BASEFLOW_sim) / P_sim,
      channel_losses = rc_hru_level - rc_sim,       # gap = channel evap + seep + ΔS
      rc_target      = rc_target,
      obj_rc         = abs(rto_q - rc_target)
    )
}

#' Compute Iteration 1 ranking (RC + BFI multisite predicate)
#'
#' Iteration 1 architecture:
#'   satisfactory_run = (obj_rc      <= threshold_obj_rc) AND
#'                      (obj_bfi_mdae <= threshold_obj_bfi)
#'
#' RC  = Q_total / P                    (single deviation, no ET clause)
#' BFI = multi-site median absolute error from Eckhardt baseflow filter
#'       applied with per-gauge (a, BFImax) constants from the BFI targets CSV
#'
#' ET is reported as a diagnostic column (rto_et) but does NOT participate in
#' the satisfactory predicate — it is the model's escape valve under the
#' closure premise (rchg_dp = 0, no deep aquifer sink).
#'
#' @param sim SWATrunR simulation object
#' @param valid_runs Character vector of run names
#' @param rc_target Numeric scalar — observed Q_total / P
#' @param threshold_obj_rc  Max acceptable |RC_obs - RC_sim|
#' @param threshold_obj_bfi Max acceptable median(|BFI_obs - BFI_sim|) across N channels
#' @param bfi_metrics Tibble from compute_bfi_multisite() (run_id + per-channel
#'                    + obj_bfi_mdae + obj_bfi_mae + n_passed_eps05)
#' @return Tibble with full diagnostic columns + satisfactory_run flag
calc_soft_metrics <- function(sim, valid_runs, rc_target,
                              threshold_obj_rc, threshold_obj_bfi,
                              bfi_metrics,
                              project_path   = NULL,
                              outlet_channel = NULL) {
  rc_df <- calc_obj_rc(sim, valid_runs, rc_target,
                       project_path   = project_path,
                       outlet_channel = outlet_channel)

  ranking <- rc_df %>%
    left_join(bfi_metrics, by = "run_id") %>%
    mutate(
      obj_total_v2 = obj_rc + obj_bfi_mdae,
      satisfactory_run =
        obj_rc       <= threshold_obj_rc &
        obj_bfi_mdae <= threshold_obj_bfi,
      run_index = as.numeric(sub("^run_0*([0-9]+).*", "\\1", run_id))
    ) %>%
    arrange(desc(satisfactory_run), obj_total_v2)

  ranking
}

#' Compute Iteration 1 screening metrics (RC + multisite BFI proxy)
#'
#' Formulation:
#'   obj_rc        = |Q_total/P − RC_obs|
#'   obj_bfi_mdae  = median over N channels of |BFI_sim_i − BFI_obs_i|
#'   obj_total     = w_RC * obj_rc + w_BFI * obj_bfi_mdae   (defaults: 1.0 / 1.0)
#'
#' ET is no longer in obj_total — under the closure premise (rchg_dp = 0),
#' ET/P is mathematically derived from Q/P, so adding both was double counting.
#' ET is captured as a diagnostic column (et_p) for post-hoc inspection.
#'
#' BFI is computed via Eckhardt forward filter at each of the N gauges in
#' bfi_metrics; see R/bfi_eckhardt.R for the implementation.
#'
#' @param sim SWATrunR simulation object
#' @param valid_runs Character vector of run names
#' @param rc_target Numeric scalar — observed Q_total / P
#' @param bfi_metrics Tibble from compute_bfi_multisite() (run_id + per-channel
#'                    + obj_bfi_mdae + obj_bfi_mae + n_passed_eps05)
#' @param w_rc Weight on obj_rc in obj_total (default 1.0)
#' @param w_bfi Weight on obj_bfi_mdae in obj_total (default 1.0)
#' @return Tibble with metrics per run
calc_iteration1_metrics <- function(sim, valid_runs, rc_target, bfi_metrics,
                                    w_rc = 1.0, w_bfi = 1.0,
                                    project_path   = NULL,
                                    outlet_channel = NULL) {
  # use channel-routed RC (calc_obj_rc) so Morris EE ranks parameters
  # against the SAME RC metric the calibration filter uses (channel outlet,
  # not HRU sum). Falls back to HRU-level if outlet output unavailable.
  rc_df <- calc_obj_rc(sim, valid_runs, rc_target,
                       project_path   = project_path,
                       outlet_channel = outlet_channel)

  # Also keep the HRU-component breakdown (surq + latq + gw_q at basin level)
  # for diagnostic / sensitivity reporting. These come from basin_wb_aa
  # (surq, latq) and basin_aqu_aa (gw_q via the `gw_q` cal_outputs key
  # which maps to basin_aqu_aa.flo in the YAML).
  hru_components <- tibble(
    run_id = valid_runs,
    surq   = as.numeric(sim$simulation$surq[nrow(sim$simulation$surq), valid_runs]),
    latq   = as.numeric(sim$simulation$latq[nrow(sim$simulation$latq), valid_runs]),
    gw_q   = as.numeric(sim$simulation$gw_q[nrow(sim$simulation$gw_q), valid_runs])
  )

  rc_df %>%
    left_join(hru_components, by = "run_id") %>%
    rename(precip = P_sim, et = ET_sim) %>%
    mutate(
      q_total = surq + latq + gw_q,             # HRU-level diagnostic
      q_p     = rc_sim,                          # channel-routed (= obj_rc denom)
      et_p    = rto_et                           # diagnostic — ET is free
    ) %>%
    left_join(bfi_metrics %>% select(run_id, obj_bfi_mdae, obj_bfi_mae,
                                      n_passed_eps05),
              by = "run_id") %>%
    mutate(
      obj_total = w_rc * obj_rc + w_bfi * obj_bfi_mdae,
      run_index = as.numeric(sub("^run_0*([0-9]+).*", "\\1", run_id))
    )
}

#' Compute water balance diagnostic (daily iterations)
#'
#' DIAGNOSTIC-ONLY for the streamflow iterations. The WB closure quantities
#' (rto_et, rto_q, rto_wyld, rto_baseflow) are reported as columns on each
#' calibration run, but they are NEVER combined with the daily NSE/KGE/PBIAS
#' filter — the multi-scale premise of the framework is that WB closure is
#' the concern of ITER1 (annual), and ITER2 and the daily iterations calibrate streamflow only.
#' These columns let you verify post-hoc whether the streamflow-calibrated
#' ensemble preserves the WB closure of ITER1.
#'
#' Column naming matches calc_soft_metrics() so the same plot_balance_scatter()
#' works for both ITER1 and the daily-iteration diagnostic outputs.
#'
#' Uses the FULL simulated streamflow Q_total = surq+latq+gw_q (i.e.,
#' basin_wb_aa.wateryld + basin_aqu_aa.flo_cha). Diagnostic columns rto_wyld
#' (surq+latq only) and rto_baseflow (gw_q only) are retained so the
#' surface/subsurface partition remains visible. Falls back to the legacy
#' wateryld-only diagnostic with a warning if the baseflow output is missing.
#'
#' @param sim SWATrunR simulation object
#' @param valid_runs Character vector of run names
#' @param targets List with et_rto, wyld_rto (target Q_obs/P)
#' @return Tibble with rto_et, rto_q, rto_wyld, rto_baseflow + error columns
calc_balance_diagnostic <- function(sim, valid_runs, targets) {
  has_baseflow <- !is.null(sim$simulation$baseflow)
  if (!has_baseflow) {
    warning("calc_balance_diagnostic: 'baseflow' output not found; falling back to ",
            "wateryld-only diagnostic (legacy). Add baseflow to cal_outputs in YAML ",
            "to enable the corrected Q_total = surq+latq+gw_q formulation.")
  }

  tibble(
    run_id   = valid_runs,
    P_sim    = as.numeric(sim$simulation$precip[nrow(sim$simulation$precip), valid_runs]),
    ET_sim   = as.numeric(sim$simulation$et[nrow(sim$simulation$et), valid_runs]),
    WYLD_sim = as.numeric(sim$simulation$wateryld[nrow(sim$simulation$wateryld), valid_runs]),
    BASEFLOW_sim = if (has_baseflow) {
      as.numeric(sim$simulation$baseflow[nrow(sim$simulation$baseflow), valid_runs])
    } else {
      rep(0, length(valid_runs))
    }
  ) %>%
    mutate(
      Q_TOTAL_sim   = WYLD_sim + BASEFLOW_sim,
      rto_et        = ET_sim       / P_sim,
      rto_wyld      = WYLD_sim     / P_sim,        # diagnostic: surq+latq only
      rto_baseflow  = BASEFLOW_sim / P_sim,        # diagnostic: gw_q
      rto_q         = Q_TOTAL_sim  / P_sim,        # full streamflow ratio
      target_et     = targets$et_rto,
      target_q      = targets$wyld_rto,            # YAML key kept for compat
      err_et_abs    = abs(rto_et   - target_et),
      err_q_abs     = abs(rto_q    - target_q),
      err_wyld_abs  = abs(rto_wyld - target_q),
      err_et_pct    = 100 * err_et_abs   / target_et,
      err_q_pct     = 100 * err_q_abs    / target_q,
      err_wyld_pct  = 100 * err_wyld_abs / target_q,
      agg_rel_error_pct = err_et_pct + err_q_pct
    )
}
