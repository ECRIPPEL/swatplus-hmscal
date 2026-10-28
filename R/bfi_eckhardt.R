# ==============================================================================
# bfi_eckhardt.R — Eckhardt (2005) digital baseflow filter + multi-site BFI
#
# Implements the HMSCal BFI multisite constraint:
#
#   For each selected gauge i:
#     bf[t] = ( (1 - BFImax_i) * a_i * bf[t-1]
#             + (1 - a_i) * BFImax_i * Q[t] )  /  (1 - a_i * BFImax_i)
#     bf[t] := clip(bf[t], 0, Q[t])               # physical constraints
#     BFI_sim_i = sum(bf) / sum(Q)               # over the calibration window
#
# WARM-UP DECISION (2026-05-04 with user):
#   No filter warm-up (`bf` is computed over the full window 1980-1989).
#   Justification: SWAT+ already includes years_skip warm-up before 1980, so the
#   simulated Q is in steady state by t=0 of the filter. The Eckhardt forward
#   recursion converges in ~30-60 days; over 10 years the transient contribution
#   to the integrated BFI is <1%. Discarding 365 days would just create
#   asymmetry with BFI_obs (which is computed over 10 full years from ANA).
#
# AGGREGATION:
#   r[i]            = |BFI_obs_i - BFI_sim_i|
#   obj_bfi_mdae    = median(r[1..N])              (continuous gradient for optimizer)
#   obj_bfi_mae     = mean(r[1..N])                (diagnostic; catches bimodal failure)
#   n_passed_eps05  = count of r[i] <= eps         (default eps = 0.05; reporting only)
#
# CHANNEL INDEXING:
#   Channels are looked up BY ID via the bfi_targets CSV — never by column order
#   in the sim matrix. The CSV is the single source of truth.
# ==============================================================================

library(tibble)
library(dplyr)

# Forward Eckhardt filter — physical constraints baked in
eckhardt_forward <- function(Q, a, BFImax) {
  n <- length(Q)
  if (n < 1) return(numeric(0))
  bf <- numeric(n)
  bf[1] <- Q[1] * BFImax
  if (is.na(bf[1])) bf[1] <- 0
  if (bf[1] > Q[1] && !is.na(Q[1])) bf[1] <- Q[1]
  if (bf[1] < 0)   bf[1] <- 0

  denom <- 1 - a * BFImax
  if (denom <= 0) {
    stop(sprintf("Eckhardt filter: 1 - a*BFImax = %g must be > 0 (a=%g, BFImax=%g)",
                 denom, a, BFImax))
  }

  for (t in 2:n) {
    if (is.na(Q[t])) {
      bf[t] <- 0
      next
    }
    bf[t] <- ((1 - BFImax) * a * bf[t - 1] + (1 - a) * BFImax * Q[t]) / denom
    if (is.na(bf[t]) || bf[t] < 0)   bf[t] <- 0
    if (bf[t] > Q[t])                bf[t] <- Q[t]
  }
  bf
}

# BFI for a single run on a single channel
# Q_daily: numeric vector aligned with `dates`
# min_days: minimum number of valid (non-NA, non-negative) days required to
#           compute BFI; below this returns NA. Production default 365 (full
#           year) per spec; smoke runs may set lower (e.g., 90).
# Returns scalar BFI (or NA if fewer than min_days valid days)
compute_bfi_one <- function(dates, Q_daily, a, BFImax,
                            period_start, period_end,
                            min_days = 365) {
  ok <- !is.na(dates) & dates >= as.Date(period_start) & dates <= as.Date(period_end)
  Q  <- Q_daily[ok]

  na_mask <- is.na(Q) | Q < 0
  if (sum(!na_mask) < min_days) return(NA_real_)

  bf <- eckhardt_forward(Q, a, BFImax)
  # Exclude days where Q is NA (or negative) from BOTH sums
  num <- sum(bf[!na_mask], na.rm = TRUE)
  den <- sum(Q[!na_mask],  na.rm = TRUE)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

# Read the bfi_targets CSV (looked up under cfg$project_path)
# Required columns: channel_id, bfi (= BFI_obs), a, bfimax_q9050
read_bfi_targets <- function(csv_path) {
  if (!file.exists(csv_path)) {
    stop("bfi_targets CSV not found: ", csv_path)
  }
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  required <- c("channel_id", "bfi", "a", "bfimax_q9050")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("bfi_targets CSV missing columns: ", paste(missing, collapse = ", "),
         "\n  File: ", csv_path)
  }
  df %>%
    transmute(
      channel_id = as.integer(channel_id),
      bfi_obs    = as.numeric(bfi),
      a          = as.numeric(a),
      bfimax     = as.numeric(bfimax_q9050),
      gauge_code = if ("gauge_code" %in% names(df)) as.character(gauge_code) else NA_character_
    )
}

# Compute BFI multi-site for ALL valid_runs.
#
# Two source options for the simulated BFI per channel:
#
#   sim_source = "swat_internal_bf"  (default, the acceptance filter)
#     BFI_sim_i = (Σ aqu_flo_cha × area over upstream FLOODPLAIN aquifers)
#                  / (Q volume at gauge channel i)
#     Tracks the SWAT+ internal aquifer-to-channel partition. Forces
#     calibration to actually activate the aquifer module rather than
#     mimicking baseflow recession via slow latq. Requires project_path
#     and extra_outputs.aquifer.enabled (so aqu_flo_cha_<i> outputs exist).
#
#   sim_source = "eckhardt"  (legacy / diagnostic only)
#     BFI_sim_i = sum(bf) / sum(Q) where bf is the Eckhardt forward filter
#     applied to simulated daily Q at gauge i, with the same (a, BFImax)
#     used to derive BFI_obs. Sees only hydrograph SHAPE, not internal
#     source — fooled by long-recession latq.
#
# sim          : SWATrunR simulation object
# output_name  : YAML key of the multi-channel daily output (e.g.,"vazao_multi_day")
# bfi_targets  : tibble from read_bfi_targets (channel_id, bfi_obs, a, bfimax)
# valid_runs   : character vector of run_NNN names
# sim_source   : "swat_internal_bf" (default) | "eckhardt"
# project_path : SWAT+ project root (required for "swat_internal_bf")
# period_start : "YYYY-MM-DD" — Eckhardt window only (ignored by SWAT internal)
# period_end   : "YYYY-MM-DD" — Eckhardt window only
# eps          : per-channel pass tolerance for n_passed_eps reporting
# min_days     : Eckhardt minimum valid days per channel
#
# Returns a tibble with columns:
#   run_id, bfi_sim_ch<id> ..., bfi_dev_ch<id> ..., obj_bfi_mdae, obj_bfi_mae,
#   n_passed_eps05, n_evaluated_ch
compute_bfi_multisite <- function(sim, output_name, bfi_targets, valid_runs,
                                  sim_source   = "swat_internal_bf",
                                  project_path = NULL,
                                  period_start = "1980-01-01",
                                  period_end   = "1989-12-31",
                                  eps          = 0.05,
                                  min_days     = 365) {

  if (!sim_source %in% c("eckhardt", "swat_internal_bf")) {
    stop("compute_bfi_multisite: unknown sim_source = '", sim_source,
         "'. Must be 'eckhardt' or 'swat_internal_bf'.")
  }

  bfi_obs_vec <- bfi_targets$bfi_obs
  names(bfi_obs_vec) <- as.character(bfi_targets$channel_id)

  bfi_mat <- matrix(NA_real_, nrow = length(valid_runs),
                    ncol = nrow(bfi_targets),
                    dimnames = list(valid_runs,
                                    as.character(bfi_targets$channel_id)))

  # ---- Source A: Eckhardt forward filter applied to simulated daily Q ------
  if (sim_source == "eckhardt") {
    per_channel <- list()
    missing_ch  <- integer(0)
    for (cid in bfi_targets$channel_id) {
      key <- sprintf("%s_%d", output_name, cid)
      df  <- sim$simulation[[key]]
      if (is.null(df)) { missing_ch <- c(missing_ch, cid); next }
      per_channel[[as.character(cid)]] <- df
    }
    if (length(missing_ch) > 0) {
      stop("compute_bfi_multisite (eckhardt): SWAT+ output not found for ch(s) ",
           paste(missing_ch, collapse = ", "),
           ".\n  Expected sim$simulation$", output_name, "_<channel_id>.")
    }

    ref_df   <- per_channel[[1]]
    date_col <- intersect(c("date", "Date"), names(ref_df))[1]
    if (is.na(date_col)) stop("eckhardt: no date column in per-channel tibbles.")
    dates <- as.Date(ref_df[[date_col]])

    for (j in seq_len(nrow(bfi_targets))) {
      cid <- bfi_targets$channel_id[j]
      a   <- bfi_targets$a[j]
      bfx <- bfi_targets$bfimax[j]
      df  <- per_channel[[as.character(cid)]]
      for (run_nm in valid_runs) {
        Q <- df[[run_nm]]
        if (is.null(Q)) next
        bfi_mat[run_nm, as.character(cid)] <-
          compute_bfi_one(dates, Q, a, bfx, period_start, period_end, min_days)
      }
    }
  }

  # ---- Source B: SWAT+ internal aquifer-to-channel partition ---------------
  # Delegates to compute_baseflow_comparison (defined in baseflow_comparison.R)
  # which uses the topology in aquifer.con + chandeg.con to aggregate
  # aqu_flo_cha across upstream FLOODPLAIN aquifers, then divides by total
  # Q volume at the gauge.
  if (sim_source == "swat_internal_bf") {
    if (is.null(project_path)) {
      stop("compute_bfi_multisite: sim_source='swat_internal_bf' requires ",
           "project_path argument (where aquifer.con + chandeg.con live).")
    }
    bf_int <- compute_baseflow_comparison(
      sim              = sim,
      project_path     = project_path,
      bfi_targets      = bfi_targets,
      valid_runs       = valid_runs,
      multi_day_output = output_name
    )
    # Pull bfi_swat_internal_chXX columns into the matrix
    for (cid in bfi_targets$channel_id) {
      col <- paste0("bfi_swat_internal_ch", cid)
      if (col %in% names(bf_int)) {
        bfi_mat[, as.character(cid)] <- bf_int[[col]][match(valid_runs, bf_int$run_id)]
      }
    }
  }

  # ---- Common: residuals, MdAE, MAE, n_passed (same for both sources) ------
  res <- tibble(run_id = valid_runs)
  for (cid in as.character(bfi_targets$channel_id)) {
    res[[paste0("bfi_sim_ch", cid)]] <- as.numeric(bfi_mat[, cid])
  }
  for (cid in as.character(bfi_targets$channel_id)) {
    res[[paste0("bfi_dev_ch", cid)]] <-
      abs(as.numeric(bfi_mat[, cid]) - bfi_obs_vec[cid])
  }

  dev_cols <- grep("^bfi_dev_ch", names(res), value = TRUE)
  dev_mat  <- as.matrix(res[, dev_cols])
  res$obj_bfi_mdae   <- apply(dev_mat, 1, function(r) median(r, na.rm = TRUE))
  res$obj_bfi_mae    <- apply(dev_mat, 1, function(r) mean(r,   na.rm = TRUE))
  res$n_passed_eps05 <- apply(dev_mat, 1, function(r) sum(r <= eps, na.rm = TRUE))
  res$n_evaluated_ch <- apply(dev_mat, 1, function(r) sum(!is.na(r)))
  res$bfi_sim_source <- sim_source

  res
}
