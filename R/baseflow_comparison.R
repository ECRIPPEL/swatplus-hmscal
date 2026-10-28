# ==============================================================================
# baseflow_comparison.R — Compare SWAT+ internal baseflow partition vs Eckhardt
#
# For each gauge channel in bfi_targets:
#   BFI_swat_internal_chid = total_aqu_flow_upstream / total_Q_at_channel
#   BFI_eckhardt_sim_chid  = (already computed in bfi_eckhardt.R)
#
# The two quantities answer different questions:
#   - BFI_eckhardt_sim : recession-shape BFI from the simulated channel
#                        hydrograph (what an external observer with only Q would
#                        infer). Same signal extracted from BFI_obs.
#   - BFI_swat_internal: actual fraction of channel volume sourced from
#                        SWAT+ aquifer module (sum of aqu_flo_cha × area_aqu
#                        across upstream FLOODPLAIN aquifers, divided by
#                        total Q volume at the gauge channel).
#
# SWAT+ aquifer topology (mirrors LSU floodplain/upland pairing):
#   - floodplain aquifer : drains directly to channel (`sdc`); flo_cha > 0
#   - upland aquifer     : drains to floodplain (`aqu`); flo_ls > 0, flo_cha = 0
#   - deep aquifer       : single regional reservoir (receives `rhg` from all);
#                          does not contribute to channel flow
# Only floodplain aquifers carry non-zero flo_cha. The upland → floodplain
# transfer (via flo_ls) is already accumulated in the floodplain's flo_cha,
# so summing flo_cha × area across upstream floodplain aquifers is the
# correct (no double-counting) channel baseflow estimate from SWAT+ internals.
#
# Wide gap between the two = SWAT+ may be reproducing the right hydrograph
# shape (Eckhardt match) via the wrong source partition (lateral interflow
# `latq` mimicking baseflow recession). Diagnostic for structural realism.
#
# REQUIRES:
#   - sim object with `aqu_flo_cha_<i>` outputs (extra_outputs.aquifer.enabled)
#   - sim object with `<bfi_multi_day_output>_<chid>` outputs (vazao_multi_day)
#   - SWAT+ project with aquifer.con + chandeg.con readable
#
# Topology cache: <project_path>/.aqu_upstream_map.rds (auto-rebuilt if stale)
# ==============================================================================

library(dplyr)
library(tibble)


# ------------------------------------------------------------------------------
# parse_swat_con — generic parser for *.con files (chandeg, rout_unit, aquifer)
# Returns list of named lists, each with id/name/gis_id/area/lat/lon/elev +
# conns = list of (obj_typ, obj_id, hyd_typ, frac) tuples.
# ------------------------------------------------------------------------------
parse_swat_con <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  lines <- readLines(path, warn = FALSE)
  if (length(lines) <= 2) stop("File too short (header-only?): ", path)
  body  <- lines[-c(1, 2)]

  out <- vector("list", length(body))
  for (k in seq_along(body)) {
    tok <- strsplit(trimws(body[k]), "\\s+")[[1]]
    if (length(tok) < 13) next
    out_tot <- as.integer(tok[13])
    conns   <- vector("list", out_tot)
    pos     <- 14
    for (j in seq_len(out_tot)) {
      if (pos + 3 > length(tok)) break
      conns[[j]] <- list(
        obj_typ = tok[pos],
        obj_id  = as.integer(tok[pos + 1]),
        hyd_typ = tok[pos + 2],
        frac    = as.numeric(tok[pos + 3])
      )
      pos <- pos + 4
    }
    out[[k]] <- list(
      id      = as.integer(tok[1]),
      name    = tok[2],
      gis_id  = as.integer(tok[3]),
      area    = as.numeric(tok[4]),
      lat     = as.numeric(tok[5]),
      lon     = as.numeric(tok[6]),
      elev    = as.numeric(tok[7]),
      out_tot = out_tot,
      conns   = conns
    )
  }
  out[!vapply(out, is.null, logical(1))]
}


# ------------------------------------------------------------------------------
# build_aqu_upstream_map — for each gauge channel, list upstream aquifers
#
# @param project_path SWAT+ TxtInOut path
# @param gauge_channels integer vector of channel ids
# @param cache TRUE to persist as .aqu_upstream_map.rds
# @return list:
#   $aqu_to_cha      : named int (aquifer_id → channel_id it drains to)
#   $aqu_area_ha     : named numeric (aquifer_id → area in ha)
#   $upstream_chans  : named list (channel_id → upstream channel ids incl. self)
#   $upstream_aqus   : named list (channel_id → aquifer ids upstream of it)
# ------------------------------------------------------------------------------
build_aqu_upstream_map <- function(project_path, gauge_channels, cache = TRUE) {
  cache_path <- file.path(project_path, ".aqu_upstream_map.rds")

  if (cache && file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(cached) &&
        all(as.character(gauge_channels) %in% names(cached$upstream_aqus))) {
      cat(sprintf("[baseflow] Loaded aquifer upstream map from cache (%d gauges).\n",
                  length(gauge_channels)))
      return(cached)
    }
  }

  cat(sprintf("[baseflow] Building aquifer upstream map from %s ...\n",
              project_path))

  # Channels — build downstream and reversed (upstream) graphs
  chans <- parse_swat_con(file.path(project_path, "chandeg.con"))
  cat(sprintf("[baseflow]   chandeg.con   : %d channels\n", length(chans)))

  chan_down <- list()
  for (ch in chans) {
    cid <- ch$id
    downs <- integer(0)
    for (cn in ch$conns) {
      if (!is.null(cn) && cn$obj_typ %in% c("sdc", "cha")) {
        downs <- c(downs, cn$obj_id)
      }
    }
    chan_down[[as.character(cid)]] <- unique(downs)
  }

  # Reversed graph: for each channel, which upstream channels feed it
  chan_up <- setNames(vector("list", length(chan_down)), names(chan_down))
  for (cid in names(chan_up)) chan_up[[cid]] <- integer(0)
  for (cid in names(chan_down)) {
    for (down in chan_down[[cid]]) {
      key <- as.character(down)
      if (!is.null(chan_up[[key]])) {
        chan_up[[key]] <- c(chan_up[[key]], as.integer(cid))
      }
    }
  }

  # BFS upstream (inclusive) from one channel
  upstream_chans_of <- function(g) {
    g <- as.integer(g)
    visited <- integer(0)
    queue   <- g
    while (length(queue) > 0) {
      v <- queue[1]; queue <- queue[-1]
      if (v %in% visited) next
      visited <- c(visited, v)
      parents <- chan_up[[as.character(v)]]
      if (!is.null(parents) && length(parents) > 0) {
        queue <- c(queue, setdiff(parents, visited))
      }
    }
    visited
  }

  # Aquifers — parse aquifer.con and split into:
  #   aqu_to_cha[i] : channel id for FLOODPLAIN aquifers (obj_typ = sdc/cha)
  #   aqu_to_aqu[i] : neighbouring floodplain id for UPLAND aquifers (obj_typ = aqu)
  # Deep aquifer (regional, 1 per basin): no Connection 1 / both fields NA.
  aqus <- parse_swat_con(file.path(project_path, "aquifer.con"))
  cat(sprintf("[baseflow]   aquifer.con   : %d aquifers\n", length(aqus)))

  aqu_to_cha  <- integer(0)
  aqu_to_aqu  <- integer(0)
  aqu_area_ha <- numeric(0)
  for (aq in aqus) {
    aid <- aq$id
    aqu_area_ha[as.character(aid)] <- aq$area
    aqu_to_cha[as.character(aid)] <- NA_integer_
    aqu_to_aqu[as.character(aid)] <- NA_integer_
    for (cn in aq$conns) {
      if (is.null(cn)) next
      if (cn$hyd_typ == "rhg") next   # skip recharge-to-deep connection
      if (cn$obj_typ %in% c("sdc", "cha")) {
        aqu_to_cha[as.character(aid)] <- cn$obj_id
      } else if (cn$obj_typ == "aqu") {
        aqu_to_aqu[as.character(aid)] <- cn$obj_id
      }
    }
  }

  n_floodplain <- sum(!is.na(aqu_to_cha))
  n_upland     <- sum(!is.na(aqu_to_aqu) & is.na(aqu_to_cha))
  n_deep       <- sum(is.na(aqu_to_cha) & is.na(aqu_to_aqu))
  cat(sprintf("[baseflow]   aquifer split : %d floodplain (→cha) + %d upland (→aqu) + %d deep (terminal)\n",
              n_floodplain, n_upland, n_deep))

  # For each floodplain aquifer, find its upstream upland chain
  # (aqu→aqu BFS — usually only 1-step deep but generic)
  upland_feeders_of <- function(fp_id) {
    visited <- integer(0)
    queue   <- as.integer(fp_id)
    seen    <- as.integer(fp_id)   # don't include floodplain itself
    feeders <- integer(0)
    while (length(queue) > 0) {
      v <- queue[1]; queue <- queue[-1]
      # who drains TO v via aqu?
      contributors <- as.integer(names(aqu_to_aqu)[
        which(aqu_to_aqu == v)
      ])
      contributors <- setdiff(contributors, seen)
      feeders <- c(feeders, contributors)
      seen    <- c(seen, contributors)
      queue   <- c(queue, contributors)
    }
    feeders
  }

  # For each gauge channel: upstream channels + aquifers (floodplain + upland)
  upstream_chans_map     <- list()
  upstream_aqus_fp_map   <- list()   # floodplain only (numerator-relevant)
  upstream_aqus_all_map  <- list()   # floodplain + upland (denominator)
  for (g in gauge_channels) {
    chans_up <- upstream_chans_of(g)
    upstream_chans_map[[as.character(g)]] <- sort(chans_up)

    fp_hits <- which(aqu_to_cha %in% chans_up)
    fp_ids  <- as.integer(names(aqu_to_cha)[fp_hits])
    upstream_aqus_fp_map[[as.character(g)]] <- sort(fp_ids)

    upland_ids <- integer(0)
    for (fp in fp_ids) upland_ids <- c(upland_ids, upland_feeders_of(fp))
    upstream_aqus_all_map[[as.character(g)]] <- sort(unique(c(fp_ids, upland_ids)))
  }

  result <- list(
    aqu_to_cha          = aqu_to_cha,
    aqu_to_aqu          = aqu_to_aqu,
    aqu_area_ha         = aqu_area_ha,
    upstream_chans      = upstream_chans_map,
    upstream_aqus_fp    = upstream_aqus_fp_map,    # only floodplain (flo_cha>0)
    upstream_aqus_all   = upstream_aqus_all_map,   # floodplain + upland (area)
    gauge_channels      = sort(unique(as.integer(gauge_channels)))
  )

  if (cache) {
    saveRDS(result, cache_path)
    cat(sprintf("[baseflow] Cached to %s\n", cache_path))
  }
  result
}


# ------------------------------------------------------------------------------
# compute_baseflow_comparison — main entry
#
# For each (run × gauge_channel):
#   total_aqu_volume_m3yr = Σ aqu_flo_cha_<i> [mm/yr] × area_<i> [ha] × 10
#   total_Q_volume_m3yr   = mean(daily_Q [m³/s] over period) × 86400 × 365
#   bfi_swat_internal     = total_aqu_volume_m3yr / total_Q_volume_m3yr
#
# @param sim SWATrunR sim object (must have aqu_flo_cha_<i> + <multi_day>_<chid>)
# @param project_path SWAT+ project root
# @param bfi_targets Tibble from read_bfi_targets (channel_id, ...)
# @param valid_runs Character vector of run names
# @param multi_day_output YAML key for daily multi-channel Q (default "vazao_multi_day")
# @return Tibble with run_id + per-channel: bfi_swat_internal_chXX, gap_internal_chXX
# ------------------------------------------------------------------------------
compute_baseflow_comparison <- function(sim, project_path, bfi_targets,
                                         valid_runs,
                                         multi_day_output = "vazao_multi_day") {

  # 1. Topology
  topo <- build_aqu_upstream_map(project_path,
                                  gauge_channels = bfi_targets$channel_id)

  # 2. Pre-extract aquifer flows (annual avg, mm) per run
  #    sim$simulation$aqu_flo_cha_<i> is a tibble with date + run_* cols
  aqu_ids <- as.integer(names(topo$aqu_to_cha))
  aqu_keys <- paste0("aqu_flo_cha_", aqu_ids)
  missing_keys <- setdiff(aqu_keys, names(sim$simulation))
  if (length(missing_keys) > 0) {
    stop("[baseflow] Missing per-aquifer outputs in sim: ",
         paste(head(missing_keys, 5), collapse = ", "), " ...\n",
         "  Enable extra_outputs.aquifer with aqu_vars including 'flo_cha'.")
  }

  # Build a matrix [run × aquifer] of annual avg flo_cha (mm/yr)
  aqu_flo_mat <- matrix(NA_real_, nrow = length(valid_runs),
                         ncol = length(aqu_ids),
                         dimnames = list(valid_runs, as.character(aqu_ids)))
  for (j in seq_along(aqu_ids)) {
    df <- sim$simulation[[aqu_keys[j]]]
    if (is.null(df)) next
    # Annual avg files are 1 row per run; if multi-row, use mean
    for (run_nm in valid_runs) {
      v <- df[[run_nm]]
      if (is.null(v)) next
      aqu_flo_mat[run_nm, as.character(aqu_ids[j])] <- mean(v, na.rm = TRUE)
    }
  }

  # 3. Pre-extract daily Q at each gauge channel
  q_mean_mat <- matrix(NA_real_, nrow = length(valid_runs),
                        ncol = nrow(bfi_targets),
                        dimnames = list(valid_runs,
                                        as.character(bfi_targets$channel_id)))
  for (j in seq_len(nrow(bfi_targets))) {
    cid <- bfi_targets$channel_id[j]
    key <- sprintf("%s_%d", multi_day_output, cid)
    df  <- sim$simulation[[key]]
    if (is.null(df)) {
      warning(sprintf("[baseflow] daily Q output missing for channel %d (%s)",
                      cid, key))
      next
    }
    for (run_nm in valid_runs) {
      v <- df[[run_nm]]
      if (is.null(v)) next
      q_mean_mat[run_nm, as.character(cid)] <- mean(v, na.rm = TRUE)
    }
  }

  # 4. For each (run × gauge), compute BFI_swat_internal AND BF_mm equivalent.
  #
  # Numerator (only floodplain aquifers carry flo_cha > 0):
  #   BF_vol_m3yr = Σ (flo_cha_i × area_i × 10)   for i ∈ upstream_aqus_fp
  #
  # Denominator for mm normalization (all upstream contributing area):
  #   A_up_total_ha = Σ area_i              for i ∈ upstream_aqus_all
  #                   (floodplain + upland of upstream LSUs;
  #                    matches basin_aqu_aa.flo_cha normalization)
  #
  # Outputs per gauge channel:
  #   bf_mm_at_X       : BF_vol / (A_up × 10)  — mm/yr depth-equivalent
  #   bfi_swat_int_X   : BF_vol / Q_vol       — dimensionless
  res <- tibble(run_id = valid_runs)

  for (j in seq_len(nrow(bfi_targets))) {
    cid     <- bfi_targets$channel_id[j]
    aqu_fp  <- topo$upstream_aqus_fp[[as.character(cid)]]
    aqu_all <- topo$upstream_aqus_all[[as.character(cid)]]
    if (length(aqu_fp) == 0) {
      warning(sprintf("[baseflow] channel %d has 0 upstream floodplain aquifers", cid))
      next
    }

    # Volume from upstream FLOODPLAIN aquifers (m³/yr) — numerator
    areas_fp  <- topo$aqu_area_ha[as.character(aqu_fp)]
    flo_sub   <- aqu_flo_mat[, as.character(aqu_fp), drop = FALSE]
    aqu_vol_m3yr <- as.numeric(flo_sub %*% areas_fp) * 10

    # Total upstream area in ha (floodplain + upland) — for mm normalization
    A_up_total_ha <- sum(topo$aqu_area_ha[as.character(aqu_all)])

    # Channel volume (m³/yr) from mean daily Q (m³/s)
    q_mean_m3s <- q_mean_mat[, as.character(cid)]
    q_vol_m3yr <- q_mean_m3s * 86400 * 365

    # mm depth equivalent of baseflow (analog of basin_aqu_aa.flo_cha)
    bf_mm <- aqu_vol_m3yr / (A_up_total_ha * 10)

    # Dimensionless BFI — guard against zero/NA Q
    bfi_int <- aqu_vol_m3yr / q_vol_m3yr
    bfi_int[!is.finite(bfi_int) | q_vol_m3yr <= 0] <- NA_real_

    res[[paste0("bf_mm_swat_ch",         cid)]] <- bf_mm
    res[[paste0("bfi_swat_internal_ch",  cid)]] <- bfi_int
    res[[paste0("aqu_vol_m3yr_ch",       cid)]] <- aqu_vol_m3yr
    res[[paste0("q_vol_m3yr_ch",         cid)]] <- q_vol_m3yr
    res[[paste0("A_up_total_ha_ch",      cid)]] <- A_up_total_ha
    res[[paste0("n_aqu_fp_ch",           cid)]] <- length(aqu_fp)
    res[[paste0("n_aqu_all_ch",          cid)]] <- length(aqu_all)
  }

  res
}


# ------------------------------------------------------------------------------
# join_baseflow_with_eckhardt — convenience joiner
#
# Adds gap_internal_eckhardt_<chid> column for each channel:
#   gap = bfi_swat_internal - bfi_eckhardt_sim
#
# Negative gap → SWAT+ internal partition UNDERESTIMATES baseflow vs Eckhardt
#                (i.e., aquifer module not generating enough; latq dominant)
# Positive gap → SWAT+ internal OVERESTIMATES baseflow vs Eckhardt
#                (i.e., aquifer too active relative to channel hydrograph shape)
# ------------------------------------------------------------------------------
join_baseflow_with_eckhardt <- function(bfi_metrics, baseflow_internal,
                                         bfi_targets) {
  out <- bfi_metrics %>% left_join(baseflow_internal, by = "run_id")

  for (cid in bfi_targets$channel_id) {
    eck_col <- paste0("bfi_sim_ch",            cid)
    int_col <- paste0("bfi_swat_internal_ch", cid)
    gap_col <- paste0("gap_internal_eckhardt_ch", cid)
    if (eck_col %in% names(out) && int_col %in% names(out)) {
      out[[gap_col]] <- out[[int_col]] - out[[eck_col]]
    }
  }
  out
}
