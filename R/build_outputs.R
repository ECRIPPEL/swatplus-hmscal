# ==============================================================================
# build_outputs.R — Construct the SWATrunR `output` list for run_swatplus()
#
# Combines two YAML blocks:
#   1) base block (cfg$cal_outputs OR cfg$morris_outputs) — flat schema with one
#      entry per (file, variable, unit). `unit` is typically a scalar.
#   2) cfg$extra_outputs (OPTIONAL) — bulk retention block matching the pattern
#      of swatplus-calibrator's `h3_outputs:`. Allows capturing per-channel,
#      per-LSU, per-HRU, per-aquifer outputs without enumerating each unit.
#
# extra_outputs schema (all fields optional):
#
#   extra_outputs:
#     enabled: true                     # master switch (default: false)
#
#     channels:                         # channel_sd_aa per reach (H_diss)
#       enabled: true
#       units: [22, 33, 79, 113, 209]   # OR "auto" → uses cha_units_file
#       cha_units_file: "data/cha_units.txt"   # optional, used when units="auto"
#       vars: ["flo_in", "flo_out", "evap", "seep"]   # default if omitted
#
#     lsu:                              # lsunit_wb_aa per LSU
#       enabled: true
#       n_lsu: 422                      # required
#       wb_vars: ["precip", "et", "wateryld", "surq_gen", "latq", "perc"]
#
#     hru:                              # hru_wb_aa per filtered HRU
#       enabled: true
#       hru_ids_file: "data/hru_ids.txt"   # required (file with one int per line)
#       wb_vars: ["precip", "et", "surq_gen", "latq", "perc", "wateryld"]
#
#     aquifer:                          # aquifer_aa per aquifer
#       enabled: true
#       n_aqu: 47                       # required
#       aqu_vars: ["flo", "rchrg", "seep", "revap", "flo_cha"]
#
# Output naming convention in sim$simulation:
#   - base block keys are preserved as-is
#   - channels   → cha_aa_<var>      (e.g., cha_aa_flo_out)
#   - lsu        → lsu_<var>         (e.g., lsu_precip)
#   - hru        → hru_<var>         (e.g., hru_et)
#   - aquifer    → aqu_<var>         (e.g., aqu_rchrg)
#
# REQUIRES: SWATrunR (define_output already loaded by the caller).
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

.read_int_file <- function(path) {
  if (!file.exists(path)) {
    stop("extra_outputs: file not found: ", path)
  }
  raw <- paste(readLines(path, warn = FALSE), collapse = " ")
  ids <- as.integer(unlist(regmatches(raw, gregexpr("[0-9]+", raw))))
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0) {
    stop("extra_outputs: no integer IDs parsed from: ", path)
  }
  ids
}

build_outputs <- function(base_block, extra_block = NULL) {
  if (is.null(base_block) || length(base_block) == 0) {
    stop("build_outputs: base block is empty (cfg$cal_outputs / cfg$morris_outputs).")
  }

  # 1) Base outputs — flat, one entry per (file, variable, unit).
  #    `unit` accepts a scalar or a vector; yaml::read_yaml returns sequences
  #    as lists, so unlist() to coerce to a numeric vector before define_output.
  outputs <- lapply(base_block, function(o) {
    unit <- as.integer(unlist(o$unit))
    define_output(file = o$file, variable = o$variable, unit = unit)
  })

  # 2) Extra outputs — bulk retention (optional)
  if (is.null(extra_block) || !isTRUE(extra_block$enabled)) {
    return(outputs)
  }

  cat("[extra_outputs] enabled — extending output list\n")

  # ---- (a) Channels (channel_sd_aa per reach) -------------------------------
  ch <- extra_block$channels
  if (!is.null(ch) && isTRUE(ch$enabled)) {
    if (identical(ch$units, "auto")) {
      if (is.null(ch$cha_units_file)) {
        stop("extra_outputs.channels.units = 'auto' requires cha_units_file.")
      }
      cha_units <- .read_int_file(ch$cha_units_file)
    } else {
      cha_units <- as.integer(unlist(ch$units))
    }
    cha_vars <- as.character(unlist(ch$vars %||%
                                      c("flo_in", "flo_out", "evap", "seep")))
    for (v in cha_vars) {
      outputs[[paste0("cha_aa_", v)]] <- define_output(
        file = "channel_sd_aa", variable = v, unit = cha_units
      )
    }
    cat(sprintf("  channels  : %d reaches x %d vars (channel_sd_aa)\n",
                length(cha_units), length(cha_vars)))
  }

  # ---- (b) LSU (lsunit_wb_aa per LSU) ---------------------------------------
  lsu <- extra_block$lsu
  if (!is.null(lsu) && isTRUE(lsu$enabled)) {
    if (is.null(lsu$n_lsu)) stop("extra_outputs.lsu.n_lsu is required.")
    n_lsu <- as.integer(lsu$n_lsu)
    wb_vars <- as.character(unlist(lsu$wb_vars %||%
                                     c("precip", "et", "wateryld",
                                       "surq_gen", "latq", "perc")))
    for (v in wb_vars) {
      outputs[[paste0("lsu_", v)]] <- define_output(
        file = "lsunit_wb_aa", variable = v, unit = seq_len(n_lsu)
      )
    }
    cat(sprintf("  lsu       : %d LSUs x %d vars (lsunit_wb_aa)\n",
                n_lsu, length(wb_vars)))
  }

  # ---- (c) HRU (hru_wb_aa per filtered HRU) ---------------------------------
  hru <- extra_block$hru
  if (!is.null(hru) && isTRUE(hru$enabled)) {
    if (is.null(hru$hru_ids_file)) {
      stop("extra_outputs.hru.hru_ids_file is required.")
    }
    hru_ids <- .read_int_file(hru$hru_ids_file)
    wb_vars <- as.character(unlist(hru$wb_vars %||%
                                     c("precip", "et", "surq_gen",
                                       "latq", "perc", "wateryld")))
    for (v in wb_vars) {
      outputs[[paste0("hru_", v)]] <- define_output(
        file = "hru_wb_aa", variable = v, unit = hru_ids
      )
    }
    cat(sprintf("  hru       : %d HRUs x %d vars (hru_wb_aa)\n",
                length(hru_ids), length(wb_vars)))
  }

  # ---- (d) Aquifer (aquifer_aa per aquifer) ---------------------------------
  aq <- extra_block$aquifer
  if (!is.null(aq) && isTRUE(aq$enabled)) {
    if (is.null(aq$n_aqu)) stop("extra_outputs.aquifer.n_aqu is required.")
    n_aqu <- as.integer(aq$n_aqu)
    aqu_vars <- as.character(unlist(aq$aqu_vars %||%
                                      c("flo", "dep_wt", "stor", "rchrg",
                                        "seep", "revap",
                                        "flo_cha", "flo_res", "flo_ls")))
    for (v in aqu_vars) {
      outputs[[paste0("aqu_", v)]] <- define_output(
        file = "aquifer_aa", variable = v, unit = seq_len(n_aqu)
      )
    }
    cat(sprintf("  aquifer   : %d aquifers x %d vars (aquifer_aa)\n",
                n_aqu, length(aqu_vars)))
  }

  # ---- (e) Aquifer monthly (aquifer_mon per aquifer) ------------------------
  # NEW 2026-06-23 — monthly per-aquifer outputs (needs print.prt monthly
  # aquifer output enabled). Same pattern as the AA block. Output keys are
  # named aqu_mon_<var> to disambiguate from aqu_<var> (annual).
  aqmon <- extra_block$aquifer_mon
  if (!is.null(aqmon) && isTRUE(aqmon$enabled)) {
    if (is.null(aqmon$n_aqu)) stop("extra_outputs.aquifer_mon.n_aqu is required.")
    n_aqu_mon <- as.integer(aqmon$n_aqu)
    aqu_mon_vars <- as.character(unlist(aqmon$aqu_vars %||%
                                          c("flo_cha")))
    for (v in aqu_mon_vars) {
      outputs[[paste0("aqu_mon_", v)]] <- define_output(
        file = "aquifer_mon", variable = v, unit = seq_len(n_aqu_mon)
      )
    }
    cat(sprintf("  aquifer_mon: %d aquifers x %d vars (aquifer_mon)\n",
                n_aqu_mon, length(aqu_mon_vars)))
  }

  outputs
}
