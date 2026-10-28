# ==============================================================================
# write_run_info.R — Write a human-readable metadata file (run_info.txt) into
# each run output folder, recording the conditions of the simulation campaign.
#
# The file is intended for offline reproducibility: when revisiting an old
# output folder months later, you can read run_info.txt to recover the iteration,
# simulation period, warmup, parameter ranges used, targets, and key result
# counts without having to load the .rds.
# ==============================================================================

#' Write run metadata to <path>/run_info.txt
#'
#' @param path        Output folder of the run (ver$path).
#' @param cfg         Full config list (yaml::read_yaml output).
#' @param run_type    "Morris screening", "CAL", or "VAL".
#' @param param_info  Tibble with columns short, par_name, min, max
#'                    (parameters actually used in the run).
#' @param summary     Named list of result counts. Recommended keys:
#'                      n_simulations    : total simulations launched
#'                      n_valid          : valid runs (no SWAT+ failure)
#'                      n_satisfactory   : runs passing primary criterion
#'                      n_accepted       : accepted subset (CAL only, optional)
#'                      top_param        : top-ranked parameter (Morris screening only)
#'                      top_param_mu     : mu_star of top parameter (Morris screening only)
#' @param extras      Optional named list of additional key/value lines to append
#'                    (free-form).
write_run_info <- function(path, cfg, run_type = c("Morris screening", "CAL", "VAL"),
                           param_info = NULL, summary = list(), extras = list()) {
  run_type <- match.arg(run_type)
  out_file <- file.path(path, "run_info.txt")

  # Resolve generic vs iteration-specific config keys
  iter      <- cfg$iteration %||% NA
  iter_lab  <- cfg$script_version %||% "(unset)"
  if (run_type == "Morris screening") {
    period_start <- cfg$morris$start_date %||% cfg$start_date
    period_end   <- cfg$morris$end_date   %||% cfg$end_date
    warmup       <- cfg$morris$warmup     %||% cfg$warmup
    n_traj       <- cfg$morris$m_traj     %||% NA
    delta        <- cfg$morris$delta      %||% NA
    tag_prefix   <- cfg$morris$tag_prefix %||% "(unset)"
  } else if (run_type == "VAL") {
    period_start <- cfg$cal$start_date %||% cfg$start_date
    period_end   <- cfg$cal$end_date   %||% cfg$end_date
    warmup       <- cfg$cal$warmup     %||% cfg$warmup
    n_traj       <- NA
    delta        <- NA
    tag_prefix   <- cfg$validation$tag_prefix %||% "VAL"
  } else {
    period_start <- cfg$cal$start_date %||% cfg$start_date
    period_end   <- cfg$cal$end_date   %||% cfg$end_date
    warmup       <- cfg$cal$warmup     %||% cfg$warmup
    n_traj       <- NA
    delta        <- NA
    tag_prefix   <- cfg$cal$tag_prefix %||% "(unset)"
  }

  L <- character()
  add <- function(x = "") {
    if (length(x) == 0L || is.null(x)) x <- ""
    L[length(L) + 1L] <<- as.character(x)[1L]
  }
  sec <- function(title) {
    add("")
    add(strrep("-", 78))
    add(title)
    add(strrep("-", 78))
  }

  # ---- HEADER ----
  add(strrep("=", 78))
  add(sprintf("  swatplus-hmscal — %s run summary", run_type))
  add(strrep("=", 78))
  add(sprintf("  generated  : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
  add(sprintf("  output dir : %s", normalizePath(path, winslash = "/", mustWork = FALSE)))
  add(sprintf("  config     : %s", iter_lab))
  add(sprintf("  iteration  : %s", as.character(iter)))
  add(sprintf("  tag_prefix : %s", tag_prefix))
  add(sprintf("  R version  : %s", R.version.string))

  # ---- SIMULATION SETTINGS ----
  sec("Simulation settings")
  add(sprintf("  project_path     : %s", cfg$project_path %||% "(not set)"))
  add(sprintf("  start_date       : %s", as.character(period_start)))
  add(sprintf("  end_date         : %s", as.character(period_end)))
  add(sprintf("  warmup (years)   : %s", as.character(warmup)))
  add(sprintf("  threads          : %s", as.character(cfg$threads %||% NA)))
  add(sprintf("  seed             : %s", as.character(cfg$seed %||% NA)))
  add(sprintf("  channel_exutorio : %s", as.character(cfg$channel_exutorio %||% NA)))
  add(sprintf("  temporal_scale   : %s", as.character(cfg$temporal_scale %||% NA)))

  # ---- TYPE-SPECIFIC SETTINGS ----
  if (run_type == "Morris screening") {
    sec("Morris screening / Morris design")
    add(sprintf("  m_traj (trajectories) : %s", as.character(n_traj)))
    add(sprintf("  delta (perturbation)  : %s", as.character(delta)))
    add(sprintf("  metrics computed      : %s",
                paste(cfg$morris_metrics %||% character(), collapse = ", ")))
    if (!is.null(cfg$morris_targets)) {
      if (!is.null(cfg$morris_targets$rc)) {
        add(sprintf("  RC target             : %s", cfg$morris_targets$rc))
      } else if (!is.null(cfg$morris_targets$q_p)) {
        add(sprintf("  WB targets (annual)   : RC=%s, ET/P=%s",
                    cfg$morris_targets$q_p, cfg$morris_targets$et_p))
      }
    }
  } else {
    sec("Calibration settings")
    add(sprintf("  cal type            : %s", cfg$cal$type %||% "(unset)"))
    add(sprintf("  n_simulations       : %s", as.character(cfg$cal$n_simulations %||% NA)))
    add(sprintf("  inheritance_source  : %s",
                cfg$cal$inheritance_source %||% "satisfactory"))
    if (!is.null(cfg$cal$threshold_nse)) {
      add(sprintf("  threshold_nse        : %s", cfg$cal$threshold_nse))
      add(sprintf("  threshold_kge        : %s", cfg$cal$threshold_kge))
      add(sprintf("  threshold_abs_pbias  : %s", cfg$cal$threshold_abs_pbias))
    }
    if (!is.null(cfg$cal$rc_target)) {
      add(sprintf("  RC target           : %s   threshold_obj_rc = %s",
                  cfg$cal$rc_target,
                  cfg$cal$threshold_obj_rc %||% NA))
      add(sprintf("  BFI proxy threshold : threshold_obj_bfi = %s   (CSV: %s)",
                  cfg$cal$threshold_obj_bfi %||% NA,
                  cfg$bfi_targets_csv %||% "(missing)"))
    } else if (!is.null(cfg$cal$targets)) {
      add(sprintf("  WB targets (soft)   : ET/P=%s, Q/P=%s, threshold_pct=%s",
                  cfg$cal$targets$et_rto, cfg$cal$targets$wyld_rto,
                  cfg$cal$threshold_pct %||% NA))
    }
    if (!is.null(cfg$diagnostic_balance)) {
      add(sprintf("  WB diagnostic       : ET/P=%s, Q/P=%s, threshold_pct=%s",
                  cfg$diagnostic_balance$et_rto, cfg$diagnostic_balance$wyld_rto,
                  cfg$diagnostic_balance$threshold_pct %||% NA))
    }
    if (!is.null(cfg$inherit)) {
      add(sprintf("  inherit from        : %s", cfg$inherit$from %||% "?"))
      if (!is.null(cfg$inherit$drop) && length(cfg$inherit$drop) > 0) {
        add(sprintf("  inherit.drop        : %s",
                    paste(cfg$inherit$drop, collapse = ", ")))
      }
      if (!is.null(cfg$inherit$override_ranges) &&
          length(cfg$inherit$override_ranges) > 0) {
        add(sprintf("  inherit.override    : %s",
                    paste(names(cfg$inherit$override_ranges), collapse = ", ")))
      }
    }
  }

  # ---- PARAMETERS ----
  if (!is.null(param_info) && nrow(param_info) > 0) {
    sec(sprintf("Parameters used (%d)", nrow(param_info)))
    add(sprintf("  %-12s %-50s %12s %12s",
                "short", "par_name", "min", "max"))
    for (i in seq_len(nrow(param_info))) {
      add(sprintf("  %-12s %-50s %12.4f %12.4f",
                  param_info$short[i],
                  substr(param_info$par_name[i], 1, 50),
                  param_info$min[i],
                  param_info$max[i]))
    }
  }

  # ---- RESULT SUMMARY ----
  if (length(summary) > 0) {
    sec("Result summary")
    for (k in names(summary)) {
      v <- summary[[k]]
      if (is.numeric(v)) v <- format(v, digits = 6)
      add(sprintf("  %-22s : %s", k, paste(v, collapse = ", ")))
    }
  }

  # ---- EXTRAS ----
  if (length(extras) > 0) {
    sec("Notes")
    for (k in names(extras)) add(sprintf("  %s : %s", k, extras[[k]]))
  }

  add("")
  add(strrep("=", 78))
  add("  end of run_info")
  add(strrep("=", 78))

  writeLines(L, out_file, useBytes = FALSE)
  invisible(out_file)
}

# Tiny helper because base R lacks a default-on-NULL operator pre-4.4
`%||%` <- function(a, b) if (!is.null(a)) a else b
