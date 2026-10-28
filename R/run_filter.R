# ==============================================================================
# run_filter.R — Invalid series filtering and valid run extraction
# ==============================================================================

library(dplyr)
library(tidyr)

#' Filter out runs with invalid time series (Inf, NA, sd==0, extreme outliers)
#'
#' @param comparativo data.frame with Date, Flow (observed), and run_* columns
#' @param valid_runs Character vector of run names to evaluate
#' @param obs_col Name of the observed column (default "Flow")
#' @param outlier_mult Multiplier of max observed for upper limit (default 10)
#' @return List with $runs_ok (filtered vector) and $runs_excluded
filter_invalid_runs <- function(comparativo, valid_runs, obs_col = "Flow",
                                outlier_mult = 10) {
  upper_limit <- max(comparativo[[obs_col]], na.rm = TRUE) * outlier_mult

  runs_ok <- valid_runs[sapply(valid_runs, function(run_nm) {
    sim_vals <- comparativo[[run_nm]]

    if (is.null(sim_vals)) return(FALSE)
    if (all(is.na(sim_vals))) return(FALSE)
    if (any(!is.finite(sim_vals))) return(FALSE)
    if (length(na.omit(sim_vals)) < 2) return(FALSE)
    if (sd(sim_vals, na.rm = TRUE) == 0) return(FALSE)
    if (max(sim_vals, na.rm = TRUE) > upper_limit) return(FALSE)

    TRUE
  })]

  runs_excluded <- setdiff(valid_runs, runs_ok)

  if (length(runs_excluded) > 0) {
    cat("Runs excluded (invalid series): ", length(runs_excluded), "\n", sep = "")
  }

  if (length(runs_ok) == 0) {
    stop("All runs were excluded by the invalid series filter.")
  }

  list(runs_ok = runs_ok, runs_excluded = runs_excluded)
}

#' Extract valid run names from a SWATrunR simulation object
#'
#' Handles both shapes that SWATrunR produces:
#'   (a) scalar-unit output: sim_output[[nm]] is a tibble with run_* columns
#'   (b) vector-unit output: sim_output[[nm]] is NULL; instead the per-unit
#'       tibbles live as sim_output[[paste0(nm, "_", id)]] (e.g.,
#'       "vazao_multi_day_5", "vazao_multi_day_19", ...). For (b) we use
#'       the FIRST matching per-unit tibble as the canonical run set.
#'
#' @param sim_output Simulation list (sim$simulation)
#' @param output_names Character vector of output names to check
#' @return Character vector of run names present in ALL outputs
get_valid_runs <- function(sim_output, output_names) {
  runs_by_output <- lapply(output_names, function(nm) {
    if (!is.null(sim_output[[nm]])) {
      grep("^run_", names(sim_output[[nm]]), value = TRUE)
    } else {
      # Vector-unit output — find first sim_output[["<nm>_<id>"]] match
      pat  <- sprintf("^%s_[0-9]+$", nm)
      keys <- grep(pat, names(sim_output), value = TRUE)
      if (length(keys) == 0) {
        warning("get_valid_runs: output '", nm,
                "' not found (no scalar key, no <name>_<id> pattern).")
        return(character(0))
      }
      grep("^run_", names(sim_output[[keys[1]]]), value = TRUE)
    }
  })

  runs_by_output <- runs_by_output[lengths(runs_by_output) > 0]
  if (length(runs_by_output) == 0) return(character(0))
  Reduce(intersect, runs_by_output)
}

#' Build observed-vs-simulated comparison data.frame
#'
#' @param sim_df Simulated data.frame with 'date' column
#' @param obs_df Observed data.frame with Date and Flow columns
#' @param runs Character vector of run names
#' @param floor_month If TRUE, floor dates to month (for monthly data)
#' @return data.frame with Date, Flow, and run columns
build_comparativo <- function(sim_df, obs_df, runs, floor_month = FALSE) {
  flow_sim <- sim_df %>%
    select(date, all_of(runs)) %>%
    mutate(Date = as.Date(date, origin = "1970-01-01")) %>%
    select(-date)

  if (floor_month) {
    flow_sim <- flow_sim %>%
      mutate(Date = lubridate::floor_date(Date, "month"))
    obs_df <- obs_df %>%
      mutate(Date = lubridate::floor_date(Date, "month"))
  }

  flow_sim <- flow_sim %>% arrange(Date)

  inner_join(obs_df, flow_sim, by = "Date") %>%
    drop_na(Flow)
}
