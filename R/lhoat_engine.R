# ==============================================================================
# lhoat_engine.R — LH-OAT trajectory generation, scaling, EE computation
# ==============================================================================

library(lhs)
library(dplyr)
library(tibble)
library(purrr)

#' Generate LH-OAT normalized trajectories in [0,1]
#'
#' @param m_traj Number of trajectories
#' @param p Number of parameters
#' @param delta Perturbation step (default 0.20)
#' @param param_names Character vector of short parameter names
#' @return data.frame with columns: param_names + traj_id + step_id + run_index
generate_lhoat_trajectories <- function(m_traj, p, delta, param_names) {
  traj_list <- vector("list", m_traj)
  run_counter <- 1

  for (i in seq_len(m_traj)) {
    x0 <- as.numeric(randomLHS(1, p))

    X <- matrix(NA_real_, nrow = p + 1, ncol = p)
    X[1, ] <- x0

    for (j in seq_len(p)) {
      xj <- X[j, ]
      xj[j] <- if (xj[j] <= (1 - delta)) xj[j] + delta else xj[j] - delta
      X[j + 1, ] <- xj
    }

    traj_df <- as.data.frame(X)
    names(traj_df) <- param_names
    traj_df$traj_id   <- i
    traj_df$step_id   <- 0:p
    traj_df$run_index <- run_counter:(run_counter + p)

    run_counter <- run_counter + p + 1
    traj_list[[i]] <- traj_df
  }

  bind_rows(traj_list)
}

#' Scale parameters from [0,1] to [min,max]
#'
#' @param design_norm Normalized data.frame (output of generate_lhoat_trajectories)
#' @param param_info Tibble with columns: short, par_name, min, max
#' @return List with $vals (scaled values) and $swat_params (SWATrunR-named tibble)
scale_parameters <- function(design_norm, param_info) {
  vals <- design_norm

  for (k in seq_len(nrow(param_info))) {
    nm <- param_info$short[k]
    vals[[nm]] <- param_info$min[k] + vals[[nm]] * (param_info$max[k] - param_info$min[k])
  }

  swat_params <- vals %>%
    dplyr::transmute(
      !!!setNames(
        lapply(param_info$short, function(x) as.name(x)),
        param_info$par_name
      )
    )

  list(vals = vals, swat_params = swat_params)
}

#' Compute Elementary Effects for a single trajectory
#'
#' @param df_traj data.frame for one trajectory (filtered by traj_id)
#' @param metric_col Name of the metric column
#' @param param_names Character vector of short parameter names
#' @return Tibble with EE per step
calc_ee <- function(df_traj, metric_col, param_names) {
  df_traj <- df_traj %>% arrange(step_id)
  out <- vector("list", nrow(df_traj) - 1)

  for (i in 2:nrow(df_traj)) {
    prev <- df_traj[i - 1, ]
    curr <- df_traj[i, ]

    diffs <- sapply(param_names, function(v) !isTRUE(all.equal(prev[[v]], curr[[v]])))
    changed <- param_names[diffs]

    if (length(changed) != 1) {
      out[[i - 1]] <- tibble(
        traj_id   = curr$traj_id,
        step_id   = curr$step_id,
        parameter = NA_character_,
        metric    = metric_col,
        delta_x   = NA_real_,
        delta_y   = NA_real_,
        EE        = NA_real_
      )
    } else {
      delta_x <- as.numeric(curr[[changed]] - prev[[changed]])
      delta_y <- as.numeric(curr[[metric_col]] - prev[[metric_col]])

      out[[i - 1]] <- tibble(
        traj_id   = curr$traj_id,
        step_id   = curr$step_id,
        parameter = changed,
        metric    = metric_col,
        delta_x   = delta_x,
        delta_y   = delta_y,
        EE        = ifelse(is.na(delta_x) || delta_x == 0, NA_real_, delta_y / delta_x)
      )
    }
  }

  bind_rows(out)
}

#' Compute EE across all metrics and trajectories
#'
#' @param metrics_df data.frame with metrics + parameters + traj_id + step_id
#' @param metric_names Character vector of metric names (e.g., c("NSE_mon", "obj_total"))
#' @param param_names Character vector of short parameter names
#' @return Tibble with all EE values
calc_all_ee <- function(metrics_df, metric_names, param_names) {
  map_dfr(
    metric_names,
    ~ map_dfr(
      split(metrics_df, metrics_df$traj_id),
      calc_ee,
      metric_col  = .x,
      param_names = param_names
    )
  )
}

#' Compute sensitivity indices (mu_star, sigma) from Elementary Effects
#'
#' @param ee_all Tibble with EE values (output of calc_all_ee)
#' @return List with:
#'   - $sensitivity : full per-metric × per-parameter table (saved as
#'                    morris_sensitivity_<TAG>.csv for analysis)
#'   - $ranking     : obj_total-only ranking (used by classify_morris and
#'                    by the run_cal.R narrowing logic; mu_med statistics
#'                    are only meaningful within one metric)
compute_sensitivity <- function(ee_all) {
  sensitivity <- ee_all %>%
    filter(!is.na(parameter)) %>%
    group_by(parameter, metric) %>%
    summarise(
      mu_star = mean(abs(EE), na.rm = TRUE),
      sigma   = sd(EE, na.rm = TRUE),
      n_ee    = sum(!is.na(EE)),
      .groups = "drop"
    ) %>%
    arrange(metric, desc(mu_star))

  ranking <- sensitivity %>%
    filter(metric == "obj_total") %>%
    arrange(desc(mu_star))

  list(sensitivity = sensitivity, ranking = ranking)
}
