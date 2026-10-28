# ==============================================================================
# param_inherit.R — Automatic discovery of inherited parameter ranges
# ==============================================================================

library(dplyr)

# Default-on-NULL helper (in case loaded standalone)
if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a)) a else b

#' Find the inherited ranges CSV from a previous iteration
#'
#' Looks under <base_path>/<from_tag_prefix>_v<N>/ for the file
#' inherited_ranges_<from_tag_prefix>_v<N>.csv. Picks the highest version
#' that contains the file. Filename is fully determined by the folder TAG —
#' there is no need (and no way) to override it from the YAML.
#'
#' @param base_path SWAT+ project base directory (e.g., "C:/path/to/swatplus_project")
#' @param from_tag_prefix Folder TAG prefix of the source iteration (e.g., "ITER1")
#' @return List with $path (CSV path), $folder, $tag, $data (tibble)
find_inherited_ranges <- function(base_path, from_tag_prefix) {
  folders <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
  pattern <- paste0("[/\\\\]", from_tag_prefix, "_v[0-9]+$")
  folders <- folders[grepl(pattern, folders)]

  if (length(folders) == 0) {
    stop("No folder ", from_tag_prefix, "_v* found in: ", base_path)
  }

  versions <- suppressWarnings(
    as.numeric(sub(paste0(".*", from_tag_prefix, "_v"), "", basename(folders)))
  )
  folders_sorted <- folders[order(versions, decreasing = TRUE, na.last = NA)]

  for (folder_i in folders_sorted) {
    tag       <- basename(folder_i)
    candidate <- file.path(folder_i, paste0("inherited_ranges_", tag, ".csv"))
    if (file.exists(candidate)) {
      data <- read.csv(candidate, stringsAsFactors = FALSE)

      stopifnot(!any(is.na(data$new_min)))
      stopifnot(!any(is.na(data$new_max)))
      stopifnot(all(data$new_min < data$new_max))

      cat("Inherited ranges from: ", folder_i, "\n", sep = "")
      cat("File: ", candidate, "\n", sep = "")
      cat("Parameters: ", nrow(data), "\n", sep = "")

      return(list(path = candidate, folder = folder_i, tag = tag, data = data))
    }
  }

  stop("File inherited_ranges_<", from_tag_prefix, "_v*>.csv not found in any ",
       from_tag_prefix, "_v* folder. ",
       "Check whether the source CAL run completed and produced an ",
       "inherited_ranges file (it is skipped only if cal.inheritance_source = ",
       "'satisfactory' and 0 satisfactory runs were found — the producer ",
       "stops loud in that case).")
}

#' Find the latest Morris screening folder (param_info + morris_ranking) for a given tag
#'
#' @param base_path SWAT+ project base directory
#' @param morris_tag_prefix Morris screening folder TAG prefix (e.g., "MORRIS1", "MORRIS2", "MORRIS3")
#' @return List with $folder, $tag, $param_info, $ranking, $metrics_file
find_morris_results <- function(base_path, morris_tag_prefix) {
  folders <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
  pattern <- paste0("[/\\\\]", morris_tag_prefix, "_v[0-9]+$")
  folders <- folders[grepl(pattern, folders)]

  if (length(folders) == 0) {
    stop("No folder ", morris_tag_prefix, "_v* found in: ", base_path)
  }

  versions <- suppressWarnings(
    as.numeric(sub(paste0(".*", morris_tag_prefix, "_v"), "", basename(folders)))
  )
  folders_sorted <- folders[order(versions, decreasing = TRUE, na.last = NA)]

  for (folder_i in folders_sorted) {
    tag       <- basename(folder_i)
    f_param   <- file.path(folder_i, paste0("param_info_",     tag, ".csv"))
    f_ranking <- file.path(folder_i, paste0("morris_ranking_", tag, ".csv"))
    f_metrics <- file.path(folder_i, paste0("morris_metrics_", tag, ".csv"))

    if (file.exists(f_param) && file.exists(f_ranking)) {
      param_info <- read.csv(f_param,   stringsAsFactors = FALSE)
      ranking    <- read.csv(f_ranking, stringsAsFactors = FALSE)

      stopifnot(all(c("short", "par_name", "min", "max") %in% names(param_info)))
      stopifnot(all(param_info$min < param_info$max))

      cat("Morris screening results from: ", folder_i, "\n", sep = "")

      return(list(
        folder       = folder_i,
        tag          = tag,
        param_info   = param_info,
        ranking      = ranking,
        metrics_file = if (file.exists(f_metrics)) f_metrics else NA_character_
      ))
    }
  }

  stop("No param_info_<TAG>.csv + morris_ranking_<TAG>.csv pair found in any ",
       morris_tag_prefix, "_v* folder")
}

#' Compute auto-drop list from a Morris ranking
#'
#' Returns the set of parameter shorts whose Morris mu*(obj_total) is at or
#' below the threshold, excluding any in `protect`. Designed for use in the
#' CAL stage to skip parameters that the Morris screening stage flagged as inert (μ*≈0)
#' — typically due to (a) physical insensitivity in the catchment regime,
#' (b) silent failures of SWATrunR's parameter passthrough, or (c) SWAT+
#' initialization-order issues where post-init calibration writes are
#' clobbered by precomputed caches (e.g., channel rating curves built
#' before proc_cal applies calibration.cal in SWAT+ rev 61.0.2).
#'
#' Note: the Morris screening Morris figure (`morris_screening_<TAG>.tif`) is generated
#' BEFORE this filter, intentionally — the figure is the EVIDENCE that
#' justifies the drop and belongs in any reproducibility audit.
#'
#' @param ranking_morris Tibble with columns `parameter` and `mu_star` (and
#'             optionally `metric`). If `metric` exists, only `obj_total`
#'             rows are used.
#' @param threshold Numeric: drop params with `mu_star <= threshold`.
#'             Set to 0 to drop only literal zeros; set higher (e.g. 0.001)
#'             to also drop physically irrelevant params.
#' @param protect Optional character vector of shorts to NEVER auto-drop
#'             (safety net for params known to be important even if Morris
#'             estimates them as low due to noise from small m_traj).
#' @return List with `$drop` (character vector of shorts to drop) and
#'             `$detail` (tibble of dropped params with their mu* values
#'             and a `zero_flag` column for literal zeros).
compute_auto_drop <- function(ranking_morris,
                              threshold = NULL,
                              protect   = character(0)) {
  if (is.null(threshold)) {
    return(list(drop = character(0),
                detail = tibble::tibble(parameter = character(),
                                         mu_star = numeric(),
                                         zero_flag = logical())))
  }

  ranking <- tibble::as_tibble(ranking_morris)
  if ("metric" %in% names(ranking)) {
    ranking <- ranking %>% dplyr::filter(metric == "obj_total")
  }

  detail <- ranking %>%
    dplyr::filter(mu_star <= threshold,
                  !parameter %in% protect) %>%
    dplyr::transmute(parameter,
                     mu_star,
                     zero_flag = (mu_star == 0))

  list(drop = detail$parameter, detail = detail)
}

#' Build param_info from inherited ranges + parameter map + drop/override/extras
#'
#' Order of operations: load inherited → drop (from inherited AND extras) →
#' override → add remaining extras.
#'
#' @param ranges_data Tibble of inherited ranges (columns: parameter, new_min, new_max)
#' @param param_map Named vector (short -> SWATrunR par_name)
#' @param extra_params Optional tibble of extra parameters (with $short, $par_name, $min, $max)
#' @param drop Optional character vector of `short` names to remove (applies to BOTH
#'             inherited and extras)
#' @param override_ranges Optional named list: list(short = list(min=, max=), ...)
#' @param strict_drop Logical: if TRUE, stop() when drop references unknown shorts;
#'             if FALSE, silently skip absent shorts (used for auto-drop where
#'             params may already be absent due to upstream drops)
#' @return Tibble param_info ready for use
build_param_info <- function(ranges_data, param_map, extra_params = NULL,
                             drop = NULL, override_ranges = NULL,
                             strict_drop = TRUE) {
  param_info <- ranges_data %>%
    mutate(
      par_name = unname(param_map[parameter]),
      short    = parameter,
      min      = new_min,
      max      = new_max
    ) %>%
    filter(!is.na(par_name)) %>%
    select(short, par_name, min, max,
           any_of(c("grupo", "q_low", "q_high", "post_min", "post_max", "mediana")))

  # ---- (1) DROP — remove from inherited AND from extras ---------------------
  if (!is.null(drop) && length(drop) > 0) {
    inh_shorts   <- param_info$short
    extra_shorts <- if (!is.null(extra_params)) extra_params$short else character()
    all_shorts   <- c(inh_shorts, extra_shorts)
    missing_drop <- setdiff(drop, all_shorts)

    if (length(missing_drop) > 0) {
      if (strict_drop) {
        stop("drop references parameters not in inherited set or extras: ",
             paste(missing_drop, collapse = ", "),
             "\nAvailable shorts: ",
             paste(all_shorts, collapse = ", "))
      } else {
        cat("[drop] note: ", length(missing_drop),
            " param(s) already absent (no-op): ",
            paste(missing_drop, collapse = ", "), "\n", sep = "")
      }
    }

    drop_present <- intersect(drop, all_shorts)
    param_info <- param_info %>% filter(!short %in% drop_present)
    if (!is.null(extra_params)) {
      extra_params <- extra_params %>% filter(!short %in% drop_present)
    }
    if (length(drop_present) > 0) {
      cat("[drop] removed (", length(drop_present), "): ",
          paste(drop_present, collapse = ", "), "\n", sep = "")
    }
  }

  # ---- (2) OVERRIDE — replace min/max for inherited and/or extra params -----
  if (!is.null(override_ranges) && length(override_ranges) > 0) {
    inh_shorts   <- param_info$short
    extra_shorts <- if (!is.null(extra_params)) extra_params$short else character()

    for (sh in names(override_ranges)) {
      if (!sh %in% c(inh_shorts, extra_shorts)) {
        stop("inherit.override_ranges references parameter '", sh,
             "' not in inherited set or extra_params. ",
             "Available shorts: ",
             paste(c(inh_shorts, extra_shorts), collapse = ", "))
      }
      ov <- override_ranges[[sh]]
      if (is.null(ov$min) || is.null(ov$max) || ov$min >= ov$max) {
        stop("inherit.override_ranges['", sh, "'] requires numeric min < max ",
             "(got min=", ov$min %||% "NULL", ", max=", ov$max %||% "NULL", ").")
      }
      if (sh %in% inh_shorts) {
        param_info$min[param_info$short == sh] <- ov$min
        param_info$max[param_info$short == sh] <- ov$max
      }
      if (sh %in% extra_shorts) {
        extra_params$min[extra_params$short == sh] <- ov$min
        extra_params$max[extra_params$short == sh] <- ov$max
      }
      cat("[inherit.override_ranges] ", sh, " -> [",
          ov$min, ", ", ov$max, "]\n", sep = "")
    }
  }

  # ---- (3) ADD — append extra_params (already overridden if applicable) -----
  if (!is.null(extra_params) && nrow(extra_params) > 0) {
    extra_params <- extra_params %>% filter(!short %in% param_info$short)
    if (nrow(extra_params) > 0) {
      param_info <- bind_rows(param_info, extra_params)
      cat("[extra_params] added (", nrow(extra_params), "): ",
          paste(extra_params$short, collapse = ", "), "\n", sep = "")
    }
  }

  stopifnot(length(unique(param_info$short)) == nrow(param_info))
  stopifnot(all(param_info$min < param_info$max))

  param_info
}
