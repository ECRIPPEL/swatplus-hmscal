# ==============================================================================
# check_range.R - Parameter range validation against cal_parms.cal (Iteration 1)
# ==============================================================================

library(dplyr)
library(tibble)

#' Extract change type from SWATrunR par_name
#' @param par_name SWATrunR string (e.g., "cn2::cn2.hru | change = pctchg")
#' @return String with change type (e.g., "pctchg", "absval")
extract_change_type <- function(par_name) {
  trimws(sub(".*change =\\s*", "", par_name))
}

#' Read and parse cal_parms.cal
#' @param path Path to cal_parms.cal
#' @return data.frame with name, obj_typ, abs_min, abs_max
read_cal_parms <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[-c(1, 2, 3)]

  parse_line <- function(x) {
    tok <- strsplit(trimws(x), "\\s+")[[1]]
    if (length(tok) < 4) return(NULL)
    data.frame(
      name    = tok[1],
      obj_typ = tok[2],
      abs_min = suppressWarnings(as.numeric(tok[3])),
      abs_max = suppressWarnings(as.numeric(tok[4])),
      stringsAsFactors = FALSE
    )
  }

  cal_tbl <- do.call(rbind, lapply(lines, parse_line))
  cal_tbl[!is.na(cal_tbl$abs_min) & !is.na(cal_tbl$abs_max), ]
}

#' Read soils.sol and return data.frame with layer-level variables
#' @param path Path to soils.sol
#' @return data.frame with soil variables per layer
read_soils_sol <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[-c(1, 2)]
  out <- list()
  i <- 1

  while (i <= length(lines)) {
    tk <- strsplit(trimws(lines[i]), "\\s+")[[1]]
    if (length(tk) >= 7) {
      soil_name <- tk[1]
      n_layers  <- as.integer(tk[2])
      if (is.na(n_layers) || n_layers < 1) { i <- i + 1; next }
      for (j in seq_len(n_layers)) {
        if ((i + j) > length(lines)) break
        x <- strsplit(trimws(lines[i + j]), "\\s+")[[1]]
        if (length(x) < 14) next
        out[[length(out) + 1]] <- data.frame(
          soil_name = soil_name,
          dp = as.numeric(x[1]), bd = as.numeric(x[2]), awc = as.numeric(x[3]),
          soil_k = as.numeric(x[4]), carbon = as.numeric(x[5]), clay = as.numeric(x[6]),
          silt = as.numeric(x[7]), sand = as.numeric(x[8]), rock = as.numeric(x[9]),
          alb = as.numeric(x[10]), usle_k = as.numeric(x[11]), ec = as.numeric(x[12]),
          caco3 = as.numeric(x[13]), ph = as.numeric(x[14]),
          stringsAsFactors = FALSE
        )
      }
      i <- i + n_layers + 1
    } else {
      i <- i + 1
    }
  }
  do.call(rbind, out)
}

#' Read generic SWAT+ table file (hydrology.hyd, topography.hyd)
#' @param path File path
#' @param col_names Character vector of column names (including "name")
#' @return data.frame
read_swat_table <- function(path, col_names) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  read.table(
    text = paste(lines[-c(1, 2)], collapse = "\n"),
    header = FALSE, stringsAsFactors = FALSE, col.names = col_names
  )
}

#' Get range (min, max) of a variable from SWAT+ input files
#' @param source_file Source file name (e.g., "soils.sol")
#' @param source_var Variable name
#' @param soil_tbl data.frame from read_soils_sol
#' @param hyd_tbl data.frame from hydrology.hyd
#' @param topo_tbl data.frame from topography.hyd
#' @return Numeric vector c(min, max)
get_input_range <- function(source_file, source_var, soil_tbl, hyd_tbl, topo_tbl) {
  if (source_file == "soils.sol") {
    vals <- soil_tbl[[source_var]]
  } else if (source_file == "hydrology.hyd") {
    vals <- hyd_tbl[[source_var]]
  } else if (source_file == "topography.hyd") {
    vals <- topo_tbl[[source_var]]
  } else {
    return(c(NA_real_, NA_real_))
  }
  if (is.null(vals)) return(c(NA_real_, NA_real_))
  vals <- suppressWarnings(as.numeric(vals))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(c(NA_real_, NA_real_))
  c(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE))
}

#' Run the full parameter range check against cal_parms.cal
#'
#' @param param_info Tibble with short, par_name, min, max
#' @param project_path SWAT+ project path
#' @param check_map List of mappings (short, source_file, source_var, cal_name, cal_obj_typ)
#' @return Tibble with check result per parameter. Stops with error if any out of range.
check_parameter_range <- function(param_info, project_path, check_map) {
  cat("\n[CHECK PARAMETER RANGE]\n")

  cal_tbl  <- read_cal_parms(file.path(project_path, "cal_parms.cal"))
  soil_tbl <- read_soils_sol(file.path(project_path, "soils.sol"))

  hyd_tbl <- read_swat_table(
    file.path(project_path, "hydrology.hyd"),
    c("name","lat_ttime","lat_sed","can_max","esco","epco","orgn_enrich","orgp_enrich",
      "cn3_swf","bio_mix","perco","lat_orgn","lat_orgp","pet_co","latq_co")
  )

  topo_tbl <- read_swat_table(
    file.path(project_path, "topography.hyd"),
    c("name","slp","slp_len","lat_len","dist_cha","depos")
  )

  check_map_df <- bind_rows(check_map)
  check_list <- vector("list", nrow(param_info))

  for (i in seq_len(nrow(param_info))) {
    short_i  <- param_info$short[i]
    rmin     <- param_info$min[i]
    rmax     <- param_info$max[i]
    change_i <- extract_change_type(param_info$par_name[i])

    if (short_i == "cn2") {
      check_list[[i]] <- tibble(parameter = short_i, result = "SKIPPED_CN2")
      next
    }

    map_row <- filter(check_map_df, short == short_i)
    if (nrow(map_row) == 0) {
      check_list[[i]] <- tibble(parameter = short_i, result = "SKIPPED_NO_MAPPING")
      next
    }

    cal_row <- filter(cal_tbl, name == map_row$cal_name[1], obj_typ == map_row$cal_obj_typ[1])
    if (nrow(cal_row) == 0) {
      check_list[[i]] <- tibble(parameter = short_i, result = "NOT_FOUND_IN_CAL")
      next
    }

    cal_min <- cal_row$abs_min[1]
    cal_max <- cal_row$abs_max[1]

    if (change_i == "absval") {
      final_min <- rmin
      final_max <- rmax
    } else if (change_i == "pctchg") {
      rng <- get_input_range(map_row$source_file[1], map_row$source_var[1], soil_tbl, hyd_tbl, topo_tbl)
      candidates <- c(rng[1] * (1 + rmin/100), rng[1] * (1 + rmax/100),
                       rng[2] * (1 + rmin/100), rng[2] * (1 + rmax/100))
      final_min <- min(candidates, na.rm = TRUE)
      final_max <- max(candidates, na.rm = TRUE)
    } else if (change_i == "abschg") {
      rng <- get_input_range(map_row$source_file[1], map_row$source_var[1], soil_tbl, hyd_tbl, topo_tbl)
      candidates <- c(rng[1] + rmin, rng[1] + rmax, rng[2] + rmin, rng[2] + rmax)
      final_min <- min(candidates, na.rm = TRUE)
      final_max <- max(candidates, na.rm = TRUE)
    } else {
      final_min <- NA_real_
      final_max <- NA_real_
    }

    result_i <- if (is.na(final_min) || is.na(final_max)) {
      "CHECK_NA"
    } else if (final_min >= cal_min && final_max <= cal_max) {
      "OKAY"
    } else if (final_min < cal_min && final_max > cal_max) {
      "OUT_BELOW_AND_ABOVE"
    } else if (final_min < cal_min) {
      "OUT_BELOW_MIN"
    } else if (final_max > cal_max) {
      "OUT_ABOVE_MAX"
    } else {
      "CHECK"
    }

    check_list[[i]] <- tibble(
      parameter = short_i, finalmin = final_min, finalmax = final_max,
      mincalib = cal_min, maxcalib = cal_max, result = result_i
    )
  }

  check_result <- bind_rows(check_list)
  print(check_result)

  failures <- filter(check_result, result %in% c("OUT_BELOW_MIN","OUT_ABOVE_MAX","OUT_BELOW_AND_ABOVE","CHECK_NA"))
  if (nrow(failures) > 0) {
    cat("\nPARAMETERS OUT OF RANGE:\n")
    print(failures)
    stop("ABORTING: some parameters exceed the cal_parms.cal bounds")
  } else {
    cat("[CHECK PARAMETER RANGE] OK\n")
  }

  check_result
}
