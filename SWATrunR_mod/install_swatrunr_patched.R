# =============================================================================
# install_swatrunr_patched.R
#
# One-command install of the SWATrunR build that HMSCal needs.
#
#   source("SWATrunR_mod/install_swatrunr_patched.R")
#   install_swatrunr_patched()
#
# It downloads SWATrunR 0.9.4 at the pinned upstream commit, applies the
# read_output_i() header guard, installs it, and verifies the result.
#
# Requirements: R (tested on 4.3.3) and an internet connection.
# NOT required: git, RTools, or any compiler. SWATrunR is a pure R package
# (no src/ directory), so nothing is compiled.
#
# See README.md in this folder for what the patch does and why.
# =============================================================================

SWATRUNR_COMMIT <- "82e6353001480027303b230adf57bf5c7ba7057e"   # v0.9.4
READR_SAFE      <- "2.1.5"   # last 2.1.x; 2.2.0 removed read_table2()

# --- the patch, as a line-level edit -----------------------------------------
# Matched on the original source lines rather than applied as a diff, so it
# works regardless of CRLF/LF and needs no git or patch utility.

.anchor_old <- c(
  "read_output_i <- function(output_i, col_names_i, thread_path, date_cols) {",
  "  fread(thread_path%//%output_i$file_full[1], skip = 3) %>%",
  "    as_tibble(.) %>%",
  "    .[,1:length(col_names_i)] %>%"
)

.anchor_new <- c(
  "read_output_i <- function(output_i, col_names_i, thread_path, date_cols) {",
  "  raw_data <- fread(thread_path%//%output_i$file_full[1], skip = 3) %>%",
  "    as_tibble(.)",
  "",
  "  # Guard against SWAT+ writer defects where the header row lists more column",
  "  # names than the data rows actually contain (observed in rev. 2026.61.0.2.61",
  "  # for lsunit_wb_aa.txt: 51 names in header, 49 numeric columns per row --",
  "  # plant_cov and mgt_ops are listed but never written). Without this guard",
  "  # `[, 1:length(col_names_i)]` aborts with \"Locations N and N+1 don't exist\"",
  "  # and the whole run_swatplus() call returns NULL after every thread finished.",
  "  n_data <- ncol(raw_data)",
  "  if (length(col_names_i) > n_data) {",
  "    dropped <- col_names_i[(n_data + 1):length(col_names_i)]",
  "    warning(sprintf(",
  "      \"[SWATrunR] '%s': header has %d cols but data has %d. Truncating header (dropping: %s).\",",
  "      output_i$file_full[1], length(col_names_i), n_data,",
  "      paste(dropped, collapse = \", \")), call. = FALSE)",
  "    col_names_i <- col_names_i[seq_len(n_data)]",
  "  }",
  "",
  "  raw_data %>%",
  "    .[,1:length(col_names_i)] %>%"
)

.find_block <- function(lines, block) {
  n <- length(block)
  hits <- which(trimws(lines) == trimws(block[1]))
  for (h in hits) {
    if (h + n - 1 <= length(lines) &&
        all(trimws(lines[h:(h + n - 1)]) == trimws(block))) {
      return(h)
    }
  }
  NA_integer_
}

install_swatrunr_patched <- function(commit = SWATRUNR_COMMIT,
                                     lib    = .libPaths()[1],
                                     pin_readr = TRUE) {

  msg <- function(...) cat("[swatrunr-patch]", ..., "\n")

  # --- 1. readr must stay below 2.2.0 ----------------------------------------
  if (pin_readr) {
    has_readr <- requireNamespace("readr", quietly = TRUE)
    v <- if (has_readr) utils::packageVersion("readr") else NULL
    if (!has_readr || v >= package_version("2.2.0")) {
      msg("readr", if (has_readr) as.character(v) else "absent",
          "-> installing", READR_SAFE, "(2.2.0 removed read_table2)")
      if (!requireNamespace("remotes", quietly = TRUE)) {
        utils::install.packages("remotes")
      }
      remotes::install_version("readr", version = READR_SAFE,
                               upgrade = "never", lib = lib)
    } else {
      msg("readr", as.character(v), "is below 2.2.0, keeping it")
    }
  }

  # --- 2. fetch the pinned upstream source -----------------------------------
  tmp <- file.path(tempdir(), paste0("swatrunr-", substr(commit, 1, 7)))
  unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)
  tarball <- file.path(tmp, "src.tar.gz")
  url <- sprintf("https://github.com/chrisschuerz/SWATrunR/archive/%s.tar.gz", commit)

  msg("downloading SWATrunR at", substr(commit, 1, 7))
  utils::download.file(url, tarball, mode = "wb", quiet = TRUE)
  utils::untar(tarball, exdir = tmp)

  srcdir <- list.dirs(tmp, recursive = FALSE)
  srcdir <- srcdir[grepl("SWATrunR", basename(srcdir))]
  if (length(srcdir) != 1) stop("could not locate the extracted source tree")

  # --- 3. apply the guard ----------------------------------------------------
  target <- file.path(srcdir, "R", "swatplus_handle_output.R")
  if (!file.exists(target)) stop("missing ", target)
  lines <- readLines(target, warn = FALSE)

  if (any(grepl("Truncating header", lines, fixed = TRUE))) {
    msg("source already carries the guard, skipping patch")
  } else {
    at <- .find_block(lines, .anchor_old)
    if (is.na(at)) {
      stop("could not find read_output_i() in its expected form. ",
           "The pinned commit may be wrong, or upstream changed. ",
           "Apply 01-read_output_i-header-guard.patch by hand instead.")
    }
    lines <- append(lines[-(at:(at + length(.anchor_old) - 1))],
                    .anchor_new, after = at - 1)
    writeLines(lines, target)
    msg("guard inserted at line", at, "of R/swatplus_handle_output.R")
  }

  # the patched file must still parse
  tryCatch(parse(target), error = function(e)
    stop("patched file does not parse: ", conditionMessage(e)))

  # --- 4. install ------------------------------------------------------------
  msg("installing into", lib)
  utils::install.packages(srcdir, repos = NULL, type = "source", lib = lib)

  # --- 5. verify the INSTALLED package, not the source -----------------------
  pkg <- file.path(lib, "SWATrunR")
  if (!dir.exists(pkg)) stop("install failed: ", pkg, " not found")
  e <- new.env()
  lazyLoad(file.path(pkg, "R", "SWATrunR"), envir = e)
  ok <- any(grepl("Truncating header",
                  paste(deparse(e$read_output_i), collapse = " "), fixed = TRUE))
  if (!ok) stop("installed package does NOT contain the guard")

  msg("OK. SWATrunR", as.character(utils::packageVersion("SWATrunR", lib.loc = lib)),
      "installed with the header guard.")
  msg("Restart R before running the pipeline.")
  invisible(TRUE)
}
