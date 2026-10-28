# ==============================================================================
# versioning.R — Auto-versioned output folder management
# ==============================================================================

#' Create an auto-versioned folder (e.g., MORRIS1_v1, v2, v3...)
#'
#' @param base_path Base directory (e.g., "C:/path/to/swatplus_project")
#' @param prefix Version prefix (e.g., "MORRIS1_v")
#' @return List with $tag (e.g., "MORRIS1_v3") and $path (full path)
create_versioned_folder <- function(base_path, prefix) {
  folders <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)

  versions <- suppressWarnings(
    as.numeric(gsub(prefix, "", folders[grep(paste0("^", prefix), folders)]))
  )

  next_v <- if (length(versions) == 0 || all(is.na(versions))) {
    1
  } else {
    max(versions, na.rm = TRUE) + 1
  }

  tag  <- paste0(prefix, next_v)
  path <- file.path(base_path, tag)

  if (!dir.exists(path)) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
  }

  list(tag = tag, path = path)
}
