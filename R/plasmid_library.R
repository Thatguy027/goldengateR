# Plasmid library: bundled GenBank/FASTA files accessible by name.
#
# Plasmids live in two places (searched in order):
#   1. inst/extdata/plasmids/  — bundled with the package, installed alongside it
#   2. options("goldengateR.plasmid_dir")  — optional user-configured directory,
#      useful for assembled plasmids or personal lab collections
#
# Bundled plasmids take precedence: if a name exists in both places the
# package copy wins. This keeps pYTK parts stable and predictable.
#
# Files can be .gb, .gbk, .genbank, .fa, .fna, or .fasta. The plasmid "name"
# is the filename without the extension. So a file named pYTK001.gb is
# loaded with load_plasmid("pYTK001").

#' Get the search path for plasmid lookups
#'
#' Returns a character vector of directories that will be searched for
#' plasmid files, in order. The bundled package library is always searched
#' first; the user-configured directory (if set) is appended after.
#'
#' @return Character vector of directory paths.
#' @export
plasmid_search_path <- function() {
  paths <- character()

  # Bundled package plasmids first
  installed_dir <- tryCatch(
    system.file("extdata", "plasmids", package = "goldengateR"),
    error = function(e) ""
  )
  if (nzchar(installed_dir) && dir.exists(installed_dir)) {
    paths <- c(paths, installed_dir)
  } else {
    # Development fallback: try to locate inst/extdata/plasmids relative to
    # this very source file. Works whether sourced via source() or
    # devtools::load_all().
    dev_candidates <- character()

    # Strategy 1: this file is in <pkg>/R/plasmid_library.R, so go up two.
    # We can't reliably get __file__, but we can search the call stack for
    # a sourced filename.
    for (i in seq_len(sys.nframe())) {
      sf <- sys.frame(i)
      ofile <- attr(sf, "srcfile")
      if (!is.null(ofile) && !is.null(ofile$filename)) {
        candidate <- normalizePath(
          file.path(dirname(ofile$filename), "..", "inst", "extdata", "plasmids"),
          mustWork = FALSE
        )
        if (dir.exists(candidate)) dev_candidates <- c(dev_candidates, candidate)
      }
    }

    # Strategy 2: package-dir option (set by users in dev workflow)
    pkg_dir_opt <- getOption("goldengateR.pkg_dir", default = NULL)
    if (!is.null(pkg_dir_opt)) {
      candidate <- file.path(path.expand(pkg_dir_opt), "inst", "extdata", "plasmids")
      if (dir.exists(candidate)) dev_candidates <- c(dev_candidates, candidate)
    }

    # Strategy 3: walk up from getwd() looking for inst/extdata/plasmids
    wd <- getwd()
    for (up in 0:3) {
      candidate <- file.path(wd, paste(rep("..", up), collapse = "/"),
                             "inst", "extdata", "plasmids")
      candidate <- normalizePath(candidate, mustWork = FALSE)
      if (dir.exists(candidate)) dev_candidates <- c(dev_candidates, candidate)
    }

    paths <- c(paths, unique(dev_candidates))
  }

  # User-configured directory appended last (supplemental, not override)
  user_dir <- getOption("goldengateR.plasmid_dir", default = NULL)
  if (!is.null(user_dir) && nzchar(user_dir)) {
    user_dir <- path.expand(user_dir)
    if (dir.exists(user_dir)) paths <- c(paths, user_dir)
  }

  unique(paths)
}

#' List bundled and user-library plasmids
#'
#' Returns a data.frame of all plasmid files visible on the search path.
#' Columns: name (filename without extension), file (basename), path
#' (full path), source ("user" if from user-configured dir, "bundled" if
#' from the package).
#'
#' @param pattern Optional regex to filter names (case-insensitive).
#' @return A data.frame.
#' @export
list_plasmids <- function(pattern = NULL) {
  paths <- plasmid_search_path()
  if (!length(paths)) {
    return(data.frame(name = character(), file = character(),
                      path = character(), source = character(),
                      stringsAsFactors = FALSE))
  }

  user_dir <- getOption("goldengateR.plasmid_dir", default = NULL)
  if (!is.null(user_dir)) user_dir <- path.expand(user_dir)

  rows <- list()
  for (p in paths) {
    files <- list.files(p, pattern = "\\.(gb|gbk|genbank|fa|fna|fasta)$",
                        full.names = TRUE, ignore.case = TRUE)
    if (!length(files)) next
    src <- if (!is.null(user_dir) && normalizePath(p, mustWork = FALSE) ==
                                     normalizePath(user_dir, mustWork = FALSE)) "user" else "bundled"
    for (f in files) {
      name <- tools::file_path_sans_ext(basename(f))
      rows[[length(rows) + 1]] <- data.frame(
        name = name, file = basename(f), path = f, source = src,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(data.frame(name = character(), file = character(),
                      path = character(), source = character(),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  # Earlier paths take precedence: drop duplicates by name keeping the first
  out <- out[!duplicated(out$name), , drop = FALSE]
  if (!is.null(pattern)) {
    out <- out[grepl(pattern, out$name, ignore.case = TRUE), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Load a plasmid by name from the library
#'
#' Looks up a plasmid by its filename (without extension) in the
#' \code{plasmid_search_path()}. If multiple matches exist across
#' directories, the first one (highest-priority) wins.
#'
#' @param name Plasmid name, e.g. "pYTK001". Case-insensitive.
#' @param topology Default topology if loading from FASTA (ignored for GenBank).
#' @return A ggr_part.
#' @export
#' @examples
#' \dontrun{
#'   vec <- load_plasmid("pYTK001")
#'   list_plasmids(pattern = "^pYTK0")
#' }
load_plasmid <- function(name, topology = "circular") {
  available <- list_plasmids()
  if (!nrow(available)) {
    stop("No plasmids found on search path. ",
         "Either install bundled plasmids or set ",
         "options(goldengateR.plasmid_dir = '/path/to/dir').")
  }
  match_idx <- which(tolower(available$name) == tolower(name))
  if (!length(match_idx)) {
    stop("Plasmid '", name, "' not found. Available plasmids: ",
         paste(head(available$name, 20), collapse = ", "),
         if (nrow(available) > 20) sprintf(" ... (%d total)", nrow(available)) else "")
  }
  read_part(available$path[match_idx[1]], topology = topology)
}

#' Set the user plasmid directory
#'
#' Convenience wrapper for \code{options(goldengateR.plasmid_dir = dir)}.
#' Once set, plasmids in this directory are findable by \code{load_plasmid()}.
#' The bundled package library is always searched first; this directory is
#' searched second, so it is additive rather than an override.
#'
#' @param dir Path to a directory containing GenBank/FASTA plasmid files,
#'   or NULL to clear.
#' @export
set_plasmid_dir <- function(dir) {
  if (is.null(dir)) {
    options(goldengateR.plasmid_dir = NULL)
    message("Cleared user plasmid directory.")
    return(invisible())
  }
  dir <- path.expand(dir)
  if (!dir.exists(dir)) stop("Directory does not exist: ", dir)
  options(goldengateR.plasmid_dir = dir)
  n <- nrow(list_plasmids())
  message(sprintf("User plasmid directory set to %s (%d plasmid files visible).",
                  dir, n))
  invisible(dir)
}
