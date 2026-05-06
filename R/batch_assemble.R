PART_COLS <- c("vec", "1", "2", "3", "3a", "3b", "234",
               "4", "4a", "4b", "5", "6", "7", "8", "8a", "8b", "678")

.load_parts_from_row <- function(row) {
  names_to_load <- character()
  for (col in PART_COLS) {
    val <- row[[col]]
    if (!is.null(val) && !is.na(val) && nzchar(trimws(as.character(val)))) {
      names_to_load <- c(names_to_load, trimws(as.character(val)))
    }
  }
  if (!length(names_to_load)) stop("No parts specified")
  lapply(names_to_load, load_plasmid)
}

#' Interactively assemble a plasmid from the console
#'
#' Prompts for an output plasmid name and a space-separated list of part
#' names, loads each via \code{load_plasmid()}, and runs
#' \code{golden_gate_assemble()}.
#'
#' @param enzyme "BsaI" or "BsmBI".
#' @param output_dir Directory to write the output GenBank file.
#' @return Invisibly, the path to the written GenBank file.
#' @export
assemble_interactive <- function(enzyme = c("BsaI", "BsmBI"), output_dir = ".") {
  enzyme <- match.arg(enzyme)

  pname <- trimws(readline("Output plasmid name: "))
  if (!nzchar(pname)) stop("No plasmid name supplied.")

  parts_raw <- trimws(readline("Parts (space-separated): "))
  if (!nzchar(parts_raw)) stop("No parts supplied.")
  part_names <- strsplit(parts_raw, "\\s+")[[1]]

  parts <- lapply(part_names, load_plasmid)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  outfile <- file.path(output_dir, paste0(pname, ".gb"))
  golden_gate_assemble(parts, enzyme = enzyme, name = pname, output = outfile)
  invisible(outfile)
}

#' Assemble multiple plasmids from a manifest table
#'
#' Reads a CSV or TSV where each row describes one assembly. A required
#' \code{name} column provides the output plasmid name. Part columns
#' (\code{vec}, \code{1}, \code{2}, \code{3}, \code{3a}, \code{3b},
#' \code{234}, \code{4}, \code{4a}, \code{4b}, \code{5}, \code{6},
#' \code{7}, \code{8}, \code{8a}, \code{8b}, \code{678}) hold plasmid
#' names resolvable via \code{load_plasmid()}. Empty cells are skipped.
#' Rows with a blank \code{name} are skipped. Tab-separated files
#' (\code{.tsv}) are detected by extension; everything else is treated as
#' comma-separated.
#'
#' @param file Path to the manifest CSV or TSV.
#' @param output_dir Directory for output GenBank files (created if absent).
#' @param enzyme "BsaI" or "BsmBI".
#' @return Invisibly, a data.frame with columns \code{name}, \code{status},
#'   \code{output}, \code{error} — one row per assembly attempted.
#' @export
#' @examples
#' \dontrun{
#'   results <- assemble_from_table("assemblies.tsv", output_dir = "./plasmids")
#'   print(results)
#' }
assemble_from_table <- function(file, output_dir = ".", enzyme = c("BsaI", "BsmBI")) {
  enzyme <- match.arg(enzyme)
  if (!file.exists(file)) stop("File not found: ", file)

  sep <- if (tolower(tools::file_ext(file)) == "tsv") "\t" else ","
  tbl <- read.csv(file, sep = sep, stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = c("", "NA"))

  if (!"name" %in% colnames(tbl)) {
    stop("Manifest must have a 'name' column for the output plasmid name.")
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  results <- vector("list", nrow(tbl))
  for (i in seq_len(nrow(tbl))) {
    row   <- as.list(tbl[i, , drop = FALSE])
    pname <- trimws(as.character(row[["name"]]))

    if (!nzchar(pname) || is.na(pname)) {
      results[[i]] <- data.frame(name = paste0("row_", i), status = "SKIPPED",
                                 output = NA_character_, error = "blank name",
                                 stringsAsFactors = FALSE)
      next
    }

    outfile <- file.path(output_dir, paste0(pname, ".gb"))
    tryCatch({
      parts <- .load_parts_from_row(row)
      golden_gate_assemble(parts, enzyme = enzyme, name = pname, output = outfile)
      message(sprintf("[%d/%d] OK      %s  ->  %s", i, nrow(tbl), pname, outfile))
      results[[i]] <- data.frame(name = pname, status = "OK",
                                 output = outfile, error = NA_character_,
                                 stringsAsFactors = FALSE)
    }, error = function(e) {
      message(sprintf("[%d/%d] FAILED  %s  —  %s", i, nrow(tbl), pname, conditionMessage(e)))
      results[[i]] <<- data.frame(name = pname, status = "FAILED",
                                  output = NA_character_, error = conditionMessage(e),
                                  stringsAsFactors = FALSE)
    })
  }

  out <- do.call(rbind, results)
  n_ok     <- sum(out$status == "OK")
  n_failed <- sum(out$status == "FAILED")
  message(sprintf("\nDone: %d OK, %d failed, %d skipped",
                  n_ok, n_failed, nrow(out) - n_ok - n_failed))
  invisible(out)
}
