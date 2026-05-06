#!/usr/bin/env Rscript
# build_entry_plasmids.R
#
# Reads gg_parts_summary.tsv and ligates each PCR product into a fresh
# copy of a single hard-coded Type 3 destination vector, producing one
# entry plasmid per row. Each output is written as a GenBank file and a
# summary table is written at the end.
#
# Workflow per row:
#   1. Read the destination vector once (hard-coded path below)
#   2. Digest + gel-purify the BACKBONE fragment (largest band, contains
#      CamR + ColE1; the GFP dropout is filtered/discarded)
#   3. Build a linear PCR product ggr_part from the `pcr_product` column
#   4. Digest + gel-purify the released CDS fragment (largest band; the
#      two primer-tail edge fragments fall off because they have no
#      overhang on their outer end)
#   5. ligate(backbone, insert) -> circular Type 3 entry plasmid
#   6. Write GenBank, record outcome
#
# Usage:
#   Rscript build_entry_plasmids.R                          # uses defaults
#   Rscript build_entry_plasmids.R input.tsv output_dir     # override paths
#
# Required TSV columns: gene_name, species, accession, fwd_primer,
#   rev_primer, codon_optimized_sequence, pcr_product

# =============================================================================
# Configuration  --- edit these to match your setup
# =============================================================================
GOLDENGATE_PKG_DIR <- "~/goldengateR"        # directory containing R/, DESCRIPTION, etc.
DESTINATION_PLASMID <- "pYTK_T3_dest_demo"   # name of plasmid in the bundled library
                                             # (or in your set_plasmid_dir() directory)
USER_PLASMID_DIR   <- NULL                   # optional: set to a path with your
                                             # personal lab plasmid collection
ENZYME             <- "BsaI"

# =============================================================================
# Setup
# =============================================================================

# Allow CLI overrides:  Rscript build_entry_plasmids.R [tsv] [outdir]
args <- commandArgs(trailingOnly = TRUE)
tsv_path   <- if (length(args) >= 1) args[1] else "gg_parts_summary.tsv"
output_dir <- if (length(args) >= 2) args[2] else "entry_plasmids"

GOLDENGATE_PKG_DIR <- path.expand(GOLDENGATE_PKG_DIR)

# Load the package (sourcing R/ files; swap for library(goldengateR) if installed)
pkg_R <- file.path(GOLDENGATE_PKG_DIR, "R")
if (!dir.exists(pkg_R)) {
  stop("goldengateR R/ directory not found at: ", pkg_R,
       "\n  Edit GOLDENGATE_PKG_DIR at the top of this script.")
}
for (f in list.files(pkg_R, pattern = "\\.R$", full.names = TRUE)) source(f)

# Tell the dev-mode plasmid search where the package lives, so bundled
# plasmids work regardless of the current working directory.
options(goldengateR.pkg_dir = GOLDENGATE_PKG_DIR)

# Optionally point at a user plasmid collection
if (!is.null(USER_PLASMID_DIR) && nzchar(USER_PLASMID_DIR)) {
  set_plasmid_dir(USER_PLASMID_DIR)
}

if (!file.exists(tsv_path)) {
  stop("Parts TSV not found at: ", tsv_path)
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Read inputs
# =============================================================================
parts_df <- read.table(tsv_path, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "",
                       comment.char = "", check.names = FALSE)

required_cols <- c("gene_name", "species", "accession",
                   "fwd_primer", "rev_primer",
                   "codon_optimized_sequence", "pcr_product")
missing_cols <- setdiff(required_cols, names(parts_df))
if (length(missing_cols)) {
  stop("TSV is missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat(sprintf("Loaded %d parts from %s\n", nrow(parts_df), tsv_path))

# Read and digest the destination vector ONCE — every clone uses the same
# gel-purified backbone fragment, no point repeating the work N times.
dest_vec <- tryCatch(
  load_plasmid(DESTINATION_PLASMID),
  error = function(e) {
    stop("Could not load destination plasmid '", DESTINATION_PLASMID, "':\n  ",
         conditionMessage(e), "\n",
         "  Available plasmids: ",
         paste(list_plasmids()$name, collapse = ", "))
  }
)
cat(sprintf("Destination plasmid: %s (%d bp %s)\n",
            dest_vec$name, nchar(dest_vec$sequence), dest_vec$topology))

backbone <- digest_and_select(dest_vec, ENZYME, select = "largest")
cat(sprintf("Backbone fragment: %d bp, overhangs L=%s R=%s\n\n",
            nchar(backbone$sequence),
            attr(backbone, "left_oh"),
            attr(backbone, "right_oh")))

# =============================================================================
# Per-row ligation
# =============================================================================
results <- data.frame(
  gene_name      = character(),
  species        = character(),
  accession      = character(),
  status         = character(),
  product_bp     = integer(),
  output_file    = character(),
  message        = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(parts_df))) {
  row <- parts_df[i, ]
  gene <- row$gene_name
  cat(sprintf("[%d/%d] %s ... ", i, nrow(parts_df), gene))

  # Sanitize gene name for filesystem
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", gene)
  out_file <- file.path(output_dir, sprintf("pYTK_T3_%s.gb", safe_name))

  result <- tryCatch({
    pcr_seq <- toupper(gsub("\\s+", "", row$pcr_product))
    if (!nchar(pcr_seq)) stop("empty pcr_product")
    if (!grepl("^[ACGTN]+$", pcr_seq)) stop("pcr_product has non-ACGTN characters")

    amp <- new_part(sequence = pcr_seq,
                    topology = "linear",
                    name = paste0(safe_name, "_amplicon"))

    # The middle fragment (released CDS) is the LARGEST after digestion
    # because the primer-tail edge fragments are tiny.
    insert <- digest_and_select(amp, ENZYME, select = "largest")

    # Capture warnings without aborting. Common case: the assembled product
    # contains 2 internal BsaI sites — this is EXPECTED for a Type 3 entry
    # plasmid (those are the cassette flanks that release the CDS in the
    # next Golden Gate round) and is not a failure.
    captured_warnings <- character()
    asm <- withCallingHandlers(
      ligate(
        fragments = list(backbone, insert),
        enzyme    = ENZYME,
        name      = sprintf("pYTK_T3_%s", safe_name),
        output    = out_file
      ),
      warning = function(w) {
        captured_warnings <<- c(captured_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    status <- if (length(captured_warnings)) "OK_WITH_WARNING" else "OK"
    msg <- if (length(captured_warnings)) {
      paste0("ligation OK with note: ",
             paste(captured_warnings, collapse = " | "))
    } else {
      sprintf("ligation OK (insert L=%s R=%s)",
              attr(insert, "left_oh"), attr(insert, "right_oh"))
    }
    list(status = status,
         product_bp = nchar(asm$sequence),
         output_file = out_file,
         message = msg)
  }, error = function(e) {
    list(status = "FAIL",
         product_bp = NA_integer_,
         output_file = NA_character_,
         message = conditionMessage(e))
  })

  cat(result$status, "—", result$message, "\n")

  results <- rbind(results, data.frame(
    gene_name   = gene,
    species     = row$species,
    accession   = row$accession,
    status      = result$status,
    product_bp  = result$product_bp,
    output_file = result$output_file %||% NA_character_,
    message     = result$message,
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# Summary
# =============================================================================
summary_path <- file.path(output_dir, "ligation_summary.tsv")
write.table(results, summary_path, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")

n_ok      <- sum(results$status == "OK")
n_ok_warn <- sum(results$status == "OK_WITH_WARNING")
n_fail    <- sum(results$status == "FAIL")

cat("\n=== Summary ===\n")
cat(sprintf("  OK:                %d\n", n_ok))
cat(sprintf("  OK with note:      %d\n", n_ok_warn))
cat(sprintf("  FAIL:              %d\n", n_fail))
cat(sprintf("\nGenBank files: %s\n", output_dir))
cat(sprintf("Summary TSV:   %s\n", summary_path))

if (n_fail > 0) {
  cat("\nFailures:\n")
  failed <- results[results$status == "FAIL", c("gene_name", "message")]
  for (i in seq_len(nrow(failed))) {
    cat(sprintf("  %s: %s\n", failed$gene_name[i], failed$message[i]))
  }
  quit(status = 1)
}
