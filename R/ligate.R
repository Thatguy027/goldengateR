# Gel purification + direct ligation
#
# These functions extend the package to support the workflow:
#   1. Digest a plasmid with BsaI/BsmBI
#   2. Run on a gel, pick a band by size
#   3. Mix that gel-purified fragment with a PCR product
#      (also digested) and ligate
#
# The crucial difference from one-pot Golden Gate (golden_gate_assemble)
# is that gel-purified fragments are PROTECTED from re-digestion: even if
# they contain internal recognition sites, those sites won't be cut because
# the enzyme isn't present in the ligation step. So we expose the
# digestion result without the "filter fragments containing sites" step
# that one-pot assembly applies.

#' Digest a part and return a single fragment selected by gel-band size
#'
#' Simulates a digestion-and-gel-purification step. The part is digested
#' with the chosen enzyme and the fragment matching the size selection
#' is returned as a linear \code{ggr_part} carrying its 5' overhangs as
#' attributes for later ligation.
#'
#' Unlike \code{golden_gate_assemble}, this function does NOT filter
#' fragments that contain internal recognition sites — the gel
#' purification physically separates the band from the enzyme, so internal
#' sites in the purified fragment don't matter.
#'
#' @param part A ggr_part (from \code{read_part}).
#' @param enzyme Either "BsaI" or "BsmBI".
#' @param select One of:
#'   \itemize{
#'     \item \code{"largest"} — return the largest fragment (typical for
#'       backbone purification: CamR + ori is much larger than the dropout)
#'     \item \code{"smallest"} — return the smallest (typical for
#'       Level 1 assembly: the released Type 3 part is small)
#'     \item a length-2 numeric vector \code{c(min_bp, max_bp)} — return
#'       the unique fragment whose total length (overhangs + core) falls
#'       in that range. Errors if 0 or >1 fragments match.
#'   }
#' @return A \code{ggr_part} with topology = "linear" and two new attributes
#'   \code{left_oh} and \code{right_oh} carrying the 4 nt overhangs.
#'   The sequence stored in \code{$sequence} is the top-strand including
#'   both overhangs (5' overhang bases are part of the top strand).
#' @export
digest_and_select <- function(part, enzyme = c("BsaI", "BsmBI"),
                              select = "largest") {
  enzyme <- match.arg(enzyme)
  if (!inherits(part, "ggr_part")) stop("part must be a ggr_part")

  # Use the existing digestion machinery, but bypass the
  # "filter out fragments with internal sites" step.
  fragments <- digest_part_no_filter(part, enzyme)
  fragments <- Filter(function(f) nchar(f$left_oh) == 4 && nchar(f$right_oh) == 4,
                      fragments)
  if (!length(fragments)) {
    stop("Digestion produced no fragments with 4 nt overhangs on both ends. ",
         "Check that the part contains at least 2 ", enzyme, " sites.")
  }

  sizes <- vapply(fragments, function(f) {
    nchar(f$left_oh) + nchar(f$sequence) + nchar(f$right_oh)
  }, integer(1))

  picked <- if (identical(select, "largest")) {
    which.max(sizes)
  } else if (identical(select, "smallest")) {
    which.min(sizes)
  } else if (is.numeric(select) && length(select) == 2) {
    in_band <- which(sizes >= select[1] & sizes <= select[2])
    if (!length(in_band)) {
      stop("No fragment falls in size range ", select[1], "-", select[2],
           " bp. Fragment sizes were: ", paste(sort(sizes), collapse = ", "))
    }
    if (length(in_band) > 1) {
      stop("Multiple fragments (", length(in_band), ") fall in size range ",
           select[1], "-", select[2], " bp (sizes ",
           paste(sizes[in_band], collapse = ", "),
           "). Tighten the range to select a single band.")
    }
    in_band
  } else {
    stop('select must be "largest", "smallest", or c(min_bp, max_bp)')
  }

  f <- fragments[[picked]]
  # Build a linear ggr_part. The top strand includes both overhang regions.
  full_seq <- paste0(f$left_oh, f$sequence, f$right_oh)

  # Shift features: in fragment storage, features were already rebased so
  # that position 1 = first base of left_oh, so they're already correct.
  out <- new_part(sequence = full_seq,
                  features = f$features,
                  topology = "linear",
                  name = paste0(part$name, "_", enzyme, "_",
                                if (is.character(select)) select else "band"))
  attr(out, "left_oh") <- f$left_oh
  attr(out, "right_oh") <- f$right_oh
  attr(out, "purified_size_bp") <- nchar(full_seq)
  message(sprintf("Selected fragment: %d bp, left_oh=%s, right_oh=%s (of %d candidate fragments: %s bp)",
                  nchar(full_seq), f$left_oh, f$right_oh,
                  length(fragments), paste(sort(sizes), collapse = ", ")))
  out
}

# Internal: digest without filtering fragments that contain internal sites.
# Mirrors digest_part() but drops the final Filter() call.
digest_part_no_filter <- function(part, enzyme) {
  L <- nchar(part$sequence)
  circular <- part$topology == "circular"
  sites <- find_sites(part$sequence, enzyme, circular = circular)
  if (nrow(sites) == 0) {
    return(list(new_fragment(
      left_oh = "", sequence = part$sequence, right_oh = "",
      features = part$features, source_name = part$name, offset = 1L
    )))
  }
  sites <- sites[order(sites$top_cut), , drop = FALSE]
  fragments <- list()

  if (circular) {
    n <- nrow(sites)
    for (i in seq_len(n)) {
      cur <- sites[i, ]
      nxt <- sites[if (i == n) 1 else i + 1, ]
      core_start <- cur$oh_end + 1L
      core_end <- nxt$oh_start - 1L
      core <- circular_substr(part$sequence, core_start, core_end, L)
      frag_features <- subset_features(part$features,
                                       cur$oh_start, nxt$oh_end, L,
                                       circular = TRUE)
      fragments[[length(fragments) + 1]] <- new_fragment(
        left_oh = cur$overhang, sequence = core, right_oh = nxt$overhang,
        features = frag_features, source_name = part$name,
        offset = cur$oh_start
      )
    }
  } else {
    n <- nrow(sites)
    for (i in seq_len(n - 1)) {
      cur <- sites[i, ]; nxt <- sites[i + 1, ]
      core <- substring(part$sequence, cur$oh_end + 1L, nxt$oh_start - 1L)
      frag_features <- subset_features(part$features,
                                       cur$oh_start, nxt$oh_end, L,
                                       circular = FALSE)
      fragments[[length(fragments) + 1]] <- new_fragment(
        left_oh = cur$overhang, sequence = core, right_oh = nxt$overhang,
        features = frag_features, source_name = part$name,
        offset = cur$oh_start
      )
    }
  }
  # Only require valid 4 nt overhangs on both ends — do NOT filter for
  # internal recognition sites.
  Filter(function(f) nchar(f$left_oh) == 4 && nchar(f$right_oh) == 4,
         fragments)
}

#' Construct a fragment directly from sequence and overhangs
#'
#' Useful for representing pre-formed fragments such as annealed oligo
#' duplexes, synthesized double-stranded gene blocks with designed
#' overhangs, or any other fragment that did not come from digesting a
#' plasmid in this package.
#'
#' @param left_oh Character(1), 4 nt 5' overhang on the left end.
#' @param sequence Character(1), the duplex core sequence (between overhangs).
#'   Pass an empty string for an oligo that consists only of the overhangs.
#' @param right_oh Character(1), 4 nt 5' overhang on the right end.
#' @param name Optional name (default "synthetic").
#' @return A linear \code{ggr_part} with overhangs as attributes, ready
#'   to pass to \code{ligate()}.
#' @export
make_fragment <- function(left_oh, sequence = "", right_oh, name = "synthetic") {
  if (nchar(left_oh) != 4 || nchar(right_oh) != 4) {
    stop("Both overhangs must be exactly 4 nt")
  }
  if (!grepl("^[ACGTNacgtn]*$", sequence)) {
    stop("sequence contains non-ACGTN characters")
  }
  full_seq <- paste0(toupper(left_oh), toupper(sequence), toupper(right_oh))
  out <- new_part(sequence = full_seq,
                  features = empty_features(),
                  topology = "linear",
                  name = name)
  attr(out, "left_oh") <- toupper(left_oh)
  attr(out, "right_oh") <- toupper(right_oh)
  out
}

#' Ligate pre-digested fragments into a unique circular product
#'
#' Performs sticky-end ligation of fragments that have already been
#' digested or otherwise prepared with 4 nt 5' overhangs. Unlike
#' \code{golden_gate_assemble}, this function does NOT digest its inputs
#' and does NOT filter for internal recognition sites — it assumes the
#' user has prepared each fragment intentionally (gel-purified band,
#' annealed oligo, designed gene block, etc.).
#'
#' This is the appropriate function for the workflow:
#'   1. Digest entry vector + PCR product with BsaI
#'   2. Gel-purify the entry vector backbone (large fragment)
#'   3. Mix with the digested PCR product (small fragment)
#'   4. Add T4 ligase
#'
#' represented as:
#' \preformatted{
#'   bb  <- digest_and_select(entry_vec, "BsaI", select = "largest")
#'   pcr <- digest_and_select(read_part("amplicon.fa", topology = "linear"),
#'                            "BsaI", select = "largest")
#'   asm <- ligate(list(bb, pcr), enzyme = "BsaI",
#'                 name = "pYTK_xylB", output = "out.gb")
#' }
#'
#' Multi-fragment ligations work the same way (e.g. backbone + 3 oligo
#' duplexes assembled in a single tube).
#'
#' @param fragments List of linear \code{ggr_part}s with overhangs (from
#'   \code{digest_and_select} or \code{make_fragment}).
#' @param enzyme Optional. If supplied, the assembled product is checked
#'   for residual recognition sites and a warning is issued if any are
#'   found (these would prevent re-cutting in subsequent rounds of
#'   Golden Gate but are biologically fine). Pass \code{NULL} to skip the check.
#' @param name Name for the assembled product.
#' @param output Optional path to write GenBank.
#' @return A circular \code{ggr_part}.
#' @export
ligate <- function(fragments, enzyme = NULL, name = "ligation", output = NULL) {
  if (!is.list(fragments) || length(fragments) < 2) {
    stop("ligate() requires a list of at least 2 fragments")
  }
  # Convert each ggr_part with overhang attributes into the internal
  # fragment representation expected by find_circular_assembly.
  internal <- lapply(seq_along(fragments), function(i) {
    p <- fragments[[i]]
    if (!inherits(p, "ggr_part")) {
      stop("fragments[[", i, "]] is not a ggr_part. ",
           "Use digest_and_select() or make_fragment() to create one.")
    }
    if (p$topology != "linear") {
      stop("fragments[[", i, "]] is circular. Only linear fragments can be ligated. ",
           "Did you mean to use golden_gate_assemble() for circular plasmid inputs?")
    }
    left_oh <- attr(p, "left_oh")
    right_oh <- attr(p, "right_oh")
    if (is.null(left_oh) || is.null(right_oh)) {
      stop("fragments[[", i, "]] is missing left_oh/right_oh attributes. ",
           "Use digest_and_select() or make_fragment() to attach them.")
    }
    if (nchar(left_oh) != 4 || nchar(right_oh) != 4) {
      stop("fragments[[", i, "]] has overhangs of length ",
           nchar(left_oh), "/", nchar(right_oh), "; both must be 4 nt.")
    }
    # Internal fragment storage: $sequence is the CORE (without overhangs),
    # and overhangs are stored separately. The ggr_part stores the full
    # top strand including overhangs, so we strip them.
    full <- p$sequence
    core <- substring(full, 5L, nchar(full) - 4L)
    # Adjust feature coordinates: the public ggr_part has positions
    # relative to the start of the left overhang; internal storage uses
    # the same convention, so no adjustment is needed.
    new_fragment(left_oh = left_oh, sequence = core, right_oh = right_oh,
                 features = p$features, source_name = p$name, offset = 1L)
  })

  res <- find_circular_assembly(internal)
  stitched <- stitch_assembly(res$path, res$oriented)

  asm <- new_part(sequence = stitched$sequence,
                  features = stitched$features,
                  topology = "circular",
                  name = name)

  if (!is.null(enzyme)) {
    e <- ENZYMES[[enzyme]]
    if (!is.null(e)) {
      # Check both strands of the (circular) product for residual sites
      doubled <- paste0(asm$sequence, asm$sequence)
      n_fwd <- length(gregexpr(e$site, doubled, fixed = TRUE)[[1]])
      n_fwd <- if (n_fwd == 1 && gregexpr(e$site, doubled, fixed = TRUE)[[1]][1] == -1) 0 else n_fwd
      n_rev <- length(gregexpr(revcomp(e$site), doubled, fixed = TRUE)[[1]])
      n_rev <- if (n_rev == 1 && gregexpr(revcomp(e$site), doubled, fixed = TRUE)[[1]][1] == -1) 0 else n_rev
      # Subtract the wrap duplicates (sites are counted on both copies of the doubled seq)
      n_internal <- (n_fwd + n_rev) %/% 2
      if (n_internal > 0) {
        warning(sprintf("Assembled product contains %d internal %s site(s). ",
                        n_internal, enzyme),
                "This is fine for a finished plasmid but means the product ",
                "cannot serve as a Level 0 entry for further Golden Gate ",
                "rounds with the same enzyme without re-cutting.")
      }
    }
  }

  if (!is.null(output)) {
    write_genbank(asm, output)
    message("Wrote ligation product to ", output, " (",
            nchar(asm$sequence), " bp, ", nrow(asm$features), " features)")
  }
  asm
}
