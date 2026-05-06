# Ligate fragments into a unique circular product.
#
# Each fragment has a left_oh and right_oh (4 nt, top-strand). A directed
# graph is built where there is an edge from fragment A to fragment B iff
# A$right_oh == B$left_oh. A valid Golden Gate assembly corresponds to a
# Hamiltonian cycle in this graph (each fragment used exactly once, returning
# to start). For correct designs the graph is a single cycle and the
# Hamiltonian cycle is unique.
#
# We also consider the reverse complement of each fragment as an alternative
# orientation, since cassette parts can be flipped. (For canonical YTK this
# is rarely needed because the overhang scheme is asymmetric, but it's
# correct to allow it.)

reverse_fragment <- function(f) {
  # Reverse complement: the new left overhang is revcomp of old right_oh,
  # the new right overhang is revcomp of old left_oh, and the sequence is
  # revcomp'd. Features get mirrored.
  L <- nchar(f$sequence)
  rc_seq <- revcomp(f$sequence)
  feats <- f$features
  if (nrow(feats)) {
    new_starts <- L - feats$end + 1L
    new_ends   <- L - feats$start + 1L
    feats$start <- new_starts
    feats$end <- new_ends
    feats$strand <- -feats$strand
  }
  new_fragment(
    left_oh = revcomp(f$right_oh),
    sequence = rc_seq,
    right_oh = revcomp(f$left_oh),
    features = feats,
    source_name = paste0(f$source_name, "_rc"),
    offset = f$offset
  )
}

# Find a unique Hamiltonian cycle through the fragment set. Each fragment
# can be used in either orientation, but only ONCE total.
#
# Approach: enumerate cycles starting from fragment #1 (anchor — by symmetry
# every cycle visits frag 1 in some orientation, so anchoring it eliminates
# rotational symmetry without missing solutions). For each anchor orientation,
# DFS extending by matching right_oh -> left_oh, requiring no fragment is
# reused. A cycle closes when len == n AND the last fragment's right_oh
# equals the anchor's left_oh.
#
# n is small in practice (typically 5-10 for YTK assemblies), so exhaustive
# search is fine.

find_circular_assembly <- function(fragments) {
  n <- length(fragments)
  if (n == 0) stop("No usable fragments after digestion")

  # Build oriented fragment list: each input fragment yields 2 oriented
  # versions (forward and revcomp). We track which "input index" each came
  # from so we don't reuse the same physical fragment.
  oriented <- list()
  source_idx <- integer()
  orient <- integer()
  for (i in seq_along(fragments)) {
    oriented[[length(oriented) + 1]] <- fragments[[i]]
    source_idx <- c(source_idx, i)
    orient <- c(orient, 1L)
    oriented[[length(oriented) + 1]] <- reverse_fragment(fragments[[i]])
    source_idx <- c(source_idx, i)
    orient <- c(orient, -1L)
  }

  # Anchor: try each orientation of fragment 1 as the cycle start.
  solutions <- list()
  for (anchor in which(source_idx == 1)) {
    used_sources <- rep(FALSE, n)
    used_sources[1] <- TRUE
    path <- c(anchor)
    extend(path, used_sources, oriented, source_idx, n, solutions_env = environment())
  }

  if (!length(solutions)) {
    stop("No circular assembly possible: fragments cannot be joined into a closed cycle. ",
         "Check that overhangs are compatible and no fragment is missing.")
  }

  # Dedupe: a cycle and its reverse-complement traversal are the same
  # physical product. Canonicalize by always picking the lexicographically
  # smallest rotation/reversal of the source-index sequence.
  canonical <- vapply(solutions, function(path) {
    src_seq <- source_idx[path]
    canonical_rotation(src_seq)
  }, character(1))
  unique_solutions <- solutions[!duplicated(canonical)]

  if (length(unique_solutions) > 1) {
    stop("Multiple distinct circular products possible (", length(unique_solutions),
         "): assembly is ambiguous. ",
         "This usually means two fragments share the same overhang. ",
         "Check input parts for duplicate overhangs.")
  }

  list(path = unique_solutions[[1]], oriented = oriented)
}

extend <- function(path, used_sources, oriented, source_idx, n, solutions_env) {
  current <- oriented[[path[length(path)]]]
  if (length(path) == n) {
    # Cycle closure check
    anchor <- oriented[[path[1]]]
    if (current$right_oh == anchor$left_oh) {
      solutions_env$solutions[[length(solutions_env$solutions) + 1]] <- path
    }
    return(invisible())
  }
  # Try extending with any unused source whose left_oh matches current$right_oh
  for (j in seq_along(oriented)) {
    if (used_sources[source_idx[j]]) next
    if (oriented[[j]]$left_oh != current$right_oh) next
    used_sources[source_idx[j]] <- TRUE
    extend(c(path, j), used_sources, oriented, source_idx, n, solutions_env)
    used_sources[source_idx[j]] <- FALSE
  }
}

canonical_rotation <- function(v) {
  # Return the lex-smallest rotation of v (and its reverse), as a string key.
  rotations <- function(x) {
    n <- length(x)
    vapply(seq_len(n), function(i) paste(c(x[i:n], if (i > 1) x[1:(i-1)]), collapse = ","),
           character(1))
  }
  forward <- rotations(v)
  backward <- rotations(rev(v))
  min(c(forward, backward))
}

# Stitch the oriented fragments along the cycle into a single circular
# top-strand sequence with merged features.
stitch_assembly <- function(path, oriented) {
  fragments <- oriented[path]
  # Each fragment contributes: left_oh + sequence (the 4 nt overhang plus
  # the duplex core). The right_oh is consumed by the next fragment's left_oh.
  pieces <- character(length(fragments))
  feature_chunks <- vector("list", length(fragments))
  cum_offset <- 0L
  for (i in seq_along(fragments)) {
    f <- fragments[[i]]
    piece <- paste0(f$left_oh, f$sequence)
    pieces[i] <- piece
    feats <- f$features
    if (nrow(feats)) {
      # Features were stored in coords where position 1 = start of left_oh.
      # Adjust by cum_offset.
      feats$start <- feats$start + cum_offset
      feats$end   <- feats$end   + cum_offset
      feature_chunks[[i]] <- feats
    }
    cum_offset <- cum_offset + nchar(piece)
  }
  full_seq <- paste(pieces, collapse = "")

  # Add a /label feature for each fragment showing its source
  source_features <- do.call(rbind, lapply(seq_along(fragments), function(i) {
    f <- fragments[[i]]
    start <- if (i == 1) 1L else sum(nchar(pieces[1:(i-1)])) + 1L
    end   <- sum(nchar(pieces[1:i]))
    df <- data.frame(type = "misc_feature", start = start, end = end,
                     strand = 1L, stringsAsFactors = FALSE)
    df$qualifiers <- list(list(label = f$source_name,
                               note = "fragment from goldengateR assembly"))
    df
  }))

  all_features <- if (length(feature_chunks)) {
    chunks <- Filter(function(x) !is.null(x) && nrow(x) > 0, feature_chunks)
    if (length(chunks)) do.call(rbind, chunks) else empty_features()
  } else {
    empty_features()
  }
  all_features <- rbind(all_features, source_features)
  rownames(all_features) <- NULL

  list(sequence = full_seq, features = all_features)
}

#' Assemble plasmids and PCR products by Golden Gate
#'
#' Performs an in silico Golden Gate assembly using BsaI or BsmBI. Each input
#' is digested at all enzyme sites (both strands), the resulting fragments
#' are joined by 4 nt sticky-end matching, and the unique closed circular
#' product is returned. Errors if the assembly is impossible or ambiguous.
#'
#' @param parts A list of ggr_part objects (use \code{read_part} to load
#'   GenBank or FASTA files). At least 2 parts are required for a meaningful
#'   assembly.
#' @param enzyme Either "BsaI" or "BsmBI".
#' @param name Name for the assembled product (becomes LOCUS in the GenBank
#'   output). Defaults to "assembly".
#' @param output Optional path. If supplied, the assembly is written to this
#'   file as GenBank.
#' @return A ggr_part representing the circular product.
#' @export
#' @examples
#' \dontrun{
#'   vec  <- read_part("pYTK090.gb")
#'   p1   <- read_part("pYTK008.gb")
#'   p234 <- read_part("pYTK047.gb")
#'   p5   <- read_part("pYTK073.gb")
#'   p6   <- read_part("pYTK074.gb")
#'   p7   <- read_part("pYTK086.gb")
#'   p8b  <- read_part("pYTK092.gb")
#'   asm <- golden_gate_assemble(list(vec, p1, p234, p5, p6, p7, p8b),
#'                               enzyme = "BsaI",
#'                               name = "pYTK096_rebuilt",
#'                               output = "pYTK096_rebuilt.gb")
#' }
golden_gate_assemble <- function(parts, enzyme = c("BsaI", "BsmBI"),
                                 name = "assembly", output = NULL) {
  enzyme <- match.arg(enzyme)
  if (!is.list(parts) || !length(parts)) stop("parts must be a non-empty list")
  for (i in seq_along(parts)) {
    if (!inherits(parts[[i]], "ggr_part")) {
      stop("parts[[", i, "]] is not a ggr_part. Use read_part() to load it.")
    }
  }

  # Digest all parts. Pool all valid fragments.
  all_fragments <- list()
  per_part <- integer(length(parts))
  for (i in seq_along(parts)) {
    fr <- digest_part(parts[[i]], enzyme)
    per_part[i] <- length(fr)
    all_fragments <- c(all_fragments, fr)
  }
  if (!length(all_fragments)) {
    stop("No fragments with valid 4 nt overhangs were produced. ",
         "Check that the parts contain ", enzyme, " sites.")
  }

  # Find the unique cycle.
  res <- find_circular_assembly(all_fragments)
  stitched <- stitch_assembly(res$path, res$oriented)

  asm <- new_part(sequence = stitched$sequence,
                  features = stitched$features,
                  topology = "circular",
                  name = name)

  if (!is.null(output)) {
    write_genbank(asm, output)
    message("Wrote assembly to ", output, " (", nchar(asm$sequence), " bp, ",
            nrow(asm$features), " features)")
  }
  asm
}
