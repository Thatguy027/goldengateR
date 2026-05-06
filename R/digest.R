# Type IIS enzyme definitions for BsaI and BsmBI.
#
# Both are Type IIS enzymes that cut downstream of an asymmetric recognition
# site, leaving a 4 nt 5' overhang. The recognition site itself is excised
# and ends up in a discarded fragment, which is what makes Golden Gate
# scarless and directional.
#
# BsaI:  5'-GGTCTC(N)1 / 5(N)-3'   recognition GGTCTC, top cut +1, bottom cut +5
# BsmBI: 5'-CGTCTC(N)1 / 5(N)-3'   recognition CGTCTC, same geometry
#
# Concretely, on a top strand 5'-...GGTCTC N|NNNN...-3'
#                            3'-...CCAGAG NNNNN|...-5'
# The top strand is cleaved between positions +1 and +2 after the recognition
# end; the bottom strand between +5 and +6. This leaves the 4 nt at +2..+5 as
# a 5' overhang on the downstream fragment.

ENZYMES <- list(
  BsaI  = list(name = "BsaI",  site = "GGTCTC", top_cut = 1L, bottom_cut = 5L),
  BsmBI = list(name = "BsmBI", site = "CGTCTC", top_cut = 1L, bottom_cut = 5L)
)

# A fragment after digestion. Stored as the double-stranded core with
# 5' overhangs on each end. left_oh = 4 nt 5' overhang on the left (upstream)
# end of the top strand; right_oh = 4 nt 5' overhang on the right (downstream)
# end of the bottom strand, expressed as the corresponding top-strand bases
# that the partner overhang must base-pair with.
#
# Equivalently for ligation matching: the right_oh of one fragment must equal
# the left_oh of its successor (both written as the 4 nt 5' single-stranded
# top-strand bases at the junction).

new_fragment <- function(left_oh, sequence, right_oh, features = empty_features(),
                         source_name = "", offset = 0L) {
  list(left_oh = left_oh,
       sequence = sequence,
       right_oh = right_oh,
       features = features,
       source_name = source_name,
       offset = offset)  # original 1-based position of fragment[1] in source
}

# Find all enzyme cut sites on a sequence (both strands).
# Returns a data.frame: strand (+1/-1), top_cut_pos (1-based, position AFTER which the top strand is cut)
find_sites <- function(sequence, enzyme, circular) {
  e <- ENZYMES[[enzyme]]
  if (is.null(e)) stop("Unknown enzyme: ", enzyme)

  L <- nchar(sequence)
  # For circular, scan a doubled sequence so sites spanning the origin are found,
  # then dedupe by modular position.
  scan_seq <- if (circular) paste0(sequence, sequence) else sequence

  fwd <- gregexpr(e$site, scan_seq, fixed = TRUE)[[1]]
  rev_site <- revcomp(e$site)
  rev <- gregexpr(rev_site, scan_seq, fixed = TRUE)[[1]]

  hits <- list()
  if (fwd[1] != -1) {
    for (p in fwd) {
      # p is start of recognition site (1-based). Top strand cut is between
      # p+6+top_cut-1 and p+6+top_cut (i.e. after position p+6).
      # Recognition site ends at p+5; cut is e$top_cut nt downstream of site end.
      top_cut <- p + nchar(e$site) - 1 + e$top_cut  # position of last nt before cut
      bottom_cut <- p + nchar(e$site) - 1 + e$bottom_cut
      if (circular) top_cut <- ((top_cut - 1) %% L) + 1
      if (circular) bottom_cut <- ((bottom_cut - 1) %% L) + 1
      if (!circular && (top_cut > L || bottom_cut > L)) next
      hits[[length(hits) + 1]] <- list(
        strand = 1L, recog_start = if (circular) ((p - 1) %% L) + 1 else p,
        top_cut = top_cut, bottom_cut = bottom_cut
      )
    }
  }
  if (rev[1] != -1) {
    for (p in rev) {
      # Reverse-strand site: recognition reads CCN on top from p..p+5 (since
      # rev_site is revcomp of GGTCTC = GAGACC). The enzyme is bound to the
      # bottom strand and cuts upstream of its recognition (in top-strand
      # coordinates). Top cut is bottom_cut (5) nt UPSTREAM of recog start;
      # bottom cut is top_cut (1) nt upstream.
      # i.e. top strand cut at p - e$bottom_cut (between p-bottom_cut and p-bottom_cut+1)
      #      bottom strand cut at p - e$top_cut
      top_cut <- p - e$bottom_cut - 1  # last top-strand base before cut, 0-based offset corrected
      bottom_cut <- p - e$top_cut - 1
      # convert: cut is "after position X", so X = p - bottom_cut - 1 means
      # last retained top-strand base before the cut. Let me recompute:
      # site on top: positions p..p+5. Bottom strand recognition is the
      # revcomp at the same positions. Enzyme cuts top strand at offset
      # -bottom_cut from the 5' end of its bottom-strand recognition.
      # The bottom-strand 5' end (reading 3'->5' on top) is at top position p+5.
      # Going "downstream" on the bottom strand = "upstream" on top by
      # bottom_cut positions => top-strand cut after position (p - 1 - 0) = p-1
      # ... this is getting confused. Use a cleaner formulation:
      #
      # By symmetry: a reverse-strand site's cuts are obtained by reflecting
      # the forward-strand geometry around the recognition site.
      # Forward site at p: top cut after (p + 5 + 1) = p+6, bottom cut after p+10.
      # Reverse site at p: bottom cut after (p - 1) - 1 + 1 = p-1 ON THE TOP STRAND
      # we want the TOP strand cut for fragment splitting.
      # For a reverse-strand BsaI site GAGACC on top at positions p..p+5,
      # the enzyme cuts top strand 5 nt to the LEFT of p, between p-5 and p-4.
      # And cuts bottom strand 1 nt to the left, between p-1 and p (in top coords).
      # Wait — bottom strand cut in TOP coords means the position on the top
      # strand where the bottom strand break is.
      top_cut <- p - e$bottom_cut - 1   # top strand cut after this position
      bottom_cut <- p - e$top_cut - 1   # bottom strand cut after this position (top coords)
      if (circular) {
        top_cut <- ((top_cut - 1) %% L) + 1
        bottom_cut <- ((bottom_cut - 1) %% L) + 1
      }
      if (!circular && (top_cut < 1 || bottom_cut < 1)) next
      hits[[length(hits) + 1]] <- list(
        strand = -1L, recog_start = if (circular) ((p - 1) %% L) + 1 else p,
        top_cut = top_cut, bottom_cut = bottom_cut
      )
    }
  }

  if (!length(hits)) {
    return(data.frame(strand = integer(), recog_start = integer(),
                      top_cut = integer(), bottom_cut = integer(),
                      oh_start = integer(), oh_end = integer(),
                      overhang = character()))
  }
  df <- do.call(rbind, lapply(hits, as.data.frame))
  df <- unique(df)

  # Compute overhangs. For a Type IIS cut with top_cut < bottom_cut, the
  # 4 nt overhang sits at positions (top_cut + 1)..bottom_cut on the top
  # strand. But on a circular molecule the cut may straddle the origin:
  # then top_cut is large (near L) and bottom_cut is small (near 1), with
  # the overhang wrapping. Detect this by: for circular topology, the
  # SHORTER of the two intervals (top_cut+1..bottom_cut or
  # bottom_cut+1..top_cut going around) is the 4 nt overhang. Equivalently,
  # the overhang always has length = bottom_cut_offset - top_cut_offset = 4
  # in modular arithmetic (mod L) where the cut offsets differ by 4.
  L <- nchar(sequence)
  df$oh_start <- integer(nrow(df))
  df$oh_end   <- integer(nrow(df))
  df$overhang <- character(nrow(df))
  for (i in seq_len(nrow(df))) {
    tc <- df$top_cut[i]; bc <- df$bottom_cut[i]
    # Try the non-wrapping interpretation first
    if (tc < bc && (bc - tc) == 4L) {
      oh_start <- tc + 1L
      oh_end   <- bc
      oh       <- substring(sequence, oh_start, oh_end)
    } else if (circular) {
      # Wrap case: overhang is (tc+1..L) ++ (1..bc), assuming the wrap
      # produces a 4 nt span.
      span_wrap <- ((bc - tc - 1L) %% L) + 1L  # length of (tc+1)..bc going forward mod L
      if (span_wrap == 4L) {
        oh_start <- (tc %% L) + 1L
        oh_end   <- bc
        oh <- if (oh_start <= oh_end) substring(sequence, oh_start, oh_end)
              else paste0(substring(sequence, oh_start, L),
                          substring(sequence, 1, oh_end))
      } else {
        # Try the other direction
        span_other <- ((tc - bc - 1L) %% L) + 1L
        if (span_other == 4L) {
          oh_start <- (bc %% L) + 1L
          oh_end   <- tc
          oh <- if (oh_start <= oh_end) substring(sequence, oh_start, oh_end)
                else paste0(substring(sequence, oh_start, L),
                            substring(sequence, 1, oh_end))
        } else {
          # Shouldn't happen with Type IIS enzymes that produce 4 nt overhangs
          oh_start <- NA_integer_; oh_end <- NA_integer_; oh <- NA_character_
        }
      }
    } else {
      # Linear, non-canonical: skip
      oh_start <- NA_integer_; oh_end <- NA_integer_; oh <- NA_character_
    }
    df$oh_start[i] <- oh_start
    df$oh_end[i] <- oh_end
    df$overhang[i] <- oh
  }
  # Drop any sites that didn't produce a valid 4 nt overhang
  df <- df[!is.na(df$overhang) & nchar(df$overhang) == 4L, , drop = FALSE]

  df <- df[order(df$top_cut), , drop = FALSE]
  rownames(df) <- NULL
  df
}

# Digest a single part with the enzyme and return the list of fragments.
# Each fragment has its 4 nt 5' overhangs on left and right ends.
#
# A 5' overhang from a Type IIS cut: top strand is cut at position T (after T),
# bottom strand at position B (after B), with B > T => bases T+1..B are the
# 5' overhang on the downstream fragment's left end (top strand) AND on the
# upstream fragment's right end (bottom strand, which equals revcomp of top
# strand T+1..B; expressed in top-strand bases for matching, this is the same
# 4 nt sequence T+1..B for the left overhang of the downstream fragment, and
# the same 4 nt sequence T+1..B for the right overhang of the upstream
# fragment when both are written as "the 4 nt that base-pair across the
# junction, on the top strand").

digest_part <- function(part, enzyme) {
  L <- nchar(part$sequence)
  circular <- part$topology == "circular"
  sites <- find_sites(part$sequence, enzyme, circular = circular)

  if (nrow(sites) == 0) {
    # No cuts. Return the whole part as a single fragment with no overhangs.
    # (Such fragments cannot participate in assembly and will be filtered out.)
    return(list(new_fragment(
      left_oh = "", sequence = part$sequence, right_oh = "",
      features = part$features, source_name = part$name, offset = 1L
    )))
  }

  # Each cut produces a top-strand break after position top_cut and a
  # bottom-strand break after position bottom_cut. The 4 nt between them
  # (positions oh_start..oh_end) is the overhang; find_sites has already
  # computed these columns.

  # For each cut, the "split point" on the top strand for fragment
  # boundaries is the top_cut for forward-strand sites and the top_cut for
  # reverse-strand sites too (both measured as "last top-strand base before
  # the top-strand break"). The 4 nt overhang then attaches to the LEFT end
  # of the downstream fragment (its top strand starts mid-overhang... no:
  # the downstream fragment top strand starts at oh_start; its left 5'
  # overhang IS oh_start..oh_end). The upstream fragment ends at top_cut on
  # the top strand and at bottom_cut on the bottom strand: its right end has
  # a 4 nt 5' overhang on the bottom strand, which in top-strand coords is
  # the same bases oh_start..oh_end (because the bottom strand extends
  # further by 4 nt = the overhang, base-pairing with what would be the
  # downstream top strand if not cut).
  #
  # BUT — for fragment matching during ligation, we use the TOP-strand
  # representation of the overhang bases at each junction. So both the
  # right_oh of fragment N and the left_oh of fragment N+1 are the SAME
  # 4-nt string (the top-strand sequence at the junction).

  # Order cuts by top_cut position for fragment generation.
  sites <- sites[order(sites$top_cut), , drop = FALSE]

  # Helper: the fragment "core" sequence is the double-stranded duplex
  # part. The top strand spans from the left end's oh_start (start of left
  # overhang) to the right end's oh_end (end of right overhang). I.e., the
  # top strand of the fragment INCLUDES both overhang regions; the bottom
  # strand is shorter on each end by 4 nt.
  # For ligation we represent the fragment as the full top-strand sequence
  # from oh_start_left to oh_end_right, with the left/right overhang bases
  # recorded separately. A junction matches when right_oh(A) == left_oh(B);
  # the joined fragment's top strand is just A_top + B_top_minus_left_overhang
  # ... actually because both fragments INCLUDE the overhang in their top
  # strand, joining them needs to drop one copy. Cleaner alternative:
  # represent the fragment top strand as JUST the duplex core (excluding
  # overhang), with overhangs as separate 4 nt strings.
  #
  # Switching to that representation now:
  # core_top = top strand from (left oh_end + 1) to (right oh_start - 1)
  # left_oh  = 4 nt at top-strand positions oh_start_left..oh_end_left
  # right_oh = 4 nt at top-strand positions oh_start_right..oh_end_right

  fragments <- list()

  if (circular) {
    # n cuts produce n fragments. Walk pairs of consecutive cuts (wrapping).
    n <- nrow(sites)
    for (i in seq_len(n)) {
      cur <- sites[i, ]
      nxt <- sites[if (i == n) 1 else i + 1, ]
      # left overhang of this fragment = current cut's overhang
      # right overhang = next cut's overhang
      # core spans from (cur$oh_end + 1) to (nxt$oh_start - 1), wrapping
      core_start <- cur$oh_end + 1L
      core_end <- nxt$oh_start - 1L
      core <- circular_substr(part$sequence, core_start, core_end, L)
      frag_features <- subset_features(part$features,
                                       cur$oh_start, nxt$oh_end, L, circular = TRUE)
      fragments[[length(fragments) + 1]] <- new_fragment(
        left_oh = cur$overhang,
        sequence = core,
        right_oh = nxt$overhang,
        features = frag_features,
        source_name = part$name,
        offset = cur$oh_start
      )
    }
  } else {
    # Linear: n cuts produce n+1 fragments. The first and last have no
    # overhang on their outer end (and so cannot ligate there) — these
    # are "edge" fragments and represent the unwanted ends of a PCR product.
    n <- nrow(sites)
    # Fragment 0: from 1 to sites[1]$top_cut, no left overhang, right overhang = sites[1]$overhang
    # ... but wait, the bases (sites[1]$oh_start)..(sites[1]$oh_end) are
    # the overhang on the DOWNSTREAM fragment, so the upstream fragment's
    # core ends at sites[1]$top_cut (= oh_start - 1).
    # We discard edge fragments (no overhang on outer end).
    for (i in seq_len(n - 1)) {
      cur <- sites[i, ]
      nxt <- sites[i + 1, ]
      core <- substring(part$sequence, cur$oh_end + 1L, nxt$oh_start - 1L)
      frag_features <- subset_features(part$features,
                                       cur$oh_start, nxt$oh_end, L, circular = FALSE)
      fragments[[length(fragments) + 1]] <- new_fragment(
        left_oh = cur$overhang,
        sequence = core,
        right_oh = nxt$overhang,
        features = frag_features,
        source_name = part$name,
        offset = cur$oh_start
      )
    }
  }

  # Filter out fragments that contain residual enzyme sites in their core
  # — these would be re-cut. Also filter fragments with empty overhangs.
  fragments <- Filter(function(f) {
    if (nchar(f$left_oh) != 4 || nchar(f$right_oh) != 4) return(FALSE)
    full <- paste0(f$left_oh, f$sequence, f$right_oh)
    e <- ENZYMES[[enzyme]]
    !grepl(e$site, full, fixed = TRUE) && !grepl(revcomp(e$site), full, fixed = TRUE)
  }, fragments)

  fragments
}

circular_substr <- function(seq, start, end, L) {
  # Inclusive 1-based, wrapping. start/end may exceed L or be < 1.
  start <- ((start - 1) %% L) + 1
  end <- ((end - 1) %% L) + 1
  if (start <= end) {
    substring(seq, start, end)
  } else {
    paste0(substring(seq, start, L), substring(seq, 1, end))
  }
}

# Subset features that fall within [from, to] (in original sequence coords).
# For circular sources where the fragment may wrap, this is approximate:
# we keep features fully inside the linear projection. Coordinates are
# rebased so that the fragment's left overhang starts at position 1.
subset_features <- function(features, from, to, L, circular) {
  if (!nrow(features)) return(features)
  keep_idx <- integer()
  new_starts <- integer()
  new_ends <- integer()

  if (circular && to < from) {
    # Wrapping fragment: kept range is [from..L] U [1..to]
    in_range <- (features$start >= from & features$end <= L) |
                (features$start >= 1 & features$end <= to)
    keep_idx <- which(in_range)
    for (i in keep_idx) {
      if (features$start[i] >= from) {
        new_starts <- c(new_starts, features$start[i] - from + 1L)
        new_ends   <- c(new_ends,   features$end[i]   - from + 1L)
      } else {
        new_starts <- c(new_starts, features$start[i] + (L - from + 1L))
        new_ends   <- c(new_ends,   features$end[i]   + (L - from + 1L))
      }
    }
  } else {
    in_range <- features$start >= from & features$end <= to
    keep_idx <- which(in_range)
    new_starts <- features$start[keep_idx] - from + 1L
    new_ends   <- features$end[keep_idx]   - from + 1L
  }

  if (!length(keep_idx)) return(empty_features())
  out <- features[keep_idx, , drop = FALSE]
  out$start <- new_starts
  out$end <- new_ends
  rownames(out) <- NULL
  out
}
