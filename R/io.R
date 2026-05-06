#' Read a part from GenBank or FASTA
#'
#' @param file Path to a .gb/.gbk or .fa/.fasta file.
#' @param topology For FASTA only, "circular" or "linear". GenBank files
#'   declare topology in the LOCUS line and override this argument.
#' @return a ggr_part
#' @export
read_part <- function(file, topology = "circular") {
  if (!file.exists(file)) stop("File not found: ", file)
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("gb", "gbk", "genbank")) {
    read_genbank(file)
  } else if (ext %in% c("fa", "fasta", "fna")) {
    read_fasta(file, topology = topology)
  } else {
    # sniff
    first <- readLines(file, n = 1, warn = FALSE)
    if (length(first) && startsWith(first, "LOCUS")) read_genbank(file)
    else if (length(first) && startsWith(first, ">")) read_fasta(file, topology = topology)
    else stop("Cannot determine format of ", file)
  }
}

read_fasta <- function(file, topology = "circular") {
  lines <- readLines(file, warn = FALSE)
  if (!length(lines) || !startsWith(lines[1], ">")) {
    stop("Not a FASTA file: ", file)
  }
  header <- sub("^>", "", lines[1])
  name <- strsplit(header, "\\s+")[[1]][1]
  seq <- paste(lines[-1], collapse = "")
  seq <- gsub("\\s+", "", seq)
  new_part(sequence = seq, topology = topology, name = name)
}

# --- GenBank parsing ---
#
# GenBank format reference: NCBI flat file format.
# - LOCUS line: name, length, "linear"|"circular", date
# - FEATURES section: starts at column 1 with "FEATURES",
#   feature key at columns 6-20, location at column 22+,
#   qualifiers prefixed with "/" at column 22, can wrap.
# - ORIGIN section: sequence in 6-column blocks, terminated by "//"

read_genbank <- function(file) {
  lines <- readLines(file, warn = FALSE)
  locus_line <- lines[grepl("^LOCUS", lines)][1]
  if (is.na(locus_line)) stop("No LOCUS line in ", file)

  # LOCUS name is field 2
  fields <- strsplit(trimws(locus_line), "\\s+")[[1]]
  name <- fields[2]
  topology <- if (any(grepl("circular", lines[1], ignore.case = TRUE))) "circular" else "linear"

  # Sequence: between ORIGIN and //
  origin_idx <- which(grepl("^ORIGIN", lines))
  end_idx <- which(grepl("^//", lines))
  if (!length(origin_idx) || !length(end_idx)) stop("Malformed GenBank: ", file)
  seq_lines <- lines[(origin_idx[1] + 1):(end_idx[1] - 1)]
  seq <- paste(seq_lines, collapse = "")
  seq <- gsub("[^A-Za-z]", "", seq)

  # Features: between FEATURES line and ORIGIN
  feat_idx <- which(grepl("^FEATURES", lines))
  features <- if (length(feat_idx)) {
    parse_features(lines[(feat_idx[1] + 1):(origin_idx[1] - 1)])
  } else {
    empty_features()
  }

  new_part(sequence = seq, features = features, topology = topology, name = name)
}

parse_features <- function(flines) {
  # Group lines into feature blocks. A new feature starts when columns 1-5
  # are blank and column 6 is non-blank (the feature key column).
  is_new_feature <- function(L) {
    nchar(L) >= 6 &&
      substr(L, 1, 5) == "     " &&
      substr(L, 6, 6) != " "
  }
  starts <- which(vapply(flines, is_new_feature, logical(1)))
  if (!length(starts)) return(empty_features())
  ends <- c(starts[-1] - 1, length(flines))

  feats <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    block <- flines[starts[i]:ends[i]]
    feats[[i]] <- parse_feature_block(block)
  }
  feats <- Filter(Negate(is.null), feats)
  if (!length(feats)) return(empty_features())

  df <- data.frame(
    type   = vapply(feats, `[[`, character(1), "type"),
    start  = vapply(feats, `[[`, integer(1),   "start"),
    end    = vapply(feats, `[[`, integer(1),   "end"),
    strand = vapply(feats, `[[`, integer(1),   "strand"),
    stringsAsFactors = FALSE
  )
  df$qualifiers <- lapply(feats, `[[`, "qualifiers")
  df
}

parse_feature_block <- function(block) {
  first <- block[1]
  type <- trimws(substr(first, 6, 20))
  loc <- trimws(substr(first, 22, nchar(first)))

  # Continuation lines for the location appear before the first /qualifier.
  # Qualifier continuations are joined too.
  rest <- block[-1]
  qual_starts <- grep("^\\s{21}/", rest)
  if (length(qual_starts)) {
    loc_cont <- rest[seq_len(qual_starts[1] - 1)]
    qual_lines <- rest[qual_starts[1]:length(rest)]
  } else {
    loc_cont <- rest
    qual_lines <- character(0)
  }
  loc <- paste0(loc, paste(trimws(loc_cont), collapse = ""))

  parsed_loc <- parse_location(loc)
  if (is.null(parsed_loc)) return(NULL)

  qualifiers <- parse_qualifiers(qual_lines)

  list(type = type,
       start = parsed_loc$start,
       end = parsed_loc$end,
       strand = parsed_loc$strand,
       qualifiers = qualifiers)
}

parse_location <- function(loc) {
  # Handle complement(...) and simple a..b. Skip joins (we don't support
  # split features in assemblies; they get dropped with a warning).
  strand <- 1L
  if (grepl("^complement\\(", loc)) {
    strand <- -1L
    loc <- sub("^complement\\(", "", loc)
    loc <- sub("\\)$", "", loc)
  }
  if (grepl("^join\\(", loc) || grepl("^order\\(", loc)) {
    return(NULL)  # caller drops
  }
  loc <- gsub("[<>]", "", loc)
  m <- regmatches(loc, regexec("^(\\d+)\\.\\.(\\d+)$", loc))[[1]]
  if (length(m) == 3) {
    return(list(start = as.integer(m[2]), end = as.integer(m[3]), strand = strand))
  }
  # Single position
  m <- regmatches(loc, regexec("^(\\d+)$", loc))[[1]]
  if (length(m) == 2) {
    p <- as.integer(m[2])
    return(list(start = p, end = p, strand = strand))
  }
  NULL
}

parse_qualifiers <- function(qual_lines) {
  if (!length(qual_lines)) return(character(0))
  # Join continuations: a new qualifier starts with /key=
  is_new_qual <- grepl("^\\s{21}/", qual_lines)
  group <- cumsum(is_new_qual)
  blocks <- split(qual_lines, group)
  blocks <- blocks[names(blocks) != "0"]  # any pre-qual junk

  out <- list()
  for (b in blocks) {
    joined <- paste(trimws(b), collapse = " ")
    joined <- sub("^/", "", joined)
    eq <- regexpr("=", joined, fixed = TRUE)
    if (eq > 0) {
      key <- substr(joined, 1, eq - 1)
      val <- substr(joined, eq + 1, nchar(joined))
      val <- gsub('^"|"$', "", val)
      out[[key]] <- val
    } else {
      out[[joined]] <- ""
    }
  }
  out
}

#' Write a part to GenBank
#'
#' @param part A ggr_part.
#' @param file Output path.
#' @export
write_genbank <- function(part, file) {
  stopifnot(inherits(part, "ggr_part"))
  seq <- part$sequence
  L <- nchar(seq)
  out <- character()

  # LOCUS line: name, length bp, dna, topology, division, date
  topo <- if (part$topology == "circular") "circular" else "linear"
  out <- c(out, sprintf("LOCUS       %-16s %d bp    DNA     %s SYN %s",
                        substr(part$name, 1, 16), L, topo,
                        format(Sys.Date(), "%d-%b-%Y")))
  out <- c(out, sprintf("DEFINITION  %s assembled by goldengateR", part$name))
  out <- c(out, "ACCESSION   .")
  out <- c(out, "VERSION     .")
  out <- c(out, "KEYWORDS    .")
  out <- c(out, "SOURCE      synthetic DNA construct")
  out <- c(out, "  ORGANISM  synthetic DNA construct")
  out <- c(out, "FEATURES             Location/Qualifiers")

  if (nrow(part$features)) {
    for (i in seq_len(nrow(part$features))) {
      out <- c(out, format_feature(part$features[i, ]))
    }
  }

  out <- c(out, "ORIGIN")
  out <- c(out, format_origin(seq))
  out <- c(out, "//")

  writeLines(out, file)
}

format_feature <- function(f) {
  # f is a 1-row data.frame
  loc <- if (f$start <= f$end) {
    sprintf("%d..%d", f$start, f$end)
  } else {
    # Wrapping feature - emit as join
    sprintf("join(%d..%d,1..%d)", f$start, attr(f, "seqlen") %||% f$start, f$end)
  }
  if (f$strand == -1) loc <- sprintf("complement(%s)", loc)

  lines <- sprintf("     %-15s %s", f$type, loc)
  quals <- f$qualifiers[[1]]
  if (length(quals)) {
    for (k in names(quals)) {
      v <- quals[[k]]
      if (nchar(v)) {
        # Quote unless it's a numeric-only value
        if (grepl("^[0-9.]+$", v)) {
          lines <- c(lines, sprintf("                     /%s=%s", k, v))
        } else {
          lines <- c(lines, sprintf('                     /%s="%s"', k, v))
        }
      } else {
        lines <- c(lines, sprintf("                     /%s", k))
      }
    }
  }
  lines
}

`%||%` <- function(a, b) if (is.null(a)) b else a

format_origin <- function(seq) {
  seq <- tolower(seq)
  L <- nchar(seq)
  out <- character()
  for (i in seq(1, L, by = 60)) {
    chunk <- substr(seq, i, min(i + 59, L))
    # Split into 6 blocks of 10
    blocks <- character(6)
    for (j in 0:5) {
      blocks[j + 1] <- substr(chunk, j * 10 + 1, j * 10 + 10)
    }
    blocks <- blocks[blocks != ""]
    out <- c(out, sprintf("%9d %s", i, paste(blocks, collapse = " ")))
  }
  out
}
