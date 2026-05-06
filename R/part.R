#' Part: a DNA sequence with features and topology
#'
#' A simple S3 representation of a plasmid or linear PCR product.
#' Fields:
#'   sequence: character(1), uppercase ACGTN, the top strand 5'->3'
#'   features: data.frame with columns type, start, end, strand, qualifiers
#'             (qualifiers is a list-column of named character vectors)
#'   topology: "circular" or "linear"
#'   name:     character(1), used as LOCUS in output GenBank
#'
#' For circular parts, features may wrap (start > end). For linear parts they may not.
#'
#' @keywords internal
new_part <- function(sequence, features = empty_features(),
                     topology = c("circular", "linear"), name = "part") {
  topology <- match.arg(topology)
  sequence <- toupper(sequence)
  if (!grepl("^[ACGTN]*$", sequence)) {
    stop("Sequence contains non-ACGTN characters")
  }
  structure(
    list(sequence = sequence,
         features = features,
         topology = topology,
         name = name),
    class = "ggr_part"
  )
}

empty_features <- function() {
  data.frame(
    type = character(),
    start = integer(),
    end = integer(),
    strand = integer(),
    stringsAsFactors = FALSE
  ) -> df
  df$qualifiers <- list()
  df
}

#' @export
print.ggr_part <- function(x, ...) {
  cat(sprintf("<ggr_part> %s [%s, %d bp, %d features]\n",
              x$name, x$topology, nchar(x$sequence), nrow(x$features)))
  invisible(x)
}

#' Reverse complement a DNA string
#' @keywords internal
revcomp <- function(s) {
  s <- toupper(s)
  chars <- strsplit(s, "")[[1]]
  comp <- c(A = "T", T = "A", G = "C", C = "G", N = "N")
  paste(rev(comp[chars]), collapse = "")
}
