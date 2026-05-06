\name{make_fragment}
\alias{make_fragment}
\title{Construct a fragment from sequence and overhangs}
\usage{
make_fragment(left_oh, sequence = "", right_oh, name = "synthetic")
}
\arguments{
\item{left_oh}{Character(1), 4 nt 5' overhang on the left end.}

\item{sequence}{Character(1), the duplex core sequence between the
overhangs. Pass an empty string for a fragment that consists only of
the overhangs (e.g., a short annealed-oligo linker).}

\item{right_oh}{Character(1), 4 nt 5' overhang on the right end.}

\item{name}{Optional name for the fragment.}
}
\value{
A linear \code{ggr_part} with \code{left_oh} and \code{right_oh}
attributes attached, ready to pass to \code{\link{ligate}}.
}
\description{
Builds a fragment representation directly from a sequence and a pair of
overhangs, without any digestion. Useful for annealed oligo duplexes,
synthesised gene blocks ordered with pre-formed overhangs, or any other
fragment that didn't come from digesting a plasmid in this package.
}
\examples{
\dontrun{
  oligo <- make_fragment(left_oh  = "CCAT",
                         sequence = "ATGAAAGCT",
                         right_oh = "AATG",
                         name     = "short_cds")
  asm <- ligate(list(backbone, oligo))
}
}
\seealso{\code{\link{ligate}}, \code{\link{digest_and_select}}}
