\name{golden_gate_assemble}
\alias{golden_gate_assemble}
\title{Assemble plasmids by one-pot Golden Gate}
\usage{
golden_gate_assemble(parts, enzyme = c("BsaI", "BsmBI"),
                     name = "assembly", output = NULL)
}
\arguments{
\item{parts}{A list of \code{ggr_part} objects. Use
\code{\link{read_part}} or \code{\link{load_plasmid}} to load each part.
At least 2 parts are required for a meaningful assembly.}

\item{enzyme}{Either \code{"BsaI"} (recognition GGTCTC, used for Level 0
to Level 1 assemblies) or \code{"BsmBI"} (recognition CGTCTC, used for
Level 1 to Level 2 assemblies). Both produce 4 nt 5' overhangs.}

\item{name}{Name for the assembled product. Becomes the LOCUS field if
written to GenBank.}

\item{output}{Optional path. If supplied, the assembly is written to
this file as GenBank.}
}
\value{
A circular \code{ggr_part} representing the unique assembled product.
}
\description{
Performs an in silico one-pot Golden Gate reaction: each input part is
digested at all enzyme sites on both strands, the resulting fragments
are joined by 4 nt sticky-end matching, and the unique closed circular
product is returned.
}
\details{
The function is strict. It returns the unique cycle that uses every
fragment-with-valid-overhangs exactly once, or it errors. Three error
conditions are distinguished:

\itemize{
  \item \strong{No fragments produced} — none of the parts contained the
  necessary recognition sites.
  \item \strong{No closed cycle} — fragments cannot be ligated into a
  circle (overhang mismatch somewhere). Most often a primer- or
  cassette-design bug.
  \item \strong{Multiple distinct cycles} — two fragments share an
  overhang. Almost always a duplicate-overhang bug in the part design.
}

Fragments that retain a recognition site after digestion are filtered
out. This is what makes a Golden Gate dropout work: when sites face
inward toward a dropout, the recognition sites end up in the dropout
fragment and are discarded; only the clean backbone participates in
ligation. For workflows where you have already gel-purified a fragment
with internal sites, use \code{\link{ligate}} instead.

Reverse-complement orientations of each fragment are considered during
the cycle search, so parts in either orientation will assemble correctly
provided the overhangs match.
}
\examples{
\dontrun{
  vec  <- load_plasmid("pYTK090")
  p1   <- load_plasmid("pYTK008")
  p234 <- load_plasmid("pYTK047")
  asm  <- golden_gate_assemble(list(vec, p1, p234, ...),
                               enzyme = "BsaI",
                               name = "pYTK096_rebuilt",
                               output = "pYTK096_rebuilt.gb")
}
}
\seealso{\code{\link{ligate}}, \code{\link{digest_and_select}}}
