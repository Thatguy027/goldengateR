\name{digest_and_select}
\alias{digest_and_select}
\title{Digest and gel-purify a single fragment}
\usage{
digest_and_select(part, enzyme = c("BsaI", "BsmBI"), select = "largest")
}
\arguments{
\item{part}{A \code{ggr_part} (circular plasmid or linear PCR product).}

\item{enzyme}{Either \code{"BsaI"} or \code{"BsmBI"}.}

\item{select}{Which fragment to return. One of:
\itemize{
  \item \code{"largest"} — the largest fragment (typical for backbone
  purification: CamR + ori is much larger than a dropout).
  \item \code{"smallest"} — the smallest fragment (typical for picking
  a released part to use in a downstream Level 1 assembly).
  \item A length-2 numeric vector \code{c(min_bp, max_bp)} — the unique
  fragment whose total length falls in that size band. Errors if zero
  or more than one fragment matches.
}}
}
\value{
A linear \code{ggr_part} with two attributes attached:
\code{left_oh} and \code{right_oh}, each a 4 nt character string giving
the 5' overhangs at each end. The \code{$sequence} field contains the
top strand including both overhangs.
}
\description{
Simulates a digest-and-gel-purify step: digests the input with the
chosen enzyme and returns one fragment selected by size, ready to feed
into \code{\link{ligate}}.
}
\details{
Unlike \code{\link{golden_gate_assemble}}, this function does \strong{not}
filter fragments that contain internal recognition sites. The wet-lab
gel purification physically separates the band from the enzyme, so
internal sites in the purified fragment are biologically fine — they
just won't be re-cut.

This is the function to use when you need to gel-purify a backbone
fragment from a plasmid where the released "insert" fragment carries
the recognition sites (e.g., constructing a Type 3 entry plasmid by
replacing the GFP dropout with a new CDS).
}
\examples{
\dontrun{
  # Gel-purify the CamR + ori backbone from a Type 3 destination vector
  vec      <- load_plasmid("pYTK_T3_dest")
  backbone <- digest_and_select(vec, "BsaI", select = "largest")

  # Or by gel band size:
  backbone <- digest_and_select(vec, "BsaI", select = c(2000, 4000))
}
}
\seealso{\code{\link{ligate}}, \code{\link{make_fragment}}}
