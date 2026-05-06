\name{ligate}
\alias{ligate}
\title{Ligate pre-prepared fragments into a circular product}
\usage{
ligate(fragments, enzyme = NULL, name = "ligation", output = NULL)
}
\arguments{
\item{fragments}{A list of linear \code{ggr_part}s with overhangs
attached, from \code{\link{digest_and_select}} or
\code{\link{make_fragment}}. At least 2 fragments are required.}

\item{enzyme}{Optional. If supplied (\code{"BsaI"} or \code{"BsmBI"}),
the assembled product is checked for residual recognition sites and a
warning is issued if any are found. Pass \code{NULL} to skip the check.}

\item{name}{Name for the assembled product.}

\item{output}{Optional path for GenBank output.}
}
\value{
A circular \code{ggr_part}.
}
\description{
Performs sticky-end ligation of fragments that have already been
digested or otherwise prepared with 4 nt 5' overhangs. Unlike
\code{\link{golden_gate_assemble}}, this function does \strong{not}
digest its inputs and does \strong{not} filter for internal recognition
sites.
}
\details{
This is the right function for sequential workflows:
\enumerate{
  \item Digest a vector with BsaI
  \item Run on a gel, gel-purify the backbone band
  \item Digest a PCR product with BsaI separately
  \item Mix backbone and digested PCR product
  \item Add T4 ligase
}

The warning emitted when \code{enzyme} is supplied and the product
contains internal sites is \strong{informational} — it is normal for a
finished Type 3 entry plasmid to contain 2 BsaI sites (the cassette
flanks that release the CDS in the next Golden Gate round). It only
matters if you intend to re-use the product as input for further Golden
Gate with the same enzyme without an intermediate cut.
}
\examples{
\dontrun{
  vec      <- load_plasmid("pYTK_T3_dest")
  backbone <- digest_and_select(vec, "BsaI", select = "largest")

  amp    <- read_part("xylB_amplicon.fa", topology = "linear")
  insert <- digest_and_select(amp, "BsaI", select = "largest")

  asm <- ligate(list(backbone, insert),
                enzyme = "BsaI",
                name   = "pYTK_T3_xylB",
                output = "pYTK_T3_xylB.gb")
}
}
\seealso{\code{\link{golden_gate_assemble}}, \code{\link{digest_and_select}},
\code{\link{make_fragment}}}
