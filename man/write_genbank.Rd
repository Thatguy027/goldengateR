\name{write_genbank}
\alias{write_genbank}
\title{Write a part to GenBank}
\usage{
write_genbank(part, file)
}
\arguments{
\item{part}{A \code{ggr_part} (from \code{\link{read_part}},
\code{\link{golden_gate_assemble}}, or \code{\link{ligate}}).}

\item{file}{Output file path.}
}
\value{
Invisibly returns NULL. Called for its side effect of writing the file.
}
\description{
Serialises a \code{ggr_part} as a standards-compliant GenBank flat file.
LOCUS topology, feature locations, qualifier strings, and the ORIGIN
sequence block are all written.
}
\details{
The output round-trips through \code{\link{read_part}} without sequence
loss. Features that span the origin in a circular plasmid are emitted
as \code{join()} locations.

Numeric qualifier values (e.g. \code{/codon_start=1}) are written
unquoted; everything else is quoted.
}
\examples{
\dontrun{
  asm <- golden_gate_assemble(list(vec, p1, p2), enzyme = "BsaI")
  write_genbank(asm, "construct.gb")
}
}
