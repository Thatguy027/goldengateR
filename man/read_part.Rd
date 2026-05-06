\name{read_part}
\alias{read_part}
\title{Read a part from GenBank or FASTA}
\usage{
read_part(file, topology = "circular")
}
\arguments{
\item{file}{Path to a \code{.gb}/\code{.gbk}/\code{.genbank} or
\code{.fa}/\code{.fasta}/\code{.fna} file.}

\item{topology}{For FASTA only, either \code{"circular"} or \code{"linear"}.
Ignored for GenBank, which declares its own topology in the LOCUS line.}
}
\value{
A \code{ggr_part} object: an S3 list with components \code{sequence}
(character), \code{features} (data.frame), \code{topology}, and \code{name}.
}
\description{
Loads a DNA sequence with annotations into the package's internal
\code{ggr_part} representation. Used as the entry point for any plasmid
or PCR product before downstream digestion or assembly.
}
\details{
The GenBank parser handles standard NCBI-style records including
multi-line locations, \code{complement(...)}, and multi-line
\code{/qualifier} continuations. \code{join(...)} and \code{order(...)}
locations are dropped (with no warning). The format is auto-detected
from the file extension; if the extension is unrecognised the function
sniffs the first line.

For FASTA, \strong{the topology defaults to \code{"circular"}}. PCR
products and other linear DNA must be loaded with
\code{topology = "linear"} or downstream digestion will compute incorrect
cut positions.
}
\examples{
\dontrun{
  vec <- read_part("pYTK001.gb")                          # GenBank
  amp <- read_part("xylB_amplicon.fa", topology = "linear")  # FASTA
}
}
\seealso{\code{\link{load_plasmid}}, \code{\link{write_genbank}}}
