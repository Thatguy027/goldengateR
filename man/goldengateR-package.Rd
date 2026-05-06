\name{goldengateR-package}
\alias{goldengateR}
\alias{goldengateR-package}
\docType{package}
\title{In silico Golden Gate assembly with BsaI and BsmBI}
\description{
Simulates Golden Gate cloning of plasmids and PCR products using the
Type IIS restriction enzymes BsaI and BsmBI. Reads GenBank or FASTA
inputs, digests each input at all enzyme sites, attempts to ligate the
resulting fragments into a single circular product via 4 nt sticky-end
matching, and writes the assembled product as an annotated GenBank file.
}
\section{Two assembly modes}{
\describe{
  \item{One-pot}{\code{\link{golden_gate_assemble}} digests and ligates
  all parts in a single call. The right tool when every input is a
  circular plasmid with the recognition sites in the standard
  inward-facing orientation.}

  \item{Sequential}{\code{\link{digest_and_select}} simulates a
  digest-and-gel-purify step, returning a single fragment. Combined with
  \code{\link{ligate}}, this models workflows where you cut a vector,
  gel-purify the backbone, and ligate it with a separately-prepared
  PCR product or annealed oligo.}
}
}
\section{Plasmid library}{
The package bundles example plasmids in \code{inst/extdata/plasmids/},
accessible by name via \code{\link{load_plasmid}}. See
\code{\link{set_plasmid_dir}} for using a personal collection.
}
\section{Key functions}{
\describe{
  \item{\code{\link{read_part}}, \code{\link{write_genbank}}}{file I/O}
  \item{\code{\link{load_plasmid}}, \code{\link{list_plasmids}},
  \code{\link{set_plasmid_dir}}, \code{\link{plasmid_search_path}}}{plasmid library}
  \item{\code{\link{golden_gate_assemble}}}{one-pot assembly}
  \item{\code{\link{digest_and_select}}, \code{\link{ligate}},
  \code{\link{make_fragment}}}{sequential workflow}
}
}
\section{Limitations}{
\itemize{
  \item Type IIS only — no support for blunt or sticky Type II enzymes.
  \item Returns the unique circular product or errors. Does not enumerate
  partial products.
  \item Features with \code{join(...)} or \code{order(...)} locations in
  input GenBank files are dropped.
}
}
\keyword{package}
