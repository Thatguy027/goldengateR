\name{assemble_from_table}
\alias{assemble_from_table}
\title{Assemble multiple plasmids from a manifest table}
\usage{
assemble_from_table(file, output_dir = ".", enzyme = c("BsaI", "BsmBI"))
}
\arguments{
\item{file}{Path to a CSV or TSV manifest. Files with a \code{.tsv} extension
are read as tab-separated; everything else is treated as comma-separated.}

\item{output_dir}{Directory for output GenBank files. Created if it does not exist.}

\item{enzyme}{Either \code{"BsaI"} or \code{"BsmBI"}.}
}
\value{
Invisibly, a data frame with columns \code{name}, \code{status}
(\code{"OK"}, \code{"FAILED"}, or \code{"SKIPPED"}), \code{output}
(path to the GenBank file, or \code{NA}), and \code{error} (message for
failed rows, otherwise \code{NA}).
}
\description{
Reads a CSV or TSV manifest where each row describes one Golden Gate
assembly. Parts are loaded via \code{\link{load_plasmid}} and assembled
with \code{\link{golden_gate_assemble}}. Per-row failures are reported
but do not abort the batch.
}
\details{
The manifest must have a \code{name} column providing the output plasmid
name for each row. Rows with a blank name are skipped.

Part columns (all optional; empty cells are skipped):

\code{vec}, \code{1}, \code{2}, \code{3}, \code{3a}, \code{3b},
\code{234}, \code{4}, \code{4a}, \code{4b}, \code{5}, \code{6},
\code{7}, \code{8}, \code{8a}, \code{8b}, \code{678}

Each cell holds a plasmid name resolvable via \code{\link{load_plasmid}}.
The column a part appears in does not affect assembly order, which is
determined by 4 nt sticky-end matching.

An example manifest (\code{assemblies_example.tsv}) is included at the
package root.
}
\examples{
\dontrun{
  results <- assemble_from_table("assemblies_example.tsv",
                                 output_dir = "./cassettes")

  # Chain assemblies: newly-built plasmids are immediately available
  # as inputs for the next round.
  set_plasmid_dir("./cassettes")
  results2 <- assemble_from_table("level2_assemblies.tsv",
                                  output_dir = "./cassettes")
}
}
\seealso{
\code{\link{assemble_interactive}}, \code{\link{golden_gate_assemble}}, \code{\link{set_plasmid_dir}}
}
