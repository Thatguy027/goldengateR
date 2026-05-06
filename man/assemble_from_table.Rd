\name{assemble_from_table}
\alias{assemble_from_table}
\title{Assemble multiple plasmids from a manifest table}
\usage{
assemble_from_table(file, output_dir = ".", enzyme = c("BsaI", "BsmBI"))
}
\arguments{
\item{file}{Path to a CSV or TSV manifest file. Tab-separated files
(\code{.tsv}) are detected by extension; everything else is treated as
comma-separated.}

\item{output_dir}{Directory for output GenBank files. Created if it does
not exist.}

\item{enzyme}{Either \code{"BsaI"} or \code{"BsmBI"}.}
}
\value{
Invisibly, a data frame with one row per assembly and columns:
\describe{
  \item{name}{Output plasmid name.}
  \item{status}{\code{"OK"}, \code{"FAILED"}, or \code{"SKIPPED"}.}
  \item{output}{Path to the written GenBank file, or \code{NA} on failure.}
  \item{error}{Error message for failed rows, otherwise \code{NA}.}
}
}
\description{
Reads a CSV or TSV manifest where each row describes one Golden Gate
assembly, loads the specified parts via \code{\link{load_plasmid}}, and
runs \code{\link{golden_gate_assemble}} for each row. A summary is
printed on completion.
}
\details{
\strong{Required column:}
\describe{
  \item{name}{Output plasmid name. Used as the GenBank LOCUS field and
  the output filename (\code{<name>.gb}). Rows with a blank name are
  skipped.}
}

\strong{Part columns} (all optional; empty cells are skipped):

\code{vec}, \code{1}, \code{2}, \code{3}, \code{3a}, \code{3b},
\code{234}, \code{4}, \code{4a}, \code{4b}, \code{5}, \code{6},
\code{7}, \code{8}, \code{8a}, \code{8b}, \code{678}

Each cell holds a plasmid name resolvable via \code{\link{load_plasmid}}.
The column a part appears in does not affect the assembly — join order is
determined by 4 nt sticky-end matching. Column names exist only to keep
the manifest readable.

Per-row failures are reported but do not abort the batch. The function
always attempts every row.

An example manifest (\code{assemblies_example.tsv}) is included at the
package root.
}
\examples{
\dontrun{
  # Run all assemblies in the example manifest
  results <- assemble_from_table(
    system.file("assemblies_example.tsv", package = "goldengateR"),
    output_dir = "./cassettes"
  )

  # Chain assemblies: point set_plasmid_dir at the output folder so
  # newly-built plasmids are immediately available for the next round
  set_plasmid_dir("./cassettes")
  results2 <- assemble_from_table("level2_assemblies.tsv",
                                  output_dir = "./cassettes")
}
}
\seealso{\code{\link{assemble_interactive}}, \code{\link{golden_gate_assemble}}, \code{\link{set_plasmid_dir}}}
