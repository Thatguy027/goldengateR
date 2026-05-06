\name{assemble_interactive}
\alias{assemble_interactive}
\title{Interactively assemble a plasmid from the console}
\usage{
assemble_interactive(enzyme = c("BsaI", "BsmBI"), output_dir = ".")
}
\arguments{
\item{enzyme}{Either \code{"BsaI"} or \code{"BsmBI"}.}

\item{output_dir}{Directory to write the output GenBank file. Created if
it does not exist.}
}
\value{
Invisibly, the path to the written GenBank file.
}
\description{
Prompts for an output plasmid name and a space-separated list of part
names, loads each part via \code{\link{load_plasmid}}, and runs
\code{\link{golden_gate_assemble}}.
}
\details{
Two prompts are shown at the console:

\enumerate{
  \item \strong{Output plasmid name} — used as the LOCUS name and the
  output filename (\code{<name>.gb}).
  \item \strong{Parts} — space-separated plasmid names, each resolved
  via \code{\link{load_plasmid}}. All names visible on the
  \code{\link{plasmid_search_path}} are valid, including bundled pYTK
  parts and anything in a directory set with \code{\link{set_plasmid_dir}}.
}

The assembly logic is identical to \code{\link{golden_gate_assemble}};
the order in which parts are listed does not matter — the solver derives
join order from 4 nt sticky-end matching.
}
\examples{
\dontrun{
  # Prompts:
  #   Output plasmid name: xyl1_cassette
  #   Parts (space-separated): pYTK083 pYTK002 pYTK009 pYTK_T3_XYL1_K270R l6_cyc1_type4 pYTK067 pYTK077 pYTK081
  assemble_interactive(enzyme = "BsaI", output_dir = "./cassettes")
}
}
\seealso{\code{\link{assemble_from_table}}, \code{\link{golden_gate_assemble}}, \code{\link{load_plasmid}}}
