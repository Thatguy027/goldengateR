\name{load_plasmid}
\alias{load_plasmid}
\title{Load a plasmid by name from the library}
\usage{
load_plasmid(name, topology = "circular")
}
\arguments{
\item{name}{Plasmid name (the filename without extension), case-insensitive.
For example \code{"pYTK001"} loads either \code{pYTK001.gb} or
\code{pYTK001.fa} from the search path.}

\item{topology}{Default topology if loading from FASTA. Ignored for
GenBank.}
}
\value{
A \code{ggr_part}.
}
\description{
Convenience lookup for plasmids stored in the package's bundled
\code{inst/extdata/plasmids/} directory or in a user-configured
collection (see \code{\link{set_plasmid_dir}}).
}
\details{
The search path is, in order:
\enumerate{
  \item Any directory set via \code{\link{set_plasmid_dir}} or
  \code{options(goldengateR.plasmid_dir = ...)}
  \item The bundled \code{inst/extdata/plasmids/} directory
}

Earlier directories take precedence, so a user-directory plasmid
shadows a bundled plasmid of the same name. This is the intended
mechanism for overriding bundled templates with lab-specific versions.

If the package isn't installed (i.e., development workflow with
\code{source()}-ing the R/ files), a fallback search locates the
\code{inst/extdata/plasmids/} directory by walking up from
\code{getwd()} or via the \code{goldengateR.pkg_dir} option.
}
\examples{
\dontrun{
  # See what's available
  list_plasmids()

  # Load by name
  vec <- load_plasmid("pYTK_T3_dest_demo")

  # Or use a personal collection
  set_plasmid_dir("~/lab_plasmids/")
  vec <- load_plasmid("my_custom_vector")
}
}
\seealso{\code{\link{list_plasmids}}, \code{\link{set_plasmid_dir}},
\code{\link{plasmid_search_path}}}
