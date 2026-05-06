\name{set_plasmid_dir}
\alias{set_plasmid_dir}
\title{Set the user plasmid directory}
\usage{
set_plasmid_dir(dir)
}
\arguments{
\item{dir}{Path to a directory containing GenBank or FASTA plasmid files,
or \code{NULL} to clear the setting.}
}
\value{
Invisibly returns the directory path (or \code{NULL} if cleared).
}
\description{
Convenience wrapper around \code{options(goldengateR.plasmid_dir = dir)}.
After calling this, plasmids in the supplied directory are findable via
\code{\link{load_plasmid}} and take precedence over bundled plasmids of
the same name.
}
\details{
This is the right place to point at a personal lab plasmid collection
that you don't want to commit into the package itself. The setting
persists for the current R session only; for a permanent default, set
\code{options(goldengateR.plasmid_dir = "...")} in your \code{.Rprofile}.
}
\examples{
\dontrun{
  set_plasmid_dir("~/lab_plasmids/")
  list_plasmids()

  set_plasmid_dir(NULL)   # back to bundled-only
}
}
\seealso{\code{\link{load_plasmid}}, \code{\link{list_plasmids}}}
