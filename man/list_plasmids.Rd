\name{list_plasmids}
\alias{list_plasmids}
\title{List bundled and user-library plasmids}
\usage{
list_plasmids(pattern = NULL)
}
\arguments{
\item{pattern}{Optional regex (case-insensitive) to filter plasmid names.}
}
\value{
A data.frame with columns \code{name}, \code{file}, \code{path}, and
\code{source} (\code{"user"} or \code{"bundled"}). Returns an empty
data.frame if nothing is found.
}
\description{
Returns the catalogue of plasmid files visible on the search path
(\code{\link{plasmid_search_path}}). Use this to discover what's
available before calling \code{\link{load_plasmid}}.
}
\examples{
\dontrun{
  list_plasmids()                     # everything
  list_plasmids(pattern = "pYTK0")    # just the integration vectors
  list_plasmids(pattern = "^p[a-z]")  # lowercase-p plasmids
}
}
\seealso{\code{\link{load_plasmid}}, \code{\link{plasmid_search_path}}}
