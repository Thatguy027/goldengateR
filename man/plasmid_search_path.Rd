\name{plasmid_search_path}
\alias{plasmid_search_path}
\title{Get the search path for plasmid lookups}
\usage{
plasmid_search_path()
}
\value{
Character vector of directory paths, in priority order (earliest first).
Returns an empty vector if no plasmid directory is configured and the
bundled directory cannot be located.
}
\description{
Returns the directories that \code{\link{load_plasmid}} and
\code{\link{list_plasmids}} will search. Useful for debugging when
\code{load_plasmid} can't find a plasmid you expected.
}
\examples{
\dontrun{
  plasmid_search_path()
}
}
\seealso{\code{\link{set_plasmid_dir}}}
