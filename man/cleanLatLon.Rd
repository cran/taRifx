\name{cleanLatLon}
\alias{cleanLatLon}
\title{
Standardize latitude/longitude coordinates
}
\description{
Standardize latitude/longitude coordinates
}
\usage{
cleanLatLon(vec)
}
\arguments{
  \item{vec}{
Character vector of lat/long coordinates
}
}
\value{
Numeric decimal lat/lon
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (vec) 
{
    if (!is.character(vec)) 
        stop("Input vector must be of type character")
    vec = sub("[NnEe]", "", vec)
    vec = sub("[WwSs]", "-", vec)
    vec = sub("^ +", "", vec)
    vec = sub(" +$", "", vec)
    vec = sub("^0+", "", vec)
    if (length(grep(" ", vec)) != 0) 
        stop("Your vector appears to be in decimal degrees, which are not yet implemented")
    return(as.numeric(vec))
  }
}

