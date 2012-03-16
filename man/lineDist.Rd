\name{lineDist}
\alias{lineDist}
\title{
Line distance in SpatialLinesDataFrame
}
\description{
Stores length of each line segment in a SpatialLinesDataFrame's data.frame
}
\usage{
lineDist(SLDF, varname = "distances")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SLDF}{
SpatialLinesDataFrame
}
  \item{varname}{
Character string containing name of variable to hold line distances.
}
}
\value{
	Returns a SpatialLinesDataFrame
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (SLDF, varname = "distances") 
{
    numSegments = length(SLDF@lines)
    Dists = rep(NA, numSegments)
    for (segnum in 1:numSegments) {
        Dists[segnum] <- cumDist(SLDF@lines[[segnum]]@Lines[[1]]@coords)
    }
    SLDF@data[[varname]] = Dists
    return(SLDF)
  }
}
