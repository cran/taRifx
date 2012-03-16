\name{SPDFtoPointsDF}
\alias{SPDFtoPointsDF}
\title{
Convert SpatialPointsDataFrame to a regular data.frame with the coordinates as "x" and "y"
}
\description{
Convert SpatialPointsDataFrame to a regular data.frame with the coordinates as "x" and "y"
}
\usage{
SPDFtoPointsDF(SPDF)
}
\arguments{
  \item{SPDF}{
SpatialPointsDataFrame
}
}
\value{data.frame with the coordinates as "x" and "y"
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (SPDF) 
{
    coords <- SPDF@coords
    colnames(coords) <- c("x", "y")
    return(cbind(SPDF@data, coords))
  }
}
