\name{cumDist}
\alias{cumDist}
\title{
Calculate cumulative distance along a matrix of x,y coordinates
}
\description{
Calculate cumulative distance along a matrix of x,y coordinates
}
\usage{
cumDist(coords)
}
\arguments{
  \item{coords}{
an [n,2] matrix of coordinates
}
}
\value{
	A single value that represents the distance along the path
}
\seealso{
See Also \code{\link{reshapeSLDF}}, \code{\link{SLDFtoLine}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (coords) 
{
    if (any(is.na(coords))) {
        return(NA)
    }
    if (dim(coords)[1] == 2) {
        return(simpledist(rbind(coords[1, ], coords[2, ])))
    }
    prevCoords = coords[1:dim(coords)[1] - 1, ]
    currCoords = coords[2:dim(coords)[1], ]
    distToPrev = rep(NA, dim(coords)[1])
    for (rowNum in 2:dim(coords)[1]) {
        distToPrev[rowNum] = simpledist(rbind(prevCoords[(rowNum - 
            1), ], currCoords[(rowNum - 1), ]))
    }
    return(sum(distToPrev, na.rm = TRUE))
  }
}
