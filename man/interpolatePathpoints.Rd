\name{interpolatePathpoints}
\alias{interpolatePathpoints}
\title{
Interpolate points along a path
}
\description{
Interpolate points along a path
}
\usage{
interpolatePathpoints(pathpoints, dens, tolerance.min = 1.2, tolerance.max = 50)
}
\arguments{
  \item{pathpoints}{
Path points as they currently exist
}
  \item{dens}{
inverse density and is in the units of the x and y in pathpoints (e.g. 1 point per density meters)
}
  \item{tolerance.min}{
The proportion of the density (e.g. 1.2 means we'll fill in gaps 20% greater than the density size)
}
  \item{tolerance.max}{
Max tolerance (see \code{tolerance.min} )
}
}
\value{
path with points interpolated
}
\seealso{
See Also \code{\link{reshapeSLDF}}, \code{\link{SLDFtoLine}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (pathpoints, dens, tolerance.min = 1.2, tolerance.max = 50) 
{
    pathpoints$xNext = c(pathpoints$x[2:length(pathpoints$x)], 
        NA)
    pathpoints$yNext = c(pathpoints$y[2:length(pathpoints$y)], 
        NA)
    pathpoints$distToNext = sqrt((pathpoints$x - pathpoints$xNext)^2 + 
        (pathpoints$y - pathpoints$yNext)^2)
    interpFlag = rep(0, length(pathpoints$x))
    interpFlag[pathpoints$distToNext > dens * tolerance.min & 
        pathpoints$distToNext < dens * tolerance.max] = 1
    interpolate = rep(0, length(interpFlag))
    interpolate[2:length(interpFlag)] = interpFlag[1:length(interpFlag) - 
        1]
    pathpoints$interpolateGroup = rep(1, length(pathpoints$x))
    for (rownum in 2:length(interpolate)) {
        pathpoints$interpolateGroup[rownum] = pathpoints$interpolateGroup[rownum - 
            1] + interpolate[rownum]
    }
    numIGroups = max(pathpoints$interpolateGroup)
    bypp = by(pathpoints, pathpoints$interpolateGroup, function(pp) {
        if (pp$interpolateGroup[1] != numIGroups) {
            numRows = length(pp[[1]])
            pointsToAdd = as.integer(floor(pp$distToNext[numRows]/dens))
            xIncrement = ((pp$xNext - pp$x)/(pointsToAdd + 1))[numRows]
            yIncrement = ((pp$yNext - pp$y)/(pointsToAdd + 1))[numRows]
            pp <- expandDF(pp, numRows, pointsToAdd)
            for (pointNum in 1:pointsToAdd) {
                pp[numRows + pointNum, "x"] = pp[numRows + pointNum - 
                  1, "x"] + xIncrement
                pp[numRows + pointNum, "y"] = pp[numRows + pointNum - 
                  1, "y"] + yIncrement
            }
        }
        return(pp)
    })
    returnpp = bypp[[1]]
    for (igroup in 2:numIGroups) {
        returnpp = rbind(returnpp, bypp[[igroup]])
    }
    return(returnpp)
  }
}
