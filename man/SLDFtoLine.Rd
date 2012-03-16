\name{SLDFtoLine}
\alias{SLDFtoLine}
\title{
 Convert a SpatialLinesDataFrame to a single line matrix with associated segment information
}
\description{
 Convert a SpatialLinesDataFrame to a single line matrix with associated segment information
}
\usage{
SLDFtoLine(lineDF, orderXY = FALSE, segments = TRUE)
}
\arguments{
  \item{lineDF}{
SpatialLinesDataFrame object
}
  \item{orderXY}{
%%     ~~Describe \code{orderXY} here~~
}
  \item{segments}{
%%     ~~Describe \code{segments} here~~
}
}
\value{
data.frame
}
\seealso{
See Also \code{\link{reshapeSLDF}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (lineDF, orderXY = F, segments = T) 
{
    X = c()
    Y = c()
    Segment = c()
    numsegments = length(lineDF@lines)
    for (seg in 1:numsegments) {
        xnew = lineDF@lines[[seg]]@Lines[[1]]@coords[, 1]
        ynew = lineDF@lines[[seg]]@Lines[[1]]@coords[, 2]
        X = append(X, xnew)
        Y = append(Y, ynew)
        if (segments) {
            Segment = append(Segment, rep.int(seg, length(xnew)))
        }
        else {
            Segment = append(Segment, rep.int(1, length(xnew)))
        }
    }
    if (orderXY) {
        Y = Y[order(X)]
        Segment = Segment[order(X)]
        X = X[order(X)]
    }
    line = data.frame(X, Y, Segment)
    line
  }
}
