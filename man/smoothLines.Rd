\name{smoothLines}
\alias{smoothLines}
\title{
Smooth the line segments in a SpatialLinesDataFrame
}
\description{
Smooth the line segments in a SpatialLinesDataFrame.  Takes a jittery SLDF (such as a GPS produces) and smooths it.
}
\usage{
smoothLines(lineDF)
}
\arguments{
  \item{lineDF}{
SpatialLinesDataFrame
}
}
\value{
	SpatialLinesDataFrame
}
\seealso{
See Also as \code{\link{SLDFtoLine}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (lineDF) 
{
    numsegments = length(lineDF@lines)
    for (seg in 1:numsegments) {
        plot(x, y, main = paste(seg))
        x = lineDF@lines[[seg]]@Lines[[1]]@coords[, 1][order(lineDF@lines[[seg]]@Lines[[1]]@coords[, 
            1])]
        y = lineDF@lines[[seg]]@Lines[[1]]@coords[, 2][order(lineDF@lines[[seg]]@Lines[[1]]@coords[, 
            1])]
        if (length(x) > 4) {
            fit <- smooth.Pspline(x, y, method = 3)
        }
        fit
    }
  }
}
