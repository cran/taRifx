\name{reshapeSLDF}
\alias{reshapeSLDF}
\title{
Reshape a spatialLinesDataFrame into a series of points with associated information
}
\description{
Reshape a spatialLinesDataFrame into a series of points with associated information (less efficient because all the segment data gets replicated over each point)
}
\usage{
reshapeSLDF(SLDF, shape = "long")
}
\arguments{
  \item{SLDF}{
spatialLinesDataFrame
}
  \item{shape}{
Do not change.  Must be "long".  For forward compatility.
}
}
\value{
data.frame
}
\seealso{
See Also \code{\link{SLDFtoLine}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (SLDF, shape = "long") 
{
    if (shape == "long") {
        for (seg in 1:length(SLDF@lines)) {
            numRecords = length(SLDF@lines[[seg]]@Lines[[1]]@coords[, 
                1])
            longit = SLDF@lines[[seg]]@Lines[[1]]@coords[, 1]
            latit = SLDF@lines[[seg]]@Lines[[1]]@coords[, 2]
            SLDFdata = as.data.frame(replicate(length(SLDF@data), 
                rep(NA, numRecords)))
            colnames(SLDFdata) <- colnames(SLDF@data)
            for (datacolnum in 1:length(SLDF@data)) {
                SLDFdata[[datacolnum]] = rep(SLDF@data[[datacolnum]][seg], 
                  numRecords)
            }
            if (seg == 1) {
                SLDFpoints = cbind(longit, latit, SLDFdata)
            }
            else {
                SLDFpoints = rbind(SLDFpoints, cbind(longit, 
                  latit, SLDFdata))
            }
        }
        return(SLDFpoints)
    }
  }
}
