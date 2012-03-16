\name{SPDFareas}
\alias{SPDFareas}
\title{
Return areas of polygons in a SpatialPolygonsDataFrame
}
\description{
Get the areas stored in the polygons and return them in the dataframe slot
}
\usage{
SPDFareas(SPDF, colname = "AREA")
}
\arguments{
  \item{SPDF}{
SpatialPolygonsDataFrame
}
  \item{colname}{
Name of the column in the data frame component of the SpatialPolygonsDataFrame in which to store the polygon areas
}
}
\value{
	SpatialPolygonsDataFrame
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (SPDF, colname = "AREA") 
{
    numPolys = length(SPDF@polygons)
    areas = data.frame(rep(NA, numPolys))
    colnames(areas) <- c(colname)
    for (polyNum in 1:numPolys) {
        areas[[colname]][polyNum] <- SPDF@polygons[[polyNum]]@area
    }
    returnSPDF = SpatialPolygonsDataFrame(polygons(SPDF), data = cbind(SPDF@data, 
        areas), match.ID = TRUE)
    return(returnSPDF)
  }
}
