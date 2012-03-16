\name{countPointsInPolys}
\alias{countPointsInPolys}
\title{
Count points within a polygon
}
\description{
 Overlays points on polygons and create a new polygon dataset with the count of the points in that polygon
}
\usage{
countPointsInPolys(points, polys, density = FALSE)
}
\arguments{
  \item{points}{
SpatialPoints
}
  \item{polys}{
SpatialPolygonsDataFrame
}
  \item{density}{
Return a density (point count divided by area) instead of a point count
}
}
\value{
SpatialPolygonsDataFrame
}
\seealso{
See Also as \code{\link[sp]{overlay}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (points, polys, density = FALSE) 
{
    grid_count = overlay(x = points, y = polys)
    pointcount = tapply(grid_count, grid_count, length)
    rownames(pointcount) <- as.character(as.numeric(rownames(pointcount)) - 
        1)
    pointcount.df = data.frame(pointcount, rownames = rownames(pointcount), 
        stringsAsFactors = FALSE)
    polys.df = data.frame(polys@data, rownames = rownames(polys@data), 
        stringsAsFactors = FALSE)
    DF = merge(polys.df, pointcount.df, all = TRUE, by.x = "rownames")
    DF = DF[order(as.numeric(DF$rownames)), ]
    rownames(DF) <- DF$rownames
    DF$pointcount <- as.numeric(DF$pointcount)
    returnSPDF = SpatialPolygonsDataFrame(polygons(polys), data = DF, 
        match.ID = T)
    return(returnSPDF)
  }
}
