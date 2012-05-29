\name{subsetSPDF}
\alias{subsetSPDF}
\alias{subset.SpatialPolygonsDataFrame}
\alias{subset.SpatialPointsDataFrame}
\title{
Subset SpatialPolygonsDataFrame or SpatialPointsDataFrame
}
\description{
Subset SpatialPolygonsDataFrame or SpatialPointsDataFrame
}
\usage{
subsetSPDF(x, tf, ...)
\method{subset}{SpatialPolygonsDataFrame}(x, tf, ...)
\method{subset}{SpatialPointsDataFrame}(x, tf, ...)
}
\arguments{
  \item{x}{
  	  SpatialPolygonsDataFrame or SpatialPointsDataFrame
}
  \item{tf}{
Boolean on which to subset
}
  \item{\dots}{
Arguments to pass to SpatialPolygonsDataFrame() or SpatialPointsDataFrame() when reconstructing objects
}
}
\value{SpatialPolygonsDataFrame or SpatialPointsDataFrame
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (SPDF, tf, ...) 
{
    selected_data <- subset(SPDF@data, tf)
    if (class(SPDF) == "SpatialPolygonsDataFrame") {
        SPDF_selected <- subset(SPDF@polygons, tf)
        centroids <- getSpPPolygonsLabptSlots(as.SpatialPolygons.PolygonsList(SPDF_selected))
        x <- centroids[, 1]
        y <- centroids[, 2]
        export <- SpatialPolygonsDataFrame(as.SpatialPolygons.PolygonsList(SPDF_selected), 
            data = data.frame(x = x, y = y, row.names = getSpPPolygonsIDSlots(as.SpatialPolygons.PolygonsList(SPDF_selected)), 
                selected_data), ...)
    }
    if (class(SPDF) == "SpatialPointsDataFrame") {
        export <- SpatialPointsDataFrame(coordinates(SPDF)[tf, 
            ], data = selected_data, coords.nrs = SPDF@coords.nrs, 
            proj4string = SPDF@proj4string, ...)
    }
    return(export)
  }
}
