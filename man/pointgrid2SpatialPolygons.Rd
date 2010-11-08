\name{pointgrid2SpatialPolygons}
\alias{pointgrid2SpatialPolygons}
\title{
Take a grid of regularly spaced points (such as those output by the centroids of Arc's fishnet function) and convert it to various grid data types
}
\description{
Take a grid of regularly spaced points (such as those output by the centroids of Arc's fishnet function) and convert it to various grid data types
}
\usage{
pointgrid2SpatialPolygons(df, type)
}
\arguments{
  \item{df}{
SpatialPointsDataFrame
}
  \item{type}{
"SpatialGrid" or "SpatialPolygons" or "SpatialPolygonsDataFrame"
}
}
\value{
SpatialGrid, SpatialPolygons, SpatialPolygonsDataFrame
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (df, type) 
{
    if (class(df)[[1]] != "SpatialPointsDataFrame" & class(df)[[1]] != 
        "SpatialPixelsDataFrame") {
        return(-99)
    }
    if (type == "SpatialGrid") {
        return(SpatialGrid(points2grid(df)))
    }
    df_sp = as.SpatialPolygons.GridTopology(points2grid(df))
    if (type == "SpatialPolygons") {
        return(df_sp)
    }
    if (type == "SpatialPolygonsDataFrame") {
        data = as.data.frame(rep(0, length(df_sp@plotOrder)))
        rownames(data) <- getSpPPolygonsIDSlots(df_sp)
        return(SpatialPolygonsDataFrame(df_sp, data))
    }
    return(-98)
  }
}
