\name{pointDistPairwise}
\alias{pointDistPairwise}
\title{
 Create all pairwise distances of points from a SpatialPointsDataFrame
}
\description{
 Create all pairwise distances of points from a SpatialPointsDataFrame
}
\usage{
pointDistPairwise(SPDF, names = "name")
}
\arguments{
  \item{SPDF}{
SpatialPointsDataFrame
}
  \item{names}{
variable name in the SPDF's dataframe used to label each point in the resulting matrix
}
}
\value{
matrix
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (SPDF, names = "name") 
{
    crds <- SPDF@coords
    numpoints <- dim(crds)[1]
    pw.mat <- matrix(rep(NaN, numpoints^2), ncol = numpoints)
    for (i in seq(numpoints)) {
        for (j in seq(numpoints)) {
            if (i >= j) {
                pw.mat[i, j] <- simpledist(rbind(crds[i, ], crds[j, 
                  ]))
            }
        }
    }
    colnames(pw.mat) <- SPDF@data[[names]]
    rownames(pw.mat) <- SPDF@data[[names]]
    return(pw.mat)
  }
}

