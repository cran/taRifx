\name{simpledist}
\alias{simpledist}
\alias{closestPoint}
\title{
Cartesian distance between points
}
\description{
Cartesian distance between points (simpledist), or Find closest point to a given point's coordinates (closestPoint).
}
\usage{
simpledist(points)
closestPoint(point,points)
}
\arguments{
  \item{points}{
points is a 2x2 matrix, where columns are x,y and rows are the two points
}
\item{point}{
the point from which you want to find the closest match in points
}
}
\value{
Distance, in same units as input
}
\examples{
points <- matrix(c(0,3,0,4),nrow=2)
simpledist(points)
}

