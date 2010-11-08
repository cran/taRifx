\name{bytable}
\alias{bytable}
\alias{as.matrix.by}
\title{
Produces a nice summary table by groupings
}
\description{
bytable produces a nice summary table by groupings, suitable for use with latex.table.by().
as.matrix.by converts a by object into a matrix (only tested on a 2d objects).  Use as.matrix(by(...)) .
}
\usage{
bytable(datavec, indices, ops = c(quote(mean)), ops.desc = list(mean = "Mean"), na.rm = TRUE)
\method{as.matrix}{by}(x,...)
}
\arguments{
  \item{datavec}{
  	  Vector to be analyzed
  }
  \item{indices}{
 Indices should be a list of grouping vectors, just like you would pass to -by-, but with sensible names for each vector
}
  \item{ops}{
  	  Vector of quote'd operations to perform
  }
  \item{ops.desc}{
  	  Vector of length length(ops) containing the column labels for the operations.
  }
  \item{na.rm}{
  	  Remove NAs or not
}
\item{x}{
  	  by() object
}
\item{\dots}{
	other arguments to pass to by
}
}
\value{
data.frame
}
\seealso{
latex.table.by
}
\examples{
bytable(runif(100),indices=list(rep(c('a','b'),50)))
}
