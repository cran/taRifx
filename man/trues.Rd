\name{trues}
\alias{trues}
\title{
Return vector of equal length containing all TRUEs
}
\description{
Takes a vector and returns a vector of equal length containing all trues (used for selecting all of a given vector)}
\usage{
trues(vec)
}

\arguments{
  \item{vec}{
  	  \code{vec} is any vector (or valid object for \code{length} )
}
}
\value{
Returns a vector of TRUEs of the length of the object passed to it
}

\examples{
x <- runif(100)
trues(x)
}
\keyword{ boolean }
\keyword{ TRUE }
