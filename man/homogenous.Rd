\name{homogenous}
\alias{homogenous}
\title{
Returns whether a vector is homogenous or not
}
\description{
Returns TRUE/FALSE if every element of vector is identical/not.
}
\usage{
homogenous(vec)
}
\arguments{
  \item{vec}{
  	  Vector to be compared
  }
}
\value{
TRUE if every element of a vector is identical; FALSE otherwise.
}

\seealso{
See also \code{\link{all}} \code{\link{any}}
}
\examples{
homogenous(c(rep("A",10),"A"))
homogenous(c(rep("A",10),"B"))
}
