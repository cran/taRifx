\name{last}
\alias{last}
\alias{first}
\title{
Convenience functions to return the last/first element of a vector}
\description{
Convenience functions to return the last/first element of a vector
}
\usage{
last(vec)
first(vec)
}
\arguments{
  \item{vec}{
  	  Vector of any type
  }
}
\value{
Vector of length 1 of same type as vec
}

\examples{
test <- seq(10)
first(test)
last(test)
}
