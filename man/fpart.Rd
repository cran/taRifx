\name{fpart}
\alias{fpart}
\title{
Obtain the fractional part of a numeric
}
\description{
Takes a numeric vector and returns a vector of the numbers after the decimal place
}
\usage{
fpart(vec)
}
\arguments{
  \item{vec}{
  	  A numeric vector of any length
  }
}
\value{
	Returns a vector of the same length as the input \code{vec}, containing
	only the decimal component.
}
\examples{
x <- runif(100)
fpart(x)
}
