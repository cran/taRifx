\name{shift}
\alias{shift}
\title{
Shifts a vector's elements left or right by N elements.
}
\description{
Shifts a vector's elements left or right by N elements.
}
\usage{
shift(vec, n = 1, wrap = TRUE, pad = FALSE)
}
\arguments{
  \item{vec}{
  	  A vector to be operated on
  }
  \item{n}{
  	  Number of rows to shift by
  }
  \item{wrap}{
  	  Whether to wrap elements or not
  }
  \item{pad}{
  	  Whether to pad with NAs or not.   pad does nothing unless wrap is false, in which case it specifies whether to pad with NAs
  }
}

\value{
vector of the same type as vec
}
\examples{
test <- seq(10)
shift(test)
}
