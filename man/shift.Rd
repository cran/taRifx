\name{shift}
\alias{shift}
\alias{shift.default}
\alias{shift.data.frame}
\title{
Shifts a vector's elements left or right by N elements.
}
\description{
Shifts a vector's elements left or right by N elements.
}
\usage{
   \method{shift}{default}(x, n = 1, wrap = TRUE, pad = FALSE,...)
   \method{shift}{data.frame}(x,...)
}
\arguments{
  \item{x}{
  	  A vector to be operated on
  }
  \item{n}{
  	  Number of rows to shift by (if negative, shift to right instead of left)
  }
  \item{wrap}{
  	  Whether to wrap elements or not
  }
  \item{pad}{
  	  Whether to pad with NAs or not.   pad does nothing unless wrap is false, in which case it specifies whether to pad with NAs
  }
  \item{\dots}{
	Other items to pass along
  }
}

\value{
vector of the same type as vec
}
\examples{
test <- seq(10)
shift(test)
}
