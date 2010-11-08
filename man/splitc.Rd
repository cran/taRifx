\name{splitc}
\alias{splitc}

\encoding{UTF-8}

\title{
Split data over columns
}

\description{
Split and run a function column-wise on \code{data.frame}, \code{matrix} and \code{array} or element-wise on a \code{list}.
}

\usage{
splitc(X, INDEX, FUN = NULL, ...)
}

\arguments{
  \item{X}{
A \code{data.frame}, \code{matrix}, \code{array} or a \code{list}.
}
  \item{INDEX}{
A factor of \code{length(X)} (number of columns or list elements). If not a factor, it will be coerced into one.
}
  \item{FUN}{
A function to be applied to individual subset of data. If not provided (\code{NULL}), raw data is returned.
}
  \item{\dots}{
Additional arguments to \code{FUN}.
}
}

\details{
Function splits a \code{data.frame}, \code{matrix} and \code{array} column-wize according to \code{INDEX} and \code{list} is sliced according to \code{INDEX}. Output is returned as a list of the same length as the number of levels in INDEX.
}
\value{
A list of the same length as there are factor levels in \code{INDEX}.
}
 
\author{
Roman Luštrik \email{roman.lustrik@gmail.com}
}
\note{
Simplification sensu \code{tapply} is not yet implemented.
}

\seealso{
\code{\link{tapply}}, \code{\link{by}}, \code{\link{aggregate}}, \code{\link{apply}}
}
\examples{
my.list <- list(a = runif(5), b = runif(5), c = runif(5), d = runif(5), e = runif(10),
		f = runif(10), g = runif(10), h = runif(10), i = runif(10), j = runif(10))
my.df <- as.data.frame(my.list)
my.matrix <- as.matrix(my.df)

ind <- factor(c(1,1,1,1, 2,3, 4,4,4,4))
ind2 <- factor(c(1,1,1,1, 2,3, 4,4,4,4), levels = 1:5)

splitc(X = my.df, INDEX = ind, FUN = mean) # applies mean to each
splitc(X = my.matrix, INDEX = ind2) # level 5 empty because non-existent 
splitc(X = my.list, INDEX = ind, FUN = sum)
}

\keyword{manip}