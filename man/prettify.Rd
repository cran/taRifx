\name{prettify}
\alias{prettify}
\title{Function to prettify the output of another function using a `var.labels` attribute
This is particularly useful in combination with read.dta et al.}
\usage{
  prettify(dat, expr)
}
\arguments{
  \item{dat}{A data.frame with attr `var.labels` giving
  descriptions of variables}

  \item{expr}{An expression to evaluate with pretty
  var.labels}
}
\value{
  The result of the expression, with variable names
  replaced with their labels
}
\description{
  Function to prettify the output of another function using
  a `var.labels` attribute This is particularly useful in
  combination with read.dta et al.
}
\examples{
testDF <- data.frame( a=seq(10),b=runif(10),c=rnorm(10) )
attr(testDF,"var.labels") <- c("Identifier","Important Data","Lies, Damn Lies, Statistics")
prettify( testDF, quote(str(dat)) )
}

