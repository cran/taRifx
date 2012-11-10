\name{iapply}
\alias{iapply}
\title{Iteratively (recursively) apply a function to its own output}
\usage{
  iapply(X, FUN, init, ...)
}
\arguments{
  \item{X}{a vector of first arguments to be passed in}

  \item{FUN}{a function taking a changing (x) and an
  initial argument (init)}

  \item{init}{an argument to be "worked on" by FUN with
  parameters x[1], x[2], etc.}

  \item{\dots}{Arguments passed to FUN.}
}
\value{
  the final value, of the same type as init
}
\description{
  Iteratively (recursively) apply a function to its own
  output
}
\examples{
vec <- "xy12"
mylist <- list( c("x","a"), c("y","b"), c("a","f") )
iapply( mylist , FUN=function(repvec,x) {
  gsub(repvec[1],repvec[2],x)
}, init=vec )
}

