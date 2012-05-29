\name{autoplot.microbenchmark}

\alias{autoplot.microbenchmark}

\title{autoplot method for microbenchmark objects}

\description{
  Uses ggplot2 to produce a more legible graph of microbenchmark timings
}

\usage{
  \method{autoplot}{microbenchmark}(object, ..., y_max=max(by(object$time,object$expr,uq)) * 1.05 )
}

\arguments{
  \item{object}{A microbenchmark object}
  \item{\dots}{Ignored}
  \item{y_max}{The upper limit of the y axis (defaults to 5 percent more than the maximum value)}
}

\value{
  A ggplot2 plot
}

