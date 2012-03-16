\name{as.data.frame.by}
\alias{as.data.frame.by}
\title{
Convert the results of by() to a data.frame.
}
\description{
Converts the results of by() to a data.frame if possible
}
\usage{
\method{as.data.frame}{by}(x, row.names = NULL, optional = FALSE, colnames = paste("IDX", seq(length(dim(x))), sep = ""), na.rm = TRUE, ...)
}
\arguments{
  \item{x}{
	The by object
  }
  \item{row.names}{
	Names of the rows.  If NULL, function tries guessing them
}
  \item{optional}{
	Ignored.
}
  \item{colnames}{
	Names of columns
}
  \item{na.rm}{
	Remove NAs or not.
}
  \item{\dots}{
	Pass-alongs.
}
}
\value{
	A data.frame.
}
